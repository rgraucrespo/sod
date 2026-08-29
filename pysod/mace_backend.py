"""MACE evaluation and fixed-cell relaxation on batched SOD configurations.

Ported from the validated prototype ``~/work/nv/batch_sod_relax.py``.  The
benchmark scaffolding (batch-size sweeps, timing tables, the legacy
split-and-rebatch edge path, ASE cross-validation) is deliberately left behind;
only the production routines are kept.

Importing this module pulls in torch, ASE, MACE and the NVIDIA ALCHEMI toolkit.
Callers should catch ImportError and report it as a missing optional dependency.

Two implementation choices are load-bearing and should not be "tidied":

* Edges are written straight into Batch storage by ``NeighborListHook._rebuild``.
  The alternative -- splitting a global neighbour list per graph and rebuilding
  AtomicData -- measured roughly 31x slower on its edge-handling portion.
* Force evaluation must NOT run under ``torch.inference_mode()``: MACE obtains
  conservative forces through autograd, even in eval mode.

This module calls a few private nvalchemi-toolkit APIs, listed in
``PRIVATE_APIS`` below, which is the authoritative record of them:
``check_private_api`` verifies they are all still present before a run starts,
and the documentation points here rather than repeating the names.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import torch
from ase.io import read as ase_read
from ase.io import write as ase_write

from nvalchemi.data import AtomicData, Batch
from nvalchemi.dynamics import FIRE2, FIRE2VariableCell, SizeAwareSampler
from nvalchemi.dynamics._ops.npt_nph import stress_to_cell_force
from nvalchemi.hooks import NeighborListHook
from nvalchemi.models.mace import MACEWrapper

DEFAULT_MODEL = "medium-0b2"
DEFAULT_BATCH_SIZE = 16
DEFAULT_FMAX = 0.05
DEFAULT_MAX_STEPS = 200
DEFAULT_FIRE_DT = 0.05
DEFAULT_FIRE_MAXSTEP = 0.1
DEFAULT_NEIGHBOR_SKIN = 0.3

# 1 eV/A^3 in GPa.
GPA_PER_EV_A3 = 160.21766208

# Private nvalchemi-toolkit APIs this backend depends on.  The toolkit exposes no
# public equivalent for either job -- edges are written straight into Batch
# storage, and the FIRE2 state has to be seeded and synced around a refill -- so
# an upgrade that renames or drops one of these breaks sod_mace.  Checking them
# once at startup turns that into an error naming the missing attribute, instead
# of an AttributeError raised deep inside a relaxation that has already been
# running for minutes.
VALIDATED_NVALCHEMI = "0.2.0"

PRIVATE_APIS: tuple[tuple[str, Any, str], ...] = (
    ("NeighborListHook", NeighborListHook, "_rebuild"),
    ("NeighborListHook", NeighborListHook, "_init_ref_positions"),
    ("FIRE2", FIRE2, "_ensure_state_initialized"),
    ("FIRE2", FIRE2, "_sync_state_to_batch"),
)


def installed_nvalchemi_version() -> str:
    """Installed nvalchemi-toolkit version, or ``"unknown"`` if it cannot be read."""
    try:
        from importlib.metadata import version

        return version("nvalchemi-toolkit")
    except Exception:
        return "unknown"


def check_private_api() -> None:
    """Fail early if a private toolkit API this backend relies on has gone.

    Called by :func:`load_model`, so every run pays a few microseconds for it and
    no run gets a surprise partway through.
    """
    missing = [
        f"{owner}.{attr}" for owner, obj, attr in PRIVATE_APIS if not hasattr(obj, attr)
    ]
    if not missing:
        return
    raise RuntimeError(
        "installed nvalchemi-toolkit is not compatible with sod_mace.\n"
        f"       Missing private API: {', '.join(missing)}\n"
        f"       Installed nvalchemi-toolkit {installed_nvalchemi_version()}, "
        f"validated against {VALIDATED_NVALCHEMI}.\n"
        "       Pin the validated version, or port pysod/mace_backend.py to the "
        "new API."
    )


def _neighbor_hook(
    neighbor_config: Any,
    neighbor_skin: float,
    variable_cell: bool = False,
) -> NeighborListHook:
    """Build the neighbour hook, disabling its position-only cache for moving cells.

    ALCHEMI 0.2.0 invalidates a skin-buffered neighbour list from Cartesian
    atomic displacements only.  A cell change can bring a periodic image inside
    the cutoff without moving either atom, so reusing that list during
    variable-cell relaxation is unsafe.  A zero skin makes every call rebuild
    the list.  Fixed-cell calculations retain the skin optimization.
    """
    skin = 0.0 if variable_cell else neighbor_skin
    return NeighborListHook(neighbor_config, skin=skin)


@dataclass
class Configuration:
    """One SOD configuration: its index, directory and structure."""

    index: int
    name: str
    path: Path
    atoms: Any
    data: AtomicData


@dataclass
class RelaxationResult:
    index: int
    initial_energy: float
    final_energy: float
    steps: int
    final_fmax: float
    converged: bool
    positions: np.ndarray
    # Variable-cell runs only; None under a fixed cell.
    cell: np.ndarray | None = None
    initial_volume: float | None = None
    final_cell_force: float | None = None
    final_pressure: float | None = None


def cueq_available() -> bool:
    """True when the cuEquivariance CUDA operations are importable and registered."""
    try:
        import cuequivariance_ops_torch  # noqa: F401
    except ImportError:
        return False
    return hasattr(torch.ops.cuequivariance, "fused_tensor_product")


def read_configurations(
    config_dirs: list[tuple[int, Path]], structure_name: str
) -> tuple[list[Configuration], list[int]]:
    """Read every configuration's structure file into ASE Atoms and AtomicData.

    Returns the configurations that could be read plus the indices of those that
    could not, so the caller can warn about a sparse result rather than abort --
    matching the ``n_missing`` behaviour of the sod_*_ener.sh collectors.
    """
    configurations: list[Configuration] = []
    missing: list[int] = []
    for index, directory in config_dirs:
        structure = directory / structure_name
        if not structure.is_file():
            missing.append(index)
            continue
        try:
            atoms = ase_read(structure)
        except Exception:  # unreadable or truncated structure file
            missing.append(index)
            continue
        configurations.append(
            Configuration(
                index=index,
                name=directory.name,
                path=directory,
                atoms=atoms,
                data=AtomicData.from_atoms(atoms, device="cpu", dtype=torch.float32),
            )
        )
    return configurations, missing


def load_model(
    device: torch.device,
    model: str | None = None,
    checkpoint: str | Path | None = None,
    enable_cueq: bool = True,
    need_forces: bool = False,
    need_stress: bool = False,
) -> MACEWrapper:
    """Load a MACE foundation model by name, or a local torch-serialised model.

    ``MACEWrapper.from_checkpoint`` does not accept older torch-serialised MACE
    models, so a local checkpoint goes through ``torch.load`` and is wrapped
    directly.  ``weights_only=True`` is deliberately not forced -- those files
    are whole pickled models.
    """
    check_private_api()
    if checkpoint is not None:
        path = Path(checkpoint).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"model checkpoint not found: {path}")
        try:
            raw = torch.load(path, map_location=device, weights_only=False)
        except Exception as exc:
            raise RuntimeError(
                f"could not load {path} as a torch-serialised MACE model.\n"
                "       Supply an ALCHEMI-compatible foundation model with -model instead."
            ) from exc
        wrapper = MACEWrapper(raw.to(device).float().eval())
    else:
        wrapper = MACEWrapper.from_checkpoint(
            model or DEFAULT_MODEL,
            device=device,
            dtype=torch.float32,
            enable_cueq=enable_cueq,
            compile_model=False,
        )
    outputs = {"energy"}
    if need_forces or need_stress:
        outputs.add("forces")
    if need_stress:
        outputs.add("stress")
    wrapper.model_config.active_outputs = outputs
    return wrapper.eval()


def _chunks(items: list[Any], size: int) -> list[list[Any]]:
    return [items[start : start + size] for start in range(0, len(items), size)]


def add_relaxation_fields(batch: Batch, with_stress: bool = False) -> None:
    """Add the node/system fields FIRE2 needs before the first update.

    ``with_stress`` additionally registers a per-graph stress tensor, which
    FIRE2VariableCell reads to derive the force on the cell.
    """
    samples = batch.to_data_list()
    batch.add_key(
        "velocities",
        [torch.zeros_like(sample.positions) for sample in samples],
        level="node",
        overwrite=True,
    )
    batch.add_key(
        "forces",
        [torch.zeros_like(sample.positions) for sample in samples],
        level="node",
        overwrite=True,
    )
    batch.add_key(
        "energy",
        [torch.zeros(1, 1, device=batch.device) for _ in samples],
        level="system",
        overwrite=True,
    )
    if with_stress:
        batch.add_key(
            "stress",
            [torch.zeros(1, 3, 3, device=batch.device) for _ in samples],
            level="system",
            overwrite=True,
        )


def force_max_per_graph(batch: Batch) -> torch.Tensor:
    """Maximum force norm within each graph of the batch."""
    atom_force_norms = torch.linalg.vector_norm(batch.forces, dim=1)
    result = torch.full(
        (batch.num_graphs,), float("-inf"), dtype=batch.forces.dtype, device=batch.device
    )
    result.scatter_reduce_(
        0, batch.batch_idx.long(), atom_force_norms, reduce="amax", include_self=False
    )
    return result


def _rebuild_edges(batch: Batch, hook: NeighborListHook) -> None:
    hook._rebuild(batch)
    if hook.skin > 0.0 and hook._ref_positions is None:
        hook._init_ref_positions(batch.positions)


def evaluate(
    batch: Batch, model: MACEWrapper, hook: NeighborListHook, pressure: float = 0.0
) -> torch.Tensor:
    """Rebuild edges, run MACE, copy energy/forces into the batch, return fmax.

    Under a target pressure the stress stored in the batch is the **effective**
    stress ``sigma + P I``, not the model's stress. That is deliberate:
    ``FIRE2VariableCell.pre_update`` derives the cell force from ``batch.stress``
    itself, so biasing the tensor here is the only way the target pressure
    reaches the dynamics -- passing it to the convergence test alone would leave
    the optimizer relaxing to zero pressure.  It follows from H = E + P V:

        F_cell = -dH/dh = -V (sigma + P I) (h^-1)^T

    Sign check: for sigma = 0, P > 0 and a cubic cell h = a I this gives
    F = -P a^2 I, i.e. compression shrinks the cell, as it must.
    *pressure* is in eV/A^3.
    """
    _rebuild_edges(batch, hook)
    outputs = model(batch)
    batch.energy.copy_(outputs["energy"].view_as(batch.energy))
    batch.forces.copy_(outputs["forces"].view_as(batch.forces))
    if "stress" in outputs and hasattr(batch, "stress"):
        stress = outputs["stress"].view_as(batch.stress)
        if pressure:
            eye = torch.eye(3, device=stress.device, dtype=stress.dtype)
            stress = stress + pressure * eye
        batch.stress.copy_(stress)
    return force_max_per_graph(batch)


def single_point(
    configurations: list[Configuration],
    model: MACEWrapper,
    device: torch.device,
    batch_size: int = DEFAULT_BATCH_SIZE,
    neighbor_skin: float = DEFAULT_NEIGHBOR_SKIN,
) -> dict[int, float]:
    """Batched single-point energies, keyed by SOD configuration index."""
    energies: dict[int, float] = {}
    for chunk in _chunks(configurations, batch_size):
        batch = Batch.from_data_list([item.data for item in chunk], device=device)
        hook = _neighbor_hook(model.model_config.neighbor_config, neighbor_skin)
        _rebuild_edges(batch, hook)
        outputs = model(batch)
        values = outputs["energy"].detach().cpu().flatten().tolist()
        for item, energy in zip(chunk, values, strict=True):
            energies[item.index] = float(energy)
        del batch
    return energies


def warm_up(
    configurations: list[Configuration],
    model: MACEWrapper,
    device: torch.device,
    batch_size: int,
    neighbor_skin: float,
    variable_cell: bool = False,
) -> None:
    """One throwaway evaluation so ragged-shape setup is not timed or measured."""
    chunk = configurations[:batch_size]
    if not chunk:
        return
    batch = Batch.from_data_list([item.data for item in chunk], device=device)
    if "forces" in model.model_config.active_outputs:
        add_relaxation_fields(batch)
    hook = _neighbor_hook(
        model.model_config.neighbor_config, neighbor_skin, variable_cell
    )
    _rebuild_edges(batch, hook)
    model(batch)
    del batch
    if device.type == "cuda":
        torch.cuda.empty_cache()


def cell_force(batch: Batch) -> torch.Tensor:
    """Force conjugate to the cell vectors, in eV/A, from ``batch.stress``."""
    volumes = torch.linalg.det(batch.cell).abs()
    return stress_to_cell_force(batch.stress, batch.cell, volumes)


def hydrostatic_pressure(batch: Batch, pressure: float = 0.0) -> torch.Tensor:
    """Physical pressure per graph in GPa.

    ``batch.stress`` holds the *effective* stress, already shifted by the target
    pressure (see :func:`evaluate`), so the target is added back to recover the
    pressure the structure is actually under.  *pressure* is in eV/A^3.
    """
    trace = batch.stress.diagonal(dim1=-2, dim2=-1).mean(dim=-1)
    return (-trace + pressure) * GPA_PER_EV_A3


def _convergence_measure(
    batch: Batch, variable_cell: bool
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor | None]:
    """Per-graph (measure, atomic fmax, cell-force max), all in eV/A.

    Under a variable cell a graph is converged only when the atomic forces *and*
    the force on the cell are both below the threshold, so the measure is the
    larger of the two -- the same way ASE's UnitCellFilter folds cell gradients
    into a single fmax array, and using the optimizer's own stress conversion
    rather than a second, independent criterion.
    """
    atom_fmax = force_max_per_graph(batch)
    if not variable_cell:
        return atom_fmax, atom_fmax, None
    cell_max = cell_force(batch).abs().amax(dim=(-2, -1))
    return torch.maximum(atom_fmax, cell_max), atom_fmax, cell_max


def _capture(
    batch: Batch,
    active: list[Configuration],
    initial_energies: dict[int, float],
    steps: dict[int, int],
    fmax_values: torch.Tensor,
    converged: bool,
    results: dict[int, RelaxationResult],
    initial_volumes: dict[int, float] | None = None,
    cell_force_values: torch.Tensor | None = None,
    pressure: float = 0.0,
) -> None:
    energies = batch.energy.detach().cpu().flatten().tolist()
    fmax_cpu = fmax_values.detach().cpu().tolist()
    positions = batch.positions.detach().cpu()
    offsets = batch.batch_ptr.detach().cpu().tolist()
    variable_cell = cell_force_values is not None
    if variable_cell:
        cells = batch.cell.detach().cpu().numpy()
        cell_force_cpu = cell_force_values.detach().cpu().tolist()
        pressures = hydrostatic_pressure(batch, pressure).detach().cpu().tolist()
    for local, item in enumerate(active):
        results[item.index] = RelaxationResult(
            index=item.index,
            initial_energy=initial_energies[item.index],
            final_energy=float(energies[local]),
            steps=steps[item.index],
            final_fmax=float(fmax_cpu[local]),
            converged=converged,
            positions=positions[offsets[local] : offsets[local + 1]].numpy(),
            cell=cells[local] if variable_cell else None,
            initial_volume=(initial_volumes or {}).get(item.index),
            final_cell_force=float(cell_force_cpu[local]) if variable_cell else None,
            final_pressure=float(pressures[local]) if variable_cell else None,
        )


def relax_chunk(
    chunk: list[Configuration],
    model: MACEWrapper,
    device: torch.device,
    fmax: float = DEFAULT_FMAX,
    max_steps: int = DEFAULT_MAX_STEPS,
    fire_dt: float = DEFAULT_FIRE_DT,
    fire_maxstep: float = DEFAULT_FIRE_MAXSTEP,
    neighbor_skin: float = DEFAULT_NEIGHBOR_SKIN,
    variable_cell: bool = False,
    pressure: float = 0.0,
) -> dict[int, RelaxationResult]:
    """FIRE2 relaxation of one batch, dropping graphs as they converge.

    With ``variable_cell`` the cell relaxes alongside the positions, driven by
    the model's stress, at target ``pressure`` (eV/A^3).
    """
    batch = Batch.from_data_list([item.data for item in chunk], device=device)
    add_relaxation_fields(batch, with_stress=variable_cell)
    hook = _neighbor_hook(
        model.model_config.neighbor_config, neighbor_skin, variable_cell
    )
    if variable_cell:
        optimizer = FIRE2VariableCell(model, dt=fire_dt, maxstep=fire_maxstep)
    else:
        optimizer = FIRE2(model, dt=fire_dt, maxstep=fire_maxstep)
    optimizer._ensure_state_initialized(batch)

    results: dict[int, RelaxationResult] = {}
    steps = {item.index: 0 for item in chunk}
    active = list(chunk)

    evaluate(batch, model, hook, pressure)
    measure, atom_fmax, cell_max = _convergence_measure(batch, variable_cell)
    initial_energies = {
        item.index: float(energy)
        for item, energy in zip(
            active, batch.energy.detach().cpu().flatten().tolist(), strict=True
        )
    }
    initial_volumes = {
        item.index: float(volume)
        for item, volume in zip(
            active, torch.linalg.det(batch.cell).abs().detach().cpu().tolist(), strict=True
        )
    }

    for step in range(max_steps + 1):
        if step > 0:
            optimizer.pre_update(batch)
            evaluate(batch, model, hook, pressure)
            measure, atom_fmax, cell_max = _convergence_measure(batch, variable_cell)
            for item in active:
                steps[item.index] += 1

        converged_mask = measure <= fmax
        if bool(converged_mask.any()):
            local = torch.nonzero(converged_mask, as_tuple=False).flatten().tolist()
            _capture(
                batch.index_select(local),
                [active[index] for index in local],
                initial_energies,
                steps,
                atom_fmax[converged_mask],
                True,
                results,
                initial_volumes,
                cell_max[converged_mask] if variable_cell else None,
                pressure,
            )
            remaining = torch.nonzero(~converged_mask, as_tuple=False).flatten()
            if remaining.numel() == 0:
                return results
            batch = batch.index_select(remaining)
            # Slice FIRE history for the survivors; never rebuild the optimizer,
            # which would discard velocity history for the other graphs.
            optimizer._sync_state_to_batch(remaining, 0, batch)
            active = [active[index] for index in remaining.tolist()]
            measure = measure[~converged_mask]
            atom_fmax = atom_fmax[~converged_mask]
            if variable_cell:
                cell_max = cell_max[~converged_mask]

    _, atom_fmax, cell_max = _convergence_measure(batch, variable_cell)
    _capture(
        batch,
        active,
        initial_energies,
        steps,
        atom_fmax,
        False,
        results,
        initial_volumes,
        cell_max if variable_cell else None,
        pressure,
    )
    return results


def relax(
    configurations: list[Configuration],
    model: MACEWrapper,
    device: torch.device,
    batch_size: int = DEFAULT_BATCH_SIZE,
    **kwargs: Any,
) -> dict[int, RelaxationResult]:
    """Relaxation of every configuration, in fixed-size chunks."""
    results: dict[int, RelaxationResult] = {}
    for chunk in _chunks(configurations, batch_size):
        results.update(relax_chunk(chunk, model, device, **kwargs))
        if device.type == "cuda":
            torch.cuda.empty_cache()
    return results


# ---------------------------------------------------------------------------
# Cached-refill relaxation
#
# Fixed chunks leave the batch progressively emptier as graphs converge. Refill
# keeps it full by admitting queued structures into the vacated slots, using
# ALCHEMI's native SizeAwareSampler + FIRE2.refill_check.
#
# The catch is that FIRE2.pre_update needs valid forces for a newcomer, so a
# naive refill must evaluate each new graph eagerly -- which costs more than it
# saves (measured ~7-12% slower than fixed chunks). Instead every structure's
# initial energy and forces are computed once up front, in batches, and admitted
# with the newcomer, so no extra model call is needed after a refill.
#
# Identity is carried by an explicit `source_index` system property: batch slots
# are reused, so a graph's position in the batch says nothing about which SOD
# configuration it is.
# ---------------------------------------------------------------------------


@dataclass
class InitialCache:
    """Initial energy, forces and (for variable cell) stress, per configuration."""

    energies: list[float]
    forces: list[torch.Tensor]
    stresses: list[torch.Tensor] | None = None


class CachedDataset:
    """SizeAwareSampler dataset admitting structures that carry cached E0/F0."""

    def __init__(self, configurations: list[Configuration], cache: InitialCache) -> None:
        self.configurations = configurations
        self.cache = cache

    def __len__(self) -> int:
        return len(self.configurations)

    def __getitem__(self, index: int) -> tuple[AtomicData, int]:
        data = self.configurations[index].data.model_copy(deep=True)
        data.add_node_property("velocities", torch.zeros_like(data.positions))
        data.add_node_property("forces", self.cache.forces[index].clone())
        data.add_system_property(
            "energy", torch.tensor([[self.cache.energies[index]]], dtype=torch.float32)
        )
        if self.cache.stresses is not None:
            # Without a cached stress a newly admitted graph would reach
            # FIRE2VariableCell.pre_update with a zero tensor and take one
            # spurious cell step before the next evaluation corrected it.
            data.add_system_property("stress", self.cache.stresses[index].clone())
        data.add_system_property("source_index", torch.tensor([[index]], dtype=torch.long))
        return data, index

    def get_metadata(self, index: int) -> tuple[int, int]:
        return int(self.configurations[index].data.positions.shape[0]), 0


def precompute_cache(
    configurations: list[Configuration],
    model: MACEWrapper,
    device: torch.device,
    batch_size: int = DEFAULT_BATCH_SIZE,
    neighbor_skin: float = DEFAULT_NEIGHBOR_SKIN,
    variable_cell: bool = False,
    pressure: float = 0.0,
) -> InitialCache:
    """Batched evaluation of every configuration's initial energy and forces."""
    energies = [0.0] * len(configurations)
    forces: list[torch.Tensor | None] = [None] * len(configurations)
    stresses: list[torch.Tensor | None] = [None] * len(configurations)
    for offset in range(0, len(configurations), batch_size):
        chunk = configurations[offset : offset + batch_size]
        batch = Batch.from_data_list([item.data for item in chunk], device=device)
        add_relaxation_fields(batch, with_stress=variable_cell)
        hook = _neighbor_hook(
            model.model_config.neighbor_config, neighbor_skin, variable_cell
        )
        evaluate(batch, model, hook, pressure)
        ptr = batch.batch_ptr.detach().cpu().tolist()
        chunk_forces = batch.forces.detach().cpu()
        chunk_stress = batch.stress.detach().cpu() if variable_cell else None
        for local in range(len(chunk)):
            index = offset + local
            energies[index] = float(batch.energy[local, 0].item())
            forces[index] = chunk_forces[ptr[local] : ptr[local + 1]].clone()
            if variable_cell:
                stresses[index] = chunk_stress[local : local + 1].clone()
        del batch
    if any(value is None for value in forces):
        raise RuntimeError("initial-force cache was incomplete")
    return InitialCache(
        energies,
        [value for value in forces if value is not None],
        [value for value in stresses if value is not None] if variable_cell else None,
    )


def _source_ids(batch: Batch) -> list[int]:
    """Dataset positions of the graphs currently in the batch."""
    return batch.source_index.detach().cpu().flatten().tolist()


def _initialize_records(
    batch: Batch,
    initial: dict[int, float],
    steps: dict[int, int],
    volumes: dict[int, float] | None = None,
) -> None:
    """Seed initial energy, volume and step count for graphs seen for the first time."""
    energies = batch.energy.detach().cpu().flatten().tolist()
    cell_volumes = torch.linalg.det(batch.cell).abs().detach().cpu().tolist()
    for source, energy, volume in zip(
        _source_ids(batch), energies, cell_volumes, strict=True
    ):
        if source not in initial:
            initial[source] = float(energy)
            steps[source] = 0
            if volumes is not None:
                volumes[source] = float(volume)


def _capture_by_source(
    batch: Batch,
    configurations: list[Configuration],
    initial: dict[int, float],
    steps: dict[int, int],
    fmax_values: torch.Tensor,
    converged: bool,
    results: dict[int, RelaxationResult],
    initial_volumes: dict[int, float] | None = None,
    cell_force_values: torch.Tensor | None = None,
    pressure: float = 0.0,
) -> None:
    energies = batch.energy.detach().cpu().flatten().tolist()
    fmax_cpu = fmax_values.detach().cpu().tolist()
    positions = batch.positions.detach().cpu()
    offsets = batch.batch_ptr.detach().cpu().tolist()
    variable_cell = cell_force_values is not None
    if variable_cell:
        cells = batch.cell.detach().cpu().numpy()
        cell_force_cpu = cell_force_values.detach().cpu().tolist()
        pressures = hydrostatic_pressure(batch, pressure).detach().cpu().tolist()
    for local, source in enumerate(_source_ids(batch)):
        configuration = configurations[source]
        results[configuration.index] = RelaxationResult(
            index=configuration.index,
            initial_energy=initial[source],
            final_energy=float(energies[local]),
            steps=steps[source],
            final_fmax=float(fmax_cpu[local]),
            converged=converged,
            positions=positions[offsets[local] : offsets[local + 1]].numpy(),
            cell=cells[local] if variable_cell else None,
            initial_volume=(initial_volumes or {}).get(source),
            final_cell_force=float(cell_force_cpu[local]) if variable_cell else None,
            final_pressure=float(pressures[local]) if variable_cell else None,
        )


def _record_selected(
    batch: Batch,
    selected: torch.Tensor,
    configurations: list[Configuration],
    initial: dict[int, float],
    steps: dict[int, int],
    fmax_values: torch.Tensor,
    converged: bool,
    results: dict[int, RelaxationResult],
    initial_volumes: dict[int, float] | None = None,
    cell_force_values: torch.Tensor | None = None,
    pressure: float = 0.0,
) -> None:
    part = batch.index_select(selected)
    _capture_by_source(
        part,
        configurations,
        initial,
        steps,
        fmax_values[selected],
        converged,
        results,
        initial_volumes,
        cell_force_values[selected] if cell_force_values is not None else None,
        pressure,
    )


def relax_refill(
    configurations: list[Configuration],
    model: MACEWrapper,
    device: torch.device,
    batch_size: int = DEFAULT_BATCH_SIZE,
    fmax: float = DEFAULT_FMAX,
    max_steps: int = DEFAULT_MAX_STEPS,
    fire_dt: float = DEFAULT_FIRE_DT,
    fire_maxstep: float = DEFAULT_FIRE_MAXSTEP,
    neighbor_skin: float = DEFAULT_NEIGHBOR_SKIN,
    variable_cell: bool = False,
    pressure: float = 0.0,
) -> dict[int, RelaxationResult]:
    """Relaxation keeping the batch full via cached-E0/F0 refill."""
    cache = precompute_cache(
        configurations, model, device, batch_size, neighbor_skin, variable_cell, pressure
    )
    dataset = CachedDataset(configurations, cache)
    sampler = SizeAwareSampler(dataset, max_batch_size=batch_size, shuffle=False)
    optimizer_class = FIRE2VariableCell if variable_cell else FIRE2
    optimizer = optimizer_class(
        model, dt=fire_dt, maxstep=fire_maxstep, sampler=sampler, max_batch_size=batch_size
    )
    batch = sampler.build_initial_batch().to(device)
    optimizer._ensure_state_initialized(batch)
    hook = _neighbor_hook(
        model.model_config.neighbor_config, neighbor_skin, variable_cell
    )

    initial: dict[int, float] = {}
    steps: dict[int, int] = {}
    volumes: dict[int, float] = {}
    results: dict[int, RelaxationResult] = {}
    _initialize_records(batch, initial, steps, volumes)
    measure, atom_fmax, cell_max = _convergence_measure(batch, variable_cell)

    while batch is not None:
        active = batch.status.squeeze(-1) < optimizer.exit_status
        out_of_steps = torch.tensor(
            [steps[source] >= max_steps for source in _source_ids(batch)],
            device=device,
            dtype=torch.bool,
        ) & active
        converged = (measure <= fmax) & active
        graduating = converged | out_of_steps
        if bool(graduating.any()):
            selected = torch.nonzero(graduating, as_tuple=False).flatten()
            done = torch.nonzero(converged, as_tuple=False).flatten()
            failed = torch.nonzero(out_of_steps & ~converged, as_tuple=False).flatten()
            if done.numel():
                _record_selected(
                    batch, done, configurations, initial, steps, atom_fmax, True,
                    results, volumes, cell_max, pressure,
                )
            if failed.numel():
                _record_selected(
                    batch, failed, configurations, initial, steps, atom_fmax, False,
                    results, volumes, cell_max, pressure,
                )
            batch.status[selected] = optimizer.exit_status

        if bool((batch.status.squeeze(-1) >= optimizer.exit_status).any()):
            batch = optimizer.refill_check(batch, optimizer.exit_status)
            if batch is None:
                break
            _initialize_records(batch, initial, steps, volumes)
            # Newcomers arrive with cached F0 (and S0 under a variable cell) and
            # retained graphs keep their own, so the convergence measure follows
            # from the batch without a model call.
            measure, atom_fmax, cell_max = _convergence_measure(batch, variable_cell)
            continue

        optimizer.pre_update(batch)
        for source in _source_ids(batch):
            steps[source] += 1
        evaluate(batch, model, hook, pressure)
        measure, atom_fmax, cell_max = _convergence_measure(batch, variable_cell)

    if len(results) != len(configurations):
        raise RuntimeError(
            f"refill relaxation produced {len(results)} of {len(configurations)} results"
        )
    return results


def write_relaxed(
    configuration: Configuration,
    result: RelaxationResult,
    filename: str,
    ase_format: str,
) -> Path:
    """Write the relaxed geometry beside its input, preserving PBC.

    Under a variable cell the new lattice is applied first, without scaling the
    atoms, because the relaxed positions are already Cartesian in that cell.
    """
    atoms = configuration.atoms.copy()
    if result.cell is not None:
        atoms.set_cell(result.cell, scale_atoms=False)
    atoms.set_positions(result.positions)
    path = configuration.path / filename
    if ase_format == "vasp":
        ase_write(path, atoms, format="vasp", direct=True, vasp5=True)
    else:
        ase_write(path, atoms, format=ase_format)
    return path
