#!/usr/bin/env python3
"""MACE machine-learning-potential energies and relaxation for SOD configurations.

This is the MLIP counterpart of the bin/sod_*_ener.sh collectors: it walks the
``cYY/`` configuration directories of a level and writes ``nXX/ENERGIES`` in the
same two-column ``m  E`` format, so sod_stat.sh consumes it unchanged.  It does
not read ENSEMBLE and knows nothing about degeneracies.

Like the shell collectors it supports both invocation modes:

    cd SODPROJECT      && sod_mace.py        # every nXX/ level in turn
    cd SODPROJECT/n04  && sod_mace.py        # that level only
    cd anywhere        && sod_mace.py n04    # a named level

Run it with an interpreter that has the MLIP stack installed; see pysod/README.md.
"""

from __future__ import annotations

import argparse
import os
import sys
import sysconfig
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cellparams  # noqa: E402
import sodpaths  # noqa: E402
from sodpaths import SodLayoutError  # noqa: E402

ENERGIES_FILE = "ENERGIES"
ENTHALPIES_FILE = "ENTHALPIES"
RELAXATION_TABLE = "MACE_RELAXATION.dat"
CELL_FILE = "CELL"
SETTINGS_FILE = "mace_settings.yaml"


def ensure_cu13_loader_path() -> None:
    """Expose the CUDA 13 cuEquivariance wheel libraries before torch is imported.

    The cu13 ops wheel can fail to load libraries such as libnvrtc.so.13 unless
    they are on LD_LIBRARY_PATH when the process starts, so the interpreter is
    re-executed once with the path prepended.  Called explicitly from main()
    rather than at import time, and only when it can matter.
    """
    cu13_lib = Path(sysconfig.get_paths()["purelib"]) / "nvidia" / "cu13" / "lib"
    if not cu13_lib.is_dir() or os.environ.get("NV_CUEQ_REEXEC") == "1":
        return
    current = os.environ.get("LD_LIBRARY_PATH", "")
    if str(cu13_lib) in current.split(":"):
        return
    environment = os.environ.copy()
    environment["LD_LIBRARY_PATH"] = f"{cu13_lib}:{current}" if current else str(cu13_lib)
    environment["NV_CUEQ_REEXEC"] = "1"
    os.execvpe(sys.executable, [sys.executable, *sys.argv], environment)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = build_parser()
    return parser.parse_args(argv)


# Every option defaults to SUPPRESS so the namespace records only what the user
# actually typed; the rest is filled from mace_settings.yaml, then DEFAULTS.
DEFAULTS: dict[str, object] = {
    "model": "medium-0b2",
    "checkpoint": None,
    "device": "cuda",
    "batch": 16,
    "relax": False,
    "relaxcell": False,
    "writerelaxed": True,
    "pressure": 0.0,
    "lattice": None,
    "batchmode": "fixed",
    "fmax": 0.05,
    "maxsteps": 200,
    "structure": None,
    "cueq": "auto",
    "force": False,
    "q": False,
}

CHOICES: dict[str, tuple[str, ...]] = {
    "device": ("cuda", "cpu"),
    "cueq": ("on", "off", "auto"),
    "batchmode": ("fixed", "refill"),
    "lattice": cellparams.LATTICE_SYSTEMS,
}

TYPES: dict[str, type] = {
    "model": str,
    "checkpoint": str,
    "device": str,
    "batch": int,
    "relax": bool,
    "relaxcell": bool,
    "writerelaxed": bool,
    "pressure": float,
    "lattice": str,
    "batchmode": str,
    "fmax": float,
    "maxsteps": int,
    "structure": str,
    "cueq": str,
    "force": bool,
    "q": bool,
}


BOOLEAN_WORDS = {
    "yes": True, "no": False,
    "true": True, "false": False,
    "on": True, "off": False,
}


def boolean_word(text: str) -> bool:
    """Parse a yes/no option value.

    Booleans elsewhere are ``store_true`` flags, but an option that defaults to
    true needs a value to switch it off, so it takes the same words YAML 1.1
    accepts for a boolean.
    """
    try:
        return BOOLEAN_WORDS[text.strip().lower()]
    except KeyError:
        raise argparse.ArgumentTypeError(
            f"expected one of {', '.join(BOOLEAN_WORDS)}, got {text!r}"
        ) from None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sod_mace.sh",
        allow_abbrev=False,
        description="MACE energies and fixed-cell relaxation over SOD configurations. "
        f"Options may also be set in a {SETTINGS_FILE} file; anything given on the "
        "command line wins.",
    )
    add = parser.add_argument
    add("level", nargs="?", default=argparse.SUPPRESS,
        help="Level directory (e.g. n04). Default: inferred from the working directory.")
    add("-model", default=argparse.SUPPRESS,
        help="MACE foundation model name (default: medium-0b2).")
    add("-checkpoint", default=argparse.SUPPRESS,
        help="Local torch-serialised MACE model (overrides -model).")
    add("-device", choices=list(CHOICES["device"]), default=argparse.SUPPRESS,
        help="(default: cuda)")
    add("-batch", type=int, default=argparse.SUPPRESS,
        help="Maximum structures per batch (default: 16).")
    add("-relax", action="store_true", default=argparse.SUPPRESS,
        help="FIRE2 relaxation of atomic positions (default: single point).")
    add("-relaxcell", action="store_true", default=argparse.SUPPRESS,
        help="Also relax the cell, driven by the MACE stress (needs -relax).")
    add("-writerelaxed", type=boolean_word, metavar="yes|no", default=argparse.SUPPRESS,
        help="Write the relaxed geometry into each cYY/ (default: yes). "
        "'no' keeps ENERGIES and MACE_RELAXATION.dat but writes no per-configuration "
        "structure, for runs with too many configurations to keep one file each.")
    add("-pressure", type=float, default=argparse.SUPPRESS,
        help="Target external pressure in GPa for -relaxcell (default: 0).")
    add("-lattice", choices=list(CHOICES["lattice"]), default=argparse.SUPPRESS,
        help="Lattice system; also writes nXX/CELL in the sod_*_cell.sh columns.")
    add("-batchmode", choices=list(CHOICES["batchmode"]), default=argparse.SUPPRESS,
        help="Batch scheduling during relaxation. fixed: independent chunks. "
        "refill: keep the batch full as structures converge (implies an initial "
        "energy/force cache; faster for many configurations). (default: fixed)")
    add("-fmax", type=float, default=argparse.SUPPRESS,
        help="Force convergence, eV/A (default: 0.05).")
    add("-maxsteps", type=int, default=argparse.SUPPRESS,
        help="Maximum relaxation steps (default: 200).")
    add("-structure", default=argparse.SUPPRESS,
        help="Structure filename in each cYY/ (default: from the INSOD FILER value).")
    add("-cueq", choices=list(CHOICES["cueq"]), default=argparse.SUPPRESS,
        help="(default: auto)")
    add("-force", action="store_true", default=argparse.SUPPRESS,
        help="Overwrite existing MACE result files.")
    add("-q", action="store_true", default=argparse.SUPPRESS,
        help="Suppress progress output.")
    return parser


def find_settings_file(root: Path, level: Path | None) -> Path | None:
    """Locate mace_settings.yaml: nXX/ wins over SODPROJECT, as INSQS does."""
    for directory in (level, root):
        if directory is None:
            continue
        candidate = directory / SETTINGS_FILE
        if candidate.is_file():
            return candidate
    return None


def load_settings(path: Path, parser: argparse.ArgumentParser) -> dict[str, object]:
    """Read and validate mace_settings.yaml into a settings dict."""
    try:
        import yaml
    except ImportError:
        fail(
            f"{path} found but PyYAML is not installed in {sys.executable}.\n"
            "       Install it with: pip install pyyaml"
        )
    try:
        data = yaml.safe_load(path.read_text())
    except yaml.YAMLError as exc:
        parser.error(f"could not parse {path}:\n{exc}")
    if data is None:
        return {}
    if not isinstance(data, dict):
        parser.error(f"{path} must contain a mapping of option names to values")

    unknown = sorted(set(data) - set(DEFAULTS))
    if unknown:
        parser.error(
            f"unknown option(s) in {path}: {', '.join(unknown)}. "
            f"Valid keys: {', '.join(sorted(DEFAULTS))}"
        )

    settings: dict[str, object] = {}
    for key, value in data.items():
        if value is None:
            continue
        expected = TYPES[key]
        # YAML 1.1 reads bare on/off/yes/no as booleans, so `cueq: off` arrives
        # as False rather than "off". Map it back for keys that accept on/off.
        if expected is str and isinstance(value, bool) and key in CHOICES:
            spelled = "on" if value else "off"
            if spelled in CHOICES[key]:
                value = spelled
        if expected is bool and not isinstance(value, bool):
            parser.error(f"{path}: '{key}' must be true or false, got {value!r}")
        if expected is int and (isinstance(value, bool) or not isinstance(value, int)):
            parser.error(f"{path}: '{key}' must be a whole number, got {value!r}")
        if expected is float:
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                parser.error(f"{path}: '{key}' must be a number, got {value!r}")
            value = float(value)
        if expected is str and not isinstance(value, str):
            parser.error(f"{path}: '{key}' must be text, got {value!r}")
        if key in CHOICES and value not in CHOICES[key]:
            parser.error(
                f"{path}: '{key}' must be one of {', '.join(CHOICES[key])}, got {value!r}"
            )
        settings[key] = value
    return settings


def validate(settings: dict[str, object], parser: argparse.ArgumentParser) -> None:
    if settings["batch"] < 1:
        parser.error("batch must be at least 1")
    if settings["fmax"] <= 0:
        parser.error("fmax must be positive")
    if settings["maxsteps"] < 1:
        parser.error("maxsteps must be at least 1")
    if settings["batchmode"] != "fixed" and not settings["relax"]:
        parser.error("batchmode refill only applies with relax")
    if settings["relaxcell"] and not settings["relax"]:
        parser.error("relaxcell only applies with relax")
    if not settings["writerelaxed"] and not settings["relax"]:
        parser.error("writerelaxed only applies with relax")
    if settings["pressure"] and not settings["relaxcell"]:
        parser.error("pressure only applies with relaxcell")
    if settings["lattice"] and not settings["relaxcell"]:
        # With a fixed cell every row of CELL would be identical.
        parser.error("lattice only applies with relaxcell")


def fail(message: str) -> None:
    print(f"Error: {message}", file=sys.stderr)
    raise SystemExit(1)


def resolve_levels(level_arg: str | None) -> tuple[Path, list[Path]]:
    """Return the SODPROJECT root and the level directories to process."""
    root = sodpaths.find_project_root(Path.cwd())
    if level_arg:
        candidate = Path(level_arg)
        level = candidate if candidate.is_absolute() else root / candidate
        if not level.is_dir():
            raise SodLayoutError(f"level directory not found: {level}")
        return root, [level]
    enclosing = sodpaths.find_level_name(root, Path.cwd())
    if enclosing:
        return root, [root / enclosing]
    levels = sodpaths.level_dirs(root)
    if not levels:
        raise SodLayoutError(f"no nXX/ folders found in {root}")
    return root, levels


def relaxed_filename(structure_name: str) -> str:
    """Relaxed geometry sits beside its input under a distinct name.

    POSCAR-in/CONTCAR-out in spirit: never write back over the input, so a
    second run cannot silently relax an already-relaxed structure.
    """
    suffix = Path(structure_name).suffix
    return f"relaxed{suffix}" if suffix else f"relaxed_{structure_name}"


def write_energies(path: Path, energies: dict[int, float], header: str) -> None:
    """Two-column ``m  E`` in eV, byte-compatible with the sod_*_ener.sh output."""
    lines = [f"# {header}\n"]
    lines += [f"{index}  {energies[index]:.10f}\n" for index in sorted(energies)]
    path.write_text("".join(lines))


def write_relaxation_table(path: Path, results: dict, variable_cell: bool = False) -> None:
    header = (
        "# configuration initial_energy_eV final_energy_eV relaxation_energy_eV "
        "steps final_fmax_eV_A converged"
    )
    if variable_cell:
        header += " initial_volume_A3 final_volume_A3 final_cellforce_eV_A final_pressure_GPa"
    lines = [header + "\n"]
    for index in sorted(results):
        result = results[index]
        row = (
            f"{index:6d} {result.initial_energy:18.10f} {result.final_energy:18.10f} "
            f"{result.final_energy - result.initial_energy:18.10f} "
            f"{result.steps:5d} {result.final_fmax:12.6f} {int(result.converged)}"
        )
        if variable_cell:
            volume = cellparams.cell_parameters(result.cell)[-1]
            row += (
                f" {result.initial_volume:14.4f} {volume:14.4f} "
                f"{result.final_cell_force:12.6f} {result.final_pressure:12.6f}"
            )
        lines.append(row + "\n")
    path.write_text("".join(lines))


def write_cell_file(level: Path, results: dict, system: str, say) -> bool:
    """Write nXX/CELL, matching what sod_vasp_cell.sh / sod_gulp_cell.sh produce.

    Rows are positional -- one per configuration in cYY order, with no index
    column -- so a gap would silently shift every later row onto the wrong
    configuration. Skip the file rather than emit a misaligned one.

    Returns False only when the ENSEMBLE itself could not be read, which is a
    hard error: coverage cannot be judged, so any existing CELL is left alone.
    Sparse coverage is an expected outcome and returns True after removing the
    now-stale CELL.
    """
    try:
        expected_count = sodpaths.ensemble_config_count(level / "ENSEMBLE")
    except SodLayoutError as exc:
        print(f"Error: not writing {level.name}/{CELL_FILE}: {exc}.", file=sys.stderr)
        return False
    expected = set(range(1, expected_count + 1))
    available = set(results)
    if available != expected:
        (level / CELL_FILE).unlink(missing_ok=True)
        print(
            f"Warning: not writing {level.name}/{CELL_FILE}: results cover "
            f"{len(available & expected)} of {expected_count} ENSEMBLE configurations, "
            "and CELL rows are positional.",
            file=sys.stderr,
        )
        return True
    cells = [results[index].cell for index in range(1, expected_count + 1)]
    cellparams.write_cell_file(level / CELL_FILE, cells, system)
    columns = cellparams.CELL_COLUMNS[cellparams.normalise_system(system)]
    say(f"  wrote {level.name}/{CELL_FILE} ({system}: {columns})")
    return True


# Level-summary results. Each pairs with the ENERGIES of the run that wrote it,
# so one left behind by a previous, wider run would be read as belonging to this
# one. Per-configuration relaxed structures are deliberately not in this set:
# they are inputs for later work, not a summary of this run.
SUMMARY_FILES = (ENERGIES_FILE, ENTHALPIES_FILE, RELAXATION_TABLE, CELL_FILE)


def planned_result_paths(level: Path, structure_name: str, args: argparse.Namespace) -> list[Path]:
    """Files a run may create, overwrite, or remove for one level."""
    paths = [level / ENERGIES_FILE]
    if args.relax:
        paths.append(level / RELAXATION_TABLE)
        if args.writerelaxed:
            relaxed_name = relaxed_filename(structure_name)
            paths.extend(directory / relaxed_name for _, directory in sodpaths.config_dirs(level))
        if args.lattice:
            paths.append(level / CELL_FILE)
    if args.relaxcell and args.pressure != 0.0:
        paths.append(level / ENTHALPIES_FILE)
    return paths


def stale_result_paths(level: Path, structure_name: str, args: argparse.Namespace) -> list[Path]:
    """Summary files this run will not rewrite, so a leftover would be stale.

    A rerun narrower than the one before it leaves the wider run's summaries
    behind: a single-point rerun does not touch MACE_RELAXATION.dat, a run
    without -lattice does not touch CELL, and a zero-pressure run does not
    touch ENTHALPIES.
    """
    planned = set(planned_result_paths(level, structure_name, args))
    return [level / name for name in SUMMARY_FILES if level / name not in planned]


def relative_path(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def summarise_paths(paths: list[Path], root: Path, limit: int = 8) -> str:
    """Name the paths for an error message, capped so a large level stays readable.

    A level can hold thousands of cYY/ directories, so listing every relaxed
    structure would bury the message it belongs to.
    """
    names = [relative_path(path, root) for path in paths[:limit]]
    remaining = len(paths) - len(names)
    if remaining > 0:
        names.append(f"and {remaining} more")
    return ", ".join(names)


def process_level(
    level: Path,
    structure_name: str,
    ase_format: str,
    model,
    device,
    args: argparse.Namespace,
    backend,
    say,
) -> bool:
    """Run one level. Returns False if a requested output could not be written."""
    ok = True
    config_dirs = sodpaths.config_dirs(level)
    if not config_dirs:
        fail(f"no cYY/ folders found in {level.name}/.")

    configurations, missing = backend.read_configurations(config_dirs, structure_name)
    if not configurations:
        fail(f"no {structure_name} files found in {level.name}/c*/.")
    say(f"  {len(configurations)} configurations read from {level.name}/c*/{structure_name}")

    backend.warm_up(
        configurations,
        model,
        device,
        args.batch,
        backend.DEFAULT_NEIGHBOR_SKIN,
        variable_cell=args.relaxcell,
    )

    backend_name = "cuEq" if getattr(model, "_sod_cueq", False) else "e3nn"
    source = args.checkpoint or args.model

    if args.relax:
        relaxer = backend.relax_refill if args.batchmode == "refill" else backend.relax
        results = relaxer(
            configurations,
            model,
            device,
            batch_size=args.batch,
            fmax=args.fmax,
            max_steps=args.maxsteps,
            variable_cell=args.relaxcell,
            pressure=args.pressure / backend.GPA_PER_EV_A3,
        )
        energies = {index: result.final_energy for index, result in results.items()}
        enthalpies = None
        if args.relaxcell and args.pressure != 0.0:
            pressure_ev_a3 = args.pressure / backend.GPA_PER_EV_A3
            enthalpies = {
                index: result.final_energy
                + pressure_ev_a3 * cellparams.cell_parameters(result.cell)[-1]
                for index, result in results.items()
            }
        if args.writerelaxed:
            relaxed_name = relaxed_filename(structure_name)
            for configuration in configurations:
                backend.write_relaxed(
                    configuration, results[configuration.index], relaxed_name, ase_format
                )
        write_relaxation_table(level / RELAXATION_TABLE, results, args.relaxcell)
        converged = sum(result.converged for result in results.values())
        say(f"  {converged}/{len(results)} converged (fmax {args.fmax} eV/A)")
        if args.writerelaxed:
            say(f"  wrote {level.name}/c*/{relaxed_name} and {level.name}/{RELAXATION_TABLE}")
        else:
            say(f"  wrote {level.name}/{RELAXATION_TABLE} (no relaxed structures: -writerelaxed no)")
        if args.relaxcell:
            residual = max(abs(result.final_pressure) for result in results.values())
            say(f"  residual pressure at most {residual:.4f} GPa")
        if args.lattice:
            ok = write_cell_file(level, results, args.lattice, say) and ok
        kind = "variable-cell" if args.relaxcell else "fixed-cell"
        target = f", {args.pressure} GPa" if args.relaxcell else ""
        header = (
            f"MACE {source}, {backend_name}, float32, "
            f"{kind} relaxation ({args.batchmode} batching, "
            f"fmax {args.fmax} eV/A, max {args.maxsteps} steps{target}), "
            "internal energy E"
        )
    else:
        energies = backend.single_point(configurations, model, device, batch_size=args.batch)
        enthalpies = None
        header = f"MACE {source}, {backend_name}, float32, single-point"

    write_energies(level / ENERGIES_FILE, energies, header)
    say(f"  wrote {level.name}/{ENERGIES_FILE} ({len(energies)} energies)")
    if enthalpies is not None:
        enthalpy_header = (
            f"MACE {source}, {backend_name}, float32, "
            f"variable-cell relaxation ({args.batchmode} batching, "
            f"fmax {args.fmax} eV/A, max {args.maxsteps} steps, {args.pressure} GPa), "
            "enthalpy H=E+PV"
        )
        write_energies(level / ENTHALPIES_FILE, enthalpies, enthalpy_header)
        say(f"  wrote {level.name}/{ENTHALPIES_FILE} ({len(enthalpies)} enthalpies)")
    if missing:
        print(
            f"Warning: missing structures for {len(missing)} configuration(s) in {level.name}/.",
            file=sys.stderr,
        )
    return ok


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    given = vars(parser.parse_args(argv))
    level_arg = given.pop("level", None)

    try:
        root, levels = resolve_levels(level_arg)
    except SodLayoutError as exc:
        fail(str(exc))

    # A settings file is looked up in the single level being processed, else in
    # SODPROJECT. Command-line options always win over the file.
    settings_path = find_settings_file(root, levels[0] if len(levels) == 1 else None)
    from_file = load_settings(settings_path, parser) if settings_path else {}
    settings = {**DEFAULTS, **from_file, **given}
    validate(settings, parser)
    args = argparse.Namespace(**settings)
    say = (lambda *_: None) if args.q else (lambda message: print(message))
    if settings_path:
        overridden = sorted(set(from_file) & set(given))
        note = f" (overridden on the command line: {', '.join(overridden)})" if overridden else ""
        say(f"Settings: {settings_path}{note}")

    try:
        structure_name = args.structure or sodpaths.structure_filename(
            sodpaths.read_filer(root / "INSOD")
        )
    except SodLayoutError as exc:
        fail(str(exc))
    ase_format = "vasp" if structure_name == "POSCAR" else "cif"

    # Refuse to destroy existing results -- they may be expensive DFT/MLIP data.
    # Checked for every target level up front, so a multi-level run cannot abort
    # half done. At zero pressure, a leftover ENTHALPIES file would be stale
    # because this run will not rewrite it.
    if args.force:
        for level in levels:
            for stale_path in stale_result_paths(level, structure_name, args):
                if stale_path.is_file():
                    stale_path.unlink()
                    say(f"  removed stale {relative_path(stale_path, root)}")
    else:
        candidates = [
            path
            for level in levels
            for path in (
                planned_result_paths(level, structure_name, args)
                + stale_result_paths(level, structure_name, args)
            )
        ]
        clashes = [path for path in candidates if path.is_file()]
        if clashes:
            plural = "s" if len(clashes) > 1 else ""
            fail(
                f"{len(clashes)} existing result file{plural}: "
                f"{summarise_paths(clashes, root)}.\n"
                "       Refusing to overwrite or leave stale result files. "
                "Use -force to replace."
            )

    want_cueq = args.cueq != "off" and args.device == "cuda"
    if want_cueq:
        ensure_cu13_loader_path()

    try:
        import mace_backend
    except ImportError as exc:
        fail(
            f"the MLIP stack is not available in this interpreter ({sys.executable}).\n"
            f"       {exc}\n"
            "       Install torch, ase, mace-torch and nvalchemi-toolkit, or run\n"
            "       sod_mace.py with an interpreter that has them (see pysod/README.md)."
        )

    import torch

    if args.device == "cuda" and not torch.cuda.is_available():
        fail("-device cuda requested but no CUDA device is visible. Use -device cpu.")
    device = torch.device(args.device)

    enable_cueq = want_cueq and mace_backend.cueq_available()
    if args.cueq == "on" and not enable_cueq:
        fail(
            "-cueq on requested but the cuEquivariance CUDA operations are unavailable.\n"
            "       Install cuequivariance, cuequivariance-torch and the matching ops wheel."
        )

    try:
        model = mace_backend.load_model(
            device,
            model=args.model,
            checkpoint=args.checkpoint,
            enable_cueq=enable_cueq,
            need_forces=args.relax,
            need_stress=args.relaxcell,
        )
    except (FileNotFoundError, RuntimeError) as exc:
        fail(str(exc))
    model._sod_cueq = enable_cueq

    if args.device == "cuda":
        say(f"GPU: {torch.cuda.get_device_name(device)}")
    say(
        f"Model: {args.checkpoint or args.model} "
        f"({'cuEq' if enable_cueq else 'e3nn'}, float32, cutoff "
        f"{float(model.model_config.neighbor_config.cutoff):.3f} A)"
    )

    status = 0
    for level in levels:
        say(f"{level.name}/")
        if not process_level(
            level, structure_name, ase_format, model, device, args, mace_backend, say
        ):
            status = 1
    return status


if __name__ == "__main__":
    raise SystemExit(main())
