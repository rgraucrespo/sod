# pysod

Python tools for SOD. Everything here is optional: the Fortran programs and the
`bin/sod_*.sh` wrappers build and run without any of it, and `bin/sod_run_tests.sh`
skips the Python tests when the dependencies are absent.

| File | Role |
|---|---|
| `../bin/sod_mace.sh` | **User entry point** — wrapper that locates SODPROJECT and the interpreter |
| `sod_mace.py` | CLI: MACE energies and fixed-cell relaxation over SOD configurations |
| `sodpaths.py` | SODPROJECT / level / configuration discovery (Python port of `bin/sod_common.sh`) |
| `mace_backend.py` | Model loading, ALCHEMI batching, evaluation and FIRE2 relaxation |

## Installation

Not installed by `make`. Create a dedicated conda environment — everything
below is on PyPI:

```bash
conda create -n nvalchemi python=3.12
conda activate nvalchemi

# 1. PyTorch for your CUDA runtime (cu130 shown; see pytorch.org for others,
#    or omit --index-url for a CPU-only build)
pip install torch --index-url https://download.pytorch.org/whl/cu130

# 2. NVIDIA ALCHEMI toolkit (pulls in nvalchemi-toolkit-ops), ASE and MACE
pip install nvalchemi-toolkit ase mace-torch

# 3. Optional cuEquivariance acceleration; match the ops wheel to your CUDA
#    major version (cu13 for CUDA 13, cu12 for CUDA 12)
pip install cuequivariance-torch cuequivariance-ops-torch-cu13
```

ALCHEMI requires Python 3.11-3.13. A GPU is optional: `-device cpu` works
without CUDA.

Check it:

```bash
python -c "import torch, ase, mace, nvalchemi; print(torch.cuda.is_available())"
```

On CUDA 13 the cuEquivariance ops wheel needs its libraries on
`LD_LIBRARY_PATH` before the process starts (otherwise `libnvrtc.so.13` may not
resolve). `sod_mace.py` handles this itself by re-executing once with the path
prepended — no shell setup required.

## Dependencies

The stack these tools were developed and validated against:

| Package | Version |
|---|---|
| Python | 3.12 |
| torch | 2.13.0 (CUDA 13) |
| nvalchemi-toolkit | 0.2.0 |
| nvalchemi-toolkit-ops | 0.4.1 |
| mace-torch | 0.3.15 |
| ASE | 3.29.0 |
| cuequivariance, cuequivariance-torch | 0.11.1 |

`nvalchemi-toolkit` **0.2.0 specifically**: `mace_backend.py` calls a few private
APIs of that release, which is deliberate — the public hook dispatch does not fit
the explicit evaluation loop. The authoritative list is `PRIVATE_APIS` in
`mace_backend.py`; it is not repeated here, because a copied list is the thing
that goes stale. `check_private_api()` verifies every entry before a run starts,
so an incompatible toolkit gives an error naming the missing attribute instead of
an `AttributeError` minutes into a relaxation, and `bin/sod_run_tests.sh` runs the
same check (it skips when the stack is absent).

cuEquivariance is optional and detected at runtime; without it the plain e3nn
backend is used. CUDA is optional too (`-device cpu`).

## Running

Use the wrapper `bin/sod_mace.sh`, like every other SOD tool. Tell it which
interpreter has the dependencies — no `conda activate` needed:

```bash
export SOD_PYTHON=~/miniconda3/envs/nvalchemi/bin/python
sod_mace.sh
```

Put that `export` in your `~/.bashrc` and it is a one-time setup. The wrapper
checks the stack is importable and fails with an actionable message otherwise.

`sod_mace.sh` follows the same two invocation modes as the `sod_*_ener.sh`
collectors:

```bash
cd SODPROJECT      && sod_mace.sh        # every nXX/ level in turn
cd SODPROJECT/n04  && sod_mace.sh        # that level only
sod_mace.sh n04                          # a named level
```

### Settings file

Options can be set in a `mace_settings.yaml` file instead of on the command
line. A copy in `nXX/` takes priority over one in SODPROJECT, matching the
`INSQS` convention; anything given on the command line overrides the file.

```yaml
model: medium-0b2
device: cuda
batch: 16
relax: true
batchmode: refill
fmax: 0.05
maxsteps: 200
```

Keys are exactly the option names below, without the leading dash. Unknown
keys, wrong types and invalid choices are rejected with a message naming the
file and key rather than being silently ignored. `examples/example01/FILER0_mace`
carries a fully commented example.

One YAML wrinkle worth knowing: bare `on`/`off`/`yes`/`no` are booleans in
YAML 1.1, so `cueq: off` arrives as `false`. That specific case is handled, but
quoting (`cueq: "off"`) is unambiguous.

Common options (single-dash, matching SOD's house style):

```
-model NAME        MACE foundation model            (default medium-0b2)
-checkpoint PATH   local torch-serialised MACE model
-device cuda|cpu   (default cuda)
-batch N           maximum structures per batch     (default 16)
-relax             FIRE2 relaxation of positions    (default: single point)
-relaxcell         also relax the cell              (-relax only)
-pressure P        target pressure in GPa           (default 0, -relaxcell only)
-lattice SYS       cub|tet|ort|hex|rho|mon|tri; also writes nXX/CELL
-batchmode M       fixed | refill batch scheduling  (default fixed, -relax only)
-fmax F            force convergence, eV/A          (default 0.05)
-maxsteps N        maximum relaxation steps         (default 200)
-structure NAME    structure file in each cYY/      (default: from INSOD FILER)
-cueq on|off|auto  (default auto)
-force             overwrite existing result files
-q                 quiet
```

### Output

`nXX/ENERGIES` in the same two-column `m  E` format (eV) that
`sod_vasp_ener.sh` and friends produce, so `sod_stat.sh` reads it unchanged:

```
# MACE medium-0b2, cuEq, float32, single-point
1  -372.6918334961
2  -372.7885131836
```

**Existing MACE result files are never overwritten without `-force`**. This
includes `ENERGIES`, `ENTHALPIES`, `MACE_RELAXATION.dat`, `CELL` when requested,
and `cYY/relaxed.cif` or `cYY/relaxed_POSCAR`. A rerun narrower than the one
before it also refuses to leave the wider run's summary files behind: a
single-point rerun reports the previous `MACE_RELAXATION.dat`, one without
`-lattice` reports the previous `CELL`, and `-force` removes them. Only these
level summaries are treated as stale — relaxed structures are inputs for later
work, so they are never removed on your behalf.

### Variable-cell relaxation

`-relaxcell` relaxes the cell alongside the positions, using the stress MACE
computes by strain autograd and ALCHEMI's `FIRE2VariableCell`. A graph is
converged when the maximum atomic force **and** the maximum cell force are both
below `-fmax` (both eV/Å), so the residual pressure at convergence is roughly
`fmax / a²` — about 0.1 GPa for `-fmax 0.05` on an 8.5 Å cell.

`-pressure P` (GPa) relaxes the enthalpy H = E + P V instead of the energy. It
is applied by storing the effective stress `σ + P·I` in the batch, because
`FIRE2VariableCell.pre_update` derives the cell force from `batch.stress`
itself — passing a pressure only to the convergence test would leave the
optimizer relaxing to zero pressure regardless. Validated against ASE
`FrechetCellFilter(scalar_pressure=...)`: volumes agree to 0.03–0.09 % at 0 and
5 GPa.

`ENERGIES` always stores the final internal energy `E`, matching the convention
used by the shell extractor scripts. At nonzero pressure, MACE also writes
`ENTHALPIES` with `H = E + P V`; `MACE_RELAXATION.dat` reports the final energy
and volume used to compute it. At zero pressure, no `ENTHALPIES` file is written.
If a zero-pressure run finds a stale `ENTHALPIES` file, it refuses to continue
unless `-force` is used to remove it.

Variable-cell relaxation is more expensive than fixed-cell relaxation. MACE
must compute stress, the optimization generally takes more steps, and SOD
rebuilds the neighbour list after every cell update. A position-only skin cache
is unsafe here: changing the periodic cell can bring new neighbours inside the
model cutoff even when Cartesian atomic positions barely move.

`-lattice` writes `nXX/CELL` in the same format as `sod_vasp_cell.sh`; see
`cellparams.py`, which reproduces the committed DFT `CELL` of
`example01/FILER11_vasp/n04` exactly.

Because `CELL` has no configuration-index column, MACE, VASP, and GULP write it
only when results cover every configuration `1..N` in the level's `ENSEMBLE`.
Sparse CPME calibration results remain valid as indexed `ENERGIES`, but cannot
produce an unambiguous positional `CELL` and are skipped with a warning.

### Batch scheduling during relaxation

`-batchmode` selects how structures are scheduled into batches. It changes no
physics: the relaxation is fixed-cell FIRE2 to the same `-fmax` either way.

`-batchmode fixed` (default) relaxes independent chunks of `-batch` structures;
the batch shrinks as graphs converge, so the GPU is progressively under-filled
towards the end of each chunk.

`-batchmode refill` keeps the batch full by admitting queued structures into
vacated slots, via ALCHEMI's `SizeAwareSampler` and `FIRE2.refill_check`. Every
structure's initial energy and forces are computed once up front in batches and
travel with the newcomer, so a refill costs no extra model call — a naive refill
that evaluates newcomers eagerly is *slower* than fixed chunks, which is why
only the cached variant is exposed.

Batching is what buys the speed. Fixed-cell relaxation of 710 configurations
of 64 atoms (RTX 5060 Ti, medium-0b2, cuEq, `fmax = 0.05`; 6620 optimizer steps
and 710/710 converged in every case), as a speedup over relaxing one structure
at a time:

| batch | `fixed` (s) | `fixed` speedup | `refill` (s) | `refill` speedup | refill gain |
|------:|------------:|----------------:|-------------:|-----------------:|------------:|
| 1     | 162.40      | 1.00x           | —            | —                | —           |
| 2     | 102.02      | 1.59x           | 87.57        | 1.85x            | 1.16x       |
| 4     | 63.76       | 2.55x           | 48.19        | 3.37x            | 1.32x       |
| 8     | 39.66       | 4.09x           | 29.60        | 5.49x            | **1.34x**   |
| 16    | 27.43       | 5.92x           | 22.48        | 7.23x            | 1.22x       |
| 32    | 22.31       | 7.28x           | 19.63        | 8.27x            | 1.14x       |
| 64    | 20.68       | 7.85x           | 17.74        | **9.15x**        | 1.17x       |

Both speedup columns share one reference: batch 1 at 162.40 s, which has no
batch mode of its own — with a single structure in flight there is nothing to
refill. Larger batches are faster over the whole range measured, so the default
`-batch 16` is a memory-safe choice (about 2.5 GiB here) rather than an optimum;
raise it if the card allows. Each doubling is worth roughly 1.6x up to batch 8,
then 1.45x, 1.23x and 1.08x.

The measurement deliberately uses 710 structures rather than the level's 71.
A 71-configuration relaxation takes about two seconds, short enough that the
one-off `refill` cache pass is a large fraction of it and run-to-run scatter
reached 20 %; at 710 the scatter is 0.3–1.4 %. Each of the 71 configurations is
written out ten times as an independent structure, and every point reports
exactly ten times the 71-configuration step count with 710/710 converged, which
is what shows the copies relax independently.

The two modes are numerically equivalent — identical step counts for every
configuration, and final energies within 2.4e-4 eV of each other, the float32
reproducibility limit. Prefer `refill` for large levels; `fixed` is simpler and
has no cache pass, so it stays the default. The raw JSON behind the table is in
`benchmarks/results/sod_mace_benchmark_710.json`, and `docs/mace.rst` has the
single-point and memory numbers with a figure.

With `-relax`, each configuration also gets `cYY/relaxed.cif` (or
`relaxed_POSCAR` for VASP input) beside its input, plus a summary table
`nXX/MACE_RELAXATION.dat`; `ENERGIES` then carries the relaxed energies. The
input structure files are never modified, so a rerun cannot silently relax an
already-relaxed geometry.

`-writerelaxed no` suppresses the per-configuration structures while keeping
every level summary. It is for levels too large to hold one file per
configuration: a 64-atom `relaxed.cif` is ~6.3 kB but occupies 8 kiB on a
4 kiB-block filesystem, so 10^6 configurations cost ~16 GiB and 3M inodes once
the `cYY/` directories and input structures are counted. The relaxation is
identical either way — same trajectories, step counts and energies — and the
time saved is small: ~4 min per 10^6 structures written, against ~7 h of
fixed-cell (or ~50 h of variable-cell) relaxation for the same count on an
RTX 5060 Ti. The saving is disk, inodes, and every later traversal, archive or
transfer of the tree. To keep a few relaxed geometries out of a very large run,
use `-writerelaxed no`, choose configurations from `ENERGIES`, and re-relax
those in a smaller level. Note that `-writerelaxed no` neither removes nor
checks for `relaxed.cif` files from an earlier run — nothing is written, so
there is nothing to protect — and a leftover one will sit beside an `ENERGIES`
it no longer matches.

Structure files are located from the INSOD `FILER` value: `0` reads
`configuration.cif`, `11` reads `POSCAR`. Other calculators write formats ASE
cannot read directly — point at something readable with `-structure`.

Note this tool does **not** read `ENSEMBLE` and knows nothing about
degeneracies; it is an energy/relaxation engine, exactly like the shell
collectors it parallels. Feed its `ENERGIES` to `sod_stat.sh` for
statistical-mechanical analysis.

## Tests

`bin/sod_run_tests.sh` includes three `sod_mace` cases that skip unless the
stack is importable. Point it at the right interpreter:

```bash
SOD_PYTHON=~/miniconda3/envs/nvalchemi/bin/python bin/sod_run_tests.sh
```

## Provenance

Ported from prototype work that measured a 6.41x end-to-end speedup over serial
ASE/MACE on 71 configurations of 64 atoms, dominated by batching rather than
cuEquivariance. Batched energies
reproduce the prototype to ~2e-4 eV, within float32 tolerance.

Measured again for this code (RTX 5060 Ti, medium-0b2, float32): single-point
throughput rises from 70 to 1119 structures/s between batch 1 and batch 64
(16x), and fixed-cell relaxation falls by 7.9x over the same range, 9.2x with
`-batchmode refill`. At batch 1 cuEquivariance is marginally *slower* than
e3nn; its real contribution is memory headroom — 9.5 GiB against e3nn's
24.3 GiB at batch 64 — which is what makes the large batches possible at all.
See the Performance section of `docs/mace.rst`, and `benchmarks/` for the
script that produces it.
