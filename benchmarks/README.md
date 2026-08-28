# benchmarks

Development tooling for measuring `sod_mace` performance. Nothing here is needed
to build or run SOD, and nothing in `pysod/` imports it — benchmark-only options
are deliberately kept out of the CLI so they never become user surface.

| File | Role |
|---|---|
| `bench_sod_mace.py` | Times single-point and relaxation throughput; writes JSON |
| `plot_throughput.py` | Renders `docs/_static/mace_throughput.svg` from that JSON |
| `plot_batching.py` | Renders `docs/_static/mace_batching.svg` from that JSON |
| `results/` | Committed JSON from the runs behind the numbers in `docs/mace.rst` |

![Batch size buys speed with diminishing returns, at a memory cost that stays
linear. Three panels against batch size for 710 configurations of 64 atoms:
(a) relaxation speedup over batch 1, reaching 7.9x with fixed batches and 9.2x
with refill; (b) wall-clock time to process the level; (c) peak GPU memory
against a 16 GiB card limit.](../docs/_static/mace_batching.svg)

## Running

Needs a CUDA GPU and the MLIP stack (see `pysod/README.md`):

```bash
~/miniconda3/envs/nvalchemi/bin/python benchmarks/bench_sod_mace.py
~/miniconda3/envs/nvalchemi/bin/python benchmarks/plot_batching.py
~/miniconda3/envs/nvalchemi/bin/python benchmarks/plot_throughput.py
```

Both plot scripts default to the 710-configuration run,
`results/sod_mace_benchmark_710.json`, which is the dataset behind the
Performance section of `docs/mace.rst`, so the commands above redraw the
committed figures as they stand. The earlier 71-configuration run is kept as
`results/sod_mace_benchmark.json` and can be plotted with `-data`.

`bench_sod_mace.py` measures `example01/FILER0_mace/n04`, 64 atoms per
configuration, MACE `medium-0b2`, float32. A full run takes roughly ten
minutes. `-level`, `-model` and `-repeats` change the workload;
`-skip-relax` / `-skip-single-point` run half of it.

## Methodology

Getting these numbers to mean what they appear to mean needs four things, all
handled by the script:

- Model loading and any foundation-model download are outside every timer.
- Every point is warmed up over the same batch shape **and** the first timed
  repeat is discarded. Ragged-shape setup and kernel autotuning otherwise
  dominate the first run — before this, the spread across repeats was as much
  as 8 s on a 2 s measurement; after it, under 0.4 s.
- `torch.cuda.synchronize` brackets every timer, so the numbers are not just
  kernel-launch latency.
- The reported value is the best of the remaining repeats, the convention the
  original prototype benchmarks used.

Peak memory is measured in a separate pass after `reset_peak_memory_stats`, and
it is the reason the e3nn curve stops where it does: past batch 16 e3nn exceeds
a 16 GiB card, so its timings there measure thrashing rather than throughput and
the plot omits them.

## Regenerating the figures after a change

Both plot scripts read everything they draw — batch sizes, memory ceilings,
GPU name, repeat count — from the JSON, so re-running them is enough; the
captions and the memory annotation follow the data. If you re-run on a
different GPU, update the prose in `docs/mace.rst` too, which quotes specific
numbers.
