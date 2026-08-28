#!/usr/bin/env python3
"""Benchmark sod_mace on a SOD level, and write the JSON the plot script reads.

This produces the numbers quoted in the Performance section of docs/mace.rst.
Run it with an interpreter that has the MLIP stack (see pysod/README.md); a CUDA
GPU is required, since the point of the exercise is batching on the GPU.

    python benchmarks/bench_sod_mace.py                     # defaults used for the docs
    python benchmarks/bench_sod_mace.py -level path/to/nXX -repeats 3

Methodology, chosen so the numbers mean what they appear to mean:

* Model loading and the foundation-model download are outside every timer.
* Each measurement is preceded by a warm-up over the same batch shape, and the
  first timed repeat is discarded as well -- ragged-shape setup and kernel
  autotuning otherwise dominate the first run and inflate the spread by seconds.
* ``torch.cuda.synchronize`` brackets every timer.
* The reported value is the **best** of the remaining repeats. Best-of is the
  convention the original prototype benchmarks used, and with the first repeat
  discarded the spread is small enough (< 0.4 s) that median and best agree
  closely.

Nothing here is imported by sod_mace itself; this is a development tool, kept
out of the CLI so that benchmark-only knobs never become user surface.
"""

from __future__ import annotations

import argparse
import json
import pathlib
import statistics
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "pysod"))

import torch  # noqa: E402

import mace_backend as mb  # noqa: E402

DEFAULT_LEVEL = ROOT / "examples" / "example01" / "FILER0_mace" / "n04"
DEFAULT_STRUCTURE = "configuration.cif"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="bench_sod_mace.py", allow_abbrev=False)
    parser.add_argument("-level", type=pathlib.Path, default=DEFAULT_LEVEL,
                        help="Level directory holding cYY/ (default: example01/FILER0_mace/n04).")
    parser.add_argument("-structure", default=DEFAULT_STRUCTURE)
    parser.add_argument("-model", default="medium-0b2")
    parser.add_argument("-repeats", type=int, default=4,
                        help="Timed repeats per point; the first is discarded (default 4).")
    parser.add_argument("-out", type=pathlib.Path,
                        default=ROOT / "benchmarks" / "results" / "sod_mace_benchmark.json")
    parser.add_argument("-sizes", default=None,
                        help="Comma-separated batch sizes, overriding the default "
                             "1,2,4,8,16,32,<n>. Sizes above the configuration count "
                             "are dropped.")
    parser.add_argument("-backends", default="cueq,e3nn",
                        help="Comma-separated single-point backends to time "
                             "(default cueq,e3nn). e3nn spills past batch 16 on a "
                             "16 GiB card, which is slow to measure.")
    parser.add_argument("-skip-varcell", action="store_true",
                        help="Skip the variable-cell plans, which dominate the "
                             "runtime when only fixed-cell points are wanted.")
    parser.add_argument("-skip-relax", action="store_true")
    parser.add_argument("-skip-single-point", action="store_true")
    return parser.parse_args()


def timed(device, fn):
    torch.cuda.synchronize(device)
    start = time.perf_counter()
    value = fn()
    torch.cuda.synchronize(device)
    return value, time.perf_counter() - start


def repeat(device, fn, repeats: int):
    """Run fn repeats+1 times, discard the first, return (best, median, worst, last value).

    Deliberately no ``empty_cache`` between repeats: releasing the caching
    allocator forces every run to re-acquire device memory from the driver, and
    that cost lands *inside* the timer. It is negligible for a multi-second
    relaxation but dominates a 0.1 s single-point pass -- with it, measured
    single-point throughput was halved. A warm allocator is also the realistic
    case, since real runs process many batches back to back.
    """
    times, value = [], None
    for _ in range(repeats + 1):
        value, elapsed = timed(device, fn)
        times.append(elapsed)
    times = times[1:]
    return min(times), statistics.median(times), max(times), value


def load_configurations(level: pathlib.Path, structure: str):
    dirs = sorted(
        ((int(p.name[1:]), p) for p in level.glob("c*")
         if p.is_dir() and p.name[1:].isdigit()),
        key=lambda pair: pair[0],
    )
    if not dirs:
        raise SystemExit(f"no cYY/ directories in {level}")
    configurations, _ = mb.read_configurations(dirs, structure)
    if not configurations:
        raise SystemExit(f"no {structure} files in {level}/c*/")
    return configurations


def batch_sizes_up_to(n: int) -> list[int]:
    sizes = [s for s in (1, 2, 4, 8, 16, 32) if s < n]
    return sizes + [n]


def main() -> int:
    args = parse_args()
    if not torch.cuda.is_available():
        raise SystemExit("a CUDA GPU is required for this benchmark")
    device = torch.device("cuda")

    configurations = load_configurations(args.level, args.structure)
    n = len(configurations)
    natoms = len(configurations[0].atoms)
    if args.sizes:
        sizes = sorted({int(s) for s in args.sizes.split(",") if int(s) <= n})
        if not sizes:
            raise SystemExit(f"no requested batch size is <= {n}")
    else:
        sizes = batch_sizes_up_to(n)
    results: dict = {
        "gpu": torch.cuda.get_device_name(device),
        "gpu_total_mib": torch.cuda.get_device_properties(device).total_memory / 1024**2,
        "model": args.model,
        "level": str(args.level),
        "n_configs": n,
        "n_atoms": natoms,
        "repeats": args.repeats,
        "single_point": [],
        "memory": [],
        "relaxation": [],
    }
    print(f"{n} configurations of {natoms} atoms on {results['gpu']} "
          f"({results['gpu_total_mib']:.0f} MiB)", flush=True)

    if not args.skip_single_point:
        print("\nsingle point (throughput and peak memory)", flush=True)
        wanted = [b.strip() for b in args.backends.split(",") if b.strip()]
        for cueq in (True, False):
            name = "cueq" if cueq else "e3nn"
            if name not in wanted:
                continue
            model = mb.load_model(device, model=args.model, enable_cueq=cueq)
            for size in sizes:
                torch.cuda.empty_cache()
                mb.warm_up(configurations, model, device, size, mb.DEFAULT_NEIGHBOR_SKIN)
                run = lambda: mb.single_point(configurations, model, device, batch_size=size)
                best, median, worst, _ = repeat(device, run, args.repeats)
                torch.cuda.empty_cache()
                torch.cuda.reset_peak_memory_stats(device)
                timed(device, run)
                peak = torch.cuda.max_memory_allocated(device) / 1024**2
                reserved = torch.cuda.max_memory_reserved(device) / 1024**2
                results["single_point"].append(
                    {"backend": name, "batch": size, "best_s": best, "median_s": median,
                     "worst_s": worst, "per_s": n / best}
                )
                results["memory"].append({"backend": name, "batch": size,
                                          "peak_mib": peak, "reserved_mib": reserved})
                print(f"  {name} batch {size:3d}: best {best:8.4f} s  {n/best:8.1f} struct/s"
                      f"  peak {peak:9.1f} MiB", flush=True)
            del model
            torch.cuda.empty_cache()

    if not args.skip_relax:
        print("\nrelaxation (fmax 0.05, max 200 steps)", flush=True)
        model = mb.load_model(device, model=args.model, enable_cueq=True,
                              need_forces=True, need_stress=True)
        plans = [
            ("fixed cell, fixed batches", mb.relax, {}, sizes),
            ("fixed cell, refill", mb.relax_refill, {}, [s for s in sizes if s != 1]),
            # Variable cell is measured at one batch size only; it is by far the
            # most expensive plan. Honour -sizes here too, so a sweep that does
            # not ask for 16 does not silently pay for it.
            ("variable cell, fixed batches", mb.relax, {"variable_cell": True},
             [] if args.skip_varcell else [s for s in (16,) if s in sizes]),
            ("variable cell, refill", mb.relax_refill, {"variable_cell": True},
             [] if args.skip_varcell else [s for s in (16,) if s in sizes]),
        ]
        for label, fn, extra, plan_sizes in plans:
            for size in plan_sizes:
                if size > n:
                    continue
                torch.cuda.empty_cache()
                mb.warm_up(
                    configurations,
                    model,
                    device,
                    size,
                    mb.DEFAULT_NEIGHBOR_SKIN,
                    variable_cell=bool(extra.get("variable_cell", False)),
                )
                run = lambda: fn(configurations, model, device, batch_size=size,
                                 fmax=0.05, max_steps=200, **extra)
                best, median, worst, value = repeat(device, run, args.repeats)
                # Peak memory needs its own pass with the counter reset, exactly
                # as the single-point loop does, so the two are comparable.
                torch.cuda.empty_cache()
                torch.cuda.reset_peak_memory_stats(device)
                timed(device, run)
                peak = torch.cuda.max_memory_allocated(device) / 1024**2
                reserved = torch.cuda.max_memory_reserved(device) / 1024**2
                results["relaxation"].append(
                    {"mode": label, "batch": size, "best_s": best, "median_s": median,
                     "worst_s": worst, "per_s": n / best,
                     "peak_mib": peak, "reserved_mib": reserved,
                     "converged": sum(r.converged for r in value.values()),
                     "total_steps": sum(r.steps for r in value.values())}
                )
                print(f"  {label:30s} batch {size:3d}: best {best:8.3f} s  "
                      f"{n/best:7.1f} struct/s  "
                      f"{sum(r.converged for r in value.values())}/{n} converged  "
                      f"{sum(r.steps for r in value.values())} steps", flush=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
