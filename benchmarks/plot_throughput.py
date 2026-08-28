#!/usr/bin/env python3
"""Regenerate docs/_static/mace_throughput.svg from a benchmark run.

    python benchmarks/bench_sod_mace.py
    python benchmarks/plot_throughput.py

Design notes, so the figure stays consistent if it is ever edited:

* Two categorical hues, blue and orange, checked for colour-vision deficiency
  separation (worst adjacent pair dE 24.7 protan, 33.6 normal vision).
* Identity is carried by direct labels at the line ends rather than a legend
  box -- stronger than colour alone, and one less thing to read.
* Values are labelled only at the two endpoints that carry the headline, never
  on every point.
* This is a static asset for Sphinx, so there is no hover layer.
* e3nn is plotted only where it fits in GPU memory: past that its timings
  measure thrashing, not throughput, and plotting them would imply a fair
  comparison that is not there.
"""

from __future__ import annotations

import argparse
import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FixedLocator, FuncFormatter  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#8b8a85"
GRID = "#e3e2dd"
CUEQ = "#2a78d6"
E3NN = "#eb6834"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="plot_throughput.py", allow_abbrev=False)
    parser.add_argument("-data", type=pathlib.Path,
                        default=ROOT / "benchmarks" / "results" / "sod_mace_benchmark_710.json")
    parser.add_argument("-out", type=pathlib.Path,
                        default=ROOT / "docs" / "_static" / "mace_throughput.svg")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    data = json.loads(args.data.read_text())
    rows, memory = data["single_point"], data["memory"]
    total_mib = data["gpu_total_mib"]

    def series(backend: str) -> list[tuple[int, float]]:
        return sorted(
            [(r["batch"], r["per_s"]) for r in rows if r["backend"] == backend],
            key=lambda pair: pair[0],
        )

    peak = {(m["backend"], m["batch"]): m["peak_mib"] for m in memory}
    cueq = series("cueq")
    # Only plot e3nn where it fits on the card.
    e3nn = [(b, v) for b, v in series("e3nn") if peak.get(("e3nn", b), 0) < total_mib]
    spilled = [b for b, _ in series("e3nn") if peak.get(("e3nn", b), 0) >= total_mib]

    fig, ax = plt.subplots(figsize=(8.0, 4.6))
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    fig.subplots_adjust(left=0.095, right=0.795, top=0.80, bottom=0.135)

    ax.grid(True, which="major", axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
        ax.spines[side].set_linewidth(0.8)

    for points, color, name in ((cueq, CUEQ, "cuEquivariance"), (e3nn, E3NN, "e3nn")):
        xs, ys = zip(*points)
        ax.plot(xs, ys, color=color, linewidth=2.0, marker="o", markersize=5.5,
                markerfacecolor=color, markeredgecolor=SURFACE, markeredgewidth=1.5,
                zorder=3, clip_on=False)
        ax.annotate(f"  {name}", (xs[-1], ys[-1]), color=color, fontsize=10.5,
                    fontweight="bold", va="center", ha="left", annotation_clip=False)

    top = max(v for _, v in cueq)
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_locator(FixedLocator([b for b, _ in cueq]))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v)}"))
    ax.set_xlim(0.93, max(b for b, _ in cueq) * 1.16)
    ax.set_ylim(0, top * 1.14)
    ax.tick_params(colors=INK_2, labelsize=9.5, length=0)
    ax.set_xlabel("batch size (structures per GPU batch)", color=INK_2, fontsize=10, labelpad=7)
    ax.set_ylabel("throughput (structures / s)", color=INK_2, fontsize=10, labelpad=7)

    ax.annotate(f"{cueq[0][1]:.0f}", cueq[0], textcoords="offset points", xytext=(9, 9),
                ha="left", color=INK, fontsize=10, fontweight="bold")
    ax.annotate(f"{cueq[-1][1]:,.0f}", cueq[-1], textcoords="offset points", xytext=(-6, 11),
                ha="right", color=INK, fontsize=10, fontweight="bold")

    if spilled:
        detail = " and ".join(
            f"{peak[('e3nn', b)]/1024:.1f} GiB ({b})" for b in spilled
        )
        ax.text(
            1.05, top * 1.065,
            f"Beyond batch {e3nn[-1][0]}, e3nn needs {detail}\n"
            f"on a {total_mib/1024:.0f} GiB card, so it is not plotted there.\n"
            f"cuEquivariance peaks at {max(peak[('cueq', b)] for b, _ in cueq)/1024:.1f} GiB"
            " - the headroom is what\nmakes the large batches, and the speed, possible.",
            color=MUTED, fontsize=8.8, ha="left", va="top", linespacing=1.5,
        )

    fig.text(0.095, 0.945, "Batching, not cuEquivariance alone, drives throughput",
             color=INK, fontsize=12.5, fontweight="bold", ha="left", va="top")
    fig.text(0.095, 0.875,
             f"{data['n_configs']} configurations of {data['n_atoms']} atoms (MgO:Ni) - "
             f"MACE {data['model']}, float32, {data['gpu'].replace('NVIDIA GeForce ', '')} - "
             f"best of {data['repeats']} runs",
             color=INK_2, fontsize=9, ha="left", va="top")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, format=args.out.suffix.lstrip("."), facecolor=SURFACE, dpi=150)
    print(f"wrote {args.out}")
    print(f"  cuEq {cueq[0][0]} -> {cueq[-1][0]}: {cueq[0][1]:.1f} -> {cueq[-1][1]:.1f} "
          f"struct/s ({cueq[-1][1]/cueq[0][1]:.1f}x)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
