#!/usr/bin/env python3
"""Regenerate docs/_static/mace_batching.svg from a benchmark run.

    python benchmarks/bench_sod_mace.py -level <nXX> -backends cueq
    python benchmarks/plot_batching.py

Three panels over one shared batch-size axis:

  (a) relaxation speedup against batch 1, for both batch modes;
  (b) wall-clock time for the three workloads that batch size actually changes;
  (c) peak GPU memory for the same three.

Design notes, so the figure stays consistent if it is ever edited:

* Panels (b) and (c) share a colour key, so a reader can carry identity from
  one to the other without a second legend.
* The two relaxation modes keep the blue/orange pair used in
  ``plot_throughput.py``; single point is neutral grey and dashed, because it
  is a different workload rather than a third peer in the fixed/refill
  comparison. Hue and dash both carry it, so colour alone is never load-bearing.
* Panel (a) has two curves and room, so it keeps the direct end labels used by
  the sibling figure. Panels (b) and (c) carry three series each in a third of
  the width, where end labels either collide with the curves or fall off the
  axes, so those use a compact frameless legend placed in the corner each
  panel leaves empty.
* Time and memory are logarithmic: both span more than a decade, and the
  interesting behaviour is the change in slope, not the absolute values.
* Values are labelled only where they carry the headline -- the endpoints of
  panel (a) and the card limit in panel (c).
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
FIXED = "#2a78d6"
REFILL = "#eb6834"
SINGLE = "#6f6e6a"

FIXED_MODE = "fixed cell, fixed batches"
REFILL_MODE = "fixed cell, refill"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(prog="plot_batching.py", allow_abbrev=False)
    parser.add_argument("-data", type=pathlib.Path,
                        default=ROOT / "benchmarks" / "results" / "sod_mace_benchmark_710.json")
    parser.add_argument("-out", type=pathlib.Path,
                        default=ROOT / "docs" / "_static" / "mace_batching.svg")
    return parser.parse_args()


def style(ax, sizes) -> None:
    ax.set_facecolor(SURFACE)
    ax.grid(True, which="major", axis="y", color=GRID, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
        ax.spines[side].set_linewidth(0.8)
    ax.set_xscale("log", base=2)
    ax.xaxis.set_major_locator(FixedLocator(sizes))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v)}"))
    ax.set_xlim(0.93, max(sizes) * 1.06)
    ax.tick_params(colors=INK_2, labelsize=9, length=0)
    ax.set_xlabel("batch size", color=INK_2, fontsize=9.5, labelpad=6)


def line(ax, points, color, dashed=False, label=None):
    xs, ys = zip(*points)
    ax.plot(xs, ys, color=color, linewidth=2.0, marker="o", markersize=4.5,
            linestyle="--" if dashed else "-", label=label,
            markerfacecolor=color, markeredgecolor=SURFACE, markeredgewidth=1.2,
            zorder=3, clip_on=False)
    return xs, ys


def key(ax, loc):
    leg = ax.legend(loc=loc, frameon=False, fontsize=9, handlelength=1.8,
                    borderpad=0.2, labelspacing=0.35, handletextpad=0.6)
    for text, handle in zip(leg.get_texts(), leg.legend_handles):
        text.set_color(handle.get_color())
        text.set_fontweight("bold")


def label_end(ax, xs, ys, text, color, dy=0):
    ax.annotate(f" {text}", (xs[-1], ys[-1]), color=color, fontsize=9.5,
                fontweight="bold", va="center", ha="left",
                textcoords="offset points", xytext=(4, dy), annotation_clip=False)


def main() -> int:
    args = parse_args()
    data = json.loads(args.data.read_text())
    total_mib = data["gpu_total_mib"]

    relax = {(r["mode"], r["batch"]): r for r in data["relaxation"]}
    single = {r["batch"]: r for r in data["single_point"] if r["backend"] == "cueq"}
    smem = {m["batch"]: m for m in data["memory"] if m["backend"] == "cueq"}

    def by_batch(mode, key):
        return sorted(((b, r[key]) for (m, b), r in relax.items() if m == mode),
                      key=lambda p: p[0])

    fixed_t = by_batch(FIXED_MODE, "best_s")
    refill_t = by_batch(REFILL_MODE, "best_s")
    fixed_m = by_batch(FIXED_MODE, "peak_mib")
    refill_m = by_batch(REFILL_MODE, "peak_mib")
    single_t = sorted((b, r["best_s"]) for b, r in single.items())
    single_m = sorted((b, m["peak_mib"]) for b, m in smem.items())
    sizes = [b for b, _ in fixed_t]

    base = dict(fixed_t)[1]
    fixed_s = [(b, base / t) for b, t in fixed_t]
    refill_s = [(b, base / t) for b, t in refill_t]

    fig, axes = plt.subplots(1, 3, figsize=(13.2, 4.3))
    fig.patch.set_facecolor(SURFACE)
    fig.subplots_adjust(left=0.055, right=0.985, top=0.775, bottom=0.135, wspace=0.26)
    a, b, c = axes
    for ax in axes:
        style(ax, sizes)

    # (a) relaxation speedup -------------------------------------------------
    xs, ys = line(a, fixed_s, FIXED, label="fixed")
    xs2, ys2 = line(a, refill_s, REFILL, label="refill")
    key(a, "upper left")
    a.set_ylabel("speedup vs batch 1", color=INK_2, fontsize=9.5, labelpad=6)
    a.set_ylim(0, max(ys2) * 1.16)
    a.annotate(f"{ys2[-1]:.1f}x", (xs2[-1], ys2[-1]), textcoords="offset points",
               xytext=(-2, 9), ha="right", color=INK, fontsize=10, fontweight="bold")
    a.set_title("(a)  relaxation speedup", color=INK, fontsize=10.5,
                fontweight="bold", loc="left", pad=10)

    # (b) wall-clock time ----------------------------------------------------
    line(b, fixed_t, FIXED, label="relax, fixed")
    line(b, refill_t, REFILL, label="relax, refill")
    line(b, single_t, SINGLE, dashed=True, label="single point")
    key(b, "upper right")
    b.set_yscale("log")
    b.set_ylabel(f"wall-clock time for {data['n_configs']} structures (s)",
                 color=INK_2, fontsize=9.5, labelpad=2)
    b.yaxis.set_major_formatter(FuncFormatter(
        lambda v, _: f"{v:.0f}" if v >= 1 else f"{v:g}"))
    b.set_title("(b)  time to process the level", color=INK, fontsize=10.5,
                fontweight="bold", loc="left", pad=10)

    # (c) peak memory --------------------------------------------------------
    c.axhline(total_mib, color=MUTED, linewidth=1.0, linestyle=(0, (4, 3)), zorder=1)
    c.annotate(f"card limit, {total_mib/1024:.0f} GiB", (sizes[0], total_mib),
               textcoords="offset points", xytext=(2, 5), ha="left",
               color=MUTED, fontsize=8.5)
    line(c, single_m, SINGLE, dashed=True, label="single point")
    line(c, refill_m, REFILL, label="relax, refill")
    line(c, fixed_m, FIXED, label="relax, fixed")
    key(c, "lower right")
    c.set_yscale("log")
    c.set_ylabel("peak GPU memory (MiB)", color=INK_2, fontsize=9.5, labelpad=2)
    c.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}"))
    c.set_ylim(top=total_mib * 2.1)
    c.set_title("(c)  peak GPU memory", color=INK, fontsize=10.5,
                fontweight="bold", loc="left", pad=10)

    fig.text(0.055, 0.955, "Batch size buys speed with diminishing returns, "
             "at a memory cost that stays linear",
             color=INK, fontsize=12.5, fontweight="bold", ha="left", va="top")
    fig.text(0.055, 0.895,
             f"{data['n_configs']} configurations of {data['n_atoms']} atoms (MgO:Ni) - "
             f"MACE {data['model']}, cuEquivariance, float32, "
             f"{data['gpu'].replace('NVIDIA GeForce ', '')} - "
             f"best of {data['repeats']} runs; relaxation to fmax 0.05 eV/A",
             color=INK_2, fontsize=9, ha="left", va="top")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, format=args.out.suffix.lstrip("."), facecolor=SURFACE, dpi=150)
    print(f"wrote {args.out}")
    for name, pts in (("fixed", fixed_s), ("refill", refill_s)):
        print(f"  {name:6s} speedup {pts[0][1]:.2f}x -> {pts[-1][1]:.2f}x "
              f"(batch {pts[0][0]} -> {pts[-1][0]})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
