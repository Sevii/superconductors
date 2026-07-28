#!/usr/bin/env python3
"""Plot the finite-torus source-tail and heat-decay summary tables."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="") as handle:
        return [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]


def plot_tail(rows: list[dict[str, float]], outdir: Path) -> None:
    grouped: dict[int, list[dict[str, float]]] = defaultdict(list)
    for row in rows:
        grouped[int(row["L"])].append(row)

    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    for L in sorted(grouped):
        data = sorted(grouped[L], key=lambda row: row["delta"])
        x = np.array([row["delta"] for row in data])
        y = np.array([row["mean_Hminus1_tail"] for row in data])
        positive = y > 0
        ax.plot(x[positive], y[positive], marker="o", label=f"L={L}")

    largest = max(grouped)
    largest_data = sorted(grouped[largest], key=lambda row: row["delta"])
    x_ref = np.array([row["delta"] for row in largest_data])
    y_ref_data = np.array([row["mean_Hminus1_tail"] for row in largest_data])
    positive = y_ref_data > 0
    if np.count_nonzero(positive):
        anchor = np.median(y_ref_data[positive] / x_ref[positive])
        ax.plot(x_ref, anchor * x_ref, linestyle="--", label="reference slope +1")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"spectral cutoff $\delta$")
    ax.set_ylabel(r"mean $L^{-2}\langle g,L^+\mathbf{1}_{(0,\delta]}g\rangle$")
    ax.set_title("Two-dimensional finite-torus source tail")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "finite_torus_source_tail_scaling.png", dpi=220)
    plt.close(fig)


def plot_heat(rows: list[dict[str, float]], outdir: Path) -> None:
    grouped: dict[int, list[dict[str, float]]] = defaultdict(list)
    for row in rows:
        grouped[int(row["L"])].append(row)

    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    for L in sorted(grouped):
        data = sorted(grouped[L], key=lambda row: row["time"])
        x = np.array([row["time"] for row in data])
        y = np.array([row["mean_heat_correlation"] for row in data])
        positive = y > 0
        ax.plot(x[positive], y[positive], marker="o", label=f"L={L}")

    largest = max(grouped)
    largest_data = sorted(grouped[largest], key=lambda row: row["time"])
    x_ref = np.array([row["time"] for row in largest_data])
    y_ref_data = np.array([row["mean_heat_correlation"] for row in largest_data])
    late = (x_ref >= 8) & (y_ref_data > 0)
    if np.count_nonzero(late):
        anchor = np.median(y_ref_data[late] * x_ref[late] ** 2)
        ax.plot(x_ref, anchor * x_ref ** -2, linestyle="--", label="reference slope -2")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"time $t$")
    ax.set_ylabel(r"mean $L^{-2}\langle g,e^{-tL}g\rangle$")
    ax.set_title("Two-dimensional drift heat correlation")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "finite_torus_heat_decay_scaling.png", dpi=220)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--indir", type=Path, default=Path("."))
    parser.add_argument("--outdir", type=Path, default=Path("."))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    plot_tail(
        read_rows(args.indir / "finite_torus_source_tail_summary.csv"),
        args.outdir,
    )
    plot_heat(
        read_rows(args.indir / "finite_torus_heat_decay_summary.csv"),
        args.outdir,
    )
    print("Wrote finite_torus_source_tail_scaling.png")
    print("Wrote finite_torus_heat_decay_scaling.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
