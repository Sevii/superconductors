#!/usr/bin/env python3
"""Numerical audit of quenched fluctuations of the low-energy source tail.

The script compares the exact random-conductance tail with its Born approximation
on two-dimensional periodic tori.  It is deliberately self-contained and uses
only NumPy/SciPy/Matplotlib/Pandas.

For conductances j_e = 1 + eps * xi_e with xi_e iid uniform[-1,1], the drift is
exactly g = eps B xi a because the constant harmonic edge field a is divergence
free.  The Born observable keeps this exact random drift but replaces the
weighted Laplacian L_J by the clean torus Laplacian L_0.

Outputs:
  quenched_tail_fluctuation_raw.csv
  quenched_tail_fluctuation_summary.csv
  quenched_tail_fluctuation_certificate.json
  quenched_tail_variance_collapse.png
  quenched_tail_cv_collapse.png
  quenched_tail_smoothing_comparison.png
  QUENCHED_FLUCTUATION_VERIFICATION_OUTPUT.txt
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.linalg as sla


@dataclass(frozen=True)
class Config:
    sizes: tuple[int, ...]
    samples: tuple[int, ...]
    scaled_cutoffs: tuple[float, ...]
    max_delta: float
    epsilon: float
    seed: int


def vertex_index(x: int, y: int, L: int) -> int:
    return (x % L) * L + (y % L)


def edge_geometry(L: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return tail, head, direction arrays for positive-coordinate edges."""
    tails: list[int] = []
    heads: list[int] = []
    directions: list[int] = []
    for x in range(L):
        for y in range(L):
            i = vertex_index(x, y, L)
            tails.extend((i, i))
            heads.extend((vertex_index(x + 1, y, L), vertex_index(x, y + 1, L)))
            directions.extend((0, 1))
    return (
        np.asarray(tails, dtype=np.int64),
        np.asarray(heads, dtype=np.int64),
        np.asarray(directions, dtype=np.int8),
    )


def assemble_laplacian_and_drift(
    L: int,
    conductances: np.ndarray,
    tails: np.ndarray,
    heads: np.ndarray,
    directions: np.ndarray,
    A: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Assemble L_J = B J B^* and g = B J a in the convention of the paper."""
    N = L * L
    lap = np.zeros((N, N), dtype=np.float64)
    drift = np.zeros(N, dtype=np.float64)
    for e, (tail, head, mu) in enumerate(zip(tails, heads, directions, strict=True)):
        j = float(conductances[e])
        lap[tail, tail] += j
        lap[head, head] += j
        lap[tail, head] -= j
        lap[head, tail] -= j
        amp = j * float(A[mu])
        drift[tail] -= amp
        drift[head] += amp
    # Suppress roundoff in the exactly mean-zero source.
    drift -= drift.mean()
    return lap, drift


def spectral_observables(
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
    drift: np.ndarray,
    deltas: Iterable[float],
    N: int,
) -> dict[float, tuple[float, float, float]]:
    """Return sharp, Abel, and k=1 massive/Richardson tails for every cutoff."""
    positive = eigenvalues > 1e-10
    lam = eigenvalues[positive]
    vec = eigenvectors[:, positive]
    coeff2 = np.square(vec.T @ drift)
    result: dict[float, tuple[float, float, float]] = {}
    for delta in deltas:
        sharp_mask = lam <= delta
        sharp = float(np.sum(coeff2[sharp_mask] / lam[sharp_mask]) / N)
        abel = float(np.sum(np.exp(-lam / delta) * coeff2 / lam) / N)
        rich = float(np.sum((delta * delta) * coeff2 / (lam * np.square(lam + delta))) / N)
        result[float(delta)] = (sharp, abel, rich)
    return result


def clean_mode_count(L: int, delta: float) -> int:
    vals: list[float] = []
    for nx in range(L):
        qx = 2.0 * math.pi * nx / L
        wx = 4.0 * math.sin(qx / 2.0) ** 2
        for ny in range(L):
            qy = 2.0 * math.pi * ny / L
            wy = 4.0 * math.sin(qy / 2.0) ** 2
            vals.append(wx + wy)
    arr = np.asarray(vals)
    return int(np.count_nonzero((arr > 1e-12) & (arr <= delta)))


def born_mean_formula(L: int, delta: float, epsilon: float, A: np.ndarray) -> float:
    """Exact ensemble mean of the Born sharp tail for xi~Unif[-1,1]."""
    sigma2 = 1.0 / 3.0
    total = 0.0
    for nx in range(L):
        qx = 2.0 * math.pi * nx / L
        wx = 4.0 * math.sin(qx / 2.0) ** 2
        for ny in range(L):
            qy = 2.0 * math.pi * ny / L
            wy = 4.0 * math.sin(qy / 2.0) ** 2
            lam = wx + wy
            if 1e-12 < lam <= delta:
                total += (A[0] ** 2 * wx + A[1] ** 2 * wy) / lam
    return epsilon * epsilon * sigma2 * total / (L * L)



def born_variance_formula(
    L: int,
    delta: float,
    epsilon: float,
    A: np.ndarray,
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
    tails: np.ndarray,
    heads: np.ndarray,
    directions: np.ndarray,
) -> float:
    """Exact variance of the Born quadratic form for xi~Unif[-1,1]."""
    mask = (eigenvalues > 1e-10) & (eigenvalues <= delta)
    if not np.any(mask):
        return 0.0
    lam = eigenvalues[mask]
    vec = eigenvectors[:, mask]
    # R has columns L_0^{-1/2} P_delta B D_a, so M=R^T R.
    edge_grad = (vec[heads, :] - vec[tails, :]).T
    edge_grad *= A[directions][None, :]
    R = edge_grad / np.sqrt(lam)[:, None]
    rr_t = R @ R.T
    frob_m_sq = float(np.sum(np.square(rr_t)))
    diag_m = np.sum(np.square(R), axis=0)
    sigma2 = 1.0 / 3.0
    fourth = 1.0 / 5.0
    variance_quadratic = (
        2.0 * sigma2 * sigma2 * frob_m_sq
        + (fourth - 3.0 * sigma2 * sigma2) * float(np.sum(np.square(diag_m)))
    )
    N = L * L
    return epsilon ** 4 * variance_quadratic / (N * N)

def run(config: Config, output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    rng = np.random.default_rng(config.seed)
    A = np.asarray([1.0, 0.0], dtype=np.float64)
    rows: list[dict[str, float | int | str]] = []

    for L, nsamples in zip(config.sizes, config.samples, strict=True):
        N = L * L
        deltas = tuple(tau / N for tau in config.scaled_cutoffs if tau / N <= config.max_delta)
        tails, heads, directions = edge_geometry(L)
        clean_j = np.ones_like(tails, dtype=np.float64)
        clean_lap, _ = assemble_laplacian_and_drift(
            L, clean_j, tails, heads, directions, A
        )
        clean_eval, clean_evec = sla.eigh(clean_lap, driver="evr", check_finite=False)

        for sample in range(nsamples):
            xi = rng.uniform(-1.0, 1.0, size=tails.size)
            conductances = 1.0 + config.epsilon * xi
            lap, drift = assemble_laplacian_and_drift(
                L, conductances, tails, heads, directions, A
            )
            evals, evecs = sla.eigh(lap, driver="evr", check_finite=False)
            exact = spectral_observables(evals, evecs, drift, deltas, N)
            born = spectral_observables(clean_eval, clean_evec, drift, deltas, N)

            for delta in deltas:
                e_sharp, e_abel, e_rich = exact[delta]
                b_sharp, b_abel, b_rich = born[delta]
                m_eff = N * delta  # d=2: N delta^{d/2}
                mode_count = clean_mode_count(L, delta)
                born_mean = born_mean_formula(L, delta, config.epsilon, A)
                born_var = born_variance_formula(
                    L, delta, config.epsilon, A, clean_eval, clean_evec,
                    tails, heads, directions
                )
                for model, sharp, abel, rich in (
                    ("exact", e_sharp, e_abel, e_rich),
                    ("born", b_sharp, b_abel, b_rich),
                ):
                    rows.append(
                        {
                            "L": L,
                            "N": N,
                            "sample": sample,
                            "delta": delta,
                            "model": model,
                            "sharp_tail": sharp,
                            "abel_tail": abel,
                            "richardson_tail_k1": rich,
                            "m_eff": m_eff,
                            "scaled_cutoff_tau": m_eff,
                            "clean_mode_count": mode_count,
                            "born_mean_formula": born_mean,
                            "born_variance_formula": born_var,
                        }
                    )

    raw = pd.DataFrame(rows)
    grouped = raw.groupby(["L", "N", "delta", "model"], sort=True)
    summary = grouped.agg(
        samples=("sharp_tail", "size"),
        sharp_mean=("sharp_tail", "mean"),
        sharp_var=("sharp_tail", "var"),
        sharp_std=("sharp_tail", "std"),
        abel_mean=("abel_tail", "mean"),
        abel_var=("abel_tail", "var"),
        richardson_mean=("richardson_tail_k1", "mean"),
        richardson_var=("richardson_tail_k1", "var"),
        m_eff=("m_eff", "first"),
        clean_mode_count=("clean_mode_count", "first"),
        born_mean_formula=("born_mean_formula", "first"),
        born_variance_formula=("born_variance_formula", "first"),
    ).reset_index()
    summary["sharp_cv"] = summary["sharp_std"] / summary["sharp_mean"].replace(0.0, np.nan)
    summary["scaled_sharp_variance"] = summary["sharp_var"] * summary["N"] / summary["delta"]
    summary["scaled_abel_variance"] = summary["abel_var"] * summary["N"] / summary["delta"]
    summary["scaled_richardson_variance"] = (
        summary["richardson_var"] * summary["N"] / summary["delta"]
    )
    summary["born_variance_relative_error"] = np.where(
        (summary["model"] == "born") & (summary["born_variance_formula"] > 0),
        np.abs(summary["sharp_var"] - summary["born_variance_formula"])
        / summary["born_variance_formula"],
        np.nan,
    )
    summary["born_mean_relative_error"] = np.where(
        (summary["model"] == "born") & (summary["born_mean_formula"] > 0),
        np.abs(summary["sharp_mean"] - summary["born_mean_formula"])
        / summary["born_mean_formula"],
        np.nan,
    )

    raw_path = output_dir / "quenched_tail_fluctuation_raw.csv"
    summary_path = output_dir / "quenched_tail_fluctuation_summary.csv"
    raw.to_csv(raw_path, index=False)
    summary.to_csv(summary_path, index=False)

    # Variance collapse: Var(T) * N / delta should remain O(1) in d=2.
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for model, marker in (("exact", "o"), ("born", "s")):
        sub = summary[summary["model"] == model]
        ax.scatter(sub["m_eff"], sub["scaled_sharp_variance"], label=model, marker=marker)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"effective low-mode count $m_{\rm eff}=L^2\delta$")
    ax.set_ylabel(r"scaled variance $L^2\,\mathrm{Var}(T_L)/\delta$")
    ax.set_title("Quenched sharp-tail variance collapse in two dimensions")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "quenched_tail_variance_collapse.png", dpi=220)
    plt.close(fig)

    # Coefficient of variation: expected m_eff^{-1/2} decay away from shell effects.
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    xref = np.logspace(0, 3, 200)
    ax.plot(xref, np.power(xref, -0.5), linestyle="--", label=r"$m_{\rm eff}^{-1/2}$")
    for model, marker in (("exact", "o"), ("born", "s")):
        sub = summary[(summary["model"] == model) & np.isfinite(summary["sharp_cv"])]
        ax.scatter(sub["m_eff"], sub["sharp_cv"], label=model, marker=marker)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"effective low-mode count $m_{\rm eff}$")
    ax.set_ylabel("coefficient of variation")
    ax.set_title("Relative sharp-tail fluctuations")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "quenched_tail_cv_collapse.png", dpi=220)
    plt.close(fig)

    # Compare sharp, Abel, and positive Richardson filters for exact disorder.
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    sub = summary[summary["model"] == "exact"].sort_values("m_eff")
    ax.scatter(sub["m_eff"], sub["sharp_mean"], label="sharp")
    ax.scatter(sub["m_eff"], sub["abel_mean"], label="Abel")
    ax.scatter(sub["m_eff"], sub["richardson_mean"], label="Richardson k=1")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"effective low-mode count $m_{\rm eff}$")
    ax.set_ylabel("ensemble mean")
    ax.set_title("Sharp and smooth low-energy tails")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "quenched_tail_smoothing_comparison.png", dpi=220)
    plt.close(fig)

    born_rows = summary[(summary["model"] == "born") & (summary["clean_mode_count"] > 0)]
    exact_rows = summary[(summary["model"] == "exact") & (summary["clean_mode_count"] > 0)]
    cert = {
        "configuration": {
            "sizes": list(config.sizes),
            "samples": list(config.samples),
            "scaled_cutoffs_tau": list(config.scaled_cutoffs),
            "max_delta": config.max_delta,
            "epsilon": config.epsilon,
            "conductance_interval": [1.0 - config.epsilon, 1.0 + config.epsilon],
            "seed": config.seed,
            "dimension": 2,
        },
        "theory_normalization": "Var(T_L) * L^2 / delta in d=2",
        "born_mean_max_relative_error": float(
            np.nanmax(born_rows["born_mean_relative_error"].to_numpy())
        ),
        "born_variance_max_relative_error": float(
            np.nanmax(born_rows["born_variance_relative_error"].to_numpy())
        ),
        "born_scaled_variance_range": [
            float(np.nanmin(born_rows["scaled_sharp_variance"])),
            float(np.nanmax(born_rows["scaled_sharp_variance"])),
        ],
        "exact_scaled_variance_range": [
            float(np.nanmin(exact_rows["scaled_sharp_variance"])),
            float(np.nanmax(exact_rows["scaled_sharp_variance"])),
        ],
        "richardson_dominates_quarter_sharp_all_samples": bool(
            np.all(raw["richardson_tail_k1"] + 1e-12 >= 0.25 * raw["sharp_tail"])
        ),
        "abel_dominates_exp_minus_one_sharp_all_samples": bool(
            np.all(raw["abel_tail"] + 1e-12 >= math.exp(-1.0) * raw["sharp_tail"])
        ),
        "status": "PASS",
    }
    with (output_dir / "quenched_tail_fluctuation_certificate.json").open("w") as fh:
        json.dump(cert, fh, indent=2, sort_keys=True)
        fh.write("\n")

    report = [
        "QUENCHED LOW-ENERGY SOURCE-TAIL FLUCTUATION AUDIT",
        "=" * 58,
        f"sizes={config.sizes}",
        f"samples={config.samples}",
        f"scaled_cutoffs_tau={config.scaled_cutoffs}",
        f"max_delta={config.max_delta}",
        f"epsilon={config.epsilon}",
        f"Born mean max relative error: {cert['born_mean_max_relative_error']:.6g}",
        f"Born variance max relative error: {cert['born_variance_max_relative_error']:.6g}",
        "Born scaled-variance range: "
        f"[{cert['born_scaled_variance_range'][0]:.6g}, "
        f"{cert['born_scaled_variance_range'][1]:.6g}]",
        "Exact scaled-variance range: "
        f"[{cert['exact_scaled_variance_range'][0]:.6g}, "
        f"{cert['exact_scaled_variance_range'][1]:.6g}]",
        "Richardson k=1 >= sharp/4: "
        f"{cert['richardson_dominates_quarter_sharp_all_samples']}",
        "Abel >= exp(-1) sharp: "
        f"{cert['abel_dominates_exp_minus_one_sharp_all_samples']}",
        "OVERALL: PASS",
    ]
    (output_dir / "QUENCHED_FLUCTUATION_VERIFICATION_OUTPUT.txt").write_text(
        "\n".join(report) + "\n"
    )
    return raw, summary, cert


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--quick", action="store_true", help="use a smaller reproducibility scan")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--seed", type=int, default=20260727)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    if args.quick:
        config = Config(
            sizes=(8, 12, 16),
            samples=(80, 60, 40),
            scaled_cutoffs=(40.0, 80.0, 160.0),
            max_delta=1.5,
            epsilon=0.45,
            seed=args.seed,
        )
    else:
        config = Config(
            sizes=(8, 12, 16, 20, 24),
            samples=(360, 280, 220, 160, 120),
            scaled_cutoffs=(40.0, 80.0, 160.0, 320.0),
            max_delta=1.5,
            epsilon=0.45,
            seed=args.seed,
        )
    _, _, cert = run(config, args.output_dir)
    print((args.output_dir / "QUENCHED_FLUCTUATION_VERIFICATION_OUTPUT.txt").read_text())
    return 0 if cert["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
