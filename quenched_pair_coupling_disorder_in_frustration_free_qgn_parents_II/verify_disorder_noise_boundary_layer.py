#!/usr/bin/env python3
"""Audit the fixed-cutoff disorder-noise boundary layer.

The exact Gaussian half-space calculation verifies the integrable t^{-1/2}
noise-energy singularity used as the model for the sharp random-conductance
cutoff.  A finite Bernoulli-majority calculation gives a discrete product-space
analogue under the edge-resampling semigroup.

The calculations are diagnostics/solvable models, not a proof of the random-
conductance theorem.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.special import gammaln
from scipy.stats import binom


def gaussian_energy(t: float, q: float = 0.7) -> float:
    """E|grad P_t(q 1_{Z<=0})|^2 for the Gaussian OU semigroup."""
    if t <= 0.0:
        return math.inf
    # Stable form of a_t^2 = 1/(exp(2t)-1).
    if t > 350.0:
        a2 = math.exp(-2.0 * t)
    else:
        a2 = 1.0 / math.expm1(2.0 * t)
    return q * q * a2 / (2.0 * math.pi * math.sqrt(1.0 + 2.0 * a2))


def gaussian_table(q: float = 0.7) -> tuple[pd.DataFrame, dict[str, float]]:
    ts = np.geomspace(1e-7, 10.0, 120)
    rows = []
    for t in ts:
        e = gaussian_energy(float(t), q)
        rows.append(
            {
                "t": float(t),
                "noise_energy": e,
                "sqrt_t_times_energy": math.sqrt(float(t)) * e,
                "small_time_limit": q * q / (4.0 * math.pi),
            }
        )
    integral, err = quad(lambda s: 2.0 * gaussian_energy(s, q), 0.0, np.inf, epsabs=1e-12, epsrel=1e-11, limit=300)
    exact_var = q * q / 4.0
    return pd.DataFrame(rows), {
        "q": q,
        "twice_integrated_noise_energy": float(integral),
        "quadrature_error_estimate": float(err),
        "exact_variance": exact_var,
        "absolute_closure_error": abs(float(integral) - exact_var),
        "small_time_limit": q * q / (4.0 * math.pi),
    }


def binomial_weights(M: int) -> np.ndarray:
    k = np.arange(M + 1)
    logw = gammaln(M + 1.0) - gammaln(k + 1.0) - gammaln(M - k + 1.0) - M * math.log(2.0)
    return np.exp(logw)


def majority_smoothed_values(M: int, t: float, q: float = 0.7) -> np.ndarray:
    """P_t[q 1_{#plus <= (M-1)/2}] as a function of initial #plus."""
    if M % 2 != 1:
        raise ValueError("M must be odd")
    rho = math.exp(-t)
    p_keep_plus = 0.5 * (1.0 + rho)
    p_flip_minus = 0.5 * (1.0 - rho)
    cutoff = (M - 1) // 2
    values = np.empty(M + 1, dtype=float)
    for k in range(M + 1):
        # Number of plus spins after noise is X+Y with
        # X~Bin(k,p_keep_plus), Y~Bin(M-k,p_flip_minus).
        x = np.arange(k + 1)
        pmf_x = binom.pmf(x, k, p_keep_plus)
        y = np.arange(M - k + 1)
        pmf_y = binom.pmf(y, M - k, p_flip_minus)
        pmf = np.convolve(pmf_x, pmf_y)
        values[k] = q * float(np.sum(pmf[: cutoff + 1]))
    return values


def majority_noise_energy(M: int, t: float, q: float = 0.7) -> float:
    f = majority_smoothed_values(M, t, q)
    w = binomial_weights(M)
    total = 0.0
    for k in range(M + 1):
        if k > 0:
            total += w[k] * k * (f[k] - f[k - 1]) ** 2 / 4.0
        if k < M:
            total += w[k] * (M - k) * (f[k + 1] - f[k]) ** 2 / 4.0
    return float(total)


def majority_table(q: float = 0.7, quick: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    Ms = [15, 31, 63] if quick else [15, 31, 63, 127]
    ts = np.geomspace(2e-4 if quick else 2e-5, 4.0, 55 if quick else 80)
    rows: list[dict[str, float | int]] = []
    summaries: list[dict[str, float | int]] = []
    for M in Ms:
        energies = []
        for t in ts:
            e = majority_noise_energy(M, float(t), q)
            energies.append(e)
            rows.append(
                {
                    "M": M,
                    "t": float(t),
                    "noise_energy": e,
                    "sqrt_t_times_energy": math.sqrt(float(t)) * e,
                    "energy_over_q2": e / (q * q),
                    "scaled_by_min_profile": e
                    / (q * q * min(math.sqrt(M), 1.0 / math.sqrt(float(t)))),
                }
            )
        # Numerically integrate on a denser logarithmic grid.  The omitted
        # interval [0,t_min] is bounded using the t=0 total influence; the
        # large-time tail is exponentially tiny at t=12.
        grid = np.geomspace(1e-7, 12.0, 220 if quick else 280)
        vals = np.array([majority_noise_energy(M, float(t), q) for t in grid])
        integral = 2.0 * float(np.trapezoid(vals, grid))
        # Add the very small [0,1e-7] rectangle; this is much below displayed precision.
        integral += 2.0 * grid[0] * majority_noise_energy(M, 0.0, q)
        summaries.append(
            {
                "M": M,
                "twice_integrated_noise_energy_numeric": integral,
                "exact_variance": q * q / 4.0,
                "absolute_closure_error": abs(integral - q * q / 4.0),
                "max_profile_ratio": max(
                    e / (q * q * min(math.sqrt(M), 1.0 / math.sqrt(float(t))))
                    for e, t in zip(energies, ts, strict=True)
                ),
            }
        )
    return pd.DataFrame(rows), pd.DataFrame(summaries)


def make_figure(gauss: pd.DataFrame, majority: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.7))
    ax.loglog(gauss["t"], gauss["noise_energy"], linewidth=2.0, label="Gaussian half-space")
    for M, group in majority.groupby("M"):
        ax.loglog(group["t"], group["noise_energy"], marker="o", markersize=2.6, linewidth=1.0, label=f"Bernoulli majority M={M}")
    t = np.geomspace(3e-4, 1.0, 100)
    q = 0.7
    ax.loglog(t, q * q / (4.0 * math.pi * np.sqrt(t)), linestyle="--", label=r"$q^2/(4\pi\sqrt{t})$")
    ax.set_xlabel("disorder-noise time t")
    ax.set_ylabel("Dirichlet/noise energy")
    ax.set_title("Integrable fixed-cutoff boundary layer")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    gauss, gaussian_summary = gaussian_table()
    majority, majority_summary = majority_table(quick=args.quick)
    gauss.to_csv(args.outdir / "gaussian_halfspace_noise_energy.csv", index=False)
    majority.to_csv(args.outdir / "bernoulli_majority_noise_energy.csv", index=False)
    majority_summary.to_csv(args.outdir / "bernoulli_majority_noise_energy_summary.csv", index=False)
    make_figure(gauss, majority, args.outdir / "disorder_noise_boundary_layer.png")

    checks = {
        "gaussian_variance_closure": gaussian_summary["absolute_closure_error"] < 2e-10,
        "gaussian_small_time_profile": abs(
            float(gauss.iloc[0]["sqrt_t_times_energy"]) - gaussian_summary["small_time_limit"]
        ) < 2e-5,
        "majority_variance_closure": float(majority_summary["absolute_closure_error"].max()) < (2e-3 if args.quick else 8e-4),
        "majority_min_profile_uniform": float(majority_summary["max_profile_ratio"].max()) < 0.5,
    }
    cert = {
        "quick": bool(args.quick),
        "gaussian": gaussian_summary,
        "majority": majority_summary.to_dict(orient="records"),
        "checks": checks,
        "all_checks_pass": bool(all(checks.values())),
        "interpretation": (
            "The Gaussian identity is exact. The Bernoulli-majority computation is an exact finite product-space model evaluated numerically and illustrates the same integrable boundary layer. Neither is a proof for the random-conductance spectral tail."
        ),
    }
    (args.outdir / "disorder_noise_boundary_layer_certificate.json").write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n")

    lines = [
        "FIXED-CUTOFF DISORDER-NOISE BOUNDARY-LAYER AUDIT",
        "=" * 56,
        f"quick: {args.quick}",
        "",
        "Gaussian half-space closure:",
        pd.Series(gaussian_summary).to_string(),
        "",
        "Bernoulli-majority summaries:",
        majority_summary.to_string(index=False),
        "",
        "Checks:",
    ]
    lines.extend(f"  {k}: {'PASS' if v else 'FAIL'}" for k, v in checks.items())
    lines.extend(["", f"OVERALL: {'PASS' if cert['all_checks_pass'] else 'FAIL'}"])
    text = "\n".join(lines) + "\n"
    (args.outdir / "DISORDER_NOISE_BOUNDARY_LAYER_VERIFICATION_OUTPUT.txt").write_text(text)
    print(text)
    return 0 if cert["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
