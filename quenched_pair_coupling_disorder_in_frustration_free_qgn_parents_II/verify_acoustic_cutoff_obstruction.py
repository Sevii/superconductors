#!/usr/bin/env python3
"""Audit the acoustic-eigenvalue obstruction to a uniform sharp-shell modulus.

The script has two logically separate parts.

1. An exact one-mode Gaussian model verifies that CLT-scale fluctuations of an
   acoustic eigenvalue are compatible with the desired centered variance but
   incompatible with any volume-uniform Hölder shell modulus at the moving
   cutoff scale eta_N ~ N^{-1/2}.

2. A finite-torus random-conductance diagnostic measures the same three scales
   in the actual weighted Laplacian: relative low-eigenvalue width, the
   current-participation density r*pi_L(r), and Var T_L(r).  The numerical fit
   is evidence only; it is not used in any proof.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.linalg as sla
from scipy.special import ndtr

from verify_sharp_cutoff_three_routes import laplacian_and_drift, torus_edges


def low_modes(
    lap: np.ndarray, drift: np.ndarray, count: int = 8
) -> tuple[np.ndarray, np.ndarray]:
    n = lap.shape[0]
    hi = min(n - 1, count)
    vals, vecs = sla.eigh(
        lap, subset_by_index=[0, hi], driver="evr", check_finite=False
    )
    vals = vals[1:]
    vecs = vecs[:, 1:]
    q = np.square(vecs.T @ drift) / (n * vals)
    return vals, q


def gaussian_one_mode_table() -> pd.DataFrame:
    a = 1.0
    q_star = 0.7
    alpha = 0.5
    probability = float(ndtr(a) - 0.5)
    density_zero = 1.0 / math.sqrt(2.0 * math.pi)
    rows: list[dict[str, float]] = []
    for n in (16, 64, 256, 1024, 4096, 16384):
        eta = math.exp(a / math.sqrt(n)) - 1.0
        q = q_star / n
        shell_second = q * q * probability
        centered_variance = q * q * 0.25
        participation_density = q * q * math.sqrt(n) * density_zero
        nominal_variance_scale = q_star * q_star / (n * n)
        holder_rhs = nominal_variance_scale * eta**alpha
        rows.append(
            {
                "N": n,
                "a": a,
                "alpha": alpha,
                "eta_N": eta,
                "shell_second_moment": shell_second,
                "centered_variance": centered_variance,
                "nominal_N_minus_2_scale": nominal_variance_scale,
                "holder_rhs_without_constant": holder_rhs,
                "shell_to_holder_ratio": shell_second / holder_rhs,
                "r_pi_at_center": participation_density,
                "r_pi_to_N_minus_2_ratio": participation_density
                / nominal_variance_scale,
            }
        )
    return pd.DataFrame(rows)


def random_conductance_scaling(
    rng: np.random.Generator, quick: bool
) -> tuple[pd.DataFrame, dict[str, float]]:
    schedule = (
        [(4, 3500), (5, 2600), (6, 1900), (8, 900), (10, 420)]
        if quick
        else [(4, 12000), (5, 9000), (6, 6000), (8, 2500), (10, 1100)]
    )
    rows: list[dict[str, float]] = []
    for L, samples in schedule:
        n = L * L
        tails, heads, dirs = torus_edges(L)
        all_vals: list[np.ndarray] = []
        all_q: list[np.ndarray] = []
        samples_data: list[tuple[np.ndarray, np.ndarray]] = []
        for _ in range(samples):
            j = rng.uniform(0.55, 1.45, size=tails.size)
            lap, drift = laplacian_and_drift(
                n, tails, heads, dirs, j, np.array([1.0, 0.0])
            )
            vals, q = low_modes(lap, drift)
            all_vals.append(vals)
            all_q.append(q)
            samples_data.append((vals, q))

        # The first acoustic multiplet consists of the first four positive
        # modes on the square torus.  Use its pooled median as a prescribed
        # center and its empirical logarithmic width as the CLT window.
        first_cluster = np.concatenate([vals[:4] for vals, _ in samples_data])
        r = float(np.median(first_cluster))
        log_sd = float(np.std(np.log(first_cluster)))
        bandwidth = max(log_sd / 4.0, 0.0025)
        vals_flat = np.concatenate(all_vals)
        q_flat = np.concatenate(all_q)
        u = np.log(vals_flat / r) / bandwidth
        kernel = np.exp(-0.5 * u * u) / (
            math.sqrt(2.0 * math.pi) * bandwidth
        )
        # Kernel density with respect to d(log r): this estimates r*pi_L(r).
        r_pi = float(np.sum(np.square(q_flat) * kernel) / samples)

        tails_at_r = np.array(
            [float(np.sum(q[vals <= r])) for vals, q in samples_data]
        )
        eta = math.exp(log_sd) - 1.0
        shells = np.array(
            [
                float(np.sum(q[(vals > r) & (vals <= (1.0 + eta) * r)]))
                for vals, q in samples_data
            ]
        )
        nominal = r / n  # d=2: N^{-1} r^{d/2}=N^{-1}r.
        rows.append(
            {
                "L": L,
                "N": n,
                "samples": samples,
                "cutoff_r": r,
                "relative_log_eigenvalue_width": log_sd,
                "moving_shell_eta": eta,
                "r_pi_kde": r_pi,
                "centered_tail_variance": float(np.var(tails_at_r)),
                "moving_shell_second_moment": float(np.mean(shells**2)),
                "nominal_variance_scale_Ninv_r": nominal,
                "r_pi_to_nominal_ratio": r_pi / nominal,
                "variance_to_nominal_ratio": float(np.var(tails_at_r)) / nominal,
            }
        )

    frame = pd.DataFrame(rows)
    exponents: dict[str, float] = {}
    for col in (
        "relative_log_eigenvalue_width",
        "r_pi_kde",
        "centered_tail_variance",
        "moving_shell_second_moment",
    ):
        slope = float(np.polyfit(np.log(frame["N"]), np.log(frame[col]), 1)[0])
        exponents[f"loglog_exponent_{col}_vs_N"] = slope
    return frame, exponents


def make_figure(frame: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.7, 4.5))
    ax.loglog(frame["N"], frame["r_pi_kde"], "o-", label=r"$r\,\pi_L(r)$")
    ax.loglog(
        frame["N"],
        frame["centered_tail_variance"],
        "s-",
        label=r"$\mathrm{Var}\,T_L(r)$",
    )
    n = frame["N"].to_numpy(dtype=float)
    c1 = float(frame["r_pi_kde"].iloc[0] * n[0] ** 1.5)
    c2 = float(frame["centered_tail_variance"].iloc[0] * n[0] ** 2.0)
    ax.loglog(n, c1 * n ** (-1.5), "--", label=r"reference $N^{-3/2}$")
    ax.loglog(n, c2 * n ** (-2.0), ":", label=r"reference $N^{-2}$")
    ax.set_xlabel(r"volume $N=L^2$")
    ax.set_ylabel("measured scale")
    ax.set_title("Acoustic cutoff: participation density versus centered variance")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--seed", type=int, default=20260727)
    args = parser.parse_args()
    out = Path(__file__).resolve().parent
    rng = np.random.default_rng(args.seed)

    toy = gaussian_one_mode_table()
    actual, exponents = random_conductance_scaling(rng, args.quick)
    toy.to_csv(out / "acoustic_gaussian_one_mode_obstruction.csv", index=False)
    actual.to_csv(out / "acoustic_low_mode_scaling.csv", index=False)
    make_figure(actual, out / "acoustic_cutoff_obstruction_scaling.png")

    checks = {
        "toy_holder_ratio_grows": bool(
            toy["shell_to_holder_ratio"].iloc[-1]
            > 4.0 * toy["shell_to_holder_ratio"].iloc[0]
        ),
        "toy_participation_density_has_sqrtN_excess": bool(
            toy["r_pi_to_N_minus_2_ratio"].iloc[-1]
            > 20.0 * toy["r_pi_to_N_minus_2_ratio"].iloc[0]
        ),
        "actual_relative_eigenvalue_width_near_N_minus_half": bool(
            -0.75
            < exponents["loglog_exponent_relative_log_eigenvalue_width_vs_N"]
            < -0.25
        ),
        "actual_participation_density_slower_than_N_minus_1p8": bool(
            exponents["loglog_exponent_r_pi_kde_vs_N"] > -1.8
        ),
        "actual_centered_variance_near_N_minus_2": bool(
            -2.45
            < exponents["loglog_exponent_centered_tail_variance_vs_N"]
            < -1.55
        ),
    }
    overall = all(checks.values())
    cert = {
        "seed": args.seed,
        "quick": bool(args.quick),
        "exact_gaussian_parameters": {"a": 1.0, "q_star": 0.7, "alpha": 0.5},
        "actual_model_loglog_exponents": exponents,
        "checks": checks,
        "overall_pass": overall,
        "interpretation": (
            "The Gaussian one-mode statements are exact. The random-conductance "
            "fits are diagnostics supporting, but not proving, the acoustic-CLT "
            "scaling mechanism."
        ),
    }
    (out / "acoustic_cutoff_obstruction_certificate.json").write_text(
        json.dumps(cert, indent=2, sort_keys=True) + "\n"
    )
    lines = [
        "ACOUSTIC CUTOFF / MAXIMAL-SHELL OBSTRUCTION AUDIT",
        "=" * 58,
        f"seed: {args.seed}",
        f"quick: {args.quick}",
        "",
        "Exact Gaussian one-mode model:",
        toy.to_string(index=False),
        "",
        "Random-conductance low-mode diagnostic:",
        actual.to_string(index=False),
        "",
        "Fitted log-log exponents:",
    ]
    lines.extend(f"  {k}: {v}" for k, v in exponents.items())
    lines.extend(["", "Checks:"])
    lines.extend(f"  {k}: {'PASS' if v else 'FAIL'}" for k, v in checks.items())
    lines.extend(["", f"OVERALL: {'PASS' if overall else 'FAIL'}"])
    text = "\n".join(lines) + "\n"
    (out / "ACOUSTIC_CUTOFF_OBSTRUCTION_VERIFICATION_OUTPUT.txt").write_text(text)
    print(text)
    if not overall:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
