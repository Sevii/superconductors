#!/usr/bin/env python3
"""Audit three proposed routes to centered sharp-cutoff fluctuations.

The script checks four things:
  1. the finite discrete relative Wegner count estimate on random 2D tori;
  2. the modulated-cycle obstruction to replacing drift weight by normalized
     eigenvalue count;
  3. exact exponential- and uniform-randomized cutoff identities;
  4. the exact Doob martingale variance decomposition for an integrated
     spectral statistic on a small Bernoulli ring.

The computations are diagnostics only. The mathematical statements are proved
in SHARP_CUTOFF_THREE_ROUTES.md and in the accompanying LaTeX paper.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.linalg as sla


def torus_edges(L: int, d: int = 2) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if d != 2:
        raise ValueError("torus_edges currently implements d=2")
    tails: list[int] = []
    heads: list[int] = []
    dirs: list[int] = []
    idx = lambda x, y: (x % L) * L + (y % L)
    for x in range(L):
        for y in range(L):
            i = idx(x, y)
            tails.extend((i, i))
            heads.extend((idx(x + 1, y), idx(x, y + 1)))
            dirs.extend((0, 1))
    return np.array(tails), np.array(heads), np.array(dirs)


def cycle_edges(N: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tails = np.arange(N, dtype=int)
    heads = (tails + 1) % N
    dirs = np.zeros(N, dtype=int)
    return tails, heads, dirs


def laplacian_and_drift(
    n_vertices: int,
    tails: np.ndarray,
    heads: np.ndarray,
    directions: np.ndarray,
    conductances: np.ndarray,
    A: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    lap = np.zeros((n_vertices, n_vertices), dtype=float)
    drift = np.zeros(n_vertices, dtype=float)
    for e, (x, y, mu) in enumerate(zip(tails, heads, directions, strict=True)):
        j = float(conductances[e])
        lap[x, x] += j
        lap[y, y] += j
        lap[x, y] -= j
        lap[y, x] -= j
        amp = j * float(A[mu])
        # B has -1 at the tail and +1 at the head.
        drift[x] -= amp
        drift[y] += amp
    drift -= drift.mean()
    return lap, drift


def spectral_data(lap: np.ndarray, drift: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    vals, vecs = sla.eigh(lap, driver="evr", check_finite=False)
    mask = vals > 1e-11
    vals = vals[mask]
    coeff2 = np.square(vecs[:, mask].T @ drift)
    return vals, coeff2


def sharp_tail(vals: np.ndarray, coeff2: np.ndarray, N: int, delta: float) -> float:
    mask = vals <= delta
    return float(np.sum(coeff2[mask] / vals[mask]) / N)


def shell_tail(
    vals: np.ndarray, coeff2: np.ndarray, N: int, delta: float, eta: float
) -> float:
    mask = (vals > delta) & (vals <= (1.0 + eta) * delta)
    return float(np.sum(coeff2[mask] / vals[mask]) / N)


def abel_tail(vals: np.ndarray, coeff2: np.ndarray, N: int, delta: float) -> float:
    return float(np.sum(np.exp(-vals / delta) * coeff2 / vals) / N)


def exponential_randomized_second_moment(
    vals: np.ndarray, coeff2: np.ndarray, N: int, delta: float
) -> float:
    weights = coeff2 / (N * vals)
    mx = np.maximum.outer(vals, vals)
    return float(np.sum(np.outer(weights, weights) * np.exp(-mx / delta)))


def uniform_randomized_mean(
    vals: np.ndarray, coeff2: np.ndarray, N: int, delta: float, eta: float
) -> float:
    top = (1.0 + eta) * delta
    q = np.where(
        vals <= delta,
        1.0,
        np.where(vals <= top, (top - vals) / (eta * delta), 0.0),
    )
    return float(np.sum(q * coeff2 / vals) / N)


def integrated_K(vals: np.ndarray, coeff2: np.ndarray, N: int, r: float) -> float:
    # K(r)=int_0^r T(s) ds = sum (r-lambda)_+ coeff^2/(N lambda).
    return float(np.sum(np.maximum(r - vals, 0.0) * coeff2 / vals) / N)


def relative_wegner_audit(rng: np.random.Generator) -> pd.DataFrame:
    L = 6
    N = L * L
    d = 2
    j_minus, j_plus = 0.5, 1.5
    rho_inf = 1.0 / (j_plus - j_minus)
    E = 0.85
    etas = (0.025, 0.05, 0.10, 0.20, 0.40)
    samples = 700
    tails, heads, dirs = torus_edges(L)
    A = np.array([1.0, 0.0])
    accum = {eta: {"count": [], "shell": []} for eta in etas}
    for _ in range(samples):
        j = rng.uniform(j_minus, j_plus, size=tails.size)
        lap, drift = laplacian_and_drift(N, tails, heads, dirs, j, A)
        vals, coeff2 = spectral_data(lap, drift)
        for eta in etas:
            upper = (1.0 + eta) * E
            count = int(np.count_nonzero((vals > E) & (vals <= upper)))
            accum[eta]["count"].append(count)
            accum[eta]["shell"].append(shell_tail(vals, coeff2, N, E, eta))
    rows: list[dict[str, float]] = []
    for eta in etas:
        counts = np.asarray(accum[eta]["count"], dtype=float)
        shells = np.asarray(accum[eta]["shell"], dtype=float)
        # The theorem is E count <= d N rho_inf j_plus |I|/E = d N rho_inf j_plus eta.
        theorem_bound = d * N * rho_inf * j_plus * eta
        rows.append(
            {
                "L": L,
                "N": N,
                "samples": samples,
                "E": E,
                "eta": eta,
                "mean_eigenvalue_count": float(counts.mean()),
                "wegner_theorem_bound": theorem_bound,
                "count_to_bound_ratio": float(counts.mean() / theorem_bound),
                "mean_drift_shell": float(shells.mean()),
                "drift_shell_second_moment": float(np.mean(shells**2)),
                "mean_drift_shell_over_eta": float(shells.mean() / eta),
            }
        )
    return pd.DataFrame(rows)


def modulated_cycle_audit() -> pd.DataFrame:
    eps = 0.5
    low_window_constant = 200.0
    rows: list[dict[str, float]] = []
    for N in (32, 64, 128, 256):
        tails, heads, dirs = cycle_edges(N)
        e = np.arange(N)
        j = 1.0 + eps * np.cos(2.0 * np.pi * e / N)
        lap, drift = laplacian_and_drift(N, tails, heads, dirs, j, np.array([1.0]))
        vals, coeff2 = spectral_data(lap, drift)
        cutoff = low_window_constant / (N * N)
        count = int(np.count_nonzero(vals <= cutoff))
        low_mass = sharp_tail(vals, coeff2, N, cutoff)
        total_mass = float(np.sum(coeff2 / vals) / N)
        first_moment = float(np.sum(coeff2) / N)
        harmonic = 1.0 / float(np.mean(1.0 / j))
        exact_network_mass = float(np.mean(j) - harmonic)
        rows.append(
            {
                "N": N,
                "epsilon": eps,
                "cutoff": cutoff,
                "positive_mode_count": count,
                "normalized_mode_count": count / N,
                "low_Hminus1_mass": low_mass,
                "total_Hminus1_mass_spectral": total_mass,
                "total_Hminus1_mass_network": exact_network_mass,
                "first_moment": first_moment,
                "low_mass_to_normalized_count": low_mass / max(count / N, 1e-300),
            }
        )
    return pd.DataFrame(rows)


def randomization_audit(rng: np.random.Generator) -> tuple[pd.DataFrame, dict[str, float]]:
    L = 8
    N = L * L
    tails, heads, dirs = torus_edges(L)
    A = np.array([1.0, 0.0])
    j = rng.uniform(0.5, 1.5, size=tails.size)
    lap, drift = laplacian_and_drift(N, tails, heads, dirs, j, A)
    vals, coeff2 = spectral_data(lap, drift)
    delta = 0.50
    eta = 0.25

    sharp = sharp_tail(vals, coeff2, N, delta)
    shell = shell_tail(vals, coeff2, N, delta, eta)
    abel = abel_tail(vals, coeff2, N, delta)
    exp_second = exponential_randomized_second_moment(vals, coeff2, N, delta)
    abel_2delta_sq = abel_tail(vals, coeff2, N, 2.0 * delta) ** 2
    uniform_mean = uniform_randomized_mean(vals, coeff2, N, delta, eta)
    K_diff = (
        integrated_K(vals, coeff2, N, (1.0 + eta) * delta)
        - integrated_K(vals, coeff2, N, delta)
    ) / (eta * delta)

    mc_samples = 250_000
    theta_exp = rng.exponential(scale=delta, size=mc_samples)
    # Evaluate T(theta) efficiently using sorted eigenvalues and cumulative weights.
    order = np.argsort(vals)
    lam = vals[order]
    modal = coeff2[order] / (N * lam)
    cumulative = np.cumsum(modal)
    idx = np.searchsorted(lam, theta_exp, side="right") - 1
    exp_values = np.where(idx >= 0, cumulative[np.maximum(idx, 0)], 0.0)

    theta_unif = rng.uniform(delta, (1.0 + eta) * delta, size=mc_samples)
    idx_u = np.searchsorted(lam, theta_unif, side="right") - 1
    unif_values = np.where(idx_u >= 0, cumulative[np.maximum(idx_u, 0)], 0.0)

    rows = pd.DataFrame(
        [
            {
                "identity": "exponential mean equals Abel tail",
                "exact_left": float(exp_values.mean()),
                "exact_right": abel,
                "absolute_error": abs(float(exp_values.mean()) - abel),
            },
            {
                "identity": "exponential second moment spectral formula",
                "exact_left": float(np.mean(exp_values**2)),
                "exact_right": exp_second,
                "absolute_error": abs(float(np.mean(exp_values**2)) - exp_second),
            },
            {
                "identity": "uniform mean equals triangular spectral filter",
                "exact_left": float(unif_values.mean()),
                "exact_right": uniform_mean,
                "absolute_error": abs(float(unif_values.mean()) - uniform_mean),
            },
            {
                "identity": "uniform mean equals integrated difference quotient",
                "exact_left": uniform_mean,
                "exact_right": K_diff,
                "absolute_error": abs(uniform_mean - K_diff),
            },
        ]
    )
    summary = {
        "L": float(L),
        "delta": delta,
        "eta": eta,
        "sharp_tail": sharp,
        "shell_tail": shell,
        "uniform_randomized_mean": uniform_mean,
        "uniform_mean_minus_sharp": uniform_mean - sharp,
        "uniform_sandwich_residual": shell - (uniform_mean - sharp),
        "abel_tail": abel,
        "exponential_randomized_variance_mc": float(np.var(exp_values, ddof=0)),
        "exponential_second_moment_exact": exp_second,
        "abel_2delta_squared": abel_2delta_sq,
        "second_moment_bound_margin": abel_2delta_sq - exp_second,
        "uniform_conditional_variance_mc": float(np.var(unif_values, ddof=0)),
        "uniform_range_variance_bound": shell * shell / 4.0,
    }
    return rows, summary


def bernoulli_ring_observable(bits: tuple[int, ...], delta: float, h: float) -> tuple[float, float]:
    N = len(bits)
    tails, heads, dirs = cycle_edges(N)
    j = np.where(np.asarray(bits, dtype=int) == 0, 0.7, 1.3)
    lap, drift = laplacian_and_drift(N, tails, heads, dirs, j, np.array([1.0]))
    vals, coeff2 = spectral_data(lap, drift)
    T = sharp_tail(vals, coeff2, N, delta)
    Y = (integrated_K(vals, coeff2, N, delta + h) - integrated_K(vals, coeff2, N, delta)) / h
    return T, Y


def doob_variance(values: np.ndarray, nbits: int) -> tuple[float, list[float]]:
    """Exact Doob variance decomposition for equiprobable bit strings.

    values are ordered lexicographically as itertools.product((0,1), repeat=nbits).
    """
    total_mean = float(values.mean())
    prev_map: dict[tuple[int, ...], float] = {(): total_mean}
    increments: list[float] = []
    configs = list(itertools.product((0, 1), repeat=nbits))
    for m in range(1, nbits + 1):
        groups: dict[tuple[int, ...], list[float]] = {}
        for cfg, val in zip(configs, values, strict=True):
            groups.setdefault(cfg[:m], []).append(float(val))
        curr_map = {prefix: float(np.mean(v)) for prefix, v in groups.items()}
        sq = 0.0
        for cfg in configs:
            curr = curr_map[cfg[:m]]
            prev = prev_map[cfg[: m - 1]]
            sq += (curr - prev) ** 2
        increments.append(sq / len(configs))
        prev_map = curr_map
    return float(np.var(values, ddof=0)), increments


def martingale_audit() -> tuple[pd.DataFrame, pd.DataFrame]:
    N = 8
    delta = 0.615
    hs = (0.20, 0.10, 0.05, 0.025)
    configs = list(itertools.product((0, 1), repeat=N))
    sharp_vals = np.asarray([bernoulli_ring_observable(c, delta, hs[0])[0] for c in configs])
    sharp_var, sharp_increments = doob_variance(sharp_vals, N)
    rows: list[dict[str, float | str]] = []
    increment_rows: list[dict[str, float | int | str]] = []
    rows.append(
        {
            "observable": "sharp",
            "h": 0.0,
            "variance": sharp_var,
            "sum_doob_quadratic_variation": float(sum(sharp_increments)),
            "decomposition_error": abs(sharp_var - sum(sharp_increments)),
            "L2_distance_to_sharp": 0.0,
        }
    )
    for edge, value in enumerate(sharp_increments, start=1):
        increment_rows.append({"observable": "sharp", "h": 0.0, "edge": edge, "increment_variance": value})
    for h in hs:
        vals = np.asarray([bernoulli_ring_observable(c, delta, h)[1] for c in configs])
        var, increments = doob_variance(vals, N)
        l2dist = float(np.mean((vals - sharp_vals) ** 2))
        rows.append(
            {
                "observable": "integrated_difference_quotient",
                "h": h,
                "variance": var,
                "sum_doob_quadratic_variation": float(sum(increments)),
                "decomposition_error": abs(var - sum(increments)),
                "L2_distance_to_sharp": l2dist,
            }
        )
        for edge, value in enumerate(increments, start=1):
            increment_rows.append(
                {
                    "observable": "integrated_difference_quotient",
                    "h": h,
                    "edge": edge,
                    "increment_variance": value,
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(increment_rows)


def make_plots(
    outdir: Path,
    wegner: pd.DataFrame,
    modulation: pd.DataFrame,
    martingale: pd.DataFrame,
) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    ax.plot(wegner["eta"], wegner["mean_eigenvalue_count"], "o-", label="empirical mean count")
    ax.plot(wegner["eta"], wegner["wegner_theorem_bound"], "s--", label="relative Wegner bound")
    ax.set_xlabel(r"relative shell width $\eta$")
    ax.set_ylabel("eigenvalue count")
    ax.set_title("Finite discrete relative Wegner audit")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "relative_wegner_count_audit.png", dpi=180)
    plt.close(fig)

    fig, ax1 = plt.subplots(figsize=(6.4, 4.3))
    ax1.plot(modulation["N"], modulation["low_Hminus1_mass"], "o-", label=r"low $H^{-1}$ mass")
    ax1.plot(modulation["N"], modulation["normalized_mode_count"], "s--", label="normalized mode count")
    ax1.set_xscale("log", base=2)
    ax1.set_yscale("log")
    ax1.set_xlabel("cycle size N")
    ax1.set_ylabel("mass / normalized count")
    ax1.set_title("Drift weight is not controlled by normalized eigenvalue count")
    ax1.legend()
    fig.tight_layout()
    fig.savefig(outdir / "modulated_cycle_drift_weight_obstruction.png", dpi=180)
    plt.close(fig)

    subset = martingale[martingale["observable"] == "integrated_difference_quotient"].sort_values("h")
    fig, ax = plt.subplots(figsize=(6.4, 4.3))
    ax.plot(subset["h"], subset["variance"], "o-", label="variance")
    ax.plot(subset["h"], subset["L2_distance_to_sharp"], "s--", label=r"$L^2$ distance to sharp")
    ax.axhline(float(martingale.loc[martingale["observable"] == "sharp", "variance"].iloc[0]), linestyle=":", label="sharp variance")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("integrated window h")
    ax.set_ylabel("quadratic quantity")
    ax.set_title("Integrated-spectral martingale audit")
    ax.legend()
    fig.tight_layout()
    fig.savefig(outdir / "integrated_martingale_audit.png", dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path("."))
    parser.add_argument("--seed", type=int, default=20260727)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    wegner = relative_wegner_audit(rng)
    modulation = modulated_cycle_audit()
    randomization, random_summary = randomization_audit(rng)
    martingale, martingale_increments = martingale_audit()

    wegner.to_csv(args.outdir / "relative_wegner_audit.csv", index=False)
    modulation.to_csv(args.outdir / "modulated_cycle_weight_obstruction.csv", index=False)
    randomization.to_csv(args.outdir / "randomized_cutoff_identity_audit.csv", index=False)
    martingale.to_csv(args.outdir / "integrated_martingale_audit.csv", index=False)
    martingale_increments.to_csv(args.outdir / "integrated_martingale_increments.csv", index=False)
    make_plots(args.outdir, wegner, modulation, martingale)

    checks = {
        "relative_wegner_empirical_below_bound": bool(np.all(wegner["mean_eigenvalue_count"] <= wegner["wegner_theorem_bound"])),
        "modulated_low_mass_stays_positive": bool(modulation["low_Hminus1_mass"].iloc[-1] > 0.02),
        "modulated_normalized_count_decays": bool(modulation["normalized_mode_count"].iloc[-1] < modulation["normalized_mode_count"].iloc[0]),
        "modulated_network_spectral_identity": bool(np.max(np.abs(modulation["total_Hminus1_mass_spectral"] - modulation["total_Hminus1_mass_network"])) < 1e-10),
        "randomized_identity_mc_tolerance": bool(randomization["absolute_error"].max() < 5e-4),
        "exponential_second_moment_bound": bool(random_summary["second_moment_bound_margin"] >= -1e-12),
        "uniform_sandwich": bool(random_summary["uniform_mean_minus_sharp"] >= -1e-12 and random_summary["uniform_sandwich_residual"] >= -1e-12),
        "uniform_variance_range_bound": bool(random_summary["uniform_conditional_variance_mc"] <= random_summary["uniform_range_variance_bound"] + 5e-6),
        "doob_decomposition_exact": bool(martingale["decomposition_error"].max() < 1e-14),
    }
    certificate = {
        "seed": args.seed,
        "checks": checks,
        "all_checks_pass": bool(all(checks.values())),
        "randomization_summary": random_summary,
        "max_randomized_identity_error": float(randomization["absolute_error"].max()),
        "max_doob_decomposition_error": float(martingale["decomposition_error"].max()),
        "largest_modulated_weight_to_normalized_count_ratio": float(modulation["low_mass_to_normalized_count"].max()),
    }
    with (args.outdir / "sharp_cutoff_three_routes_certificate.json").open("w") as fh:
        json.dump(certificate, fh, indent=2, sort_keys=True)

    lines = [
        "SHARP CUTOFF THREE-ROUTE AUDIT",
        "================================",
        f"seed: {args.seed}",
        "",
        "Relative Wegner count:",
        wegner.to_string(index=False),
        "",
        "Modulated-cycle obstruction:",
        modulation.to_string(index=False),
        "",
        "Randomized cutoff identities:",
        randomization.to_string(index=False),
        "",
        "Randomization summary:",
        json.dumps(random_summary, indent=2, sort_keys=True),
        "",
        "Integrated martingale:",
        martingale.to_string(index=False),
        "",
    ]
    for name, passed in checks.items():
        lines.append(f"{name}: {'PASS' if passed else 'FAIL'}")
    lines.append("")
    lines.append(f"OVERALL: {'PASS' if certificate['all_checks_pass'] else 'FAIL'}")
    (args.outdir / "SHARP_CUTOFF_THREE_ROUTES_VERIFICATION_OUTPUT.txt").write_text("\n".join(lines) + "\n")

    print("\n".join(lines[-12:]))
    return 0 if certificate["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
