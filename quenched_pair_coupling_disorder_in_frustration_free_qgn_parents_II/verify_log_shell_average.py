#!/usr/bin/env python3
"""Verify the exact logarithmic shell identity and heat-mixture no-go.

The script is a numerical audit only.  The corresponding statements are proved
in SHARP_CUTOFF_THREE_ROUTES.md and in the LaTeX paper.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

from verify_sharp_cutoff_three_routes import (
    laplacian_and_drift,
    spectral_data,
    sharp_tail,
    torus_edges,
)


def interval_sweep_integral(log_lam: np.ndarray, weights: np.ndarray, eta: float) -> float:
    """Compute int S(r,eta)^2 dr/r by an independent event sweep in log r."""
    width = math.log1p(eta)
    events: dict[float, float] = {}
    for u, w in zip(log_lam, weights, strict=True):
        events[u - width] = events.get(u - width, 0.0) + float(w)
        events[u] = events.get(u, 0.0) - float(w)
    current = 0.0
    integral = 0.0
    previous: float | None = None
    for x in sorted(events):
        if previous is not None:
            integral += current * current * (x - previous)
        current += events[x]
        previous = x
    return integral


def pair_kernel_integral(log_lam: np.ndarray, weights: np.ndarray, eta: float) -> float:
    width = math.log1p(eta)
    kernel = np.maximum(width - np.abs(log_lam[:, None] - log_lam[None, :]), 0.0)
    return float(weights @ kernel @ weights)


def localized_pair_integral(
    log_lam: np.ndarray, weights: np.ndarray, eta: float, R: float
) -> float:
    width = math.log1p(eta)
    cap = math.log(R)
    lows = log_lam - width
    highs = log_lam
    lo = np.maximum(lows[:, None], lows[None, :])
    hi = np.minimum(np.minimum(highs[:, None], highs[None, :]), cap)
    lengths = np.maximum(hi - lo, 0.0)
    return float(weights @ lengths @ weights)


def bad_log_measure(
    log_lam: np.ndarray,
    weights: np.ndarray,
    eta: float,
    lower_r: float,
    upper_r: float,
    threshold: float,
) -> float:
    """Exact logarithmic measure where the shell square exceeds a threshold."""
    width = math.log1p(eta)
    lo = math.log(lower_r)
    hi = math.log(upper_r)
    events: dict[float, float] = {}
    current = 0.0
    for u, w in zip(log_lam, weights, strict=True):
        left = float(u - width)
        right = float(u)
        if left <= lo < right:
            current += float(w)
        if lo < left < hi:
            events[left] = events.get(left, 0.0) + float(w)
        if lo < right < hi:
            events[right] = events.get(right, 0.0) - float(w)
    bad = 0.0
    previous = lo
    for x in sorted(events):
        if current * current > threshold:
            bad += x - previous
        current += events[x]
        previous = x
    if current * current > threshold:
        bad += hi - previous
    return bad


def heat_mixture_audit(rng: np.random.Generator, samples: int = 500) -> pd.DataFrame:
    rows: list[dict[str, float | int]] = []
    for sample in range(samples):
        count = int(rng.integers(2, 9))
        rates = np.exp(rng.uniform(math.log(0.15), math.log(8.0), size=count))
        probs = rng.dirichlet(np.ones(count))
        mean = float(np.sum(probs / rates))
        second = float(2.0 * np.sum(probs / rates**2))
        variance = second - mean * mean
        rows.append(
            {
                "sample": sample,
                "components": count,
                "mean_threshold": mean,
                "variance_threshold": variance,
                "variance_to_mean_squared": variance / (mean * mean),
            }
        )
    return pd.DataFrame(rows)


def gamma_threshold_table() -> pd.DataFrame:
    rows = []
    for shape in (1, 2, 4, 8, 16, 32, 64):
        rows.append(
            {
                "gamma_shape": shape,
                "coefficient_of_variation": 1.0 / math.sqrt(shape),
                "variance_to_mean_squared": 1.0 / shape,
                "compatible_with_positive_heat_mixture_bound": shape == 1,
            }
        )
    return pd.DataFrame(rows)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path("."))
    parser.add_argument("--seed", type=int, default=20260727)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    L = 8
    N = L * L
    tails, heads, dirs = torus_edges(L)
    j = rng.uniform(0.5, 1.5, size=tails.size)
    lap, drift = laplacian_and_drift(
        N, tails, heads, dirs, j, np.array([1.0, 0.0])
    )
    vals, coeff2 = spectral_data(lap, drift)
    weights = coeff2 / (N * vals)
    log_lam = np.log(vals)
    eta = 0.30
    R = 0.80

    pair_value = pair_kernel_integral(log_lam, weights, eta)
    sweep_value = interval_sweep_integral(log_lam, weights, eta)
    localized_value = localized_pair_integral(log_lam, weights, eta, R)
    localized_bound = math.log1p(eta) * sharp_tail(
        vals, coeff2, N, (1.0 + eta) * R
    ) ** 2

    # Exact pointwise maximalization check.  In log energy, the shell
    # (x,x+h] is contained in every doubled shell (y,y+2h] for
    # y in [x-h,x].
    r_point = float(vals[len(vals) // 3] / (1.0 + 0.5 * eta))
    h_point = math.log1p(eta)
    eta_double = math.expm1(2.0 * h_point)
    shell_point = sharp_tail(vals, coeff2, N, (1.0 + eta) * r_point) - sharp_tail(
        vals, coeff2, N, r_point
    )
    max_integral = localized_pair_integral(
        log_lam, weights, eta_double, r_point
    ) - localized_pair_integral(
        log_lam, weights, eta_double, r_point / (1.0 + eta)
    )
    maximalization_rhs = max_integral / h_point

    shell_table = pd.DataFrame(
        [
            {
                "L": L,
                "eta": eta,
                "R": R,
                "full_pair_kernel_integral": pair_value,
                "full_interval_sweep_integral": sweep_value,
                "full_identity_absolute_error": abs(pair_value - sweep_value),
                "localized_integral": localized_value,
                "localized_deterministic_bound": localized_bound,
                "localized_bound_margin": localized_bound - localized_value,
                "pointwise_cutoff": r_point,
                "pointwise_shell_square": shell_point * shell_point,
                "pointwise_maximalization_rhs": maximalization_rhs,
                "pointwise_maximalization_margin": maximalization_rhs - shell_point * shell_point,
            }
        ]
    )
    shell_table.to_csv(args.outdir / "logarithmic_shell_identity_audit.csv", index=False)

    beta = 0.5
    dyadic_localized = localized_pair_integral(log_lam, weights, eta, R) - localized_pair_integral(
        log_lam, weights, eta, R / 2.0
    )
    cap_tail = sharp_tail(vals, coeff2, N, (1.0 + eta) * R)
    bad_threshold = eta**beta * cap_tail**2
    bad_measure = bad_log_measure(
        log_lam, weights, eta, R / 2.0, R, bad_threshold
    )
    markov_bound = dyadic_localized / bad_threshold if bad_threshold > 0 else 0.0
    analytic_bound = math.log1p(eta) * cap_tail**2 / bad_threshold if bad_threshold > 0 else 0.0
    exceptional = pd.DataFrame(
        [
            {
                "L": L,
                "eta": eta,
                "beta": beta,
                "window_lower": R / 2.0,
                "window_upper": R,
                "bad_threshold": bad_threshold,
                "exact_bad_log_measure": bad_measure,
                "exact_markov_bound": markov_bound,
                "localized_analytic_bound": analytic_bound,
                "markov_margin": markov_bound - bad_measure,
            }
        ]
    )
    exceptional.to_csv(args.outdir / "exceptional_cutoff_set_audit.csv", index=False)

    mixtures = heat_mixture_audit(rng)
    mixtures.to_csv(args.outdir / "positive_heat_mixture_no_go_audit.csv", index=False)
    gamma = gamma_threshold_table()
    gamma.to_csv(args.outdir / "gamma_threshold_concentration_audit.csv", index=False)

    checks = {
        "log_shell_identity_exact": bool(abs(pair_value - sweep_value) < 1e-12),
        "localized_log_shell_bound": bool(localized_value <= localized_bound + 1e-14),
        "pointwise_log_shell_maximalization": bool(
            shell_point * shell_point <= maximalization_rhs + 1e-14
        ),
        "exceptional_set_markov_bound": bool(bad_measure <= markov_bound + 1e-14),
        "positive_heat_mixture_cv_bound": bool(float(
            mixtures["variance_to_mean_squared"].min()
        )
        >= 1.0 - 1e-12),
        "gamma_narrowness_requires_non_heat_mixture": bool(
            np.all(
                gamma.loc[
                    gamma["gamma_shape"] > 1,
                    "variance_to_mean_squared",
                ]
                < 1.0
            )
        ),
    }
    certificate = {
        "seed": args.seed,
        "checks": checks,
        "all_checks_pass": bool(all(checks.values())),
        "max_log_shell_identity_error": abs(pair_value - sweep_value),
        "localized_bound_margin": localized_bound - localized_value,
        "pointwise_maximalization_margin": maximalization_rhs - shell_point * shell_point,
        "exceptional_set_bad_log_measure": bad_measure,
        "exceptional_set_markov_bound": markov_bound,
        "minimum_positive_heat_mixture_variance_ratio": float(
            mixtures["variance_to_mean_squared"].min()
        ),
        "median_positive_heat_mixture_variance_ratio": float(
            mixtures["variance_to_mean_squared"].median()
        ),
    }
    with (args.outdir / "log_shell_average_certificate.json").open("w") as fh:
        json.dump(certificate, fh, indent=2, sort_keys=True)

    lines = [
        "LOGARITHMIC SHELL / RANDOMIZATION NO-GO AUDIT",
        "================================================",
        f"seed: {args.seed}",
        "",
        "Exact logarithmic shell identity:",
        shell_table.to_string(index=False),
        "",
        "Exceptional-cutoff set audit:",
        exceptional.to_string(index=False),
        "",
        "Positive heat-mixture variance ratios:",
        mixtures["variance_to_mean_squared"].describe().to_string(),
        "",
        "Gamma threshold comparison:",
        gamma.to_string(index=False),
        "",
    ]
    for name, passed in checks.items():
        lines.append(f"{name}: {'PASS' if passed else 'FAIL'}")
    lines.append("")
    lines.append(f"OVERALL: {'PASS' if certificate['all_checks_pass'] else 'FAIL'}")
    text = "\n".join(lines) + "\n"
    (args.outdir / "LOG_SHELL_AVERAGE_VERIFICATION_OUTPUT.txt").write_text(text)
    print(text)
    return 0 if certificate["all_checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
