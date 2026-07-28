#!/usr/bin/env python3
"""Verify the rank-one single-crossing and participation-coarea formulas.

The script checks exact finite-dimensional identities underlying the refined
pointwise martingale reduction for the nonlinear sharp source tail:

  * a one-bond rank-one path crosses a prescribed energy at most once;
  * Krein formulas for the crossing location, eigenvector, speed, and jump;
  * the scalar affine-source resolvent identity;
  * continuity after adding back the unique crossing jump;
  * the conditional one-bond variance decomposition;
  * Monte Carlo agreement between a direct log-energy participation density
    and the leave-one-edge coarea formula.

The Monte Carlo comparison is a diagnostic, not an input to any theorem.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.linalg as sla

from verify_sharp_cutoff_three_routes import (
    laplacian_and_drift,
    spectral_data,
    torus_edges,
)


def incidence_columns(
    n_vertices: int, tails: np.ndarray, heads: np.ndarray
) -> np.ndarray:
    b = np.zeros((n_vertices, tails.size), dtype=float)
    b[tails, np.arange(tails.size)] = -1.0
    b[heads, np.arange(heads.size)] = 1.0
    return b


def sharp_from_matrix(lap: np.ndarray, drift: np.ndarray, r: float) -> float:
    vals, coeff2 = spectral_data(lap, drift)
    mask = vals <= r
    return float(np.sum(coeff2[mask] / vals[mask]) / lap.shape[0])


def count_below(lap: np.ndarray, r: float) -> int:
    vals = sla.eigvalsh(lap, check_finite=False)
    return int(np.count_nonzero((vals > 1e-11) & (vals <= r)))


def modal_data(lap: np.ndarray, drift: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    vals, vecs = sla.eigh(lap, driver="evr", check_finite=False)
    mask = vals > 1e-11
    vals = vals[mask]
    vecs = vecs[:, mask]
    overlaps = vecs.T @ drift
    q = np.square(overlaps) / (lap.shape[0] * vals)
    return vals, q


def leave_one_crossing(
    lap_without: np.ndarray,
    drift_without: np.ndarray,
    b: np.ndarray,
    a_e: float,
    r: float,
    j_minus: float,
    j_plus: float,
) -> dict[str, float] | None:
    """Return the unique rank-one crossing coordinates, when present."""
    n = lap_without.shape[0]
    # Work on 1^perp to avoid the constant zero mode.
    qmat = sla.null_space(np.ones((1, n), dtype=float))
    a0 = qmat.T @ lap_without @ qmat
    bb = qmat.T @ b
    gg = qmat.T @ drift_without
    resolvent = sla.solve(a0 - r * np.eye(n - 1), np.eye(n - 1), assume_a="sym")
    rb = resolvent @ bb
    m = float(bb @ rb)
    m_prime = float(rb @ rb)
    h = float(bb @ (resolvent @ gg))
    if abs(m) < 1e-14 or m_prime <= 0.0:
        return None
    t_cross = -1.0 / m
    if not (j_minus < t_cross < j_plus):
        return None
    psi_red = rb / math.sqrt(m_prime)
    psi = qmat @ psi_red
    speed = m * m / m_prime
    overlap = (h - a_e) / math.sqrt(m_prime)
    q_jump = overlap * overlap / (n * r)
    return {
        "t": float(t_cross),
        "m": m,
        "m_prime": m_prime,
        "h": h,
        "speed": float(speed),
        "q": float(q_jump),
        "psi": psi,
    }


def affine_resolvent_formula(
    lap_without: np.ndarray,
    drift_without: np.ndarray,
    b: np.ndarray,
    a_e: float,
    t: float,
    z: complex,
) -> tuple[complex, complex]:
    n = lap_without.shape[0]
    qmat = sla.null_space(np.ones((1, n), dtype=float))
    a0 = qmat.T @ lap_without @ qmat
    bb = qmat.T @ b
    gg = qmat.T @ drift_without
    r0 = np.linalg.inv(a0.astype(complex) - z * np.eye(n - 1))
    m = bb @ r0 @ bb
    h = bb @ r0 @ gg
    h0 = gg @ r0 @ gg
    direct_lap = a0 + t * np.outer(bb, bb)
    direct_g = gg + t * a_e * bb
    direct = direct_g @ np.linalg.solve(
        direct_lap.astype(complex) - z * np.eye(n - 1), direct_g
    )
    formula = h0 + t * (a_e * a_e * m * t + 2.0 * a_e * h - h * h) / (1.0 + t * m)
    return complex(direct), complex(formula)


def locate_robust_case(
    rng: np.random.Generator,
    L: int,
    j_minus: float,
    j_plus: float,
) -> dict[str, object]:
    n = L * L
    tails, heads, dirs = torus_edges(L)
    bmat = incidence_columns(n, tails, heads)
    avec = np.where(dirs == 0, 1.0, 0.35)
    for _ in range(400):
        j = rng.uniform(j_minus, j_plus, size=tails.size)
        full_lap, full_drift = laplacian_and_drift(
            n, tails, heads, dirs, j, np.array([1.0, 0.35])
        )
        vals, _ = modal_data(full_lap, full_drift)
        candidates = vals[(vals > 0.25) & (vals < 2.5)]
        if candidates.size == 0:
            continue
        for r in candidates[:: max(1, candidates.size // 8)]:
            # Move slightly off the realized eigenvalue.
            r = float(r * 1.015)
            for e in rng.permutation(tails.size):
                j0 = j.copy()
                j0[e] = 0.0
                lap0, drift0 = laplacian_and_drift(
                    n, tails, heads, dirs, j0, np.array([1.0, 0.35])
                )
                info = leave_one_crossing(
                    lap0, drift0, bmat[:, e], float(avec[e]), r, j_minus, j_plus
                )
                if info is None:
                    continue
                t_cross = float(info["t"])
                if min(t_cross - j_minus, j_plus - t_cross) > 0.08:
                    return {
                        "j": j,
                        "edge": int(e),
                        "r": r,
                        "lap0": lap0,
                        "drift0": drift0,
                        "b": bmat[:, e],
                        "a_e": float(avec[e]),
                        "info": info,
                        "tails": tails,
                        "heads": heads,
                        "dirs": dirs,
                    }
    raise RuntimeError("failed to locate a robust one-edge crossing")


def exact_identity_audit(rng: np.random.Generator) -> tuple[dict[str, float], pd.DataFrame]:
    L = 4
    n = L * L
    j_minus, j_plus = 0.55, 1.45
    case = locate_robust_case(rng, L, j_minus, j_plus)
    lap0 = case["lap0"]
    drift0 = case["drift0"]
    b = case["b"]
    a_e = float(case["a_e"])
    r = float(case["r"])
    info = case["info"]
    t_star = float(info["t"])
    q_formula = float(info["q"])
    psi_formula = np.asarray(info["psi"])

    def matrices(t: float) -> tuple[np.ndarray, np.ndarray]:
        return lap0 + t * np.outer(b, b), drift0 + t * a_e * b

    grid = np.linspace(j_minus, j_plus, 401)
    counts = np.array([count_below(matrices(t)[0], r) for t in grid])
    count_drop = int(counts[0] - counts[-1])
    total_variation = int(np.sum(np.abs(np.diff(counts))))

    lap_star, drift_star = matrices(t_star)
    vals, vecs = sla.eigh(lap_star, driver="evr", check_finite=False)
    idx = int(np.argmin(np.abs(vals - r)))
    lambda_error = abs(float(vals[idx]) - r)
    psi_direct = vecs[:, idx]
    overlap_alignment = abs(float(psi_direct @ psi_formula))
    q_direct = float((psi_direct @ drift_star) ** 2 / (n * r))
    q_error = abs(q_direct - q_formula)

    eps = 2.0e-6
    lam_minus = sla.eigvalsh(matrices(t_star - eps)[0], check_finite=False)
    lam_plus = sla.eigvalsh(matrices(t_star + eps)[0], check_finite=False)
    lm = float(lam_minus[np.argmin(np.abs(lam_minus - r))])
    lp = float(lam_plus[np.argmin(np.abs(lam_plus - r))])
    fd_speed = (lp - lm) / (2.0 * eps)
    speed_error = abs(fd_speed - float(info["speed"]))

    f_minus = sharp_from_matrix(*matrices(t_star - eps), r)
    f_plus = sharp_from_matrix(*matrices(t_star + eps), r)
    jump_fd = f_minus - f_plus
    jump_error = abs(jump_fd - q_formula)

    z = complex(0.73 * r, 0.19 * r)
    direct_res, formula_res = affine_resolvent_formula(
        lap0, drift0, b, a_e, 0.83 * t_star, z
    )
    resolvent_error = abs(direct_res - formula_res)

    # The unique jump correction makes the one-edge path continuous.
    eps_cont = 5.0e-7
    left = sharp_from_matrix(*matrices(t_star - eps_cont), r)
    right = sharp_from_matrix(*matrices(t_star + eps_cont), r) + q_formula
    regularized_continuity_error = abs(left - right)

    # Conditional one-dimensional variance inequality.
    fine = np.linspace(j_minus, j_plus, 5001)
    fvals = np.array([sharp_from_matrix(*matrices(t), r) for t in fine])
    indicator = (fine >= t_star).astype(float)
    ftilde = fvals + q_formula * indicator
    var_f = float(np.var(fvals))
    var_ftilde = float(np.var(ftilde))
    p = (j_plus - t_star) / (j_plus - j_minus)
    variance_rhs = 2.0 * var_ftilde + 2.0 * q_formula**2 * p * (1.0 - p)
    variance_slack = variance_rhs - var_f

    rows = pd.DataFrame(
        {
            "t": fine,
            "sharp_tail": fvals,
            "regularized_tail": ftilde,
            "crossing_indicator": indicator,
        }
    )
    metrics = {
        "L": L,
        "N": n,
        "edge": int(case["edge"]),
        "cutoff_r": r,
        "crossing_t": t_star,
        "count_drop": count_drop,
        "count_total_variation": total_variation,
        "crossing_eigenvalue_absolute_error": lambda_error,
        "crossing_eigenvector_alignment": overlap_alignment,
        "jump_weight_formula": q_formula,
        "jump_weight_direct": q_direct,
        "jump_weight_absolute_error": q_error,
        "jump_finite_difference_absolute_error": jump_error,
        "speed_finite_difference_absolute_error": speed_error,
        "affine_resolvent_absolute_error": resolvent_error,
        "regularized_continuity_absolute_error": regularized_continuity_error,
        "conditional_variance": var_f,
        "conditional_variance_bound": variance_rhs,
        "conditional_variance_bound_slack": variance_slack,
    }
    return metrics, rows


def participation_coarea_audit(
    rng: np.random.Generator, samples: int
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Compare direct and leave-one-edge formulas for r*pi_L(r).

    The direct estimate uses a symmetric logarithmic shell of half-width eps:
      sum q_k^2 1_{|log(lambda_k/r)|<=eps}/(2 eps).
    This converges to r*pi_L(r).  The leave-one-edge formula is
      sum_e t_e rho(t_e) q_e(r)^2.
    """
    L = 5
    n = L * L
    j_minus, j_plus = 0.55, 1.45
    rho = 1.0 / (j_plus - j_minus)
    tails, heads, dirs = torus_edges(L)
    bmat = incidence_columns(n, tails, heads)
    avec = np.where(dirs == 0, 1.0, 0.35)
    Aprobe = np.array([1.0, 0.35])
    cutoffs = np.array([1.30, 2.70, 4.00, 5.40])
    eps_log = 0.08

    direct_sum = np.zeros_like(cutoffs)
    direct_sq = np.zeros_like(cutoffs)
    coarea_sum = np.zeros_like(cutoffs)
    coarea_sq = np.zeros_like(cutoffs)
    crossing_counts = np.zeros_like(cutoffs)

    for _ in range(samples):
        j = rng.uniform(j_minus, j_plus, size=tails.size)
        lap, drift = laplacian_and_drift(n, tails, heads, dirs, j, Aprobe)
        vals, q = modal_data(lap, drift)
        logvals = np.log(vals)
        for i, r in enumerate(cutoffs):
            mask = np.abs(logvals - math.log(r)) <= eps_log
            estimate = float(np.sum(q[mask] ** 2) / (2.0 * eps_log))
            direct_sum[i] += estimate
            direct_sq[i] += estimate * estimate

        # A fresh leave-one-out environment for the coarea expectation.  Since
        # the omitted edge value is not used, the remaining coordinates have
        # exactly their product marginal law.
        for i, r in enumerate(cutoffs):
            total = 0.0
            ncross = 0
            for e in range(tails.size):
                lap0 = lap - j[e] * np.outer(bmat[:, e], bmat[:, e])
                drift0 = drift - j[e] * float(avec[e]) * bmat[:, e]
                info = leave_one_crossing(
                    lap0,
                    drift0,
                    bmat[:, e],
                    float(avec[e]),
                    float(r),
                    j_minus,
                    j_plus,
                )
                if info is not None:
                    total += float(info["t"]) * rho * float(info["q"]) ** 2
                    ncross += 1
            coarea_sum[i] += total
            coarea_sq[i] += total * total
            crossing_counts[i] += ncross

    rows: list[dict[str, float]] = []
    rel_errors: list[float] = []
    for i, r in enumerate(cutoffs):
        direct = direct_sum[i] / samples
        coarea = coarea_sum[i] / samples
        direct_se = math.sqrt(max(direct_sq[i] / samples - direct**2, 0.0) / samples)
        coarea_se = math.sqrt(max(coarea_sq[i] / samples - coarea**2, 0.0) / samples)
        denom = max(abs(direct), abs(coarea), 1.0e-14)
        rel = abs(direct - coarea) / denom
        rel_errors.append(rel)
        rows.append(
            {
                "L": L,
                "N": n,
                "samples": samples,
                "cutoff_r": float(r),
                "log_half_width": eps_log,
                "direct_log_shell_r_pi": direct,
                "direct_standard_error": direct_se,
                "leave_one_out_r_pi": coarea,
                "leave_one_out_standard_error": coarea_se,
                "relative_discrepancy": rel,
                "mean_crossing_channels": crossing_counts[i] / samples,
            }
        )
    summary = {
        "samples": samples,
        "maximum_relative_density_discrepancy": max(rel_errors),
        "mean_relative_density_discrepancy": float(np.mean(rel_errors)),
    }
    return pd.DataFrame(rows), summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--seed", type=int, default=20260727)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    out = Path(__file__).resolve().parent

    exact, path_rows = exact_identity_audit(rng)
    samples = 160 if args.quick else 900
    density_rows, density_summary = participation_coarea_audit(rng, samples)

    path_rows.to_csv(out / "maximal_shell_single_crossing_path.csv", index=False)
    density_rows.to_csv(out / "maximal_shell_coarea_density.csv", index=False)

    tolerances = {
        "count_total_variation_max": 1,
        "crossing_eigenvalue_absolute_error": 2.0e-9,
        "crossing_eigenvector_alignment_min": 1.0 - 2.0e-8,
        "jump_weight_absolute_error": 2.0e-9,
        "jump_finite_difference_absolute_error": 2.5e-5,
        "speed_finite_difference_absolute_error": 2.5e-6,
        "affine_resolvent_absolute_error": 2.0e-10,
        "regularized_continuity_absolute_error": 2.5e-5,
        "conditional_variance_bound_slack_min": -2.0e-10,
    }
    checks = {
        "single_crossing": exact["count_total_variation"] <= 1,
        "count_drop_nonnegative": exact["count_drop"] in (0, 1),
        "crossing_eigenvalue": exact["crossing_eigenvalue_absolute_error"]
        <= tolerances["crossing_eigenvalue_absolute_error"],
        "crossing_eigenvector": exact["crossing_eigenvector_alignment"]
        >= tolerances["crossing_eigenvector_alignment_min"],
        "jump_formula": exact["jump_weight_absolute_error"]
        <= tolerances["jump_weight_absolute_error"],
        "jump_finite_difference": exact["jump_finite_difference_absolute_error"]
        <= tolerances["jump_finite_difference_absolute_error"],
        "crossing_speed": exact["speed_finite_difference_absolute_error"]
        <= tolerances["speed_finite_difference_absolute_error"],
        "affine_resolvent": exact["affine_resolvent_absolute_error"]
        <= tolerances["affine_resolvent_absolute_error"],
        "regularized_continuity": exact["regularized_continuity_absolute_error"]
        <= tolerances["regularized_continuity_absolute_error"],
        "conditional_variance_bound": exact["conditional_variance_bound_slack"]
        >= tolerances["conditional_variance_bound_slack_min"],
    }
    overall = all(checks.values())
    certificate = {
        "seed": args.seed,
        "quick": bool(args.quick),
        "exact_identity_metrics": exact,
        "density_diagnostic": density_summary,
        "tolerances": tolerances,
        "checks": checks,
        "overall_pass": overall,
        "note": (
            "The direct/coarea density comparison is a finite-sample, "
            "finite-bandwidth diagnostic and is not a theorem pass criterion."
        ),
    }
    (out / "maximal_shell_coarea_certificate.json").write_text(
        json.dumps(certificate, indent=2, sort_keys=True) + "\n"
    )

    lines = [
        "MAXIMAL CURRENT-WEIGHTED SHELL / COAREA AUDIT",
        "=" * 57,
        f"seed: {args.seed}",
        f"quick: {args.quick}",
        "",
        "Exact rank-one identities:",
    ]
    for key, value in exact.items():
        lines.append(f"  {key}: {value}")
    lines.extend(["", "Checks:"])
    for key, value in checks.items():
        lines.append(f"  {key}: {'PASS' if value else 'FAIL'}")
    lines.extend(
        [
            "",
            "Participation-density diagnostic:",
            density_rows.to_string(index=False),
            "",
            f"OVERALL: {'PASS' if overall else 'FAIL'}",
        ]
    )
    text = "\n".join(lines) + "\n"
    (out / "MAXIMAL_SHELL_COAREA_VERIFICATION_OUTPUT.txt").write_text(text)
    print(text)
    if not overall:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
