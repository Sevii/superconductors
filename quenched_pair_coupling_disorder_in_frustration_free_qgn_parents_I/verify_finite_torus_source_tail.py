#!/usr/bin/env python3
"""Numerical audit for the quantitative finite-torus QGN source-tail theorem.

The theorem proved in the accompanying draft is annealed and analytic.  This
script does not replace that proof.  It checks the exact finite-torus
intertwining identities and samples the predicted two-dimensional scalings

    E H_L(t)  = O((1+t)^-2),
    E T_L(delta) = O(delta),

for i.i.d. uniformly elliptic conductances, where

    H_L(t) = |Lambda_L|^-1 <g_L, exp(-t L_L) g_L>,
    T_L(delta) = |Lambda_L|^-1 <g_L, L_L^+ 1_(0,delta](L_L) g_L>.

The default run is intentionally modest so that the certificate can be
reproduced on a laptop.  Increase --samples and --sizes for a research scan.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


@dataclass(frozen=True)
class TorusSample:
    L: int
    d: int
    conductances: np.ndarray  # shape (N,d), positive-direction edges
    direction: np.ndarray  # harmonic bond direction A
    laplacian: sp.csr_matrix
    drift: np.ndarray

    @property
    def n_sites(self) -> int:
        return self.L**self.d


def index_of(coord: Sequence[int], L: int) -> int:
    out = 0
    for value in coord:
        out = out * L + int(value)
    return out


def shifted(coord: Sequence[int], mu: int, step: int, L: int) -> tuple[int, ...]:
    out = list(coord)
    out[mu] = (out[mu] + step) % L
    return tuple(out)


def make_torus_sample(
    L: int,
    d: int,
    rng: np.random.Generator,
    j_min: float,
    j_max: float,
    direction: np.ndarray,
) -> TorusSample:
    if L < 3:
        raise ValueError("L must be at least 3")
    if direction.shape != (d,):
        raise ValueError("direction must have shape (d,)")
    n = L**d
    j = rng.uniform(j_min, j_max, size=(n, d))

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    diagonal = np.zeros(n, dtype=float)
    drift = np.zeros(n, dtype=float)

    for coord_arr in np.ndindex(*(L for _ in range(d))):
        coord = tuple(int(v) for v in coord_arr)
        x = index_of(coord, L)
        for mu in range(d):
            y_coord = shifted(coord, mu, +1, L)
            y = index_of(y_coord, L)
            weight = float(j[x, mu])
            diagonal[x] += weight
            diagonal[y] += weight
            rows.extend((x, y))
            cols.extend((y, x))
            vals.extend((-weight, -weight))

            flux = weight * float(direction[mu])
            drift[x] -= flux
            drift[y] += flux

    rows.extend(range(n))
    cols.extend(range(n))
    vals.extend(diagonal.tolist())
    lap = sp.coo_matrix((vals, (rows, cols)), shape=(n, n), dtype=float).tocsr()
    lap.sum_duplicates()
    return TorusSample(
        L=L,
        d=d,
        conductances=j,
        direction=direction.copy(),
        laplacian=lap,
        drift=drift,
    )


def build_environment_orbit_generator(sample: TorusSample) -> tuple[sp.csr_matrix, np.ndarray]:
    """Build the environment generator and drift from translation differences.

    Orbit states are indexed by ``x`` and represent ``tau_x J``.  This routine
    implements ``D* j(0) D`` row by row, rather than reusing the physical
    incidence assembly in :func:`make_torus_sample`.  It also evaluates the
    local divergence datum ``D*(j(0) A)`` directly from translated coefficient
    blocks.  Agreement with the physical matrix and drift therefore audits the
    exact environment--site intertwining at the generator level.
    """
    L = sample.L
    d = sample.d
    n = sample.n_sites
    j = sample.conductances
    direction = sample.direction

    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    drift = np.zeros(n, dtype=float)

    for coord_arr in np.ndindex(*(L for _ in range(d))):
        coord = tuple(int(v) for v in coord_arr)
        x = index_of(coord, L)
        diagonal = 0.0
        for mu in range(d):
            plus_coord = shifted(coord, mu, +1, L)
            minus_coord = shifted(coord, mu, -1, L)
            xp = index_of(plus_coord, L)
            xm = index_of(minus_coord, L)
            outgoing = float(j[x, mu])
            incoming = float(j[xm, mu])

            # D_mu^* j_mu(0) D_mu evaluated at tau_x J.
            diagonal += outgoing + incoming
            rows.extend((x, x))
            cols.extend((xp, xm))
            vals.extend((-outgoing, -incoming))

            # D_mu^*(j_mu(0) A_mu) at tau_x J.
            drift[x] += (incoming - outgoing) * float(direction[mu])

        rows.append(x)
        cols.append(x)
        vals.append(diagonal)

    generator = sp.coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
    generator.sum_duplicates()
    return generator, drift


def low_spectrum(
    laplacian: sp.csr_matrix,
    max_cutoff: float,
    initial_k: int = 48,
    tolerance: float = 1e-9,
) -> tuple[np.ndarray, np.ndarray]:
    """Return all eigenpairs needed up to max_cutoff, including zero mode.

    A tiny diagonal shift makes shift-invert robust while preserving
    eigenvectors.  The shift is subtracted from the returned eigenvalues.
    """
    n = laplacian.shape[0]
    if n <= 2:
        dense = laplacian.toarray()
        return np.linalg.eigh(dense)
    k = min(max(8, initial_k), n - 1)
    eps = 1e-9
    shifted_matrix = laplacian + eps * sp.eye(n, format="csr")

    while True:
        vals, vecs = spla.eigsh(
            shifted_matrix,
            k=k,
            sigma=0.0,
            which="LM",
            tol=tolerance,
            maxiter=max(5000, 20 * n),
        )
        vals = vals - eps
        order = np.argsort(vals)
        vals = vals[order]
        vecs = vecs[:, order]
        if vals[-1] > 1.15 * max_cutoff or k >= n - 1:
            return vals, vecs
        k = min(n - 1, max(k + 16, 2 * k))


def heat_correlation(sample: TorusSample, time: float) -> float:
    evolved = spla.expm_multiply((-time) * sample.laplacian, sample.drift)
    value = float(np.dot(sample.drift, evolved) / sample.n_sites)
    # Numerical roundoff may create tiny negative values at very long times.
    return max(0.0, value)


def spectral_tail_from_pairs(
    sample: TorusSample,
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
    cutoff: float,
) -> tuple[float, float]:
    positive = eigenvalues > 1e-8
    low = positive & (eigenvalues <= cutoff * (1.0 + 1e-10))
    coeff = eigenvectors.T @ sample.drift
    source_mass = float(np.sum(np.abs(coeff[low]) ** 2) / sample.n_sites)
    hminus1_tail = float(
        np.sum(np.abs(coeff[low]) ** 2 / eigenvalues[low]) / sample.n_sites
    )
    return source_mass, hminus1_tail


def exact_total_hminus1(sample: TorusSample) -> float:
    """Compute N^-1 <g,L^+g> by a zero-mean constrained solve."""
    n = sample.n_sites
    # Add the constant-mode projector.  It does not affect a zero-mean right
    # hand side and makes the matrix strictly positive.
    ones = np.ones(n)
    projector = sp.csr_matrix(np.outer(ones, ones) / n)
    solution = spla.spsolve(sample.laplacian + projector, sample.drift)
    return float(np.dot(sample.drift, solution) / n)


def write_csv(path: Path, rows: list[dict[str, float | int]]) -> None:
    if not rows:
        raise ValueError(f"No rows for {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def fit_log_slope(xs: Iterable[float], ys: Iterable[float]) -> float:
    x = np.asarray(list(xs), dtype=float)
    y = np.asarray(list(ys), dtype=float)
    mask = (x > 0) & (y > 0)
    if np.count_nonzero(mask) < 2:
        return float("nan")
    return float(np.polyfit(np.log(x[mask]), np.log(y[mask]), 1)[0])


def run_scan(
    outdir: Path,
    seed: int,
    sizes: Sequence[int],
    samples: int,
    j_min: float,
    j_max: float,
    cutoffs: Sequence[float],
    times: Sequence[float],
) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    d = 2
    direction = np.array([1.0, 0.0])
    tail_raw: list[dict[str, float | int]] = []
    heat_raw: list[dict[str, float | int]] = []
    identity_errors: list[float] = []
    total_hminus1_values: list[float] = []

    max_cutoff = float(max(cutoffs))
    for L in sizes:
        for sample_index in range(samples):
            sample = make_torus_sample(L, d, rng, j_min, j_max, direction)
            orbit_generator, orbit_drift = build_environment_orbit_generator(sample)
            gen_error = float(
                spla.norm(orbit_generator - sample.laplacian, ord="fro")
            )
            drift_error = float(np.max(np.abs(orbit_drift - sample.drift)))
            identity_errors.append(max(gen_error, drift_error))
            if identity_errors[-1] > 1e-11:
                raise AssertionError(
                    f"environment/physical intertwining failed: {identity_errors[-1]}"
                )

            vals, vecs = low_spectrum(sample.laplacian, max_cutoff=max_cutoff)
            if vals[-1] <= max_cutoff and vals.size < sample.n_sites:
                raise AssertionError(
                    f"Low-spectrum audit did not pass cutoff: L={L}, "
                    f"largest={vals[-1]}, cutoff={max_cutoff}"
                )

            total_hminus1 = exact_total_hminus1(sample)
            total_hminus1_values.append(total_hminus1)
            for delta in cutoffs:
                mass, tail = spectral_tail_from_pairs(sample, vals, vecs, delta)
                tail_raw.append(
                    {
                        "L": L,
                        "N": sample.n_sites,
                        "sample": sample_index,
                        "delta": float(delta),
                        "source_mass_per_site": mass,
                        "Hminus1_tail_per_site": tail,
                        "tail_over_delta": tail / float(delta),
                        "tail_fraction_of_total_Hminus1": (
                            tail / total_hminus1 if total_hminus1 > 0 else 0.0
                        ),
                    }
                )

            for time in times:
                heat = heat_correlation(sample, float(time))
                heat_raw.append(
                    {
                        "L": L,
                        "N": sample.n_sites,
                        "sample": sample_index,
                        "time": float(time),
                        "heat_correlation": heat,
                        "scaled_t2_heat": heat * (1.0 + float(time)) ** 2,
                    }
                )
            print(
                f"PASS sample L={L:>2} idx={sample_index:>2}: "
                f"lambda_max_audited={vals[-1]:.4f}, "
                f"total_Hminus1/N={total_hminus1:.6g}"
            )

    write_csv(outdir / "finite_torus_source_tail_raw.csv", tail_raw)
    write_csv(outdir / "finite_torus_heat_decay_raw.csv", heat_raw)

    tail_summary: list[dict[str, float | int]] = []
    for L in sizes:
        for delta in cutoffs:
            values = np.asarray(
                [
                    float(row["Hminus1_tail_per_site"])
                    for row in tail_raw
                    if row["L"] == L and row["delta"] == float(delta)
                ]
            )
            ratios = values / float(delta)
            tail_summary.append(
                {
                    "L": L,
                    "N": L**d,
                    "samples": samples,
                    "delta": float(delta),
                    "mean_Hminus1_tail": float(np.mean(values)),
                    "std_Hminus1_tail": float(np.std(values, ddof=1))
                    if samples > 1
                    else 0.0,
                    "mean_tail_over_delta": float(np.mean(ratios)),
                    "max_tail_over_delta": float(np.max(ratios)),
                }
            )
    write_csv(outdir / "finite_torus_source_tail_summary.csv", tail_summary)

    heat_summary: list[dict[str, float | int]] = []
    for L in sizes:
        for time in times:
            values = np.asarray(
                [
                    float(row["heat_correlation"])
                    for row in heat_raw
                    if row["L"] == L and row["time"] == float(time)
                ]
            )
            heat_summary.append(
                {
                    "L": L,
                    "N": L**d,
                    "samples": samples,
                    "time": float(time),
                    "mean_heat_correlation": float(np.mean(values)),
                    "std_heat_correlation": float(np.std(values, ddof=1))
                    if samples > 1
                    else 0.0,
                    "mean_scaled_t2_heat": float(
                        np.mean(values) * (1.0 + float(time)) ** 2
                    ),
                }
            )
    write_csv(outdir / "finite_torus_heat_decay_summary.csv", heat_summary)

    largest_L = max(sizes)
    largest_tail_means = [
        float(row["mean_Hminus1_tail"])
        for row in tail_summary
        if row["L"] == largest_L
    ]
    largest_heat_means = [
        float(row["mean_heat_correlation"])
        for row in heat_summary
        if row["L"] == largest_L
    ]
    tail_slope = fit_log_slope(cutoffs, largest_tail_means)
    # Use the later half of the time grid for the decay fit.
    start = len(times) // 2
    heat_slope = fit_log_slope(times[start:], largest_heat_means[start:])

    max_ratio = max(float(row["max_tail_over_delta"]) for row in tail_summary)
    max_scaled_heat = max(float(row["mean_scaled_t2_heat"]) for row in heat_summary)
    certificate = {
        "seed": seed,
        "dimension": d,
        "direction": direction.tolist(),
        "conductance_interval": [j_min, j_max],
        "sizes": list(sizes),
        "samples_per_size": samples,
        "cutoffs": list(cutoffs),
        "times": list(times),
        "max_intertwining_error": max(identity_errors, default=0.0),
        "largest_L_tail_loglog_slope": tail_slope,
        "largest_L_heat_loglog_slope_late_times": heat_slope,
        "max_observed_tail_over_delta": max_ratio,
        "max_observed_mean_scaled_t2_heat": max_scaled_heat,
        "mean_total_Hminus1_per_site": float(np.mean(total_hminus1_values)),
        "checks": {
            "intertwining_machine_precision": max(identity_errors, default=0.0)
            < 1e-11,
            "all_tails_nonnegative": all(
                float(row["Hminus1_tail_per_site"]) >= -1e-12 for row in tail_raw
            ),
            "all_heat_correlations_nonnegative": all(
                float(row["heat_correlation"]) >= -1e-12 for row in heat_raw
            ),
            # Diagnostic, deliberately loose: the theorem has an unspecified
            # constant depending on ellipticity.  This check only detects a
            # catastrophic finite-size growth in the seeded scan.
            "tail_over_delta_stays_bounded_in_scan": max_ratio < 10.0,
            "scaled_t2_heat_stays_bounded_in_scan": max_scaled_heat < 20.0,
        },
    }
    certificate["overall"] = (
        "PASS" if all(certificate["checks"].values()) else "FAIL"
    )
    (outdir / "finite_torus_tail_certificate.json").write_text(
        json.dumps(certificate, indent=2) + "\n"
    )

    lines = [
        "Quantitative finite-torus source-tail numerical audit",
        "====================================================",
        f"seed: {seed}",
        f"dimension: {d}",
        f"conductances: Uniform[{j_min}, {j_max}]",
        f"sizes: {list(sizes)}",
        f"samples per size: {samples}",
        f"max environment/physical intertwining error: "
        f"{certificate['max_intertwining_error']:.3e}",
        f"largest-L tail log-log slope: {tail_slope:.6f} "
        "(theorem predicts exponent 1 in d=2)",
        f"largest-L late-time heat log-log slope: {heat_slope:.6f} "
        "(theorem predicts exponent -2 in d=2 before finite-size saturation)",
        f"max observed tail/delta ratio: {max_ratio:.6f}",
        f"max observed mean (1+t)^2 H_L(t): {max_scaled_heat:.6f}",
        f"OVERALL: {certificate['overall']}",
        "",
        "Interpretation: the exact theorem is proved analytically in the paper.",
        "The fitted slopes are finite-size diagnostics and are not used as proof.",
    ]
    (outdir / "FINITE_TORUS_VERIFICATION_OUTPUT.txt").write_text(
        "\n".join(lines) + "\n"
    )
    print("\n" + "\n".join(lines))
    return certificate


def parse_float_list(value: str) -> tuple[float, ...]:
    return tuple(float(part.strip()) for part in value.split(",") if part.strip())


def parse_int_list(value: str) -> tuple[int, ...]:
    return tuple(int(part.strip()) for part in value.split(",") if part.strip())


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path("."))
    parser.add_argument("--seed", type=int, default=20260727)
    parser.add_argument("--sizes", type=parse_int_list, default=(16, 24, 32, 48))
    parser.add_argument("--samples", type=int, default=8)
    parser.add_argument("--j-min", type=float, default=0.5)
    parser.add_argument("--j-max", type=float, default=2.0)
    parser.add_argument(
        "--cutoffs",
        type=parse_float_list,
        default=(0.08, 0.12, 0.18, 0.27, 0.40),
    )
    parser.add_argument(
        "--times",
        type=parse_float_list,
        default=(0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0, 64.0),
    )
    args = parser.parse_args()
    if args.j_min <= 0 or args.j_max <= args.j_min:
        raise ValueError("Require 0 < j_min < j_max")
    if args.samples < 1:
        raise ValueError("samples must be positive")
    args.outdir.mkdir(parents=True, exist_ok=True)
    cert = run_scan(
        outdir=args.outdir,
        seed=args.seed,
        sizes=args.sizes,
        samples=args.samples,
        j_min=args.j_min,
        j_max=args.j_max,
        cutoffs=args.cutoffs,
        times=args.times,
    )
    return 0 if cert["overall"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
