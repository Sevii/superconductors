#!/usr/bin/env python3
"""Verify exact identities used in the quenched Abel-tail fluctuation proof.

Checks on small random periodic tori:
  1. Abel tail equals the integrated heat correlation.
  2. Abel tail equals the parabolic cell-energy excess.
  3. The exact single-bond derivative formula agrees with centered finite differences.
  4. The samplewise envelopes T_sharp(delta) <= e A_Abel(delta) and
     T_sharp(delta) <= 4 R_1(delta) hold.

The checks are numerical diagnostics only; the paper gives the algebraic proofs.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import scipy.linalg as sla
import scipy.integrate as sint


def geometry(L: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    tails, heads, directions = [], [], []
    idx = lambda x, y: (x % L) * L + (y % L)
    for x in range(L):
        for y in range(L):
            v = idx(x, y)
            for mu, w in enumerate((idx(x + 1, y), idx(x, y + 1))):
                tails.append(v)
                heads.append(w)
                directions.append(mu)
    return np.asarray(tails), np.asarray(heads), np.asarray(directions)


def assemble(
    L: int,
    j: np.ndarray,
    tails: np.ndarray,
    heads: np.ndarray,
    directions: np.ndarray,
    A: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    N = L * L
    E = len(j)
    B = np.zeros((N, E))
    B[tails, np.arange(E)] = -1.0
    B[heads, np.arange(E)] = 1.0
    a = A[directions]
    lap = (B * j[None, :]) @ B.T
    g = B @ (j * a)
    g -= g.mean()
    return B, lap, g


def spectral_data(lap: np.ndarray, g: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lam, vec = sla.eigh(lap, driver="evr", check_finite=False)
    keep = lam > 1e-11
    lam = lam[keep]
    vec = vec[:, keep]
    coeff = vec.T @ g
    return lam, vec, coeff


def abel_tail(lap: np.ndarray, g: np.ndarray, delta: float) -> float:
    lam, _, coeff = spectral_data(lap, g)
    return float(np.sum(np.exp(-lam / delta) * coeff * coeff / lam) / lap.shape[0])


def sharp_tail(lap: np.ndarray, g: np.ndarray, delta: float) -> float:
    lam, _, coeff = spectral_data(lap, g)
    mask = lam <= delta
    return float(np.sum(coeff[mask] ** 2 / lam[mask]) / lap.shape[0])


def richardson_one(lap: np.ndarray, g: np.ndarray, delta: float) -> float:
    lam, _, coeff = spectral_data(lap, g)
    return float(np.sum(delta**2 * coeff**2 / (lam * (lam + delta) ** 2)) / lap.shape[0])


def heat_integral(lap: np.ndarray, g: np.ndarray, delta: float) -> float:
    lam, _, coeff = spectral_data(lap, g)
    tau = 1.0 / delta
    # Independent numerical quadrature of the scalar heat correlation.
    integrand = lambda t: float(np.sum(np.exp(-t * lam) * coeff * coeff) / lap.shape[0])
    val, _ = sint.quad(integrand, tau, np.inf, epsabs=2e-11, epsrel=2e-11, limit=250)
    return float(val)


def cell_energy(j: np.ndarray, a: np.ndarray, B: np.ndarray, phi: np.ndarray) -> float:
    r = a + B.T @ phi
    return float(np.dot(j * r, r) / B.shape[0])


def parabolic_excess(
    j: np.ndarray,
    a: np.ndarray,
    B: np.ndarray,
    lap: np.ndarray,
    g: np.ndarray,
    delta: float,
) -> float:
    lam, vec, coeff = spectral_data(lap, g)
    inv_g = vec @ (coeff / lam)
    t = 0.5 / delta
    exp_inv_g = vec @ (np.exp(-t * lam) * coeff / lam)
    phi_star = -inv_g
    phi_t = -(inv_g - exp_inv_g)
    return cell_energy(j, a, B, phi_t) - cell_energy(j, a, B, phi_star)


def derivative_spectral(
    B: np.ndarray,
    lap: np.ndarray,
    g: np.ndarray,
    a: np.ndarray,
    edge: int,
    delta: float,
) -> float:
    """Frechet derivative of <g, f(L)g>/N with f(x)=exp(-x/delta)/x."""
    lam, vec, coeff = spectral_data(lap, g)
    bcoord = vec.T @ B[:, edge]
    f = np.exp(-lam / delta) / lam
    fp = -np.exp(-lam / delta) * (1.0 / (delta * lam) + 1.0 / lam**2)
    # Loewner divided-difference matrix.
    li = lam[:, None]
    lj = lam[None, :]
    denom = li - lj
    dd = np.empty_like(denom)
    off = np.abs(denom) > 1e-10
    dd[off] = (f[:, None] - f[None, :])[off] / denom[off]
    dd[~off] = ((fp[:, None] + fp[None, :]) * 0.5)[~off]
    first = 2.0 * a[edge] * float(np.dot(bcoord * f, coeff))
    second = float(np.sum((coeff * bcoord)[:, None] * (coeff * bcoord)[None, :] * dd))
    return (first + second) / lap.shape[0]


def derivative_time_formula(
    B: np.ndarray,
    lap: np.ndarray,
    g: np.ndarray,
    a: np.ndarray,
    edge: int,
    delta: float,
) -> float:
    """Evaluate the exact Duhamel time formula by spectral quadrature."""
    lam, vec, coeff = spectral_data(lap, g)
    beta = (vec.T @ B[:, edge]) * coeff
    tau = 1.0 / delta
    # Linear term is explicit.
    linear = 2.0 * a[edge] * float(np.sum(beta * np.exp(-tau * lam) / lam))
    # The double-time term is the Loewner derivative written in an explicit
    # spectral form.  For i != j, integral_{s=tau}^inf integral_0^s
    # exp(-(s-r)li-r lj) dr ds = (e^-tau*lj/lj-e^-tau*li/li)/(li-lj).
    li = lam[:, None]
    lj = lam[None, :]
    denom = li - lj
    kernel = np.empty_like(denom)
    h = np.exp(-tau * lam) / lam
    off = np.abs(denom) > 1e-10
    kernel[off] = (h[None, :] - h[:, None])[off] / denom[off]
    # diagonal integral: int_tau^inf s exp(-s lambda) ds
    diagonal = np.exp(-tau * lam) * (tau / lam + 1.0 / lam**2)
    kernel[~off] = ((diagonal[:, None] + diagonal[None, :]) * 0.5)[~off]
    quadratic = float(np.sum(beta[:, None] * beta[None, :] * kernel))
    return (linear - quadratic) / lap.shape[0]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--seed", type=int, default=20260727)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(args.seed)
    L = 6
    A = np.asarray([1.0, -0.35])
    tails, heads, directions = geometry(L)
    a = A[directions]
    j = rng.uniform(0.55, 1.45, size=len(tails))
    B, lap, g = assemble(L, j, tails, heads, directions, A)
    delta = 0.72

    abel = abel_tail(lap, g, delta)
    heat = heat_integral(lap, g, delta)
    excess = parabolic_excess(j, a, B, lap, g, delta)
    sharp = sharp_tail(lap, g, delta)
    rich = richardson_one(lap, g, delta)

    derivative_rows = []
    h = 2.5e-6
    for edge in rng.choice(len(j), size=8, replace=False):
        jp = j.copy(); jp[edge] += h
        jm = j.copy(); jm[edge] -= h
        _, lp, gp = assemble(L, jp, tails, heads, directions, A)
        _, lm, gm = assemble(L, jm, tails, heads, directions, A)
        fd = (abel_tail(lp, gp, delta) - abel_tail(lm, gm, delta)) / (2.0 * h)
        exact = derivative_spectral(B, lap, g, a, int(edge), delta)
        time = derivative_time_formula(B, lap, g, a, int(edge), delta)
        derivative_rows.append(
            {
                "edge": int(edge),
                "finite_difference": fd,
                "spectral_derivative": exact,
                "duhamel_derivative": time,
                "fd_abs_error": abs(fd - exact),
                "duhamel_abs_error": abs(time - exact),
            }
        )

    cert = {
        "configuration": {"L": L, "delta": delta, "seed": args.seed, "A": A.tolist()},
        "abel_tail": abel,
        "integrated_heat_tail": heat,
        "parabolic_energy_excess": excess,
        "sharp_tail": sharp,
        "richardson_k1_tail": rich,
        "abel_heat_abs_error": abs(abel - heat),
        "abel_energy_excess_abs_error": abs(abel - excess),
        "max_finite_difference_derivative_error": max(r["fd_abs_error"] for r in derivative_rows),
        "max_duhamel_derivative_error": max(r["duhamel_abs_error"] for r in derivative_rows),
        "abel_envelope_holds": bool(sharp <= math.e * abel + 1e-11),
        "richardson_envelope_holds": bool(sharp <= 4.0 * rich + 1e-11),
        "derivatives": derivative_rows,
    }
    thresholds = {
        "abel_heat_abs_error": 2e-9,
        "abel_energy_excess_abs_error": 2e-10,
        "max_finite_difference_derivative_error": 2e-7,
        "max_duhamel_derivative_error": 2e-10,
    }
    cert["thresholds"] = thresholds
    cert["status"] = "PASS" if (
        all(cert[k] <= v for k, v in thresholds.items())
        and cert["abel_envelope_holds"]
        and cert["richardson_envelope_holds"]
    ) else "FAIL"

    json_path = args.output_dir / "abel_tail_identity_certificate.json"
    json_path.write_text(json.dumps(cert, indent=2, sort_keys=True) + "\n")
    report = [
        "ABEL-TAIL IDENTITY AND SENSITIVITY AUDIT",
        "=" * 49,
        f"L={L}, delta={delta}, seed={args.seed}",
        f"Abel vs heat absolute error: {cert['abel_heat_abs_error']:.6e}",
        f"Abel vs energy-excess absolute error: {cert['abel_energy_excess_abs_error']:.6e}",
        "Maximum finite-difference derivative error: "
        f"{cert['max_finite_difference_derivative_error']:.6e}",
        "Maximum Duhamel/spectral derivative error: "
        f"{cert['max_duhamel_derivative_error']:.6e}",
        f"Sharp <= e Abel: {cert['abel_envelope_holds']}",
        f"Sharp <= 4 Richardson-k1: {cert['richardson_envelope_holds']}",
        f"OVERALL: {cert['status']}",
    ]
    (args.output_dir / "ABEL_TAIL_IDENTITY_VERIFICATION_OUTPUT.txt").write_text("\n".join(report) + "\n")
    print("\n".join(report))
    return 0 if cert["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
