#!/usr/bin/env python3
"""Deterministic checks for the abstract zero-frequency response theorems.

The script verifies:
  1. the positive-square target-space finite-frequency identity;
  2. the static projected-source formula;
  3. the quantitative H^{-1} tail bounds;
  4. the exact soft-target order-of-limits counterexample;
  5. transverse diamagnetic/paramagnetic cancellation for a clean free gas.

These finite-dimensional and numerical checks support, but do not replace, the proofs.
"""
from __future__ import annotations

import platform
import sys

import numpy as np
import scipy
import scipy.integrate
import scipy.linalg

TOL = 2.0e-10


def opnorm(a: np.ndarray) -> float:
    return float(scipy.linalg.norm(a, 2))


def herm(a: np.ndarray) -> np.ndarray:
    return (a + a.conj().T) / 2.0


def random_unitary(rng: np.random.Generator, n: int) -> np.ndarray:
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, r = np.linalg.qr(z)
    phases = np.diag(r)
    phases = np.where(np.abs(phases) > 0, phases / np.abs(phases), 1.0)
    return q @ np.diag(phases.conj())


def psd_function(a: np.ndarray, fn, cutoff: float = 1.0e-12) -> np.ndarray:
    vals, vecs = scipy.linalg.eigh(herm(a))
    out = np.array([fn(float(x)) if x > cutoff else fn(0.0) for x in vals])
    return (vecs * out) @ vecs.conj().T


def check(name: str, error: float, tol: float = TOL) -> None:
    status = "PASS" if error <= tol else "FAIL"
    print(f"{name:<72s} {status}  error={error:.3e}  tol={tol:.1e}")
    if error > tol:
        raise AssertionError(name)


def report(name: str, value: float) -> None:
    print(f"{name:<72s} VALUE {value:.12g}")


def check_positive_square_identity() -> None:
    rng = np.random.default_rng(2026072601)
    max_dynamic = 0.0
    max_static = 0.0
    max_defect_psd = 0.0

    for trial in range(12):
        n = 7
        m = 9
        rank = 4
        ut = random_unitary(rng, m)
        vs = random_unitary(rng, n)
        singular = np.array([0.45, 0.8, 1.25, 2.1]) * (1.0 + 0.03 * trial)
        sigma = np.zeros((m, n), dtype=complex)
        sigma[:rank, :rank] = np.diag(singular)
        d = ut @ sigma @ vs.conj().T
        dp = rng.normal(size=(m, n)) + 1j * rng.normal(size=(m, n))
        dp /= np.sqrt(2.0 * m)

        p = vs[:, rank:] @ vs[:, rank:].conj().T
        q = np.eye(n) - p
        h = herm(d.conj().T @ d / 2.0)
        ell = herm(d @ d.conj().T / 2.0)
        jop = herm((dp.conj().T @ d + d.conj().T @ dp) / 2.0)
        top = herm(dp.conj().T @ dp)
        s = dp @ p

        zeta = 0.17 + 0.09 * trial
        hvals, hvecs = scipy.linalg.eigh(h)
        h_filter = np.where(hvals > 1.0e-12, hvals / (hvals * hvals + zeta * zeta), 0.0)
        hreg = (hvecs * h_filter) @ hvecs.conj().T
        k_state = p @ top @ p - 2.0 * p @ jop @ q @ hreg @ q @ jop @ p

        lvals, lvecs = scipy.linalg.eigh(ell)
        l_filter = zeta * zeta / (lvals * lvals + zeta * zeta)
        freg = (lvecs * l_filter) @ lvecs.conj().T
        k_target = s.conj().T @ freg @ s
        max_dynamic = max(max_dynamic, opnorm(k_state - k_target))

        hplus_filter = np.where(hvals > 1.0e-12, 1.0 / hvals, 0.0)
        hplus = (hvecs * hplus_filter) @ hvecs.conj().T
        c_state = p @ top @ p - 2.0 * p @ jop @ q @ hplus @ q @ jop @ p
        pker_l = lvecs[:, lvals < 1.0e-12] @ lvecs[:, lvals < 1.0e-12].conj().T
        c_target = s.conj().T @ pker_l @ s
        max_static = max(max_static, opnorm(c_state - c_target))

        eigmin = float(np.min(scipy.linalg.eigvalsh(herm(k_state - c_state))))
        max_defect_psd = max(max_defect_psd, max(0.0, -eigmin))

    check("positive-square finite-frequency target identity", max_dynamic, 5.0e-10)
    check("positive-square static projected-source identity", max_static, 5.0e-10)
    check("finite-frequency defect is positive semidefinite", max_defect_psd, 5.0e-10)


def check_hminus1_bounds() -> None:
    rng = np.random.default_rng(2026072602)
    max_lower_violation = 0.0
    max_upper_violation = 0.0

    dim = 5
    energies = np.array([0.018, 0.061, 0.19, 0.73, 2.4, 7.0])
    weights = []
    for _ in energies:
        a = rng.normal(size=(3, dim)) + 1j * rng.normal(size=(3, dim))
        w = a.conj().T @ a
        weights.append(w / (8.0 * opnorm(w)))

    m_total = sum(w / e for w, e in zip(weights, energies))
    M = opnorm(m_total)

    for zeta, delta in [(0.025, 0.2), (0.08, 0.4), (0.23, 1.0)]:
        defect = sum(
            2.0 * zeta * zeta / (e * (e * e + zeta * zeta)) * w
            for w, e in zip(weights, energies)
        )
        low_zeta = sum((w / e for w, e in zip(weights, energies) if e <= zeta), np.zeros((dim, dim), complex))
        low_delta = sum((w / e for w, e in zip(weights, energies) if e <= delta), np.zeros((dim, dim), complex))
        lower_violation = max(0.0, opnorm(low_zeta) - opnorm(defect))
        upper_rhs = 2.0 * opnorm(low_delta) + 2.0 * (zeta * zeta / (delta * delta)) * M
        upper_violation = max(0.0, opnorm(defect) - upper_rhs)
        max_lower_violation = max(max_lower_violation, lower_violation)
        max_upper_violation = max(max_upper_violation, upper_violation)

    check("quantitative H^{-1} lower-tail bound", max_lower_violation)
    check("quantitative H^{-1} upper-tail bound", max_upper_violation)
    report("total current H^{-1} norm in certificate", M)


def check_soft_target_counterexample() -> None:
    alpha = 0.37
    c0 = 0.21
    zeta = 0.12
    deltas = [0.8, 0.25, 0.06, 0.009]
    max_formula_error = 0.0

    for L, delta in enumerate(deltas, start=3):
        volume = float(L**3)
        # Target basis is (r,z), state basis is (g,e).
        d0 = np.array([[0.0, np.sqrt(2.0 * delta)], [0.0, 0.0]], dtype=complex)
        dp = np.array(
            [[np.sqrt(2.0 * alpha * volume), 0.0], [np.sqrt(c0 * volume), 0.0]],
            dtype=complex,
        )
        p = np.diag([1.0, 0.0])
        s = dp @ p
        ell = d0 @ d0.conj().T / 2.0
        vals, vecs = scipy.linalg.eigh(ell)
        fz = (vecs * (zeta * zeta / (vals * vals + zeta * zeta))) @ vecs.conj().T
        pker = vecs[:, vals < 1.0e-12] @ vecs[:, vals < 1.0e-12].conj().T
        dynamic = float(np.real((s.conj().T @ fz @ s)[0, 0]) / volume)
        static = float(np.real((s.conj().T @ pker @ s)[0, 0]) / volume)
        dynamic_exact = c0 + 2.0 * alpha * zeta * zeta / (delta * delta + zeta * zeta)
        max_formula_error = max(max_formula_error, abs(dynamic - dynamic_exact), abs(static - c0))

    check("soft-target counterexample formulas", max_formula_error)
    thermodynamic_first = c0 + 2.0 * alpha
    adiabatic_first = c0
    report("soft-target thermodynamic-first limit", thermodynamic_first)
    report("soft-target adiabatic-first limit", adiabatic_first)
    check("soft-target limits differ by 2 alpha", abs((thermodynamic_first - adiabatic_first) - 2.0 * alpha))


def check_free_gas_cancellation() -> None:
    # Two-dimensional continuum gas with m=1.  A smooth Fermi function avoids
    # distributional derivatives; integration by parts is exact as kmax -> infinity.
    mass = 1.0
    beta = 9.0
    mu = 1.15

    def fermi(e: float) -> float:
        x = beta * (e - mu)
        if x > 45.0:
            return 0.0
        if x < -45.0:
            return 1.0
        return 1.0 / (np.exp(x) + 1.0)

    def fermi_prime(e: float) -> float:
        f = fermi(e)
        return -beta * f * (1.0 - f)

    kmax = 9.0
    density, _ = scipy.integrate.quad(
        lambda k: k * fermi(k * k / (2.0 * mass)) / (2.0 * np.pi),
        0.0,
        kmax,
        epsabs=1.0e-12,
        epsrel=1.0e-12,
        limit=300,
    )
    param, _ = scipy.integrate.quad(
        lambda k: k * (k * k / (2.0 * mass * mass))
        * fermi_prime(k * k / (2.0 * mass))
        / (2.0 * np.pi),
        0.0,
        kmax,
        epsabs=1.0e-12,
        epsrel=1.0e-12,
        limit=300,
    )
    residual = density / mass + param
    check("clean-gas finite-beta transverse cancellation", abs(residual), 2.0e-10)
    report("clean-gas Kohn/Drude density n/m", density / mass)
    report("clean-gas finite-beta paramagnetic term", param)


def main() -> int:
    print("ZERO-FREQUENCY DYNAMICAL RESPONSE VERIFIER")
    print("python", platform.python_version())
    print("numpy", np.__version__)
    print("scipy", scipy.__version__)
    print("-" * 108)
    check_positive_square_identity()
    check_hminus1_bounds()
    check_soft_target_counterexample()
    check_free_gas_cancellation()
    print("-" * 108)
    print("OVERALL: PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"OVERALL: FAIL ({type(exc).__name__}: {exc})", file=sys.stderr)
        raise
