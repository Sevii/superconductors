#!/usr/bin/env python3
"""Deterministic certificates for infrared-stable Meissner transfer.

Checks random target-metric identities, the current H^{-1} defect bound,
local row-ideal factorization, and a 2D infrared example separating anisotropic
and point-group-scalar target metrics.
"""
from __future__ import annotations
import math
import sys
from dataclasses import dataclass
import numpy as np
import scipy
import scipy.linalg as la

TOL = 2.5e-9
RNG = np.random.default_rng(20260726)


def hermitian(a: np.ndarray) -> np.ndarray:
    return (a + a.conj().T) / 2.0


def sqrt_pos(a: np.ndarray) -> np.ndarray:
    w, v = la.eigh(hermitian(a))
    if np.min(w) <= 0:
        raise ValueError("matrix is not positive definite")
    return (v * np.sqrt(w)) @ v.conj().T


def psd_min_eig(a: np.ndarray) -> float:
    return float(np.min(la.eigvalsh(hermitian(a))))


def opnorm(a: np.ndarray) -> float:
    return float(la.norm(a, 2))


def random_unitary(n: int) -> np.ndarray:
    z = RNG.normal(size=(n, n)) + 1j * RNG.normal(size=(n, n))
    q, r = la.qr(z)
    diag = np.diag(r)
    phase = np.where(np.abs(diag) > 0, diag / np.abs(diag), 1.0)
    return q @ np.diag(phase.conj())


@dataclass
class RandomCheckResult:
    delta: float
    static_lower_margin: float
    dynamic_lower_margin: float
    dynamic_upper_margin: float
    defect_bound_margin: float
    target_state_formula_error: float


def random_target_metric_check(delta: float, eta: float = 0.08, zeta: float = 0.031) -> RandomCheckResult:
    n_phys, n_target, rank = 7, 10, 5
    ground_dim = n_phys - rank
    u, v = random_unitary(n_target), random_unitary(n_phys)
    singular = np.array([delta, 0.37, 0.71, 1.12, 1.53])
    d = u[:, :rank] @ np.diag(singular) @ v[:, :rank].conj().T

    left_null = u[:, rank:]
    coeff = RNG.normal(size=(n_target - rank, ground_dim)) + 1j * RNG.normal(
        size=(n_target - rank, ground_dim)
    )
    s = left_null @ coeff
    s /= math.sqrt(opnorm(hermitian(s.conj().T @ s)))
    c0 = hermitian(s.conj().T @ s)

    x = hermitian(RNG.normal(size=(n_target, n_target)) + 1j * RNG.normal(size=(n_target, n_target)))
    x /= opnorm(x)
    g = np.eye(n_target) + 2.0 * eta * x
    eig_g = la.eigvalsh(g)
    m, M = float(np.min(eig_g)), float(np.max(eig_g))
    if m <= 0:
        raise AssertionError("target metric lost positivity")

    gs = sqrt_pos(g)
    dt, st = gs @ d, gs @ s
    ell = hermitian(0.5 * dt @ dt.conj().T)
    p_range = dt @ la.pinv(dt, rtol=1e-12)
    c = hermitian(st.conj().T @ (np.eye(n_target) - p_range) @ st)
    reg = zeta**2 * la.inv(ell @ ell + zeta**2 * np.eye(n_target))
    k_target = hermitian(st.conj().T @ reg @ st)

    h = hermitian(0.5 * dt.conj().T @ dt)
    jp = 0.5 * dt.conj().T @ st
    stress = hermitian(st.conj().T @ st)
    hreg = h @ la.inv(h @ h + zeta**2 * np.eye(n_phys))
    k_state = hermitian(stress - 2.0 * jp.conj().T @ hreg @ jp)

    rhs_defect = (4.0 * eta**2 / m) * hermitian((x @ s).conj().T @ (x @ s))
    defect = hermitian(k_target - c)
    return RandomCheckResult(
        delta,
        psd_min_eig(c - m * c0),
        psd_min_eig(k_target - c),
        psd_min_eig(M * c0 - k_target),
        psd_min_eig(rhs_defect - defect),
        opnorm(k_target - k_state),
    )


def local_row_ideal_check() -> tuple[float, float]:
    n_phys, n_target, rank = 8, 11, 6
    u, v = random_unitary(n_target), random_unitary(n_phys)
    singular = np.array([0.11, 0.29, 0.52, 0.83, 1.21, 1.67])
    d = u[:, :rank] @ np.diag(singular) @ v[:, :rank].conj().T
    dplus = la.pinv(d, rtol=1e-13)
    q = hermitian(dplus @ d)
    raw = hermitian(RNG.normal(size=(n_phys, n_phys)) + 1j * RNG.normal(size=(n_phys, n_phys)))
    r = hermitian(q @ raw @ q)
    x = hermitian(dplus.conj().T @ r @ dplus)
    residual = opnorm(d.conj().T @ x @ d - r)
    ratio = opnorm(x) / (opnorm(dplus) ** 2 * opnorm(r))
    return residual, ratio


def cycle_symbols(L: int) -> tuple[complex, complex]:
    qx = 2.0 * math.pi / L
    return np.exp(1j * qx) - 1.0, 0.0j


def torus_infrared_check(eta: float = 0.1) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    x_swap = np.array([[0.0, 1.0], [1.0, 0.0]])
    t = np.array([0.0, 1.0], dtype=complex)
    for L in [8, 16, 32, 64, 128, 256]:
        dx, dy = cycle_symbols(L)
        dvec = np.array([dx, dy], dtype=complex)
        g = np.eye(2) + 2.0 * eta * x_swap
        energy = 0.5 * float(np.real(np.vdot(dvec, g @ dvec)))
        jamp = 0.5 * np.vdot(dvec, g @ t)
        stress = float(np.real(np.vdot(t, g @ t)))
        static = stress - (2.0 * abs(jamp) ** 2 / energy if energy > 0 else 0.0)
        zeta = 1.0 / L
        dynamic = stress - 2.0 * abs(jamp) ** 2 * energy / (energy**2 + zeta**2)
        hminus1 = abs(jamp) ** 2 / energy if energy > 0 else 0.0

        g_iso = (1.0 + 2.0 * eta) * np.eye(2)
        energy_iso = 0.5 * float(np.real(np.vdot(dvec, g_iso @ dvec)))
        j_iso = 0.5 * np.vdot(dvec, g_iso @ t)
        stress_iso = float(np.real(np.vdot(t, g_iso @ t)))
        static_iso = stress_iso - (2.0 * abs(j_iso) ** 2 / energy_iso if energy_iso > 0 else 0.0)
        rows.append(
            dict(
                L=float(L),
                coclosed_error=abs(np.vdot(dvec, t)),
                energy=energy,
                static_aniso=static,
                dynamic_aniso=dynamic,
                drude_minus_meissner=stress - static,
                hminus1_weight=hminus1,
                static_iso=static_iso,
                isotropic_source=abs(j_iso),
            )
        )
    return rows


def main() -> int:
    print("MICROSCOPIC MEISSNER H^-1 TRANSFER CERTIFICATE")
    print("numpy", np.__version__)
    print("scipy", scipy.__version__)
    print()
    ok = True

    print("[1] Random target-metric inequalities with a closing singular value")
    for delta in [1e-1, 1e-2, 1e-3, 1e-4, 1e-6]:
        r = random_target_metric_check(delta)
        print(
            f"delta={r.delta:8.1e} static-margin={r.static_lower_margin:+.3e} "
            f"K-C={r.dynamic_lower_margin:+.3e} upper-margin={r.dynamic_upper_margin:+.3e} "
            f"defect-margin={r.defect_bound_margin:+.3e} target/state err={r.target_state_formula_error:.3e}"
        )
        ok &= min(r.static_lower_margin, r.dynamic_lower_margin, r.dynamic_upper_margin, r.defect_bound_margin) > -TOL
        ok &= r.target_state_formula_error < 5e-8

    print()
    print("[2] Local two-sided row-ideal factorization")
    residual, ratio = local_row_ideal_check()
    print(f"||D^* X D-r|| = {residual:.3e}")
    print(f"||X||/(||D^+||^2 ||r||) = {ratio:.6f}")
    ok &= residual < 5e-8 and ratio <= 1.0 + 5e-8

    print()
    print("[3] 2D transverse infrared audit")
    print("Expected anisotropic mismatch = 4 eta^2 = 0.04; isotropic mismatch = 0")
    rows = torus_infrared_check()
    for row in rows:
        print(
            f"L={int(row['L']):3d} E={row['energy']:.3e} C_aniso={row['static_aniso']:.10f} "
            f"K(zeta=1/L)={row['dynamic_aniso']:.10f} D-C={row['drude_minus_meissner']:.10f} "
            f"H^-1={row['hminus1_weight']:.10f} C_iso={row['static_iso']:.10f} "
            f"|J_iso|={row['isotropic_source']:.3e}"
        )
        ok &= row["coclosed_error"] < 1e-12
        ok &= abs(row["static_aniso"] - 0.96) < 2e-10
        ok &= abs(row["drude_minus_meissner"] - 0.04) < 2e-10
        ok &= abs(row["hminus1_weight"] - 0.02) < 2e-10
        ok &= abs(row["static_iso"] - 1.2) < 2e-10
        ok &= row["isotropic_source"] < 1e-12
    ok &= abs(rows[-1]["dynamic_aniso"] - 1.0) < 5e-4

    print()
    print("OVERALL:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
