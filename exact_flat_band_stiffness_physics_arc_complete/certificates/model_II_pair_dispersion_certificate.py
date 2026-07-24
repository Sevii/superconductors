#!/usr/bin/env python3
"""Certificate for the Model-II pair-dispersion lemma and stiffness theorem.

The script checks four independent statements for the projected PSD-Hubbard
Model II of Gao-Han-Khalaf:

1. the exact one-pair formula
       E_pair(Q) = |U|/4 * (1 - |F_N(Q)|),
2. the flat-connection identity E_1(A) = E_pair(2A),
3. the exact finite-grid mass m_pair^{-1} = |U| xi^2/8 I for N_x,N_y >= 3,
4. the finite-size many-pair curvature
       C_n = n(V-n)/(V-1) * C_1.

The direct Q-block matrix and the many-body fixed-number implementation are
constructed independently.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import eigvalsh
from scipy.special import j0

ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import ghk_models_ed_search as ghk  # independent fixed-number implementation

Array = NDArray[np.complex128]


def momentum_grid(nx: int, ny: int) -> NDArray[np.float64]:
    return np.array(
        [(2 * np.pi * ix / nx, 2 * np.pi * iy / ny)
         for iy in range(ny) for ix in range(nx)],
        dtype=float,
    )


def alpha(k: NDArray[np.float64], xi: float) -> NDArray[np.float64]:
    return xi * (np.cos(k[:, 0]) + np.cos(k[:, 1]))


def lower_spinor(k: NDArray[np.float64], xi: float) -> Array:
    """A smooth lower-band gauge for H_II=-t(sigma_x sin alpha+sigma_y cos alpha)."""
    a = alpha(k, xi)
    return np.column_stack((np.exp(0.5j * a), 1j * np.exp(-0.5j * a))) / np.sqrt(2)


def form_factor(nx: int, ny: int, q: NDArray[np.float64], xi: float) -> complex:
    k = momentum_grid(nx, ny)
    ap = alpha(k + q[None, :] / 2, xi)
    am = alpha(k - q[None, :] / 2, xi)
    return complex(np.mean(np.exp(1j * (ap - am))))


def form_factor_thermodynamic(q: NDArray[np.float64], xi: float) -> float:
    return float(j0(2 * xi * np.sin(q[0] / 2)) * j0(2 * xi * np.sin(q[1] / 2)))


def pair_energy_formula(
    nx: int, ny: int, q: NDArray[np.float64], xi: float, U: float = 1.0
) -> float:
    return float(abs(U) * (1 - abs(form_factor(nx, ny, q, xi))) / 4)


def pair_q_block(
    nx: int, ny: int, q: NDArray[np.float64], xi: float, U: float = 1.0
) -> Array:
    """Direct one-up/one-down Q block of the projected PSD-Hubbard Hamiltonian."""
    k = momentum_grid(nx, ny)
    up_plus = lower_spinor(k + q[None, :] / 2, xi)
    up_minus = lower_spinor(k - q[None, :] / 2, xi)
    g = np.conjugate(up_plus) * up_minus  # columns are orbital pair form factors
    v = nx * ny
    kernel = sum(np.outer(g[:, a], np.conjugate(g[:, a])) for a in range(2)) / v
    h = abs(U) * (0.5 * np.eye(v, dtype=complex) - kernel)
    return (h + h.conj().T) / 2


def five_point_second(f, h: float = 2e-4) -> float:
    return float((-f(2*h) + 16*f(h) - 30*f(0.0) + 16*f(-h) - f(-2*h)) / (12*h*h))


def check_pair_formula() -> float:
    worst = 0.0
    for nx, ny in ((3, 3), (4, 5)):
        for xi in (0.3, 1.0, 1.4):
            for q in (np.array([0.274, -0.102]), np.array([0.04, 0.06])):
                direct = float(eigvalsh(pair_q_block(nx, ny, q, xi), subset_by_index=[0, 0])[0])
                formula = pair_energy_formula(nx, ny, q, xi)
                worst = max(worst, abs(direct - formula))
    assert worst < 2e-13, worst
    return worst


def check_twist_identity() -> float:
    worst = 0.0
    for nx, ny in ((3, 3), (4, 5)):
        for xi in (0.3, 1.0, 1.4):
            for avec in (np.array([0.137, -0.051]), np.array([0.02, 0.03])):
                h = ghk.HubbardSector(nx, ny, 1, avec, ghk.model_II_d, {"xi": xi}).dense()
                e_twist = float(eigvalsh(h, subset_by_index=[0, 0])[0])
                e_pair = pair_energy_formula(nx, ny, 2 * avec, xi)
                worst = max(worst, abs(e_twist - e_pair))
    assert worst < 2e-12, worst
    return worst


def check_mass() -> float:
    worst = 0.0
    for nx, ny in ((3, 3), (4, 5), (6, 7)):
        for xi in (0.3, 0.8, 1.4):
            target = xi * xi / 8
            for axis in (0, 1):
                def e(qval: float) -> float:
                    q = np.zeros(2)
                    q[axis] = qval
                    return pair_energy_formula(nx, ny, q, xi)
                measured = five_point_second(e)
                worst = max(worst, abs(measured - target))
    assert worst < 2e-7, worst
    return worst


def check_bessel_limit() -> float:
    worst = 0.0
    xi = 1.0
    for q in (np.array([0.2, -0.15]), np.array([0.05, 0.03])):
        finite = form_factor(160, 160, q, xi)
        infinite = form_factor_thermodynamic(q, xi)
        worst = max(worst, abs(finite - infinite))
    assert worst < 2e-14, worst
    return worst


def check_many_body_curvature() -> float:
    worst = 0.0
    nx = ny = 3
    v = nx * ny
    xi = 1.0
    c1_target = xi * xi / 2  # C_1 = 4 m_pair^{-1}
    for n in (1, 2, 3):
        row = ghk.curvature(nx, ny, n, ghk.model_II_d, {"xi": xi}, h=5e-4)
        rho = n * (v - n) / (v - 1)
        target = rho * c1_target
        worst = max(worst, abs(row["curvature"] - target))
        assert row["deg"] == 1
        assert row["gap"] > 0
    assert worst < 5e-7, worst
    return worst


def main() -> None:
    checks = {
        "direct Q-block vs exact pair formula": check_pair_formula(),
        "electronic twist vs Q=2A": check_twist_identity(),
        "finite-grid pair-mass Hessian": check_mass(),
        "thermodynamic Bessel form factor": check_bessel_limit(),
        "many-body finite-size curvature": check_many_body_curvature(),
    }
    print("Model-II pair-dispersion certificate")
    for name, err in checks.items():
        print(f"  {name}: max error = {err:.3e}")
    print("  exact mass tensor: m_pair^{-1} = |U| xi^2/8 * I_2")
    print("  exact finite-size stiffness: kappa_n = V/(V-1) nu(1-nu) |U|xi^2/8 * I_2")
    print("PASS")


if __name__ == "__main__":
    main()
