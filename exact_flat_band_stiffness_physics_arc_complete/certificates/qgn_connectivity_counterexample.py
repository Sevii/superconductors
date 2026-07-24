#!/usr/bin/env python3
"""Fixed-local QGN counterexample candidate based on projector reducibility.

This script studies a one-dimensional, two-orbital, rank-one flat band per spin
with Bloch spinors

    u_up(k)   = (cos k,  sin k)^T,
    u_down(k) = (cos k, -sin k)^T = u_up^*(-k).

The parent h_sigma(k)=I-P_sigma(k) has exact flat eigenvalues 0 and 1 and
Fourier range two, independent of system size.  Spin conservation and time
reversal put the model in the uniform-pairing p-p QGN class.  The interaction
is the projected positive-semidefinite Hubbard convention

    H(A) = (U/2) sum_{R,a} [nbar_{R a up}(A)-nbar_{R a down}(A)]^2,

with the electronic flat connection inserted as k -> k+A.

For every even M=2L, the projector contains only even real-space displacement,
so the band Hilbert space splits into even- and odd-cell components.  At n=L
pairs, filling one component completely and leaving the other empty gives an
exact zero-energy state for every A.  Hence the exact ground-energy curvature
and stiffness vanish for all even M.

The one-pair branch is nevertheless

    E_pair(Q) = (U/2) sin^2 Q                      (M >= 6),

so m_pair^{-1}=U.  At nu=1/2 and N_flat=2, the conjectured stiffness is U/4,
while the exact value is zero: R_M(1/2)=0 along the full even-M sequence.

The script independently builds the many-body Hamiltonian in a fixed-number
bit basis and checks finite-size curvatures for M=5,6,8.  M=5 is a useful
connected control: gcd(2,5)=1 and the filling law is recovered numerically.
The even tori expose the missing connectivity/irreducibility hypothesis.
"""
from __future__ import annotations

import argparse
import csv
import functools
import itertools
import math
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Iterable

import numpy as np
from scipy.linalg import eigh


@functools.lru_cache(None)
def masks(M: int, n: int) -> tuple[int, ...]:
    return tuple(sum(1 << i for i in comb) for comb in itertools.combinations(range(M), n))


@functools.lru_cache(None)
def bilinears(M: int, n: int) -> np.ndarray:
    """Return matrices of d_i^dagger d_j in the fixed-n real-space basis."""
    basis = masks(M, n)
    index = {mask: i for i, mask in enumerate(basis)}
    dim = len(basis)
    out = np.zeros((M, M, dim, dim), dtype=np.complex128)
    for i in range(M):
        for j in range(M):
            mat = out[i, j]
            for col, mask in enumerate(basis):
                if not ((mask >> j) & 1):
                    continue
                sign = -1 if (mask & ((1 << j) - 1)).bit_count() % 2 else 1
                after = mask ^ (1 << j)
                if (after >> i) & 1:
                    continue
                if (after & ((1 << i) - 1)).bit_count() % 2:
                    sign *= -1
                final = after | (1 << i)
                mat[index[final], col] = sign
    return out


def projected_annihilator_vector(M: int, R: int, orbital: int, A: float) -> np.ndarray:
    """Coefficients v_j in cbar_{R,orb}(A)=sum_j v_j d_j.

    This is the exact Fourier transform of u(k+A)=(cos(k+A),sin(k+A)).
    The down-spin second component differs by an overall minus sign, which
    drops out of its density operator.
    """
    v = np.zeros(M, dtype=np.complex128)
    plus = (R + 1) % M
    minus = (R - 1) % M
    if orbital == 0:
        v[plus] += 0.5 * np.exp(1j * A)
        v[minus] += 0.5 * np.exp(-1j * A)
    elif orbital == 1:
        v[plus] += np.exp(1j * A) / (2j)
        v[minus] -= np.exp(-1j * A) / (2j)
    else:
        raise ValueError("orbital must be 0 or 1")
    return v


def local_density(M: int, n: int, R: int, orbital: int, A: float) -> np.ndarray:
    v = projected_annihilator_vector(M, R, orbital, A)
    one_body = np.outer(np.conjugate(v), v)
    B = bilinears(M, n)
    mat = np.einsum("ij,ijab->ab", one_body, B, optimize=True)
    return (mat + mat.conj().T) / 2.0


@functools.lru_cache(None)
def identity_fixed_n(M: int, n: int) -> np.ndarray:
    return np.eye(len(masks(M, n)), dtype=np.complex128)


def hamiltonian(M: int, n: int, A: float, U: float = 1.0) -> np.ndarray:
    """Dense PSD projected-Hubbard Hamiltonian in N_up=N_down=n."""
    I = identity_fixed_n(M, n)
    dim = len(I)
    H = np.zeros((dim * dim, dim * dim), dtype=np.complex128)
    for R in range(M):
        for orbital in (0, 1):
            density = local_density(M, n, R, orbital, A)
            diff = np.kron(density, I) - np.kron(I, density)
            H += (U / 2.0) * (diff @ diff)
    return (H + H.conj().T) / 2.0


def low_spectrum(M: int, n: int, A: float, U: float = 1.0, count: int = 8) -> np.ndarray:
    H = hamiltonian(M, n, A, U=U)
    high = min(count - 1, H.shape[0] - 1)
    return eigh(H, eigvals_only=True, subset_by_index=[0, high], driver="evr")


def five_point_curvature(M: int, n: int, U: float = 1.0, h: float = 7.5e-4) -> tuple[float, float, int, float]:
    spec0 = low_spectrum(M, n, 0.0, U=U)
    e0 = float(spec0[0])
    deg = int(np.sum(np.abs(spec0 - e0) < 1e-9))
    gap = next((float(x - e0) for x in spec0[1:] if x - e0 > 1e-9), float("nan"))
    energies = {}
    for s in (-2, -1, 1, 2):
        energies[s] = float(low_spectrum(M, n, s * h, U=U, count=1)[0])
    c = (
        -energies[2]
        + 16.0 * energies[1]
        - 30.0 * e0
        + 16.0 * energies[-1]
        - energies[-2]
    ) / (12.0 * h * h)
    return e0, gap, deg, float(c)


def pair_energy_exact(Q: float, U: float = 1.0) -> float:
    return 0.5 * U * math.sin(Q) ** 2


def inverse_pair_mass_exact(U: float = 1.0) -> float:
    return U


def projector_up(k: np.ndarray) -> np.ndarray:
    """Gauge-invariant rank-one projector P_up(k)."""
    c = np.cos(k)
    s = np.sin(k)
    P = np.empty((len(k), 2, 2), dtype=np.float64)
    P[:, 0, 0] = c * c
    P[:, 0, 1] = c * s
    P[:, 1, 0] = c * s
    P[:, 1, 1] = s * s
    return P


def structural_diagnostics(M: int, A: float = 0.217) -> dict[str, float]:
    k = 2.0 * np.pi * np.arange(M) / M
    P0 = projector_up(k)
    PA = projector_up(k + A)
    idempotency = float(np.max(np.abs(np.einsum("kab,kbc->kac", P0, P0) - P0)))
    upc0 = np.mean(np.diagonal(P0, axis1=1, axis2=2), axis=0)
    upcA = np.mean(np.diagonal(PA, axis1=1, axis2=2), axis=0)
    return {
        "idempotency_error": idempotency,
        "upc0_orb0": float(upc0[0]),
        "upc0_orb1": float(upc0[1]),
        "upcA_orb0": float(upcA[0]),
        "upcA_orb1": float(upcA[1]),
        "period_pi_error": float(np.max(np.abs(projector_up(k + np.pi) - P0))),
    }


def component_polarized_residual(M: int, A: float, U: float = 1.0) -> float:
    """Maximum ||(nbar_up-nbar_down)|Phi_even>|| over PSD channels.

    It is enough to check each linear factor rather than forming the full
    C(M,M/2)^2 Hamiltonian.  For the product state with the even component
    full for both spins, the up and down density matrices are identical.
    """
    del U  # the null-vector check is independent of the positive coefficient
    if M % 2:
        raise ValueError("component-polarized state requires even M")
    n = M // 2
    basis = masks(M, n)
    full_even = sum(1 << j for j in range(0, M, 2))
    idx = basis.index(full_even)
    dim = len(basis)
    e = np.zeros(dim, dtype=np.complex128)
    e[idx] = 1.0
    worst = 0.0
    for R in range(M):
        for orbital in (0, 1):
            density = local_density(M, n, R, orbital, A)
            w = density @ e
            diff_state = np.kron(w, e) - np.kron(e, w)
            worst = max(worst, float(np.linalg.norm(diff_state)))
    return worst


@dataclass
class EDRow:
    M: int
    n: int
    nu: float
    connected_mod_M: bool
    E0: float
    gap_above_sampled_ground_space: float
    sampled_degeneracy: int
    electronic_twist_curvature: float
    rho_full: float
    predicted_from_one_pair: float
    ratio: float
    defect: float


def evaluate_ed_cases(cases: Iterable[tuple[int, int]], U: float, h: float) -> list[EDRow]:
    rows: list[EDRow] = []
    c1_by_M: dict[int, float] = {}
    raw: list[tuple[int, int, float, float, int, float]] = []
    for M, n in cases:
        e0, gap, deg, curvature = five_point_curvature(M, n, U=U, h=h)
        raw.append((M, n, e0, gap, deg, curvature))
        if n == 1:
            c1_by_M[M] = curvature
    for M, n, e0, gap, deg, curvature in raw:
        if M not in c1_by_M:
            _, _, _, c1 = five_point_curvature(M, 1, U=U, h=h)
            c1_by_M[M] = c1
        rho = n * (M - n) / (M - 1)
        predicted = rho * c1_by_M[M]
        rows.append(
            EDRow(
                M=M,
                n=n,
                nu=n / M,
                connected_mod_M=math.gcd(2, M) == 1,
                E0=e0,
                gap_above_sampled_ground_space=gap,
                sampled_degeneracy=deg,
                electronic_twist_curvature=curvature,
                rho_full=rho,
                predicted_from_one_pair=predicted,
                ratio=curvature / predicted if abs(predicted) > 1e-13 else float("nan"),
                defect=curvature - predicted,
            )
        )
    return rows


def print_rows(rows: list[EDRow]) -> None:
    print("\nExact-diagonalization checks (raw electronic-twist curvature):")
    print(
        f"{'M':>3} {'n':>3} {'nu':>7} {'connected':>10} {'E0':>13} {'gap':>12} "
        f"{'deg*':>5} {'Epp':>14} {'pred':>14} {'R':>12} {'defect':>14}"
    )
    for r in rows:
        print(
            f"{r.M:3d} {r.n:3d} {r.nu:7.4f} {str(r.connected_mod_M):>10} "
            f"{r.E0:13.6g} {r.gap_above_sampled_ground_space:12.6g} "
            f"{r.sampled_degeneracy:5d} {r.electronic_twist_curvature:14.9g} "
            f"{r.predicted_from_one_pair:14.9g} {r.ratio:12.9g} {r.defect:14.9g}"
        )
    print("*deg is the degeneracy visible in the requested low-spectrum window.")


def write_csv(rows: list[EDRow], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(asdict(rows[0]).keys()))
        writer.writeheader()
        for row in rows:
            writer.writerow(asdict(row))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--U", type=float, default=1.0)
    parser.add_argument("--step", type=float, default=7.5e-4)
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path(__file__).resolve().parent / "qgn_connectivity_counterexample_results.csv",
    )
    args = parser.parse_args()
    if args.U <= 0:
        raise SystemExit("U must be positive")

    diagnostics = structural_diagnostics(8)
    print("Structural diagnostics on M=8:")
    for key, value in diagnostics.items():
        print(f"  {key}: {value:.15g}")

    # M=5 is the connected control.  M=6 and M=8 are the reducible even tori.
    cases = [(5, 1), (5, 2), (6, 1), (6, 2), (6, 3), (8, 1), (8, 2)]
    rows = evaluate_ed_cases(cases, U=args.U, h=args.step)
    print_rows(rows)
    write_csv(rows, args.csv)

    print("\nExact one-pair statement for every M>=6:")
    print(f"  E_pair(Q) = (U/2) sin^2 Q; m_pair^(-1) = {inverse_pair_mass_exact(args.U):.12g}")
    print("At nu=1/2, N_flat=2:")
    print(f"  conjectured D_s = U/4 = {args.U/4:.12g}")
    print("  exact D_s = 0 on every even M, hence R_M(1/2)=0 and R_M-1=-1.")

    for M in (6, 8):
        for A in (0.0, 0.137, -0.241):
            residual = component_polarized_residual(M, A, U=args.U)
            print(f"  residual ||H(A)|Phi_even>|| for M={M}, A={A:+.3f}: {residual:.3e}")
            if residual > 2e-12:
                raise AssertionError((M, A, residual))

    # Headline regression assertions.
    m6 = {r.n: r for r in rows if r.M == 6}
    m8 = {r.n: r for r in rows if r.M == 8}
    m5 = {r.n: r for r in rows if r.M == 5}
    assert abs(m5[2].ratio - 1.0) < 2e-6
    assert abs(m6[2].ratio - 5.0 / 8.0) < 3e-5
    assert abs(m6[3].ratio) < 3e-5
    assert abs(m8[2].ratio - 7.0 / 9.0) < 3e-5
    assert abs(m6[1].electronic_twist_curvature - 4.0 * args.U) < 3e-5
    assert abs(m8[1].electronic_twist_curvature - 4.0 * args.U) < 3e-5

    print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()
