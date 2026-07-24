#!/usr/bin/env python3
"""Independent Phase-I/Phase-II verifier for the resonant winding example.

This code deliberately uses a different assembly route from the exact SymPy
certificate: it constructs the projected local densities in a fixed-number bit
basis and forms

    H_cross(A) = -U sum_{R,a} nbar_{R a up}(A) nbar_{R a down}(A)

or, optionally, the positive-semidefinite QGN/Hubbard convention

    H_psd(A) = (U/2) sum_{R,a} [nbar_{R a up}(A)-nbar_{R a down}(A)]^2.

The two-orbital flat-band spinor is

    u_up(k) = [cos(pi/4 + c sin(r k)),
               exp(i b sin k) sin(pi/4 + c sin(r k))]^T,

with u_down(k)=u_up(-k)^*.  The twist is inserted before projection by
k -> k+A.

Default checks:
  1. M=4, r=4: exact period-M resonance and filling-law violation.
  2. M=8, r=4: the same hidden harmonic becomes alternating under twist;
     scalar factorization fails and the filling law is restored numerically.
  3. M=8, r=8: period-M resonance returns, with the analytic defect predicted
     by the general scalar-deformation formula.

The script is intended for small-system verification, not large-scale ED.
"""

from __future__ import annotations

import argparse
import csv
import itertools
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.linalg import eigh


@dataclass(frozen=True)
class Case:
    label: str
    M: int
    harmonic: int
    c: float
    pair_numbers: tuple[int, ...]


@lru_cache(maxsize=None)
def masks(M: int, n: int) -> tuple[int, ...]:
    return tuple(
        sum(1 << i for i in combination)
        for combination in itertools.combinations(range(M), n)
    )


@lru_cache(maxsize=None)
def bilinears(M: int, n: int) -> tuple[np.ndarray, ...]:
    """Dense matrices of d_p^dagger d_q in the fixed-n spin sector."""
    basis = masks(M, n)
    index = {mask: i for i, mask in enumerate(basis)}
    dimension = len(basis)
    result: list[np.ndarray] = []

    for p in range(M):
        for q in range(M):
            matrix = np.zeros((dimension, dimension), dtype=np.complex128)
            for column, mask in enumerate(basis):
                if not ((mask >> q) & 1):
                    continue

                sign = -1 if (mask & ((1 << q) - 1)).bit_count() % 2 else 1
                after_annihilation = mask ^ (1 << q)
                if (after_annihilation >> p) & 1:
                    continue

                if (after_annihilation & ((1 << p) - 1)).bit_count() % 2:
                    sign *= -1
                final_mask = after_annihilation | (1 << p)
                matrix[index[final_mask], column] = sign
            result.append(matrix)

    return tuple(result)


def band_vectors(
    M: int,
    twist: float,
    c: float,
    harmonic: int,
    b: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    momenta = 2.0 * np.pi * np.arange(M) / M
    shifted = momenta + twist

    theta_up = np.pi / 4.0 + c * np.sin(harmonic * shifted)
    up = np.column_stack(
        (
            np.cos(theta_up),
            np.exp(1j * b * np.sin(shifted)) * np.sin(theta_up),
        )
    )

    # u_down(k+A)=u_up(-k-A)^*.  Because sin is odd, the phase shown below
    # has the same sign as the spin-up phase, while the mixing-angle
    # deformation changes sign.
    theta_down = np.pi / 4.0 - c * np.sin(harmonic * shifted)
    down = np.column_stack(
        (
            np.cos(theta_down),
            np.exp(1j * b * np.sin(shifted)) * np.sin(theta_down),
        )
    )
    return momenta, up, down


def local_projected_density(
    M: int,
    n: int,
    momenta: np.ndarray,
    spinors: np.ndarray,
    cell: int,
    orbital: int,
) -> np.ndarray:
    """Matrix of a projected local density nbar_{R,orbital}."""
    one_body = bilinears(M, n)
    dimension = len(masks(M, n))
    density = np.zeros((dimension, dimension), dtype=np.complex128)
    phase = np.exp(1j * momenta * cell)

    for p in range(M):
        for q in range(M):
            coefficient = (
                np.conjugate(spinors[p, orbital])
                * spinors[q, orbital]
                * np.conjugate(phase[p])
                * phase[q]
                / M
            )
            density += coefficient * one_body[p * M + q]
    return 0.5 * (density + density.conjugate().T)


def projected_hamiltonian(
    M: int,
    n: int,
    twist: float,
    c: float,
    harmonic: int,
    b: float = np.pi / 2.0,
    U: float = 1.0,
    convention: str = "cross",
) -> np.ndarray:
    momenta, up, down = band_vectors(M, twist, c, harmonic, b)
    dimension = len(masks(M, n))
    identity = np.eye(dimension, dtype=np.complex128)
    hamiltonian = np.zeros((dimension * dimension, dimension * dimension), dtype=np.complex128)

    for cell in range(M):
        for orbital in range(2):
            n_up = local_projected_density(M, n, momenta, up, cell, orbital)
            n_down = local_projected_density(M, n, momenta, down, cell, orbital)

            if convention == "cross":
                hamiltonian -= U * np.kron(n_up, n_down)
            elif convention == "psd":
                difference = np.kron(n_up, identity) - np.kron(identity, n_down)
                hamiltonian += (U / 2.0) * (difference @ difference)
            else:
                raise ValueError(f"Unknown Hamiltonian convention: {convention}")

    return 0.5 * (hamiltonian + hamiltonian.conjugate().T)


def lowest_two(matrix: np.ndarray) -> tuple[float, float]:
    values = eigh(matrix, eigvals_only=True, subset_by_index=[0, 1])
    return float(values[0]), float(values[1])


def energy_curvature(
    M: int,
    n: int,
    c: float,
    harmonic: int,
    b: float,
    U: float,
    convention: str,
    step: float,
) -> tuple[float, float, float]:
    energies: list[float] = []
    for twist in (-2.0 * step, -step, 0.0, step, 2.0 * step):
        hamiltonian = projected_hamiltonian(
            M=M,
            n=n,
            twist=twist,
            c=c,
            harmonic=harmonic,
            b=b,
            U=U,
            convention=convention,
        )
        ground, first_excited = lowest_two(hamiltonian)
        energies.append(ground)
        if twist == 0.0:
            zero_energy = ground
            gap = first_excited - ground

    em2, em1, e0, ep1, ep2 = energies
    curvature = (-ep2 + 16.0 * ep1 - 30.0 * e0 + 16.0 * em1 - em2) / (
        12.0 * step * step
    )
    return zero_energy, gap, float(curvature)


def filling_multiplier(M: int, n: int) -> float:
    """Multiplier rho_n with E_n'' = rho_n E_1'' if the filling law holds."""
    return n * (M - n) / (M - 1)


def factorization_diagnostic(
    M: int,
    n: int,
    c: float,
    harmonic: int,
    twist: float,
    b: float,
    U: float,
    convention: str,
) -> tuple[complex, float, float]:
    h_c = projected_hamiltonian(M, n, twist, c, harmonic, b, U, convention)
    h_0 = projected_hamiltonian(M, n, twist, 0.0, harmonic, b, U, convention)

    if convention == "cross" and harmonic % M == 0:
        scalar: complex = np.cos(2.0 * c * np.sin(harmonic * twist)) ** 2
    else:
        scalar = np.vdot(h_0, h_c) / np.vdot(h_0, h_0)

    difference = h_c - scalar * h_0
    max_element = float(np.max(np.abs(difference)))
    relative_frobenius = float(np.linalg.norm(difference) / np.linalg.norm(h_c))
    return scalar, max_element, relative_frobenius


def resonant_raw_defect(M: int, n: int, c: float, harmonic: int, U: float) -> float:
    """Extra raw-curvature defect from a scalar period-M resonance.

    This is the defect relative to the c=0 model:
        Delta_n(c)-Delta_n(0)
        = 4 U c^2 r^2 n(n-1)/(M-1)
    for the cross-attraction convention and r divisible by M.
    """
    if harmonic % M != 0:
        raise ValueError("The analytic resonance formula requires harmonic divisible by M")
    return 4.0 * U * c * c * harmonic * harmonic * n * (n - 1) / (M - 1)


def run_case(
    case: Case,
    b: float,
    U: float,
    convention: str,
    step: float,
) -> list[dict[str, float | int | str]]:
    print(f"\n[{case.label}] M={case.M}, r={case.harmonic}, c={case.c:g}, convention={convention}")
    rows: list[dict[str, float | int | str]] = []
    one_pair_curvature: float | None = None

    for n in case.pair_numbers:
        energy, gap, curvature = energy_curvature(
            case.M, n, case.c, case.harmonic, b, U, convention, step
        )
        if n == 1:
            one_pair_curvature = curvature
        assert one_pair_curvature is not None

        rho = filling_multiplier(case.M, n)
        prediction = rho * one_pair_curvature
        defect = curvature - prediction
        reduced_ratio = curvature / rho

        print(
            f"  n={n}: E(0)={energy:+.12f}, gap={gap:.12f}, "
            f"E''={curvature:.12f}, rho_n E_1''={prediction:.12f}, "
            f"defect={defect:+.3e}, E''/rho_n={reduced_ratio:.12f}"
        )
        rows.append(
            {
                "case": case.label,
                "M": case.M,
                "harmonic": case.harmonic,
                "c": case.c,
                "n": n,
                "energy_zero": energy,
                "gap": gap,
                "curvature": curvature,
                "rho_n": rho,
                "predicted_curvature": prediction,
                "defect": defect,
                "reduced_ratio": reduced_ratio,
            }
        )

    scalar, max_element, relative = factorization_diagnostic(
        case.M,
        min(2, max(case.pair_numbers)),
        case.c,
        case.harmonic,
        twist=0.137,
        b=b,
        U=U,
        convention=convention,
    )
    print(
        "  factorization diagnostic at A=0.137: "
        f"best/analytic scalar={scalar}, max|difference|={max_element:.3e}, "
        f"relative Frobenius={relative:.3e}"
    )

    if convention == "cross" and case.harmonic % case.M == 0:
        for n in case.pair_numbers:
            if n >= 2:
                analytic = resonant_raw_defect(case.M, n, case.c, case.harmonic, U)
                numerical = float(rows[n - 1]["defect"])
                print(
                    f"  resonance formula, n={n}: predicted defect={analytic:.12f}, "
                    f"numerical={numerical:.12f}"
                )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", type=float, default=7.0e-4, help="finite-difference step")
    parser.add_argument("--b", type=float, default=float(np.pi / 2.0), help="phase amplitude")
    parser.add_argument("--U", type=float, default=1.0, help="attraction strength")
    parser.add_argument(
        "--convention",
        choices=("cross", "psd"),
        default="cross",
        help="projected interaction convention",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("phase12_results.csv"),
        help="output CSV path",
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="also evaluate n=3 in the M=8 cases (substantially slower)",
    )
    args = parser.parse_args()

    m8_pairs = (1, 2, 3) if args.full else (1, 2)
    cases = (
        Case("M4_r4_resonant", 4, 4, 1.0 / 8.0, (1, 2, 3)),
        Case("M8_r4_nonresonant_control", 8, 4, 1.0 / 8.0, m8_pairs),
        Case("M8_r8_resonant", 8, 8, 1.0 / 8.0, m8_pairs),
    )

    all_rows: list[dict[str, float | int | str]] = []
    for case in cases:
        all_rows.extend(run_case(case, args.b, args.U, args.convention, args.step))

    args.csv.parent.mkdir(parents=True, exist_ok=True)
    with args.csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(all_rows[0].keys()))
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nWrote {args.csv}")


if __name__ == "__main__":
    main()
