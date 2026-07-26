#!/usr/bin/env python3
"""Exact finite-size verifier for a covariantly twisted two-block null-row model.

The model is a one-dimensional periodic ring with two spinful pairing channels
(a=0,1), V cells, and four fermion modes per cell.  It has:

  * onsite seniority projectors q_{a,x}=(n_up-n_down)^2;
  * nearest-neighbor fermionic site-swap projectors Pi^-_{a,x};
  * a local cross-block row K_x=Pi^-_{12,x} B_x, where B_x is a crossed
    spin-preserving interblock hopping quadrature.

Every bond operator is Peierls twisted by conjugating with the charge-position
operator on the lifted right endpoint of that oriented bond.  This twists every
monomial in the full many-body row, including the swap projector.

At zero twist, the exact even-particle kernel is the product-Dicke/product-AGP
span, one state for every allowed block composition.  In the one-pair sector,
in the ordered composition basis (|1,0>, |0,1>), the exact curvature matrix is

    [[2 J1 + 4 g^2,  -4 g^2],
     [ -4 g^2,      2 J2 + 4 g^2]].

This script audits the kernel, evaluates the degenerate least-squares curvature,
and checks the curvature against direct finite differences.

Dependencies: numpy, scipy
"""

from __future__ import annotations

import argparse
import itertools
import math
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp


def fixed_weight_basis(number_modes: int, number_particles: int) -> list[int]:
    """Return bitstrings with the requested particle number."""
    if not 0 <= number_particles <= number_modes:
        return []
    return [
        state
        for state in range(1 << number_modes)
        if state.bit_count() == number_particles
    ]


def annihilate(state: int, mode: int) -> tuple[int, int] | None:
    if not ((state >> mode) & 1):
        return None
    parity = (state & ((1 << mode) - 1)).bit_count() & 1
    return state ^ (1 << mode), (-1 if parity else 1)


def create(state: int, mode: int) -> tuple[int, int] | None:
    if (state >> mode) & 1:
        return None
    parity = (state & ((1 << mode) - 1)).bit_count() & 1
    return state | (1 << mode), (-1 if parity else 1)


def mode_index(volume: int, block: int, site: int, spin: int) -> int:
    return ((block * volume + site) * 2 + spin)


def cdagger_c_matrix(basis: Sequence[int], dst: int, src: int) -> sp.csr_matrix:
    index = {state: row for row, state in enumerate(basis)}
    rows: list[int] = []
    cols: list[int] = []
    data: list[complex] = []

    for col, state in enumerate(basis):
        first = annihilate(state, src)
        if first is None:
            continue
        state_1, sign_1 = first
        second = create(state_1, dst)
        if second is None:
            continue
        state_2, sign_2 = second
        rows.append(index[state_2])
        cols.append(col)
        data.append(sign_1 * sign_2)

    dim = len(basis)
    return sp.csr_matrix((data, (rows, cols)), shape=(dim, dim), dtype=complex)


def induced_mode_permutation(
    basis: Sequence[int], permutation: Sequence[int]
) -> sp.csr_matrix:
    """Second-quantized unitary induced by old_mode -> new_mode."""
    index = {state: row for row, state in enumerate(basis)}
    rows: list[int] = []
    cols: list[int] = []
    data: list[complex] = []

    for col, state in enumerate(basis):
        occupied = [mode for mode in range(len(permutation)) if (state >> mode) & 1]
        mapped = [permutation[mode] for mode in occupied]
        inversions = sum(
            mapped[left] > mapped[right]
            for left in range(len(mapped))
            for right in range(left + 1, len(mapped))
        )
        sign = -1 if inversions & 1 else 1
        new_state = 0
        for mode in mapped:
            new_state |= 1 << mode
        rows.append(index[new_state])
        cols.append(col)
        data.append(sign)

    dim = len(basis)
    return sp.csr_matrix((data, (rows, cols)), shape=(dim, dim), dtype=complex)


@dataclass(frozen=True)
class ModelOperators:
    volume: int
    pairs: int
    basis: list[int]
    onsite_rows: list[list[sp.csr_matrix]]
    swap_projectors: list[list[sp.csr_matrix]]
    cross_rows: list[sp.csr_matrix]
    lifted_right_charge: list[np.ndarray]



def build_operators(volume: int, pairs: int) -> ModelOperators:
    if volume < 2:
        raise ValueError("volume must be at least 2")
    if not 0 <= pairs <= 2 * volume:
        raise ValueError("pairs must satisfy 0 <= pairs <= 2*volume")

    number_modes = 4 * volume
    basis = fixed_weight_basis(number_modes, 2 * pairs)
    dim = len(basis)
    identity = sp.eye(dim, format="csr", dtype=complex)

    onsite_rows: list[list[sp.csr_matrix]] = [[None] * volume for _ in range(2)]  # type: ignore[list-item]
    for block in range(2):
        for site in range(volume):
            up = mode_index(volume, block, site, 0)
            down = mode_index(volume, block, site, 1)
            diagonal = np.asarray(
                [
                    (((state >> up) & 1) - ((state >> down) & 1)) ** 2
                    for state in basis
                ],
                dtype=float,
            )
            onsite_rows[block][site] = sp.diags(diagonal, format="csr", dtype=complex)

    swap_projectors: list[list[sp.csr_matrix]] = [[None] * volume for _ in range(2)]  # type: ignore[list-item]
    cross_rows: list[sp.csr_matrix] = []
    lifted_right_charge: list[np.ndarray] = []

    for left in range(volume):
        right = (left + 1) % volume
        swaps: list[sp.csr_matrix] = []

        for block in range(2):
            permutation = list(range(number_modes))
            for spin in range(2):
                mode_left = mode_index(volume, block, left, spin)
                mode_right = mode_index(volume, block, right, spin)
                permutation[mode_left], permutation[mode_right] = (
                    permutation[mode_right],
                    permutation[mode_left],
                )
            swap = induced_mode_permutation(basis, permutation)
            swaps.append(swap)
            swap_projectors[block][left] = 0.5 * (identity - swap)

        simultaneous_swap = swaps[0] @ swaps[1]
        simultaneous_antisymmetrizer = 0.5 * (identity - simultaneous_swap)

        directed = sp.csr_matrix((dim, dim), dtype=complex)
        for spin in range(2):
            directed += cdagger_c_matrix(
                basis,
                mode_index(volume, 0, right, spin),
                mode_index(volume, 1, left, spin),
            )
            directed += cdagger_c_matrix(
                basis,
                mode_index(volume, 0, left, spin),
                mode_index(volume, 1, right, spin),
            )
        bridge_quadrature = directed + directed.getH()

        commutator = simultaneous_antisymmetrizer @ bridge_quadrature - bridge_quadrature @ simultaneous_antisymmetrizer
        if sp.linalg.norm(commutator) > 1.0e-12:
            raise RuntimeError("local bridge failed the simultaneous-swap commutation check")

        cross_rows.append(simultaneous_antisymmetrizer @ bridge_quadrature)

        # Covering-lattice lift: the physical right endpoint has coordinate +1
        # relative to the left endpoint, including the periodic boundary bond.
        right_charge = np.asarray(
            [
                sum(
                    (state >> mode_index(volume, block, right, spin)) & 1
                    for block in range(2)
                    for spin in range(2)
                )
                for state in basis
            ],
            dtype=float,
        )
        lifted_right_charge.append(right_charge)

    return ModelOperators(
        volume=volume,
        pairs=pairs,
        basis=basis,
        onsite_rows=onsite_rows,
        swap_projectors=swap_projectors,
        cross_rows=cross_rows,
        lifted_right_charge=lifted_right_charge,
    )


def gauge_conjugate(operator: sp.csr_matrix, charge: np.ndarray, vector_potential: float) -> sp.csr_matrix:
    phase = np.exp(-1j * vector_potential * charge)
    return sp.diags(phase, format="csr") @ operator @ sp.diags(np.conjugate(phase), format="csr")


def factor_rows(
    operators: ModelOperators,
    vector_potential: float,
    onsite_strength: float,
    swap_strengths: tuple[float, float],
    cross_strength: float,
) -> list[sp.csr_matrix]:
    rows: list[sp.csr_matrix] = []
    volume = operators.volume

    for block in range(2):
        for site in range(volume):
            rows.append(math.sqrt(onsite_strength) * operators.onsite_rows[block][site])

    for block in range(2):
        for left in range(volume):
            right = (left + 1) % volume
            block_right_charge = np.asarray(
                [
                    sum(
                        (state >> mode_index(volume, block, right, spin)) & 1
                        for spin in range(2)
                    )
                    for state in operators.basis
                ],
                dtype=float,
            )
            twisted = gauge_conjugate(
                operators.swap_projectors[block][left],
                block_right_charge,
                vector_potential,
            )
            rows.append(math.sqrt(swap_strengths[block]) * twisted)

    for left in range(volume):
        twisted = gauge_conjugate(
            operators.cross_rows[left],
            operators.lifted_right_charge[left],
            vector_potential,
        )
        rows.append(cross_strength * twisted)

    return rows


def factor_rows_and_derivatives(
    operators: ModelOperators,
    onsite_strength: float,
    swap_strengths: tuple[float, float],
    cross_strength: float,
) -> tuple[list[sp.csr_matrix], list[sp.csr_matrix]]:
    rows_zero: list[sp.csr_matrix] = []
    rows_prime: list[sp.csr_matrix] = []
    volume = operators.volume
    dim = len(operators.basis)
    zero = sp.csr_matrix((dim, dim), dtype=complex)

    for block in range(2):
        for site in range(volume):
            row = math.sqrt(onsite_strength) * operators.onsite_rows[block][site]
            rows_zero.append(row)
            rows_prime.append(zero)

    for block in range(2):
        for left in range(volume):
            right = (left + 1) % volume
            charge = np.asarray(
                [
                    sum(
                        (state >> mode_index(volume, block, right, spin)) & 1
                        for spin in range(2)
                    )
                    for state in operators.basis
                ],
                dtype=float,
            )
            generator = sp.diags(charge, format="csr", dtype=complex)
            row = math.sqrt(swap_strengths[block]) * operators.swap_projectors[block][left]
            derivative = -1j * (generator @ row - row @ generator)
            rows_zero.append(row)
            rows_prime.append(derivative)

    for left in range(volume):
        generator = sp.diags(
            operators.lifted_right_charge[left], format="csr", dtype=complex
        )
        row = cross_strength * operators.cross_rows[left]
        derivative = -1j * (generator @ row - row @ generator)
        rows_zero.append(row)
        rows_prime.append(derivative)

    return rows_zero, rows_prime


def hamiltonian(
    operators: ModelOperators,
    vector_potential: float,
    onsite_strength: float,
    swap_strengths: tuple[float, float],
    cross_strength: float,
) -> sp.csr_matrix:
    dim = len(operators.basis)
    result = sp.csr_matrix((dim, dim), dtype=complex)
    for row in factor_rows(
        operators,
        vector_potential,
        onsite_strength,
        swap_strengths,
        cross_strength,
    ):
        result += 0.5 * (row.getH() @ row)
    return result


def product_agp_basis(operators: ModelOperators) -> tuple[np.ndarray, list[tuple[int, int]]]:
    volume = operators.volume
    pairs = operators.pairs
    basis_index = {state: row for row, state in enumerate(operators.basis)}
    vectors: list[np.ndarray] = []
    compositions: list[tuple[int, int]] = []

    lower = max(0, pairs - volume)
    upper = min(volume, pairs)
    for pairs_block_0 in range(upper, lower - 1, -1):
        pairs_block_1 = pairs - pairs_block_0
        vector = np.zeros(len(operators.basis), dtype=complex)

        for occupied_0 in itertools.combinations(range(volume), pairs_block_0):
            for occupied_1 in itertools.combinations(range(volume), pairs_block_1):
                state = 0
                for site in occupied_0:
                    state |= 1 << mode_index(volume, 0, site, 0)
                    state |= 1 << mode_index(volume, 0, site, 1)
                for site in occupied_1:
                    state |= 1 << mode_index(volume, 1, site, 0)
                    state |= 1 << mode_index(volume, 1, site, 1)
                vector[basis_index[state]] += 1.0

        norm = la.norm(vector)
        if norm == 0:
            raise RuntimeError("constructed a zero product-AGP vector")
        vectors.append(vector / norm)
        compositions.append((pairs_block_0, pairs_block_1))

    return np.column_stack(vectors), compositions


def least_squares_curvature(
    operators: ModelOperators,
    onsite_strength: float,
    swap_strengths: tuple[float, float],
    cross_strength: float,
) -> tuple[np.ndarray, list[tuple[int, int]], float]:
    zero_basis, compositions = product_agp_basis(operators)
    rows_zero, rows_prime = factor_rows_and_derivatives(
        operators,
        onsite_strength,
        swap_strengths,
        cross_strength,
    )

    null_map = np.vstack([row.toarray() for row in rows_zero])
    source = np.vstack([row @ zero_basis for row in rows_prime])

    # Minimize ||D chi + source||.  The residual is the projection onto ker D^dagger.
    correction, *_ = la.lstsq(
        null_map,
        -source,
        cond=1.0e-11,
        lapack_driver="gelsy",
    )
    residual = null_map @ correction + source
    curvature = residual.conj().T @ residual
    curvature = 0.5 * (curvature + curvature.conj().T)
    return curvature, compositions, float(la.norm(residual))


def expected_kernel_dimension(volume: int, pairs: int) -> int:
    return min(volume, pairs) - max(0, pairs - volume) + 1


def count_zero_eigenvalues(matrix: sp.csr_matrix, tolerance: float = 1.0e-9) -> tuple[int, np.ndarray]:
    eigenvalues = la.eigvalsh(matrix.toarray())
    return int(np.count_nonzero(np.abs(eigenvalues) < tolerance)), eigenvalues


def audit(
    onsite_strength: float,
    swap_strengths: tuple[float, float],
    cross_strength: float,
    tolerance: float,
) -> None:
    print("two-block covariant many-body null-row verifier")
    print("composition basis order: descending block-1 occupation; one-pair basis = (|1,0>, |0,1>)")
    print()
    print("=== Exact zero-manifold audit ===")
    audit_cases = [(2, pairs) for pairs in range(5)] + [(3, pairs) for pairs in range(7)]
    for volume, pairs in audit_cases:
        operators = build_operators(volume, pairs)
        zero_basis, compositions = product_agp_basis(operators)
        h_zero = hamiltonian(
            operators,
            0.0,
            onsite_strength,
            swap_strengths,
            cross_strength,
        )
        kernel_dimension, eigenvalues = count_zero_eigenvalues(h_zero, tolerance)
        expected = expected_kernel_dimension(volume, pairs)
        residual = la.norm(h_zero @ zero_basis)
        gap = eigenvalues[expected] if expected < len(eigenvalues) else math.inf
        status = "PASS" if kernel_dimension == expected and residual < tolerance else "FAIL"
        print(
            f"{status}  V={volume:2d}, n={pairs:2d}: "
            f"kernel={kernel_dimension}, expected={expected}, "
            f"AGP residual={residual:.3e}, next gap={gap:.6g}, "
            f"compositions={compositions}"
        )
        if status == "FAIL":
            raise AssertionError("zero-manifold audit failed")

    print("\n=== One-pair composition curvature ===")
    j_1, j_2 = swap_strengths
    expected_matrix = np.asarray(
        [
            [2.0 * j_1 + 4.0 * cross_strength**2, -4.0 * cross_strength**2],
            [-4.0 * cross_strength**2, 2.0 * j_2 + 4.0 * cross_strength**2],
        ],
        dtype=float,
    )

    for volume in range(2, 6):
        operators = build_operators(volume, 1)
        curvature, compositions, _ = least_squares_curvature(
            operators,
            onsite_strength,
            swap_strengths,
            cross_strength,
        )
        error = float(np.max(np.abs(curvature - expected_matrix)))
        status = "PASS" if error < tolerance else "FAIL"
        print(
            f"{status}  V={volume:2d}, basis={compositions}, "
            f"max formula error={error:.3e}"
        )
        print(np.real_if_close(curvature))
        if status == "FAIL":
            raise AssertionError("one-pair curvature formula failed")

    expected_eigenvalues = np.asarray(
        [
            j_1 + j_2 + 4.0 * cross_strength**2
            - math.sqrt((j_1 - j_2) ** 2 + 16.0 * cross_strength**4),
            j_1 + j_2 + 4.0 * cross_strength**2
            + math.sqrt((j_1 - j_2) ** 2 + 16.0 * cross_strength**4),
        ]
    )
    print(f"Exact one-pair curvature eigenvalues: {expected_eigenvalues}")

    print("\n=== Higher-composition spot check ===")
    operators_n2 = build_operators(3, 2)
    curvature_n2, compositions_n2, _ = least_squares_curvature(
        operators_n2,
        onsite_strength,
        swap_strengths,
        cross_strength,
    )
    off_diagonal = curvature_n2 - np.diag(np.diag(curvature_n2))
    off_diagonal_norm = float(la.norm(off_diagonal))
    status = "PASS" if off_diagonal_norm > tolerance else "FAIL"
    print(f"{status}  V=3, n=2, basis={compositions_n2}")
    print(np.real_if_close(curvature_n2))
    print(f"      off-diagonal Frobenius norm = {off_diagonal_norm:.6g}")
    if status == "FAIL":
        raise AssertionError("higher-composition response is diagonal")

    print("\n=== Direct finite-difference branch check ===")
    operators = build_operators(3, 1)
    curvature, _, _ = least_squares_curvature(
        operators,
        onsite_strength,
        swap_strengths,
        cross_strength,
    )
    target = la.eigvalsh(curvature)
    step = 1.0e-3
    h_step = hamiltonian(
        operators,
        step,
        onsite_strength,
        swap_strengths,
        cross_strength,
    )
    low_energies = la.eigvalsh(h_step.toarray(), subset_by_index=[0, 1])
    finite_difference = 2.0 * low_energies / step**2
    error = float(np.max(np.abs(finite_difference - target)))
    status = "PASS" if error < 5.0e-5 else "FAIL"
    print(f"{status}  least-squares eigenvalues = {target}")
    print(f"      finite-difference values   = {finite_difference}")
    print(f"      maximum discrepancy        = {error:.3e}")
    if status == "FAIL":
        raise AssertionError("finite-difference curvature check failed")

    print("\nAll retained checks passed.")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--onsite", type=float, default=2.0, help="onsite seniority strength U")
    parser.add_argument("--j1", type=float, default=1.0, help="block-1 swap strength")
    parser.add_argument("--j2", type=float, default=1.3, help="block-2 swap strength")
    parser.add_argument("--g", type=float, default=0.7, help="cross null-row strength")
    parser.add_argument("--tolerance", type=float, default=1.0e-8)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.onsite <= 0 or args.j1 <= 0 or args.j2 <= 0:
        raise ValueError("onsite, j1, and j2 must be positive")
    if args.g == 0:
        print("Warning: g=0 removes composition hopping by construction.")
    audit(
        onsite_strength=args.onsite,
        swap_strengths=(args.j1, args.j2),
        cross_strength=args.g,
        tolerance=args.tolerance,
    )


if __name__ == "__main__":
    main()
