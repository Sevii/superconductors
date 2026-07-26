#!/usr/bin/env python3
"""LEGACY FIXED-KERNEL VERIFIER -- retained for audit only.

This program holds the multi-site bridge kernel fixed while twisting the active
projector.  It therefore reproduces the algebraic source used by the superseded
draft, but it is not a physical Peierls-covariance check.  Use
``check_bridge_covariance.py`` for the corrected comparison: when the bridge
kernel is shifted as B(k+A) together with P(A), the projected zero row and its
derivative cancel.

The Fock-space routines below remain useful for checking the conditional Jacobi
coefficients once a legitimate microscopic transfer source has been supplied.
"""

from __future__ import annotations

import argparse
import itertools
import math
from collections import defaultdict
from typing import Dict, Iterable, Tuple

import numpy as np

Array = np.ndarray
Vector = Dict[int, complex]


def haar_filter(z: complex) -> Array:
    """2 x 2 Laurent-polynomial paraunitary Haar filter."""
    return 0.5 * np.array([[1.0 + z, 1.0 - z],
                           [1.0 - z, 1.0 + z]], dtype=complex)


D_PHASE = np.diag([1.0, 1.0j])


def frame(kx: float, ky: float, block: int) -> Array:
    """Columns are (active, bright); block 2 is the C4-rotated filter."""
    hx = haar_filter(np.exp(1.0j * kx))
    hy = haar_filter(np.exp(1.0j * ky))
    if block == 1:
        return hx @ D_PHASE @ hy
    if block == 2:
        return hy @ D_PHASE @ hx
    raise ValueError("block must be 1 or 2")


def active(kx: float, ky: float, block: int) -> Array:
    return frame(kx, ky, block)[:, 0]


def bright(kx: float, ky: float, block: int) -> Array:
    return frame(kx, ky, block)[:, 1]


def momentum_grid(nx: int, ny: int) -> list[tuple[float, float]]:
    return [
        (2.0 * np.pi * ix / nx, 2.0 * np.pi * iy / ny)
        for ix in range(nx)
        for iy in range(ny)
    ]


def numerical_metric(block: int, grid: int = 100, step: float = 1.0e-5) -> Array:
    total = np.zeros((2, 2), dtype=float)
    for ix in range(grid):
        kx = 2.0 * np.pi * ix / grid
        for iy in range(grid):
            ky = 2.0 * np.pi * iy / grid
            u = active(kx, ky, block)
            derivatives = []
            for dx, dy in ((step, 0.0), (0.0, step)):
                up = active(kx + dx, ky + dy, block)
                um = active(kx - dx, ky - dy, block)
                derivatives.append((up - um) / (2.0 * step))
            for i in range(2):
                for j in range(2):
                    value = (
                        np.vdot(derivatives[i], derivatives[j])
                        - np.vdot(derivatives[i], u) * np.vdot(u, derivatives[j])
                    )
                    total[i, j] += float(np.real(value))
    return total / (grid * grid)


EXPECTED_METRIC = {
    1: np.diag([1.0 / 8.0, 1.0 / 4.0]),
    2: np.diag([1.0 / 4.0, 1.0 / 8.0]),
}


def projected_orbital(
    nx: int,
    ny: int,
    cell: tuple[int, int],
    microscopic_orbital: int,
    vector_potential: tuple[float, float],
    spin: str,
    block: int,
) -> Array:
    """Rank-one projected microscopic orbital on the abstract active band."""
    ks = momentum_grid(nx, ny)
    volume = len(ks)
    vector = np.zeros(volume, dtype=complex)
    ax, ay = vector_potential
    rx, ry = cell
    for index, (kx, ky) in enumerate(ks):
        if spin == "up":
            u = active(kx + ax, ky + ay, block)
        elif spin == "down":
            # Time-reversed sector: u_down(k; A) = u_up^*(-k-A).
            u = np.conjugate(active(-kx - ax, -ky - ay, block))
        else:
            raise ValueError("spin must be 'up' or 'down'")
        vector[index] = (
            np.exp(-1.0j * (kx * rx + ky * ry))
            * np.conjugate(u[microscopic_orbital])
            / math.sqrt(volume)
        )
    return np.outer(vector, np.conjugate(vector))


def one_pair_hamiltonian(
    nx: int,
    ny: int,
    vector_potential: tuple[float, float],
    block: int,
    coupling: float = 1.0,
) -> Array:
    """One-up/one-down projected density-difference positive-square model."""
    volume = nx * ny
    identity = np.eye(volume)
    hamiltonian = np.zeros((volume * volume, volume * volume), dtype=complex)
    for rx in range(nx):
        for ry in range(ny):
            for mu in range(2):
                p_up = projected_orbital(
                    nx, ny, (rx, ry), mu, vector_potential, "up", block
                )
                p_down = projected_orbital(
                    nx, ny, (rx, ry), mu, vector_potential, "down", block
                )
                difference = np.kron(p_up, identity) - np.kron(identity, p_down)
                hamiltonian += 0.5 * coupling * (difference @ difference)
    return hamiltonian


def one_pair_mass(
    nx: int,
    ny: int,
    block: int,
    direction: int,
    step: float = 1.0e-3,
) -> float:
    """Return d^2 E/dQ_i^2; the electronic twist obeys Q=2A."""
    def ground_energy(a: float) -> float:
        vector_potential = [0.0, 0.0]
        vector_potential[direction] = a
        eigenvalues = np.linalg.eigvalsh(
            one_pair_hamiltonian(nx, ny, tuple(vector_potential), block)
        )
        return float(eigenvalues[0])

    e_zero = ground_energy(0.0)
    e_plus = ground_energy(step)
    e_minus = ground_energy(-step)
    curvature_a = (e_plus + e_minus - 2.0 * e_zero) / (step * step)
    return curvature_a / 4.0


def bridge_active_matrix(
    kx: float,
    ky: float,
    ay: float,
    quadrature: str,
) -> Array:
    """Project the fixed local bridge into the two twisted active bands."""
    u1_a = active(kx, ky + ay, 1)
    u2_a = active(kx, ky + ay, 2)
    u1_0 = active(kx, ky, 1)
    u2_0 = active(kx, ky, 2)
    r1_0 = bright(kx, ky, 1)

    # M = |bright_1(0)><active_2(0)| in microscopic orbital space.
    m12 = np.outer(r1_0, np.conjugate(u2_0))
    coefficient_12 = np.vdot(u1_a, m12 @ u2_a)
    coefficient_21 = np.conjugate(coefficient_12)

    if quadrature == "X":
        return np.array([[0.0, coefficient_12],
                         [coefficient_21, 0.0]], dtype=complex)
    if quadrature == "Y":
        return np.array([[0.0, -1.0j * coefficient_12],
                         [1.0j * coefficient_21, 0.0]], dtype=complex)
    raise ValueError("quadrature must be X or Y")


def bridge_derivative_error(step: float = 1.0e-6, grid: int = 17) -> tuple[float, float]:
    expected_x = 0.5 * np.array([[0.0, 1.0j], [-1.0j, 0.0]], dtype=complex)
    expected_y = 0.5 * np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex)
    max_error_x = 0.0
    max_error_y = 0.0
    for ix in range(grid):
        kx = 2.0 * np.pi * ix / grid
        for iy in range(grid):
            ky = 2.0 * np.pi * iy / grid
            derivative_x = (
                bridge_active_matrix(kx, ky, step, "X")
                - bridge_active_matrix(kx, ky, -step, "X")
            ) / (2.0 * step)
            derivative_y = (
                bridge_active_matrix(kx, ky, step, "Y")
                - bridge_active_matrix(kx, ky, -step, "Y")
            ) / (2.0 * step)
            max_error_x = max(max_error_x, float(np.max(np.abs(derivative_x - expected_x))))
            max_error_y = max(max_error_y, float(np.max(np.abs(derivative_y - expected_y))))
    return max_error_x, max_error_y


# --- Direct spinful Fock-space verification of the local Jacobi source. ---

def _create(state: int, orbital: int) -> Tuple[int, int] | None:
    if (state >> orbital) & 1:
        return None
    sign = -1 if (state & ((1 << orbital) - 1)).bit_count() % 2 else 1
    return state | (1 << orbital), sign


def _annihilate(state: int, orbital: int) -> Tuple[int, int] | None:
    if not ((state >> orbital) & 1):
        return None
    sign = -1 if (state & ((1 << orbital) - 1)).bit_count() % 2 else 1
    return state & ~(1 << orbital), sign


def _cdag_c(vector: Vector, create_orbital: int, annihilate_orbital: int) -> Vector:
    output: defaultdict[int, complex] = defaultdict(complex)
    for state, amplitude in vector.items():
        first = _annihilate(state, annihilate_orbital)
        if first is None:
            continue
        state_1, sign_1 = first
        second = _create(state_1, create_orbital)
        if second is None:
            continue
        state_2, sign_2 = second
        output[state_2] += amplitude * sign_1 * sign_2
    return dict(output)


def _inner(left: Vector, right: Vector) -> complex:
    return sum(
        np.conjugate(amplitude) * right.get(state, 0.0)
        for state, amplitude in left.items()
    )


def _normalize(vector: Vector) -> Vector:
    norm = math.sqrt(float(np.real(_inner(vector, vector))))
    if norm == 0.0:
        raise ValueError("Cannot normalize a zero vector")
    return {state: amplitude / norm for state, amplitude in vector.items()}


def _block_agp(capacity: int, block: int, pair_number: int) -> Vector:
    output: defaultdict[int, complex] = defaultdict(complex)
    for occupied_levels in itertools.combinations(range(capacity), pair_number):
        state = 0
        amplitude = 1
        for level in occupied_levels:
            up = ((block * capacity + level) * 2)
            down = up + 1
            result_down = _create(state, down)
            assert result_down is not None
            state, sign_down = result_down
            result_up = _create(state, up)
            assert result_up is not None
            state, sign_up = result_up
            amplitude *= sign_down * sign_up
        output[state] += amplitude
    return _normalize(dict(output))


def product_agp(capacity: int, n1: int, n2: int) -> Vector:
    left = _block_agp(capacity, 0, n1)
    right = _block_agp(capacity, 1, n2)
    output: defaultdict[int, complex] = defaultdict(complex)
    for state_left, amplitude_left in left.items():
        for state_right, amplitude_right in right.items():
            output[state_left | state_right] += amplitude_left * amplitude_right
    return _normalize(dict(output))


def _linear_combination(terms: Iterable[Tuple[complex, Vector]]) -> Vector:
    output: defaultdict[int, complex] = defaultdict(complex)
    for coefficient, vector in terms:
        for state, amplitude in vector.items():
            output[state] += coefficient * amplitude
    return dict(output)


def _apply_x_level(capacity: int, level: int, vector: Vector, adjoint: bool = False) -> Vector:
    output: defaultdict[int, complex] = defaultdict(complex)
    for spin in (0, 1):
        block_1 = ((0 * capacity + level) * 2) + spin
        block_2 = ((1 * capacity + level) * 2) + spin
        create_orbital, annihilate_orbital = (
            (block_2, block_1) if adjoint else (block_1, block_2)
        )
        term = _cdag_c(vector, create_orbital, annihilate_orbital)
        for state, amplitude in term.items():
            output[state] += amplitude
    return dict(output)


def _apply_bx_level(capacity: int, level: int, vector: Vector) -> Vector:
    return _linear_combination((
        (1.0, _apply_x_level(capacity, level, vector)),
        (1.0, _apply_x_level(capacity, level, vector, True)),
    ))


def _apply_by_level(capacity: int, level: int, vector: Vector) -> Vector:
    return _linear_combination((
        (-1.0j, _apply_x_level(capacity, level, vector)),
        (1.0j, _apply_x_level(capacity, level, vector, True)),
    ))


def local_bridge_gram_fock(
    capacity: int,
    total_pairs: int,
    microscopic_gx: float,
    microscopic_gy: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Gram matrix of dS_X=-g_X B_Y/2 and dS_Y=g_Y B_X/2."""
    r_min = max(0, total_pairs - capacity)
    r_max = min(capacity, total_pairs)
    r_values = np.arange(r_min, r_max + 1)
    states = [product_agp(capacity, int(r), total_pairs - int(r)) for r in r_values]
    gram = np.zeros((len(states), len(states)), dtype=complex)

    for level in range(capacity):
        bx_images = [_apply_bx_level(capacity, level, state) for state in states]
        by_images = [_apply_by_level(capacity, level, state) for state in states]
        for row in range(len(states)):
            for column in range(len(states)):
                gram[row, column] += (
                    (microscopic_gy * microscopic_gy / 4.0)
                    * _inner(bx_images[row], bx_images[column])
                    + (microscopic_gx * microscopic_gx / 4.0)
                    * _inner(by_images[row], by_images[column])
                )
    return r_values, gram


def local_bridge_gram_formula(
    capacity: int,
    total_pairs: int,
    microscopic_gx: float,
    microscopic_gy: float,
) -> tuple[np.ndarray, np.ndarray]:
    r_min = max(0, total_pairs - capacity)
    r_max = min(capacity, total_pairs)
    r_values = np.arange(r_min, r_max + 1)
    matrix = np.zeros((len(r_values), len(r_values)), dtype=float)
    diagonal_scale = 0.5 * (microscopic_gx**2 + microscopic_gy**2)
    hopping_scale = 0.5 * (microscopic_gy**2 - microscopic_gx**2)
    for index, r_value in enumerate(r_values):
        r = int(r_value)
        numerator = (
            (total_pairs - r) * (capacity - r)
            + r * (capacity - total_pairs + r)
        )
        matrix[index, index] = diagonal_scale * numerator / capacity
        if index + 1 < len(r_values):
            root = math.sqrt(
                (total_pairs - r)
                * (capacity - r)
                * (r + 1)
                * (capacity - total_pairs + r + 1)
            )
            matrix[index, index + 1] = hopping_scale * root / capacity
            matrix[index + 1, index] = matrix[index, index + 1]
    return r_values, matrix


def rho(capacity: int, pair_number: int) -> float:
    return pair_number * (capacity - pair_number) / (capacity - 1)


def total_y_jacobi(
    capacity: int,
    total_pairs: int,
    interaction_u: float,
    microscopic_gx: float,
    microscopic_gy: float,
) -> Array:
    """Base QGN curvature plus the local material bridge along y."""
    r_values, matrix = local_bridge_gram_formula(
        capacity, total_pairs, microscopic_gx, microscopic_gy
    )
    for index, r_value in enumerate(r_values):
        r = int(r_value)
        # m_1,yy^{-1}=U*K_11*(1/4)=U/2, K_11=2
        # m_2,yy^{-1}=U*K_22*(1/8)=5U/8, K_22=5
        matrix[index, index] += (
            2.0 * interaction_u * rho(capacity, r)
            + 2.5 * interaction_u * rho(capacity, total_pairs - r)
        )
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--torus", type=int, default=3, help="N for the N x N one-pair check")
    parser.add_argument("--fock-v", type=int, default=4, help="Capacity for the Fock-space Jacobi check")
    parser.add_argument("--fock-n", type=int, default=None, help="Total pair number for the Fock check")
    parser.add_argument("--gx", type=float, default=0.8, help="Microscopic X-quadrature coupling")
    parser.add_argument("--gy", type=float, default=0.3, help="Microscopic Y-quadrature coupling")
    parser.add_argument("--u", type=float, default=1.0)
    args = parser.parse_args()

    maximum_unitarity_error = 0.0
    for block in (1, 2):
        for kx, ky in momentum_grid(17, 19):
            matrix = frame(kx, ky, block)
            maximum_unitarity_error = max(
                maximum_unitarity_error,
                float(np.max(np.abs(matrix.conj().T @ matrix - np.eye(2)))),
            )
    print(f"maximum paraunitarity error: {maximum_unitarity_error:.3e}")

    for block in (1, 2):
        metric = numerical_metric(block, grid=60)
        metric_error = float(np.max(np.abs(metric - EXPECTED_METRIC[block])))
        print(f"block {block} average metric:\n{metric}")
        print(f"block {block} metric error: {metric_error:.3e}")
        masses = np.array([
            one_pair_mass(args.torus, args.torus, block, 0),
            one_pair_mass(args.torus, args.torus, block, 1),
        ])
        print(f"block {block} one-pair inverse masses at unit coupling: {masses}")
        print(
            "block " + str(block) + " mass/metric error: "
            + f"{float(np.max(np.abs(masses - np.diag(EXPECTED_METRIC[block])))):.3e}"
        )

    error_x, error_y = bridge_derivative_error()
    print(f"maximum bridge dB_X/dA_y error: {error_x:.3e}")
    print(f"maximum bridge dB_Y/dA_y error: {error_y:.3e}")

    total_pairs = args.fock_v if args.fock_n is None else args.fock_n
    r_fock, gram_fock = local_bridge_gram_fock(
        args.fock_v, total_pairs, args.gx, args.gy
    )
    r_formula, gram_formula = local_bridge_gram_formula(
        args.fock_v, total_pairs, args.gx, args.gy
    )
    if not np.array_equal(r_fock, r_formula):
        raise RuntimeError("Composition grids disagree")
    gram_error = float(np.max(np.abs(gram_fock - gram_formula)))
    print(f"composition values: {r_fock.tolist()}")
    print(f"maximum Fock/Jacobi bridge error: {gram_error:.3e}")

    jacobi = total_y_jacobi(
        args.fock_v, total_pairs, args.u, args.gx, args.gy
    )
    eigenvalues = np.linalg.eigvalsh(jacobi)
    print(f"total y-curvature Jacobi eigenvalues: {eigenvalues}")

    tolerance = 2.0e-8
    if maximum_unitarity_error > tolerance or error_x > tolerance or error_y > tolerance:
        raise SystemExit("Single-particle verification failed")
    if gram_error > tolerance:
        raise SystemExit("Many-body Jacobi verification failed")
    print("All checks passed.")


if __name__ == "__main__":
    main()
