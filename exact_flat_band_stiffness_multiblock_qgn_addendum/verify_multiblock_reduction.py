#!/usr/bin/env python3
"""Numerical certificates for the multi-block AGP source decomposition.

The tests use canonical skew blocks and exact finite Fock-space vectors.
They verify:
  1. the complete diagonal + intrablock + transfer norm decomposition;
  2. the transfer coefficient tau_{a<-b};
  3. the nearest-neighbor composition-space matrix formula;
  4. the two-block counterexample with positive diagonal responses and
     zero total-sector curvature.

Only NumPy is required.
"""

from __future__ import annotations

import itertools
import math
from collections import defaultdict
from typing import Dict, Iterable, List, Sequence

import numpy as np

State = Dict[int, complex]


def annihilate(state: State, mode: int) -> State:
    out: dict[int, complex] = defaultdict(complex)
    lower = (1 << mode) - 1
    for mask, amp in state.items():
        if (mask >> mode) & 1:
            sign = -1 if (mask & lower).bit_count() % 2 else 1
            out[mask ^ (1 << mode)] += sign * amp
    return dict(out)


def create(state: State, mode: int) -> State:
    out: dict[int, complex] = defaultdict(complex)
    lower = (1 << mode) - 1
    for mask, amp in state.items():
        if not ((mask >> mode) & 1):
            sign = -1 if (mask & lower).bit_count() % 2 else 1
            out[mask | (1 << mode)] += sign * amp
    return dict(out)


def apply_one_body(state: State, matrix: np.ndarray) -> State:
    out: dict[int, complex] = defaultdict(complex)
    rows, cols = np.nonzero(np.abs(matrix) > 1e-14)
    for i, j in zip(rows.tolist(), cols.tolist()):
        tmp = create(annihilate(state, j), i)
        for mask, amp in tmp.items():
            out[mask] += matrix[i, j] * amp
    return dict(out)


def inner(left: State, right: State) -> complex:
    return sum(np.conjugate(amp) * right.get(mask, 0.0) for mask, amp in left.items())


def norm_sq(state: State) -> float:
    return float(np.real(inner(state, state)))


def offsets(block_sizes: Sequence[int]) -> List[int]:
    result: List[int] = []
    cursor = 0
    for size in block_sizes:
        result.append(cursor)
        cursor += 2 * size
    return result


def product_agp(block_sizes: Sequence[int], occupations: Sequence[int]) -> State:
    """Normalized product of canonical block AGPs."""
    if len(block_sizes) != len(occupations):
        raise ValueError("block_sizes and occupations must have equal length")
    offs = offsets(block_sizes)
    choices = [
        list(itertools.combinations(range(size), occupation))
        for size, occupation in zip(block_sizes, occupations)
    ]
    normalization = math.sqrt(
        math.prod(math.comb(size, occupation) for size, occupation in zip(block_sizes, occupations))
    )
    state: State = {}
    for selected_blocks in itertools.product(*choices):
        mask = 0
        for off, selected_pairs in zip(offs, selected_blocks):
            for pair in selected_pairs:
                mask |= 1 << (off + 2 * pair)
                mask |= 1 << (off + 2 * pair + 1)
        state[mask] = 1.0 / normalization
    return state


def embedded_block(
    block_sizes: Sequence[int], matrix: np.ndarray, target: int, source: int
) -> np.ndarray:
    total_modes = 2 * sum(block_sizes)
    offs = offsets(block_sizes)
    out = np.zeros((total_modes, total_modes), dtype=complex)
    target_slice = slice(offs[target], offs[target] + 2 * block_sizes[target])
    source_slice = slice(offs[source], offs[source] + 2 * block_sizes[source])
    out[target_slice, source_slice] = matrix
    return out


def seed_transfer(
    block_sizes: Sequence[int], full_matrix: np.ndarray, target: int, source: int
) -> State:
    offs = offsets(block_sizes)
    target_slice = slice(offs[target], offs[target] + 2 * block_sizes[target])
    source_slice = slice(offs[source], offs[source] + 2 * block_sizes[source])
    block = full_matrix[target_slice, source_slice]
    occupation = [0] * len(block_sizes)
    occupation[source] = 1
    return apply_one_body(
        product_agp(block_sizes, occupation),
        embedded_block(block_sizes, block, target, source),
    )


def test_complete_source_norm(rng: np.random.Generator) -> float:
    block_sizes = [2, 3]
    occupations = [1, 2]
    total_modes = 2 * sum(block_sizes)
    raw = rng.normal(size=(total_modes, total_modes)) + 1j * rng.normal(
        size=(total_modes, total_modes)
    )
    matrix = (raw + raw.conj().T) / 2
    state = product_agp(block_sizes, occupations)
    actual = norm_sq(apply_one_body(state, matrix))

    offs = offsets(block_sizes)
    betas: List[complex] = []
    diagonal_seed_norms: List[float] = []
    for a, size in enumerate(block_sizes):
        sl = slice(offs[a], offs[a] + 2 * size)
        block = matrix[sl, sl]
        beta = np.trace(block) / size
        betas.append(beta)
        occupation = [0] * len(block_sizes)
        occupation[a] = 1
        pair_state = product_agp(block_sizes, occupation)
        diagonal_operator = embedded_block(block_sizes, block, a, a)
        acted = apply_one_body(pair_state, diagonal_operator)
        seed = dict(acted)
        for mask, amp in pair_state.items():
            seed[mask] = seed.get(mask, 0.0) - beta * amp
        diagonal_seed_norms.append(norm_sq(seed))

    predicted = abs(sum(n * beta for n, beta in zip(occupations, betas))) ** 2
    for n, size, seed_norm in zip(occupations, block_sizes, diagonal_seed_norms):
        rho = n * (size - n) / (size - 1) if size > 1 else 0.0
        predicted += rho * seed_norm

    for target in range(len(block_sizes)):
        for source in range(len(block_sizes)):
            if target == source:
                continue
            seed = seed_transfer(block_sizes, matrix, target, source)
            tau = occupations[source] * (
                block_sizes[target] - occupations[target]
            ) / block_sizes[target]
            predicted += tau * norm_sq(seed)

    return abs(actual - predicted)


def test_transfer_polarization(rng: np.random.Generator) -> float:
    block_sizes = [3, 4]
    occupations = [2, 3]
    target, source = 0, 1
    x = rng.normal(size=(6, 8)) + 1j * rng.normal(size=(6, 8))
    y = rng.normal(size=(6, 8)) + 1j * rng.normal(size=(6, 8))
    bx = embedded_block(block_sizes, x, target, source)
    by = embedded_block(block_sizes, y, target, source)
    state = product_agp(block_sizes, occupations)
    lhs = inner(apply_one_body(state, by), apply_one_body(state, bx))

    one_pair = [0, 1]
    seed_state = product_agp(block_sizes, one_pair)
    seed_inner = inner(
        apply_one_body(seed_state, by), apply_one_body(seed_state, bx)
    )
    tau = occupations[source] * (
        block_sizes[target] - occupations[target]
    ) / block_sizes[target]
    return abs(lhs - tau * seed_inner)


def test_composition_matrix(rng: np.random.Generator) -> float:
    block_sizes = [2, 2]
    total_modes = 8
    x = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    matrix = np.zeros((total_modes, total_modes), dtype=complex)
    matrix[:4, 4:] = x
    matrix[4:, :4] = x.conj().T

    compositions = [[0, 2], [1, 1], [2, 0]]
    states = [product_agp(block_sizes, comp) for comp in compositions]
    sources = [apply_one_body(state, matrix) for state in states]
    actual = np.array(
        [[inner(sources[i], sources[j]) for j in range(3)] for i in range(3)]
    )

    z_0_from_1 = seed_transfer(block_sizes, matrix, 0, 1)
    z_1_from_0 = seed_transfer(block_sizes, matrix, 1, 0)
    k01 = norm_sq(z_0_from_1)
    k10 = norm_sq(z_1_from_0)
    omega = inner(z_1_from_0, z_0_from_1)

    predicted = np.zeros((3, 3), dtype=complex)
    for index, comp in enumerate(compositions):
        n0, n1 = comp
        tau_0_from_1 = n1 * (2 - n0) / 2
        tau_1_from_0 = n0 * (2 - n1) / 2
        predicted[index, index] = tau_0_from_1 * k01 + tau_1_from_0 * k10

    # [0,2] -> [1,1] and [1,1] -> [2,0]
    for lower in (0, 1):
        n0, n1 = compositions[lower]
        coefficient = math.sqrt(
            n1 * (2 - n0) * (n0 + 1) * (2 - n1 + 1) / 4
        )
        predicted[lower + 1, lower] = coefficient * omega
        predicted[lower, lower + 1] = np.conjugate(predicted[lower + 1, lower])

    return float(np.max(np.abs(actual - predicted)))


def test_minimal_counterexample() -> tuple[float, np.ndarray]:
    block_sizes = [1, 1]
    matrix = np.zeros((4, 4), dtype=complex)
    matrix[:2, 2:] = np.eye(2)
    matrix[2:, :2] = np.eye(2)
    pair_a = product_agp(block_sizes, [1, 0])
    pair_b = product_agp(block_sizes, [0, 1])
    source_a = apply_one_body(pair_a, matrix)
    source_b = apply_one_body(pair_b, matrix)
    equality_error = math.sqrt(norm_sq({
        mask: source_a.get(mask, 0.0) - source_b.get(mask, 0.0)
        for mask in set(source_a) | set(source_b)
    }))
    q_matrix = np.array(
        [
            [inner(source_a, source_a), inner(source_a, source_b)],
            [inner(source_b, source_a), inner(source_b, source_b)],
        ]
    )
    eigenvalues = np.linalg.eigvalsh(q_matrix)
    return equality_error, eigenvalues


def main() -> None:
    rng = np.random.default_rng(20260724)
    errors: List[float] = []
    for _ in range(20):
        errors.append(test_complete_source_norm(rng))
        errors.append(test_transfer_polarization(rng))
        errors.append(test_composition_matrix(rng))

    equality_error, eigenvalues = test_minimal_counterexample()
    errors.append(equality_error)

    maximum_error = max(errors)
    print(f"maximum identity error: {maximum_error:.3e}")
    print("minimal-counterexample eigenvalues:", eigenvalues)
    if maximum_error > 5e-11:
        raise SystemExit("FAIL: an identity exceeded tolerance")
    if abs(eigenvalues[0]) > 5e-12 or abs(eigenvalues[1] - 4.0) > 5e-12:
        raise SystemExit("FAIL: unexpected counterexample spectrum")
    print("PASS")


if __name__ == "__main__":
    main()
