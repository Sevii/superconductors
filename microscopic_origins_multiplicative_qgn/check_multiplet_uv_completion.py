#!/usr/bin/env python3
"""Exact local audit for a multiplet UV origin of the QGN escape row.

The physical bond contains two blocks a=1,2, two endpoints L,R, and spin.
It constructs

    Pi12 = (1-W1 W2)/2,
    B = sum_sigma(c^†_{1L,s} c_{2R,s}+c^†_{1R,s} c_{2L,s}+h.c.),
    K = Pi12 B = B Pi12.

The script verifies:
  1. the exact normal-order identity for B^2;
  2. the seniority-zero hard-core identities for B^2 and K^2;
  3. the rank-three projector decomposition of K^2;
  4. the passive odd-multiplet UV completion, whose zero-energy Schur
     complement is H_base - (t^2/Delta_a) K^2;
  5. the compensated two-multiplet UV completion, whose zero-energy Schur
     complement is H_base + alpha K^2;
  6. exact kernel preservation, positivity, Peierls covariance, and spectral
     convergence of the UV cluster to its effective Hamiltonian.

Dependencies: numpy, scipy.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence

import numpy as np
import scipy.linalg as la

TOL = 2.0e-10


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


def mode(block: int, endpoint: int, spin: int) -> int:
    return ((block * 2 + endpoint) * 2 + spin)


@dataclass
class BondOperators:
    basis: list[int]
    identity: np.ndarray
    W1: np.ndarray
    W2: np.ndarray
    Pi1: np.ndarray
    Pi2: np.ndarray
    Pi12: np.ndarray
    B: np.ndarray
    K: np.ndarray
    q_rows: list[np.ndarray]
    right_charge: np.ndarray


def cdagger_c(basis: Sequence[int], dst: int, src: int) -> np.ndarray:
    index = {state: row for row, state in enumerate(basis)}
    M = np.zeros((len(basis), len(basis)), dtype=complex)
    for col, state in enumerate(basis):
        first = annihilate(state, src)
        if first is None:
            continue
        state_1, sign_1 = first
        second = create(state_1, dst)
        if second is None:
            continue
        state_2, sign_2 = second
        M[index[state_2], col] = sign_1 * sign_2
    return M


def induced_permutation(basis: Sequence[int], permutation: Sequence[int]) -> np.ndarray:
    index = {state: row for row, state in enumerate(basis)}
    M = np.zeros((len(basis), len(basis)), dtype=complex)
    for col, state in enumerate(basis):
        occupied = [m for m in range(len(permutation)) if (state >> m) & 1]
        mapped = [permutation[m] for m in occupied]
        inversions = sum(
            mapped[i] > mapped[j]
            for i in range(len(mapped))
            for j in range(i + 1, len(mapped))
        )
        sign = -1 if inversions & 1 else 1
        new_state = 0
        for m in mapped:
            new_state |= 1 << m
        M[index[new_state], col] = sign
    return M


def build_full_bond() -> BondOperators:
    basis = list(range(1 << 8))
    I = np.eye(len(basis), dtype=complex)
    swaps: list[np.ndarray] = []
    for block in range(2):
        perm = list(range(8))
        for spin in range(2):
            left = mode(block, 0, spin)
            right = mode(block, 1, spin)
            perm[left], perm[right] = perm[right], perm[left]
        swaps.append(induced_permutation(basis, perm))
    W1, W2 = swaps
    Pi1 = 0.5 * (I - W1)
    Pi2 = 0.5 * (I - W2)
    Pi12 = 0.5 * (I - W1 @ W2)

    B = np.zeros_like(I)
    for spin in range(2):
        B += cdagger_c(basis, mode(0, 0, spin), mode(1, 1, spin))
        B += cdagger_c(basis, mode(1, 1, spin), mode(0, 0, spin))
        B += cdagger_c(basis, mode(0, 1, spin), mode(1, 0, spin))
        B += cdagger_c(basis, mode(1, 0, spin), mode(0, 1, spin))
    K = Pi12 @ B

    q_rows: list[np.ndarray] = []
    for block in range(2):
        for endpoint in range(2):
            up = mode(block, endpoint, 0)
            down = mode(block, endpoint, 1)
            diagonal = np.asarray(
                [
                    (((state >> up) & 1) - ((state >> down) & 1)) ** 2
                    for state in basis
                ],
                dtype=float,
            )
            q_rows.append(np.diag(diagonal).astype(complex))

    right_charge = np.asarray(
        [
            sum(
                (state >> mode(block, 1, spin)) & 1
                for block in range(2)
                for spin in range(2)
            )
            for state in basis
        ],
        dtype=float,
    )
    return BondOperators(
        basis=basis,
        identity=I,
        W1=W1,
        W2=W2,
        Pi1=Pi1,
        Pi2=Pi2,
        Pi12=Pi12,
        B=B,
        K=K,
        q_rows=q_rows,
        right_charge=right_charge,
    )


def pair_create(basis: Sequence[int], up: int, down: int) -> np.ndarray:
    index = {state: row for row, state in enumerate(basis)}
    M = np.zeros((len(basis), len(basis)), dtype=complex)
    for col, state in enumerate(basis):
        first = create(state, down)
        if first is None:
            continue
        state_1, sign_1 = first
        second = create(state_1, up)
        if second is None:
            continue
        state_2, sign_2 = second
        M[index[state_2], col] = sign_1 * sign_2
    return M


def edge_square_normal_order(ops: BondOperators, i: tuple[int, int], j: tuple[int, int]) -> np.ndarray:
    """Normal-ordered expression for [sum_sigma(c_i^†c_j+h.c.)]^2."""
    basis = ops.basis
    result = np.zeros_like(ops.identity)
    for mi, mj in zip(i, j):
        ni = cdagger_c(basis, mi, mi)
        nj = cdagger_c(basis, mj, mj)
        result += ni + nj - 2.0 * ni @ nj

    Pi_dag = pair_create(basis, i[0], i[1])
    Pj_dag = pair_create(basis, j[0], j[1])
    result += 2.0 * (Pi_dag @ Pj_dag.conj().T + Pj_dag @ Pi_dag.conj().T)

    Si_plus = cdagger_c(basis, i[0], i[1])
    Sj_plus = cdagger_c(basis, j[0], j[1])
    result -= 2.0 * (
        Si_plus @ Sj_plus.conj().T + Si_plus.conj().T @ Sj_plus
    )
    return result


def hard_core_operators() -> tuple[list[np.ndarray], list[np.ndarray], list[np.ndarray]]:
    dim = 16
    lowers: list[np.ndarray] = []
    for orbital in range(4):
        M = np.zeros((dim, dim), dtype=complex)
        for state in range(dim):
            if (state >> orbital) & 1:
                M[state ^ (1 << orbital), state] = 1.0
        lowers.append(M)
    raises = [M.conj().T for M in lowers]
    numbers = [raises[i] @ lowers[i] for i in range(4)]
    return lowers, raises, numbers


def seniority_zero_indices(ops: BondOperators) -> list[int]:
    result: list[int] = []
    index = {state: row for row, state in enumerate(ops.basis)}
    for config in range(16):
        state = 0
        for orbital in range(4):
            if (config >> orbital) & 1:
                block, endpoint = divmod(orbital, 2)
                state |= 1 << mode(block, endpoint, 0)
                state |= 1 << mode(block, endpoint, 1)
        result.append(index[state])
    return result


def hard_core_formula(ops: BondOperators) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    lowers, raises, numbers = hard_core_operators()
    I = np.eye(16, dtype=complex)
    D = np.zeros((16, 16), dtype=complex)
    for i, j in ((0, 3), (1, 2)):
        D += (
            numbers[i]
            + numbers[j]
            - 2.0 * numbers[i] @ numbers[j]
            + raises[i] @ lowers[j]
            + raises[j] @ lowers[i]
        )

    swaps: list[np.ndarray] = []
    for left, right in ((0, 1), (2, 3)):
        swaps.append(
            I
            - numbers[left]
            - numbers[right]
            + 2.0 * numbers[left] @ numbers[right]
            + raises[left] @ lowers[right]
            + raises[right] @ lowers[left]
        )
    W1, W2 = swaps
    K2_formula = (I - W1 @ W2) @ D

    empty = np.eye(16, dtype=complex)
    full = np.eye(16, dtype=complex)
    for n in numbers:
        empty = empty @ (I - n)
        full = full @ n
    A = lowers[0] - lowers[1] - lowers[2] + lowers[3]
    C = lowers[0] @ lowers[2] - lowers[1] @ lowers[3]
    rank_three = A.conj().T @ empty @ A + 2.0 * C.conj().T @ empty @ C + A @ full @ A.conj().T
    return 2.0 * D, K2_formula, rank_three


def base_hamiltonian(ops: BondOperators, onsite: float, j1: float, j2: float) -> np.ndarray:
    H = sum((onsite / 2.0) * q for q in ops.q_rows)
    H += (j1 / 2.0) * ops.Pi1 + (j2 / 2.0) * ops.Pi2
    return 0.5 * (H + H.conj().T)


def gauge_conjugate(operator: np.ndarray, charge: np.ndarray, A: float) -> np.ndarray:
    phase = np.exp(-1j * A * charge)
    U = np.diag(phase)
    return U @ operator @ U.conj().T


def projector_basis(projector: np.ndarray, threshold: float = 0.5) -> np.ndarray:
    values, vectors = la.eigh(projector)
    return vectors[:, values > threshold]


def count_zero(matrix: np.ndarray, tolerance: float = 1.0e-8) -> int:
    return int(np.count_nonzero(np.abs(la.eigvalsh(matrix)) < tolerance))


def fixed_particle_indices(ops: BondOperators, particles: int) -> list[int]:
    return [i for i, state in enumerate(ops.basis) if state.bit_count() == particles]


def principal(matrix: np.ndarray, indices: Sequence[int]) -> np.ndarray:
    return matrix[np.ix_(indices, indices)]


def audit() -> None:
    ops = build_full_bond()
    I = ops.identity
    B2 = ops.B @ ops.B
    K2 = ops.K @ ops.K

    print("multiplet UV completion audit")
    print("physical mode order: (1L up/down, 1R up/down, 2L up/down, 2R up/down)")
    print()

    print("=== Algebraic identities ===")
    comm = la.norm(ops.Pi12 @ ops.B - ops.B @ ops.Pi12)
    pi_identity = la.norm(ops.Pi12 - (ops.Pi1 + ops.Pi2 - 2.0 * ops.Pi1 @ ops.Pi2))
    print(f"[Pi12,B] residual                         = {comm:.3e}")
    print(f"Pi12=Pi1+Pi2-2 Pi1 Pi2 residual          = {pi_identity:.3e}")

    e1 = ((mode(0, 0, 0), mode(0, 0, 1)), (mode(1, 1, 0), mode(1, 1, 1)))
    e2 = ((mode(0, 1, 0), mode(0, 1, 1)), (mode(1, 0, 0), mode(1, 0, 1)))
    T1 = np.zeros_like(I)
    T2 = np.zeros_like(I)
    for mi, mj in zip(*e1):
        T1 += cdagger_c(ops.basis, mi, mj) + cdagger_c(ops.basis, mj, mi)
    for mi, mj in zip(*e2):
        T2 += cdagger_c(ops.basis, mi, mj) + cdagger_c(ops.basis, mj, mi)
    normal_order = (
        edge_square_normal_order(ops, *e1)
        + edge_square_normal_order(ops, *e2)
        + 2.0 * T1 @ T2
    )
    print(f"full-Fock normal-order B^2 residual       = {la.norm(B2-normal_order):.3e}")

    sz = seniority_zero_indices(ops)
    B2_sz = principal(B2, sz)
    K2_sz = principal(K2, sz)
    B2_formula, K2_formula, rank_three = hard_core_formula(ops)
    print(f"seniority-zero B^2 formula residual       = {la.norm(B2_sz-B2_formula):.3e}")
    print(f"seniority-zero K^2 formula residual       = {la.norm(K2_sz-K2_formula):.3e}")
    print(f"rank-three decomposition residual         = {la.norm(K2_sz-rank_three):.3e}")
    print(f"spec(K^2|seniority zero)                  = {np.round(la.eigvalsh(K2_sz),10)}")

    # Compression of the ordinary B^2 term to the one-pair product-AGP basis.
    one_pair_positions = [1, 2, 4, 8]  # hard-core configurations 1L,1R,2L,2R
    B2_one = B2_sz[np.ix_(one_pair_positions, one_pair_positions)]
    agp1 = np.asarray([1.0, 1.0, 0.0, 0.0]) / math.sqrt(2.0)
    agp2 = np.asarray([0.0, 0.0, 1.0, 1.0]) / math.sqrt(2.0)
    Z = np.column_stack([agp1, agp2])
    unwanted = Z.conj().T @ B2_one @ Z
    print("P_Z B^2 P_Z in basis (|1,0>,|0,1>):")
    print(np.real_if_close(unwanted))

    odd_bound = np.min(la.eigvalsh(9.0 * ops.Pi12 - K2))
    print(f"min spec(9 Pi12-K^2)                      = {odd_bound:.3e}")
    if min(comm, 1.0) > TOL or abs(pi_identity) > TOL or la.norm(B2-normal_order) > TOL:
        raise AssertionError("full-Fock algebraic identity failed")
    if max(la.norm(B2_sz-B2_formula), la.norm(K2_sz-K2_formula), la.norm(K2_sz-rank_three)) > TOL:
        raise AssertionError("hard-core identity failed")
    if odd_bound < -TOL:
        raise AssertionError("K^2 <= 9 Pi12 failed")

    onsite, j1, j2 = 2.0, 1.0, 1.3
    Hbase = base_hamiltonian(ops, onsite, j1, j2)
    expected_kernel = 9  # (V+1)^2 for a two-site bond and two blocks
    print("\n=== Base parent ===")
    print(f"base kernel dimension                       = {count_zero(Hbase)} (expected {expected_kernel})")
    print(f"base next eigenvalue                        = {la.eigvalsh(Hbase)[expected_kernel]:.12f}")

    print("\n=== Passive odd-multiplet route ===")
    Delta_a = 20.0
    t = 0.8
    gamma = t * t / Delta_a
    sufficient = min(j1, j2) / 18.0
    print(f"gamma=t^2/Delta_a                          = {gamma:.12f}")
    print(f"uniform full-Fock sufficient bound          = {sufficient:.12f}")
    Hpass = Hbase - gamma * K2
    print(f"min spec(H_base-gamma K^2)                  = {la.eigvalsh(Hpass)[0]:.3e}")
    print(f"effective kernel dimension                  = {count_zero(Hpass)}")

    Qa_basis = projector_basis(ops.Pi12)
    K_rect = Qa_basis.conj().T @ ops.K
    Hhigh = Delta_a * np.eye(Qa_basis.shape[1], dtype=complex)
    Huv_pass = np.block(
        [
            [Hbase, t * K_rect.conj().T],
            [t * K_rect, Hhigh],
        ]
    )
    Schur_pass = Hbase - t * t * K_rect.conj().T @ la.solve(Hhigh, K_rect)
    print(f"passive zero-energy Schur residual          = {la.norm(Schur_pass-Hpass):.3e}")
    print(f"passive UV min eigenvalue                   = {la.eigvalsh(Huv_pass)[0]:.3e}")
    print(f"passive UV kernel dimension                 = {count_zero(Huv_pass)}")
    if gamma > sufficient + 1e-14:
        raise AssertionError("chosen passive parameters violate the stated sufficient bound")
    if la.eigvalsh(Hpass)[0] < -TOL or la.eigvalsh(Huv_pass)[0] < -TOL:
        raise AssertionError("passive UV completion is not positive")
    if count_zero(Huv_pass) != expected_kernel:
        raise AssertionError("passive UV kernel changed")

    Mcomp = np.asarray([[1.0, -1.0], [-1.0, 1.0]])
    Qpass = 2.0 * np.diag([j1, j2]) - 8.0 * gamma * Mcomp
    print("passive one-pair composition curvature:")
    print(Qpass)
    print(f"passive curvature eigenvalues               = {la.eigvalsh(Qpass)}")

    print("\n=== Compensated two-multiplet route ===")
    Delta_s, Delta_a2 = 12.0, 20.0
    beta = t * t / Delta_s
    alpha = t * t * (1.0 / Delta_s - 1.0 / Delta_a2)
    Qs = I - ops.Pi12
    Hdelta = Delta_s * Qs + Delta_a2 * ops.Pi12
    Hlow = Hbase + beta * B2
    Huv_pos = np.block([[Hlow, t * ops.B], [t * ops.B, Hdelta]])
    Hpos = Hbase + alpha * K2
    Schur_pos = Hlow - t * t * ops.B @ la.solve(Hdelta, ops.B)
    print(f"beta=t^2/Delta_s                           = {beta:.12f}")
    print(f"alpha=t^2(1/Delta_s-1/Delta_a)             = {alpha:.12f}")
    print(f"positive zero-energy Schur residual         = {la.norm(Schur_pos-Hpos):.3e}")
    print(f"compensated UV min eigenvalue               = {la.eigvalsh(Huv_pos)[0]:.3e}")
    print(f"compensated UV kernel dimension             = {count_zero(Huv_pos)}")
    if la.norm(Schur_pos-Hpos) > TOL or la.eigvalsh(Huv_pos)[0] < -TOL:
        raise AssertionError("compensated UV completion failed")
    if count_zero(Huv_pos) != expected_kernel:
        raise AssertionError("compensated UV kernel changed")

    Qpos = 2.0 * np.diag([j1, j2]) + 8.0 * alpha * Mcomp
    print("compensated one-pair composition curvature:")
    print(Qpos)
    print(f"compensated curvature eigenvalues            = {la.eigvalsh(Qpos)}")

    print("\n=== Peierls covariance of the Schur complement ===")
    A = 0.173
    B_A = gauge_conjugate(ops.B, ops.right_charge, A)
    Pi_A = gauge_conjugate(ops.Pi12, ops.right_charge, A)
    Pi1_A = gauge_conjugate(ops.Pi1, ops.right_charge, A)
    Pi2_A = gauge_conjugate(ops.Pi2, ops.right_charge, A)
    Hbase_A = sum((onsite / 2.0) * q for q in ops.q_rows) + (j1 / 2.0) * Pi1_A + (j2 / 2.0) * Pi2_A
    Hdelta_A = Delta_s * (I - Pi_A) + Delta_a2 * Pi_A
    Hlow_A = Hbase_A + beta * (B_A @ B_A)
    Schur_A = Hlow_A - t * t * B_A @ la.solve(Hdelta_A, B_A)
    K_A = Pi_A @ B_A
    target_A = Hbase_A + alpha * (K_A @ K_A)
    print(f"twisted compensated Schur residual           = {la.norm(Schur_A-target_A):.3e}")
    if la.norm(Schur_A-target_A) > TOL:
        raise AssertionError("Peierls-covariant Schur reduction failed")

    print("\n=== Isolated positive-penalty spectral scaling ===")
    # Work in the two-particle sector with H_base set to zero.  The nonzero
    # low multiplet has B^2 eigenvalue z=4, so its effective energy is 4 alpha.
    indices = fixed_particle_indices(ops, 2)
    B_N = principal(ops.B, indices)
    Pi_N = principal(ops.Pi12, indices)
    I_N = np.eye(len(indices), dtype=complex)
    for t_scale in (1.6, 0.8, 0.4, 0.2, 0.1):
        beta_s = t_scale * t_scale / Delta_s
        alpha_s = t_scale * t_scale * (1.0 / Delta_s - 1.0 / Delta_a2)
        Hd = Delta_s * (I_N - Pi_N) + Delta_a2 * Pi_N
        Huv = np.block(
            [
                [beta_s * (B_N @ B_N), t_scale * B_N],
                [t_scale * B_N, Hd],
            ]
        )
        low = la.eigvalsh(Huv)
        positive_low = low[(low > 1.0e-9) & (low < Delta_s / 2.0)]
        exact_low = float(positive_low[0])
        target = 4.0 * alpha_s
        error = target - exact_low
        print(
            f"t={t_scale:4.2f}: E_UV={exact_low:.12f}, 4alpha={target:.12f}, "
            f"error={error:.3e}, error/t^4={error/t_scale**4:.12f}"
        )
        if error < -TOL:
            raise AssertionError("UV low energy exceeded the second-order target")

    print("\nAll microscopic multiplet checks passed.")


if __name__ == "__main__":
    audit()
