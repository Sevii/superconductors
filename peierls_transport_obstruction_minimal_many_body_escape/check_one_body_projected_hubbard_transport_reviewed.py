#!/usr/bin/env python3
"""Audit the covariant-transport obstruction for two-block projected-Hubbard factors.

The zero-twist one-spin factors are the repaired L_A=L_B=2 frame model from
`multiblock_qgn_addendum_reviewed(3).tex`.  They generate u(2)_A + u(2)_B
and have an exact two-dimensional one-pair zero manifold spanned by the two
block AGPs.

Three first-source families are compared:

1. common transport: b_lam = i [G, x_lam] with the same Hermitian G for every
   factor.  This is the infinitesimal form of a fully covariant transport of
   the complete factor family.  Its least-squares source must vanish.

2. common transport plus block-preserving intrinsic derivatives.  The common
   transport is removed and the response remains diagonal in composition.

3. a deliberately non-covariant, factor-label-dependent interblock source.
   This need not be a common commutator and produces a nonzero off-diagonal
   composition-space curvature, illustrating exactly what a physical model
   would have to realize.

Dependencies: numpy, scipy
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import scipy.linalg as la

TOL = 1.0e-10


def fixed_weight_basis(modes: int, particles: int) -> list[int]:
    return [state for state in range(1 << modes) if state.bit_count() == particles]


def annihilate(state: int, mode: int):
    if not ((state >> mode) & 1):
        return None
    parity = (state & ((1 << mode) - 1)).bit_count() & 1
    return state ^ (1 << mode), (-1 if parity else 1)


def create(state: int, mode: int):
    if (state >> mode) & 1:
        return None
    parity = (state & ((1 << mode) - 1)).bit_count() & 1
    return state | (1 << mode), (-1 if parity else 1)


def cdagger_c(basis: list[int], dst: int, src: int) -> np.ndarray:
    index = {state: row for row, state in enumerate(basis)}
    out = np.zeros((len(basis), len(basis)), dtype=complex)
    for col, state in enumerate(basis):
        step1 = annihilate(state, src)
        if step1 is None:
            continue
        state1, sign1 = step1
        step2 = create(state1, dst)
        if step2 is None:
            continue
        state2, sign2 = step2
        out[index[state2], col] += sign1 * sign2
    return out


def one_body_second_quantized(basis: list[int], matrix: np.ndarray, spin: int) -> np.ndarray:
    """Second quantize a one-spin matrix.

    One-spin indices are 0..3.  Spin-up modes are 0..3 and spin-down modes
    are 4..7.
    """
    offset = 0 if spin == 0 else matrix.shape[0]
    out = np.zeros((len(basis), len(basis)), dtype=complex)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if abs(matrix[i, j]) > 0:
                out += matrix[i, j] * cdagger_c(basis, offset + i, offset + j)
    return out


def qgn_factor(basis: list[int], x: np.ndarray) -> np.ndarray:
    return one_body_second_quantized(basis, x, 0) - one_body_second_quantized(basis, np.conjugate(x), 1)


def normalized(v: np.ndarray) -> np.ndarray:
    n = la.norm(v)
    if n < TOL:
        raise ValueError("zero vector")
    return v / n


def rank_one(weight: float, vector: np.ndarray) -> np.ndarray:
    vector = normalized(np.asarray(vector, dtype=complex))
    return weight * np.outer(vector, np.conjugate(vector))


def repaired_factors() -> list[np.ndarray]:
    e1 = np.array([1.0, 0.0])
    e2 = np.array([0.0, 1.0])

    def f(theta: float) -> np.ndarray:
        return np.array([np.cos(theta), np.sin(theta)])

    def g(theta: float) -> np.ndarray:
        return np.array([-np.sin(theta), np.cos(theta)])

    p_a = [
        rank_one(1 / 3, e1),
        rank_one(1 / 3, e2),
        rank_one(2 / 3, f(np.pi / 4)),
        rank_one(2 / 3, g(np.pi / 4)),
    ]
    p_b = [
        rank_one(2 / 3, e1),
        rank_one(2 / 3, e2),
        rank_one(1 / 3, f(np.pi / 6)),
        rank_one(1 / 3, g(np.pi / 6)),
    ]
    return [la.block_diag(a, b) for a, b in zip(p_a, p_b, strict=True)]


def pair_state(basis: list[int], block: int) -> np.ndarray:
    """One normalized onsite AGP pair in a capacity-two block."""
    index = {state: row for row, state in enumerate(basis)}
    out = np.zeros(len(basis), dtype=complex)
    start = 2 * block
    for orbital in range(start, start + 2):
        state = (1 << orbital) | (1 << (4 + orbital))
        # With the chosen mode ordering, c_up^dagger c_down^dagger |0>
        # differs by a fixed sign from ascending-order bit basis.  The same
        # sign occurs for every orbital and is immaterial after normalization.
        out[index[state]] = 1.0
    return normalized(out)


def target_projector_ker_d_dagger(factors: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    d = np.vstack(factors)
    # ker(D^dagger) is the orthogonal complement of ran(D).
    u, s, _ = la.svd(d, full_matrices=True)
    rank = int(np.sum(s > TOL))
    ker = u[:, rank:]
    pi = ker @ np.conjugate(ker.T)
    return d, pi


def physical_twist_source(basis: list[int], b: np.ndarray) -> np.ndarray:
    """TR-covariant electromagnetic source: up b plus down conjugate(b)."""
    return one_body_second_quantized(basis, b, 0) + one_body_second_quantized(basis, np.conjugate(b), 1)


def source_column_from_ops(source_ops: list[np.ndarray], state: np.ndarray) -> np.ndarray:
    return np.concatenate([op @ state for op in source_ops])


def curvature_on_zero_manifold(
    basis: list[int],
    factors: list[np.ndarray],
    source_ops: list[np.ndarray],
    zero_states: list[np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    fock_factors = [qgn_factor(basis, x) for x in factors]
    _, pi = target_projector_ker_d_dagger(fock_factors)
    t_columns = [pi @ source_column_from_ops(source_ops, psi) for psi in zero_states]
    t = np.column_stack(t_columns)
    q = np.conjugate(t.T) @ t
    return q, t


def solve_common_generator(factors: list[np.ndarray], b_factors: list[np.ndarray]) -> tuple[np.ndarray, float]:
    """Least-squares solve b_lam^{AB}=i(G_AB x_B - x_A G_AB)."""
    la_dim = 2
    lb_dim = 2
    rows = []
    rhs = []
    for x, b in zip(factors, b_factors, strict=True):
        xa = x[:la_dim, :la_dim]
        xb = x[la_dim:, la_dim:]
        bab = b[:la_dim, la_dim:]
        # vec(G xb - xa G) = (xb^T kron I - I kron xa) vec(G)
        mat = np.kron(xb.T, np.eye(la_dim)) - np.kron(np.eye(lb_dim), xa)
        rows.append(1j * mat)
        rhs.append(bab.reshape(-1, order="F"))
    design = np.vstack(rows)
    target = np.concatenate(rhs)
    vec_g, *_ = la.lstsq(design, target)
    g_ab = vec_g.reshape((la_dim, lb_dim), order="F")
    residual = la.norm(design @ vec_g - target)
    return g_ab, float(residual)


def main() -> None:
    rng = np.random.default_rng(20260725)
    factors = repaired_factors()
    basis = fixed_weight_basis(8, 2)
    zero_states = [pair_state(basis, 0), pair_state(basis, 1)]

    fock_factors = [qgn_factor(basis, x) for x in factors]
    h0 = 0.5 * sum(s @ s for s in fock_factors)
    evals = la.eigvalsh(h0)
    nullity = int(np.sum(np.abs(evals) < TOL))

    # Random Hermitian common generator with only AB/BA blocks.
    z = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
    g = np.block([[np.zeros((2, 2)), z], [np.conjugate(z.T), np.zeros((2, 2))]])
    common_one_spin = [1j * (g @ x - x @ g) for x in factors]
    g_fock = one_body_second_quantized(basis, g, 0) + one_body_second_quantized(basis, np.conjugate(g), 1)
    common_ops = [1j * (g_fock @ s - s @ g_fock) for s in fock_factors]

    q_common, t_common = curvature_on_zero_manifold(basis, factors, common_ops, zero_states)
    _, common_fit_residual = solve_common_generator(factors, common_one_spin)

    # Add block-preserving intrinsic derivatives.
    intrinsic: list[np.ndarray] = []
    for _ in factors:
        ha = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        hb = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        ha = (ha + np.conjugate(ha.T)) / 2
        hb = (hb + np.conjugate(hb.T)) / 2
        intrinsic.append(la.block_diag(ha, hb))
    intrinsic_ops = [physical_twist_source(basis, d) for d in intrinsic]
    transported_plus_intrinsic_ops = [c + d for c, d in zip(common_ops, intrinsic_ops, strict=True)]
    q_covariant, _ = curvature_on_zero_manifold(
        basis, factors, transported_plus_intrinsic_ops, zero_states
    )

    # Deliberately non-covariant label-dependent interblock source.
    noncovariant = [np.zeros((4, 4), dtype=complex) for _ in factors]
    for lam in (0, 2):
        z_lam = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2))
        b = np.block(
            [
                [np.zeros((2, 2)), z_lam],
                [np.conjugate(z_lam.T), np.zeros((2, 2))],
            ]
        )
        noncovariant[lam] = b
    noncovariant_ops = [physical_twist_source(basis, b) for b in noncovariant]
    q_noncov, _ = curvature_on_zero_manifold(basis, factors, noncovariant_ops, zero_states)
    _, noncov_fit_residual = solve_common_generator(factors, noncovariant)

    print("two-block projected-Hubbard transport audit")
    print("composition basis order: (|1,0>, |0,1>) = (block A pair, block B pair)")
    print("random seed: 20260725")
    print(f"one-pair H0 nullity: {nullity} (expected 2)")
    print(f"first positive H0 eigenvalue: {evals[nullity]:.12g}")
    print()
    print("common covariant transport")
    print(f"  common-Sylvester residual: {common_fit_residual:.3e}")
    print(f"  ||T||: {la.norm(t_common):.3e}")
    print("  Q =")
    print(np.array2string(q_common, precision=12, suppress_small=True))
    print()
    print("common transport + block-preserving intrinsic source")
    print("  Q =")
    print(np.array2string(q_covariant, precision=12, suppress_small=True))
    print(f"  |Q_12|: {abs(q_covariant[0, 1]):.3e}")
    print()
    print("factor-label-dependent non-covariant source")
    print(f"  common-Sylvester residual: {noncov_fit_residual:.3e}")
    print("  Q =")
    print(np.array2string(q_noncov, precision=12, suppress_small=True))
    print(f"  |Q_12|: {abs(q_noncov[0, 1]):.3e}")

    assert nullity == 2
    assert common_fit_residual < 1e-9
    assert la.norm(t_common) < 1e-9
    assert abs(q_covariant[0, 1]) < 1e-9
    assert noncov_fit_residual > 1e-3
    assert abs(q_noncov[0, 1]) > 1e-6


if __name__ == "__main__":
    main()
