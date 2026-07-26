#!/usr/bin/env python3
"""Certificates for the global translationally invariant downfolding theorem.

Scope: both global audits use active particle number two (one pair).  They test
translation invariance and the predicted overlap/scaling powers; they are not
all-filling or volume-uniform numerical certificates.

Checks:
1. A three-cell overlapping periodic lattice has an exact zero-energy Feshbach
   map F(t_B)=t_B^2 H_2+O(t_B^4), with F(+t_B)=F(-t_B).  This verifies exact
   second-order additivity and that the first overlap correction is even and
   fourth order.
2. A two-cell periodic lattice is diagonalized exactly.  The canonical
   Schrieffer--Wolff Hamiltonian, obtained by Löwdin orthonormalization of the
   projected low eigenspace, agrees on seniority zero with the weighted target
       lambda^6 alpha_bar sum_x K_x^2
   up to O(lambda^8) after the order-four and order-six B^2 matching.
3. The extracted second-order active operator is translation invariant.

Dependencies: numpy, scipy.
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Sequence

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla

TOL = 2.0e-10


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


def fixed_basis(nmodes: int, npart: int) -> list[int]:
    return [sum(1 << i for i in occ) for occ in itertools.combinations(range(nmodes), npart)]


def cdagger_c(basis: Sequence[int], dst: int, src: int) -> sp.csr_matrix:
    index = {state: row for row, state in enumerate(basis)}
    rows, cols, vals = [], [], []
    for col, state in enumerate(basis):
        first = annihilate(state, src)
        if first is None:
            continue
        state1, sign1 = first
        second = create(state1, dst)
        if second is None:
            continue
        state2, sign2 = second
        rows.append(index[state2])
        cols.append(col)
        vals.append(sign1 * sign2)
    return sp.csr_matrix((vals, (rows, cols)), shape=(len(basis), len(basis)), dtype=complex)


def induced_permutation(basis: Sequence[int], permutation: Sequence[int]) -> np.ndarray:
    index = {state: row for row, state in enumerate(basis)}
    out = np.zeros((len(basis), len(basis)), dtype=complex)
    for col, state in enumerate(basis):
        occupied = [m for m in range(len(permutation)) if (state >> m) & 1]
        mapped = [permutation[m] for m in occupied]
        inversions = sum(
            mapped[i] > mapped[j]
            for i in range(len(mapped))
            for j in range(i + 1, len(mapped))
        )
        new_state = 0
        for m in mapped:
            new_state |= 1 << m
        out[index[new_state], col] = -1.0 if inversions & 1 else 1.0
    return out


def mode(volume: int, block: int, cell: int, spin: int) -> int:
    return 2 * (block * volume + cell) + spin


@dataclass
class ActiveData:
    basis: list[int]
    eye: sp.csr_matrix
    B: list[sp.csr_matrix]
    Cs: list[sp.csr_matrix]
    Ca: list[sp.csr_matrix]


def active_operators(volume: int, npart: int = 2) -> ActiveData:
    basis = fixed_basis(4 * volume, npart)
    dim = len(basis)
    eye = sp.eye(dim, format="csr", dtype=complex)
    Bs, Css, Cas = [], [], []

    for x in range(volume):
        y = (x + 1) % volume
        B = sp.csr_matrix((dim, dim), dtype=complex)
        for spin in (0, 1):
            for dst_block, dst_cell, src_block, src_cell in (
                (0, y, 1, x),
                (0, x, 1, y),
            ):
                hop = cdagger_c(
                    basis,
                    mode(volume, dst_block, dst_cell, spin),
                    mode(volume, src_block, src_cell, spin),
                )
                B = B + hop + hop.getH()

        currents = []
        for block in (0, 1):
            J = sp.csr_matrix((dim, dim), dtype=complex)
            for spin in (0, 1):
                hop = cdagger_c(
                    basis,
                    mode(volume, block, x, spin),
                    mode(volume, block, y, spin),
                )
                J = J + 1j * hop - 1j * hop.getH()
            currents.append(J)

        Bs.append(B.tocsr())
        Css.append((currents[0] + currents[1]).tocsr())
        Cas.append((currents[0] - currents[1]).tocsr())

    return ActiveData(basis, eye, Bs, Css, Cas)


def local_control_matrix(i: int, j: int) -> sp.csr_matrix:
    out = np.zeros((4, 4), dtype=complex)
    out[i, j] = 1.0
    return sp.csr_matrix(out)


def embed_control(local: sp.csr_matrix, site: int, volume: int) -> sp.csr_matrix:
    out = sp.csr_matrix([[1.0 + 0.0j]])
    eye4 = sp.eye(4, format="csr", dtype=complex)
    for x in range(volume):
        out = sp.kron(out, local if x == site else eye4, format="csr")
    return out


def control_ops():
    E = {(i, j): local_control_matrix(i, j) for i in range(4) for j in range(4)}
    return E


def translation_matrix(active: ActiveData, volume: int) -> np.ndarray:
    perm = list(range(4 * volume))
    for block in (0, 1):
        for x in range(volume):
            for spin in (0, 1):
                perm[mode(volume, block, x, spin)] = mode(
                    volume, block, (x + 1) % volume, spin
                )
    return induced_permutation(active.basis, perm)


def resummed_global_feshbach(volume: int, t_bridge: float):
    active = active_operators(volume)
    dim_active = active.eye.shape[0]
    dim_control = 4**volume
    E = control_ops()

    delta_m, delta_s, delta_a = 10.0, 8.0, 8.0
    t_s, t_a = 1.0, 0.5

    shape = (dim_active * dim_control, dim_active * dim_control)
    H_control = sp.csr_matrix(shape, dtype=complex)
    V_bridge = sp.csr_matrix(shape, dtype=complex)

    for x in range(volume):
        e_m = embed_control(E[1, 1], x, volume)
        e_s = embed_control(E[2, 2], x, volume)
        e_a = embed_control(E[3, 3], x, volume)
        e_ms = embed_control(E[1, 2] + E[2, 1], x, volume)
        e_ma = embed_control(E[1, 3] + E[3, 1], x, volume)
        e_gm = embed_control(E[0, 1] + E[1, 0], x, volume)

        H_control += sp.kron(
            active.eye,
            delta_m * e_m + delta_s * e_s + delta_a * e_a,
            format="csr",
        )
        H_control += sp.kron(t_s * active.Cs[x], e_ms, format="csr")
        H_control += sp.kron(t_a * active.Ca[x], e_ma, format="csr")
        V_bridge += sp.kron(active.B[x], e_gm, format="csr")

    H = H_control + t_bridge * V_bridge
    p_indices = np.asarray([i * dim_control for i in range(dim_active)], dtype=int)
    mask = np.ones(dim_active * dim_control, dtype=bool)
    mask[p_indices] = False
    q_indices = np.where(mask)[0]

    H_qq = H[q_indices, :][:, q_indices].tocsc()
    H_qp = H[q_indices, :][:, p_indices].tocsc()
    H_pq = H[p_indices, :][:, q_indices].tocsc()
    lu = spla.splu(H_qq)
    feshbach = -(H_pq @ lu.solve(H_qp.toarray()))

    second = np.zeros((dim_active, dim_active), dtype=complex)
    for B, Cs, Ca in zip(active.B, active.Cs, active.Ca):
        D = (
            delta_m * active.eye
            - (t_s**2 / delta_s) * (Cs @ Cs)
            - (t_a**2 / delta_a) * (Ca @ Ca)
        ).toarray()
        Bd = B.toarray()
        second += -Bd @ la.solve(D, Bd, assume_a="her")

    Hc_qq = H_control[q_indices, :][:, q_indices]
    gap = float(
        spla.eigsh(Hc_qq, k=1, which="SA", return_eigenvectors=False, tol=1.0e-10)[0]
    )
    trans = translation_matrix(active, volume)
    translation_error = float(la.norm(trans @ second - second @ trans))
    return feshbach, second, gap, translation_error


def active_sw_data(volume: int = 2):
    active = active_operators(volume)
    dim = active.eye.shape[0]
    K = []
    seniority_zero = []

    for x in range(volume):
        y = (x + 1) % volume
        perm = list(range(4 * volume))
        for block in (0, 1):
            for spin in (0, 1):
                left = mode(volume, block, x, spin)
                right = mode(volume, block, y, spin)
                perm[left], perm[right] = perm[right], perm[left]
        W = induced_permutation(active.basis, perm)
        Q_a = 0.5 * (np.eye(dim) - W)
        K.append(Q_a @ active.B[x].toarray())

    for state in active.basis:
        ok = True
        for block in (0, 1):
            for x in range(volume):
                up = (state >> mode(volume, block, x, 0)) & 1
                down = (state >> mode(volume, block, x, 1)) & 1
                if up != down:
                    ok = False
        seniority_zero.append(ok)

    return active, K, np.asarray(seniority_zero, dtype=bool)


def exact_weighted_sw(lam: float, volume: int = 2):
    active, K, seniority_zero = active_sw_data(volume)
    dim_active = active.eye.shape[0]
    dim_control = 4**volume
    E = control_ops()

    delta_m, delta_s, delta_a = 10.0, 8.0, 8.0
    s, a, b = 1.0, 0.5, 0.7
    t_s, t_a, t_b = lam * s, lam * a, lam**2 * b
    beta4 = b**2 / delta_m
    beta6 = 4.0 * b**2 * s**2 / (delta_m**2 * delta_s)
    beta = lam**4 * beta4 + lam**6 * beta6

    shape = (dim_active * dim_control, dim_active * dim_control)
    H = sp.csr_matrix(shape, dtype=complex)
    identity_control = sp.eye(dim_control, format="csr", dtype=complex)

    for x in range(volume):
        e_m = embed_control(E[1, 1], x, volume)
        e_s = embed_control(E[2, 2], x, volume)
        e_a = embed_control(E[3, 3], x, volume)
        e_ms = embed_control(E[1, 2] + E[2, 1], x, volume)
        e_ma = embed_control(E[1, 3] + E[3, 1], x, volume)
        e_gm = embed_control(E[0, 1] + E[1, 0], x, volume)

        H += sp.kron(
            active.eye,
            delta_m * e_m + delta_s * e_s + delta_a * e_a,
            format="csr",
        )
        H += sp.kron(t_s * active.Cs[x], e_ms, format="csr")
        H += sp.kron(t_a * active.Ca[x], e_ma, format="csr")
        H += sp.kron(t_b * active.B[x], e_gm, format="csr")
        H += sp.kron(beta * (active.B[x] @ active.B[x]), identity_control, format="csr")

    dense = H.toarray()
    eigenvalues, eigenvectors = la.eigh(dense, subset_by_index=[0, dim_active - 1])

    p_indices = np.asarray([i * dim_control for i in range(dim_active)], dtype=int)
    projected = eigenvectors[p_indices, :]
    gram = projected.conj().T @ projected
    gram_values, gram_vectors = la.eigh(gram)
    gram_inverse_sqrt = gram_vectors @ np.diag(1.0 / np.sqrt(gram_values)) @ gram_vectors.conj().T
    low_unitary = projected @ gram_inverse_sqrt
    exact_heff = low_unitary @ np.diag(eigenvalues) @ low_unitary.conj().T

    alpha_bar = 4.0 * b**2 / delta_m**2 * (s**2 / delta_s - a**2 / delta_a)
    target = lam**6 * alpha_bar * sum(
        (operator @ operator for operator in K),
        start=np.zeros((dim_active, dim_active), dtype=complex),
    )

    sz = np.where(seniority_zero)[0]
    error = float(la.norm((exact_heff - target)[np.ix_(sz, sz)], 2))
    return error, alpha_bar, len(sz)


def main() -> None:
    print("GLOBAL TRANSLATION-INVARIANT LATTICE DOWNFOLDING CERTIFICATE")
    print("===========================================================")

    print("\n1. Overlapping three-cell Feshbach check")
    ratios = []
    translation_error = None
    gap = None
    for t_bridge in (0.20, 0.14, 0.10, 0.07, 0.05):
        exact, second, gap, translation_error = resummed_global_feshbach(3, t_bridge)
        error = float(la.norm(exact - t_bridge**2 * second, 2))
        ratio = error / t_bridge**4
        ratios.append(ratio)
        print(
            f"t_B={t_bridge:5.3f}  error={error:.12e}  "
            f"error/t_B^4={ratio:.12e}"
        )

    plus, _, _, _ = resummed_global_feshbach(3, 0.08)
    minus, _, _, _ = resummed_global_feshbach(3, -0.08)
    even_error = float(la.norm(plus - minus, 2))
    print(f"one-control high gap: {gap:.12e}")
    print(f"second-order translation commutator: {translation_error:.12e}")
    print(f"F(+t_B)-F(-t_B): {even_error:.12e}")

    print("\n2. Weighted exact Schrieffer--Wolff check")
    weighted_ratios = []
    alpha_bar = None
    sz_dim = None
    for lam in (0.34, 0.30, 0.26, 0.22, 0.18, 0.15):
        error, alpha_bar, sz_dim = exact_weighted_sw(lam)
        ratio = error / lam**8
        weighted_ratios.append(ratio)
        print(
            f"lambda={lam:5.3f}  error={error:.12e}  "
            f"error/lambda^8={ratio:.12e}"
        )
    print(f"alpha_bar: {alpha_bar:.12e}")
    print(f"seniority-zero test dimension: {sz_dim}")

    checks = [
        max(abs(ratios[-1] - ratios[-2]), abs(ratios[-2] - ratios[-3])) < 5.0e-5,
        even_error < TOL,
        translation_error < TOL,
        gap > 0.0,
        max(abs(weighted_ratios[-1] - weighted_ratios[-2]), abs(weighted_ratios[-2] - weighted_ratios[-3])) < 1.0e-5,
        alpha_bar > 0.0,
    ]
    print("\nOVERALL:", "PASS" if all(checks) else "FAIL")
    if not all(checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
