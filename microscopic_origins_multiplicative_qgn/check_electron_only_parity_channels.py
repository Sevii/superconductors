#!/usr/bin/env python3
"""Electron-only microscopic channels for the multiplicative null-row program.

Checks four exact statements.

1. Fluxed attractive-Hubbard dimer: ordinary one-electron hopping couples the
   low empty/doublon manifold to only one bond-pair parity.  At phi=pi/2 the
   bright pair is the odd bond-parity state.  The second-order effective
   operator is a rank-one parity projector.

2. Four-site attractive-Hubbard plaquette: ordinary hopping, a seniority
   penalty, and density repulsion generate an isolated even/odd simultaneous
   bond-parity doublet.  Its splitting is fourth order and agrees with an
   exact path-sum coefficient.

3. Crossed molecular shells: the crossed one-body transfer B is exactly the
   occupation imbalance N_+-N_- of bonding/antibonding cluster orbitals.

4. A conventional constant-interaction multiorbital Hamiltonian with
   repulsions U (same shell) and V (opposite shell) contains the positive
   channel ((U-V)/4) B^2, plus a function only of the total cluster charge.

Dependencies: numpy, scipy.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from itertools import combinations
from typing import Sequence

import numpy as np
import scipy.linalg as la

TOL = 5e-11


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


def cdagger_c(basis: Sequence[int], dst: int, src: int) -> np.ndarray:
    index = {state: row for row, state in enumerate(basis)}
    M = np.zeros((len(basis), len(basis)), dtype=complex)
    for col, state in enumerate(basis):
        first = annihilate(state, src)
        if first is None:
            continue
        state1, sign1 = first
        second = create(state1, dst)
        if second is None:
            continue
        state2, sign2 = second
        if state2 in index:
            M[index[state2], col] += sign1 * sign2
    return M


def ket_from_modes(basis: Sequence[int], modes: Sequence[int]) -> np.ndarray:
    state = sum(1 << m for m in modes)
    v = np.zeros(len(basis), dtype=complex)
    v[{state: i for i, state in enumerate(basis)}[state]] = 1.0
    return v


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
        new_state = 0
        for m in mapped:
            new_state |= 1 << m
        M[index[new_state], col] = -1.0 if inversions & 1 else 1.0
    return M


def dimer_test() -> dict[str, float]:
    # Modes L up, L down, R up, R down.
    basis = [s for s in range(16) if s.bit_count() == 2]
    DL = ket_from_modes(basis, [0, 1])
    DR = ket_from_modes(basis, [2, 3])
    Lup_Rdn = ket_from_modes(basis, [0, 3])
    Ldn_Rup = ket_from_modes(basis, [1, 2])
    singlet = (Lup_Rdn - Ldn_Rup) / math.sqrt(2.0)
    P = np.column_stack([DL, DR])

    max_bright_error = 0.0
    max_sw_error = 0.0
    t = 0.73
    gap = 4.2
    for phi in (0.0, math.pi / 7.0, math.pi / 4.0, math.pi / 2.0):
        T = t * (
            np.exp(1j * phi) * (cdagger_c(basis, 0, 2) + cdagger_c(basis, 1, 3))
            + np.exp(-1j * phi) * (cdagger_c(basis, 2, 0) + cdagger_c(basis, 3, 1))
        )
        pair_phi = (np.exp(1j * phi) * DL + np.exp(-1j * phi) * DR) / math.sqrt(2.0)
        target = 2.0 * t * np.outer(singlet, pair_phi.conj())
        err = la.norm(T @ P @ P.conj().T - target)
        max_bright_error = max(max_bright_error, float(err))

        sw_direct = -(P.conj().T @ T.conj().T @ np.outer(singlet, singlet.conj()) @ T @ P) / gap
        coeff = P.conj().T @ pair_phi
        sw_formula = -(4.0 * t * t / gap) * np.outer(coeff, coeff.conj())
        max_sw_error = max(max_sw_error, float(la.norm(sw_direct - sw_formula)))

    odd = (DL - DR) / math.sqrt(2.0)
    phi = math.pi / 2.0
    pair_phi = (np.exp(1j * phi) * DL + np.exp(-1j * phi) * DR) / math.sqrt(2.0)
    odd_overlap = abs(np.vdot(odd, pair_phi))
    return {
        "bright_identity_error": max_bright_error,
        "sw_projector_error": max_sw_error,
        "odd_overlap_at_pi_over_2": float(odd_overlap),
    }


def site_mode(site: int, spin: int) -> int:
    return 2 * site + spin


@dataclass
class PlaquetteResult:
    coefficient_exact: float
    coefficient_feshbach: float
    coefficient_identity_error: float
    table: list[tuple[float, float, float, float, float, float]]


def plaquette_hamiltonian(t: float, Delta: float, V: float):
    # Four sites ordered 1L,1R,2L,2R; N=4 electrons (two local pairs).
    basis = [s for s in range(1 << 8) if s.bit_count() == 4]
    index = {s: i for i, s in enumerate(basis)}
    H = np.zeros((len(basis), len(basis)), dtype=complex)
    for i, state in enumerate(basis):
        q_count = 0
        n = []
        for site in range(4):
            up = (state >> site_mode(site, 0)) & 1
            dn = (state >> site_mode(site, 1)) & 1
            q_count += (up - dn) ** 2
            n.append(up + dn)
        # Cycle density repulsion leaves only the two diagonal pair configurations low.
        H[i, i] = (Delta / 2.0) * q_count + (V / 4.0) * (n[0] * n[1] + n[1] * n[2] + n[2] * n[3] + n[3] * n[0])
    for left, right in ((0, 1), (2, 3)):
        for spin in (0, 1):
            H -= t * (
                cdagger_c(basis, site_mode(left, spin), site_mode(right, spin))
                + cdagger_c(basis, site_mode(right, spin), site_mode(left, spin))
            )

    perm = list(range(8))
    for left, right in ((0, 1), (2, 3)):
        for spin in (0, 1):
            a, b = site_mode(left, spin), site_mode(right, spin)
            perm[a], perm[b] = perm[b], perm[a]
    W = induced_permutation(basis, perm)

    def pair_config(sites: Sequence[int]) -> np.ndarray:
        modes = []
        for site in sites:
            modes.extend([site_mode(site, 0), site_mode(site, 1)])
        return ket_from_modes(basis, modes)

    A = pair_config([0, 2])
    Bstate = pair_config([1, 3])
    return basis, H, W, A, Bstate


def plaquette_test() -> PlaquetteResult:
    Delta = 10.0
    V = 3.0
    # Exact t^4 coefficient for E_odd-E_even from the zero-energy Feshbach path sum.
    coefficient = 256.0 * (2.0 * Delta + 3.0 * V) / (
        V * (2.0 * Delta + V) * (4.0 * Delta + 3.0 * V) ** 2
    )
    # Independent fourth-order zero-energy Feshbach path sum.
    basis0, H0, _W0, A0, B0 = plaquette_hamiltonian(0.0, Delta, V)
    # Unit hopping operator: H(t)=H0+t*Tunit.
    _basis1, H1, _W1, _A1, _B1 = plaquette_hamiltonian(1.0, Delta, V)
    Tunit = H1 - H0
    ia = int(np.argmax(np.abs(A0)))
    ib = int(np.argmax(np.abs(B0)))
    qinds = [i for i in range(len(basis0)) if i not in (ia, ib)]
    HQ = H0[np.ix_(qinds, qinds)]
    R = la.inv(HQ)
    TQP = Tunit[np.ix_(qinds, [ia, ib])]
    TQQ = Tunit[np.ix_(qinds, qinds)]
    F4 = -(TQP.conj().T @ R @ TQQ @ R @ TQQ @ R @ TQP)
    coefficient_feshbach = float(2.0 * abs(F4[0, 1]))
    identity_error = abs(coefficient_feshbach - coefficient)

    table = []
    for t in (0.30, 0.20, 0.15, 0.10, 0.07, 0.05, 0.03, 0.02):
        basis, H, W, A, Bstate = plaquette_hamiltonian(t, Delta, V)
        vals, vecs = la.eigh(H)
        low_basis = np.column_stack([A, Bstate])
        weights = np.sum(np.abs(low_basis.conj().T @ vecs) ** 2, axis=0)
        candidates = np.argsort(weights)[-2:]
        states = []
        for idx in candidates:
            parity = float(np.real(vecs[:, idx].conj() @ W @ vecs[:, idx]))
            states.append((float(vals[idx]), parity, float(weights[idx])))
        even = max(states, key=lambda item: item[1])
        odd = min(states, key=lambda item: item[1])
        split = odd[0] - even[0]
        table.append((t, even[0], odd[0], split / t**4, even[2], odd[2]))
    return PlaquetteResult(
        coefficient_exact=coefficient,
        coefficient_feshbach=coefficient_feshbach,
        coefficient_identity_error=identity_error,
        table=table,
    )


def full_fock_cdagger_c(n_modes: int, dst: int, src: int) -> np.ndarray:
    basis = list(range(1 << n_modes))
    return cdagger_c(basis, dst, src)


def shell_channel_test() -> dict[str, float]:
    # Original modes: four sites x spin, ordered 1L,1R,2L,2R.
    n_modes = 8
    basis = list(range(1 << n_modes))
    dim = len(basis)
    I = np.eye(dim, dtype=complex)

    # Crossed edges e=(1L,2R), f=(1R,2L).
    B = np.zeros((dim, dim), dtype=complex)
    N_plus = np.zeros_like(B)
    N_minus = np.zeros_like(B)
    for i, j in ((0, 3), (1, 2)):
        for spin in (0, 1):
            mi, mj = site_mode(i, spin), site_mode(j, spin)
            B += cdagger_c(basis, mi, mj) + cdagger_c(basis, mj, mi)
            ni = cdagger_c(basis, mi, mi)
            nj = cdagger_c(basis, mj, mj)
            tij = cdagger_c(basis, mi, mj) + cdagger_c(basis, mj, mi)
            N_plus += 0.5 * (ni + nj + tij)
            N_minus += 0.5 * (ni + nj - tij)
    N = N_plus + N_minus

    imbalance_error = la.norm(B - (N_plus - N_minus))
    number_error = la.norm(N - np.diag([state.bit_count() for state in basis]))

    U = 2.0
    V = 1.25
    H_shell = 0.5 * U * (N_plus @ (N_plus - I) + N_minus @ (N_minus - I)) + V * (N_plus @ N_minus)
    beta = 0.25 * (U - V)
    formula = 0.25 * (U + V) * (N @ N) - 0.5 * U * N + beta * (B @ B)
    decomposition_error = la.norm(H_shell - formula)

    fixed_sector_error = 0.0
    for charge in range(n_modes + 1):
        inds = [i for i, state in enumerate(basis) if state.bit_count() == charge]
        if not inds:
            continue
        block = np.ix_(inds, inds)
        constant = 0.25 * (U + V) * charge**2 - 0.5 * U * charge
        fixed_sector_error = max(
            fixed_sector_error,
            float(la.norm(H_shell[block] - constant * np.eye(len(inds)) - beta * (B @ B)[block])),
        )

    # A concrete matching example: beta=t_c^2/Delta_s with all repulsions positive.
    tc = 0.8
    Delta_s = 12.0
    beta_target = tc**2 / Delta_s
    U_match = 2.0
    V_match = U_match - 4.0 * beta_target
    return {
        "imbalance_identity_error": float(imbalance_error),
        "number_identity_error": float(number_error),
        "shell_decomposition_error": float(decomposition_error),
        "fixed_charge_reduction_error": float(fixed_sector_error),
        "beta_example": beta,
        "beta_target": beta_target,
        "U_match": U_match,
        "V_match": V_match,
    }



def bilayer_one_pair_test() -> dict[str, object]:
    """Ordinary vertical hopping couples an active pair dimer to a remote parity dimer.

    The remote dimer's even/odd pair states are split by ordinary remote
    one-electron hopping.  The active parity splitting first appears at fourth
    order in the vertical hopping and is computed both from the exact Feshbach
    path sum and full exact diagonalization.
    """
    basis = [s for s in range(1 << 8) if s.bit_count() == 2]
    index = {state: i for i, state in enumerate(basis)}

    def pair(site: int) -> np.ndarray:
        return ket_from_modes(basis, [site_mode(site, 0), site_mode(site, 1)])

    AL, AR, RL, RR = [pair(site) for site in range(4)]
    P = np.column_stack([AL, AR])

    perm = list(range(8))
    for left, right in ((0, 1), (2, 3)):
        for spin in (0, 1):
            a, b = site_mode(left, spin), site_mode(right, spin)
            perm[a], perm[b] = perm[b], perm[a]
    W = induced_permutation(basis, perm)

    Delta = 10.0
    eps = 4.0
    r = 0.8

    def build(tau: float) -> np.ndarray:
        H = np.zeros((len(basis), len(basis)), dtype=complex)
        for i, state in enumerate(basis):
            energy = 0.0
            for site in range(4):
                up = (state >> site_mode(site, 0)) & 1
                dn = (state >> site_mode(site, 1)) & 1
                energy += (Delta / 2.0) * (up - dn) ** 2
                if site >= 2:
                    energy += eps * (up + dn)
            H[i, i] = energy
        for spin in (0, 1):
            H -= r * (
                cdagger_c(basis, site_mode(2, spin), site_mode(3, spin))
                + cdagger_c(basis, site_mode(3, spin), site_mode(2, spin))
            )
        for active, remote in ((0, 2), (1, 3)):
            for spin in (0, 1):
                H -= tau * (
                    cdagger_c(basis, site_mode(active, spin), site_mode(remote, spin))
                    + cdagger_c(basis, site_mode(remote, spin), site_mode(active, spin))
                )
        return H

    H0 = build(0.0)
    vals0, vecs0 = la.eigh(H0)
    remote_even = (RL + RR) / math.sqrt(2.0)
    remote_odd = (RL - RR) / math.sqrt(2.0)
    weights_even = np.abs(remote_even.conj() @ vecs0) ** 2
    weights_odd = np.abs(remote_odd.conj() @ vecs0) ** 2
    idx_even = int(np.argmax(weights_even))
    idx_odd = int(np.argmax(weights_odd))
    Delta_s = float(vals0[idx_even])
    Delta_a = float(vals0[idx_odd])

    # Vertical hopping operator with unit coefficient.
    Vop = np.zeros_like(H0)
    for active, remote in ((0, 2), (1, 3)):
        for spin in (0, 1):
            Vop -= (
                cdagger_c(basis, site_mode(active, spin), site_mode(remote, spin))
                + cdagger_c(basis, site_mode(remote, spin), site_mode(active, spin))
            )

    pinds = [int(np.argmax(AL)), int(np.argmax(AR))]
    qinds = [i for i in range(len(basis)) if i not in pinds]
    HQ = H0[np.ix_(qinds, qinds)]
    R = la.inv(HQ)
    VQP = Vop[np.ix_(qinds, pinds)]
    VQQ = Vop[np.ix_(qinds, qinds)]
    F4 = -(VQP.conj().T @ R @ VQQ @ R @ VQQ @ R @ VQP)
    split_coefficient = float(2.0 * abs(F4[0, 1]))

    table = []
    for tau in (0.50, 0.35, 0.25, 0.15, 0.10, 0.07, 0.05):
        vals, vecs = la.eigh(build(tau))
        low_weights = np.sum(np.abs(P.conj().T @ vecs) ** 2, axis=0)
        candidates = np.argsort(low_weights)[-2:]
        states = []
        for idx in candidates:
            parity = float(np.real(vecs[:, idx].conj() @ W @ vecs[:, idx]))
            states.append((float(vals[idx]), parity, float(low_weights[idx])))
        even = max(states, key=lambda item: item[1])
        odd = min(states, key=lambda item: item[1])
        split = odd[0] - even[0]
        table.append((tau, split / tau**4, even[2], odd[2]))

    return {
        "remote_even_energy": Delta_s,
        "remote_odd_energy": Delta_a,
        "remote_odd_minus_even": Delta_a - Delta_s,
        "feshbach_tau4_split_coefficient": split_coefficient,
        "table": table,
    }


def main() -> None:
    dimer = dimer_test()
    plaquette = plaquette_test()
    shell = shell_channel_test()
    bilayer = bilayer_one_pair_test()

    print("ELECTRON-ONLY PARITY / B^2 MICROSCOPIC CERTIFICATE")
    print("===================================================")
    print("\n1. Fluxed attractive-Hubbard dimer")
    for key, value in dimer.items():
        print(f"{key}: {value:.12e}")

    print("\n2. Four-site parity plaquette")
    print(f"closed fourth-order split coefficient: {plaquette.coefficient_exact:.12f}")
    print(f"independent Feshbach coefficient: {plaquette.coefficient_feshbach:.12f}")
    print(f"coefficient identity error: {plaquette.coefficient_identity_error:.12e}")
    print("t       E_even          E_odd           split/t^4       overlap_even  overlap_odd")
    for row in plaquette.table:
        print(f"{row[0]:.3f}  {row[1]: .12e}  {row[2]: .12e}  {row[3]: .12f}  {row[4]:.9f}  {row[5]:.9f}")
    last_error = abs(plaquette.table[-1][3] - plaquette.coefficient_exact)
    print(f"smallest-t coefficient error: {last_error:.6e}")

    print("\n3. One-pair active--remote bilayer proof of principle")
    for key in ("remote_even_energy", "remote_odd_energy", "remote_odd_minus_even", "feshbach_tau4_split_coefficient"):
        print(f"{key}: {bilayer[key]:.12e}")
    print("tau     split/tau^4       overlap_even  overlap_odd")
    for row in bilayer["table"]:
        print(f"{row[0]:.3f}   {row[1]:.12e}   {row[2]:.9f}  {row[3]:.9f}")
    bilayer_error = abs(bilayer["table"][-1][1] - bilayer["feshbach_tau4_split_coefficient"])
    print(f"smallest-t bilayer coefficient error: {bilayer_error:.6e}")

    print("\n4. Crossed molecular-shell Coulomb channel")
    for key, value in shell.items():
        print(f"{key}: {value:.12e}")

    checks = [
        dimer["bright_identity_error"] < TOL,
        dimer["sw_projector_error"] < TOL,
        abs(dimer["odd_overlap_at_pi_over_2"] - 1.0) < TOL,
        shell["imbalance_identity_error"] < TOL,
        shell["number_identity_error"] < TOL,
        shell["shell_decomposition_error"] < TOL,
        shell["fixed_charge_reduction_error"] < TOL,
        plaquette.coefficient_identity_error < TOL,
        last_error < 2e-5,
        shell["V_match"] > 0.0,
        bilayer["remote_odd_minus_even"] > 0.0,
        bilayer_error < 3e-7,
    ]
    print("\nOVERALL:", "PASS" if all(checks) else "FAIL")
    if not all(checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
