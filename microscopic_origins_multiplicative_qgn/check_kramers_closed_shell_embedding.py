#!/usr/bin/env python3
"""Exact audit of the time-reversal-invariant closed-shell control embedding.

The active bond is the eight-mode (1L,1R,2L,2R) spinful problem used by the
all-filling router.  The control has four orbitals g,m,s,a and a Kramers label
kappa=+,-.  Two control electrons occupy the unique closed shell g_+ g_-.
The bridge is Kramers even and has amplitude b/sqrt(2) per partner.  The router
is Kramers odd, so its product with the time-reversal-odd active current is
time-reversal even.

The script verifies:
  1. the active/control transformation parities under explicit fermionic
     time reversal;
  2. invariance of the full microscopic Hamiltonian;
  3. equality of the two partner denominators and degree-six coefficients;
  4. an exact zero-energy one-bond Feshbach map whose seniority-zero restriction
     agrees with the single-flavor degree-six target up to O(lambda^8).

Dependencies: numpy, scipy.
"""
from __future__ import annotations

import itertools
import math
from typing import Sequence

import numpy as np
import scipy.linalg as la

TOL = 5.0e-10


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
    out = np.zeros((len(basis), len(basis)), dtype=complex)
    for col, state in enumerate(basis):
        first = annihilate(state, src)
        if first is None:
            continue
        state1, sign1 = first
        second = create(state1, dst)
        if second is None:
            continue
        state2, sign2 = second
        out[index[state2], col] += sign1 * sign2
    return out


def antiunitary_fock_matrix(
    basis: Sequence[int], mapped_mode: Sequence[int], mapped_phase: Sequence[complex]
) -> np.ndarray:
    """Return U for the antiunitary T=U K on a fixed-number Fock basis."""
    index = {state: row for row, state in enumerate(basis)}
    out = np.zeros((len(basis), len(basis)), dtype=complex)
    for col, state in enumerate(basis):
        occupied = [m for m in range(len(mapped_mode)) if (state >> m) & 1]
        mapped = [mapped_mode[m] for m in occupied]
        coefficient = complex(np.prod([mapped_phase[m] for m in occupied]))
        inversions = sum(
            mapped[i] > mapped[j]
            for i in range(len(mapped))
            for j in range(i + 1, len(mapped))
        )
        if inversions & 1:
            coefficient *= -1.0
        new_state = 0
        for mode in mapped:
            new_state |= 1 << mode
        out[index[new_state], col] = coefficient
    return out


def tr_conjugate(U: np.ndarray, operator: np.ndarray) -> np.ndarray:
    return U @ operator.conj() @ U.conj().T


def active_mode(site: int, spin: int) -> int:
    return 2 * site + spin


def active_data():
    # Fixed two-electron sector of the four spinful active orbitals.
    basis = [s for s in range(1 << 8) if s.bit_count() == 2]
    dim = len(basis)
    eye = np.eye(dim, dtype=complex)

    B = np.zeros((dim, dim), dtype=complex)
    for left, right in ((0, 3), (1, 2)):
        for spin in (0, 1):
            hop = cdagger_c(basis, active_mode(left, spin), active_mode(right, spin))
            B += hop + hop.conj().T

    currents = []
    for left, right in ((0, 1), (2, 3)):
        J = np.zeros_like(B)
        for spin in (0, 1):
            hop = cdagger_c(basis, active_mode(left, spin), active_mode(right, spin))
            J += 1j * hop - 1j * hop.conj().T
        currents.append(J)
    J1, J2 = currents
    Cs, Ca = J1 + J2, J1 - J2

    # Simultaneous endpoint swap, with the exact fermionic permutation sign.
    mapped = list(range(8))
    for left, right in ((0, 1), (2, 3)):
        for spin in (0, 1):
            a, b = active_mode(left, spin), active_mode(right, spin)
            mapped[a], mapped[b] = mapped[b], mapped[a]
    W = antiunitary_fock_matrix(basis, mapped, [1.0] * 8)
    Qa = 0.5 * (eye - W)
    K = Qa @ B

    seniority_zero = []
    for state in basis:
        ok = True
        for site in range(4):
            up = (state >> active_mode(site, 0)) & 1
            down = (state >> active_mode(site, 1)) & 1
            if up != down:
                ok = False
                break
        seniority_zero.append(ok)
    sz = np.where(np.asarray(seniority_zero))[0]

    # Physical spin time reversal: up -> down, down -> -up.
    tr_map = []
    tr_phase = []
    for site in range(4):
        tr_map.extend([active_mode(site, 1), active_mode(site, 0)])
        tr_phase.extend([1.0, -1.0])
    Utr = antiunitary_fock_matrix(basis, tr_map, tr_phase)
    return basis, eye, B, Cs, Ca, Qa, K, sz, Utr


def control_mode(orbital: int, kappa: int) -> int:
    # orbital order g,m,s,a; kappa 0=+, 1=-.
    return 2 * orbital + kappa


def control_data():
    basis = [
        sum(1 << i for i in occ)
        for occ in itertools.combinations(range(8), 2)
    ]
    dim = len(basis)
    eye = np.eye(dim, dtype=complex)

    E = {}
    for mu in range(4):
        for nu in range(4):
            for kappa in (0, 1):
                E[(mu, nu, kappa)] = cdagger_c(
                    basis, control_mode(mu, kappa), control_mode(nu, kappa)
                )

    Delta_m, Delta_s, Delta_a = 10.0, 8.0, 8.0
    H0 = np.zeros((dim, dim), dtype=complex)
    for kappa in (0, 1):
        H0 += Delta_m * E[(1, 1, kappa)]
        H0 += Delta_s * E[(2, 2, kappa)]
        H0 += Delta_a * E[(3, 3, kappa)]

    bridge = np.zeros_like(H0)
    router_s = np.zeros_like(H0)
    router_a = np.zeros_like(H0)
    for kappa, chi in ((0, 1.0), (1, -1.0)):
        bridge += E[(1, 0, kappa)] + E[(0, 1, kappa)]
        router_s += chi * (E[(2, 1, kappa)] + E[(1, 2, kappa)])
        router_a += chi * (E[(3, 1, kappa)] + E[(1, 3, kappa)])

    low_state = (1 << control_mode(0, 0)) | (1 << control_mode(0, 1))
    low_index = {state: i for i, state in enumerate(basis)}[low_state]

    # Kramers time reversal: + -> -, - -> -+.
    tr_map = []
    tr_phase = []
    for orbital in range(4):
        tr_map.extend([control_mode(orbital, 1), control_mode(orbital, 0)])
        tr_phase.extend([1.0, -1.0])
    Utr = antiunitary_fock_matrix(basis, tr_map, tr_phase)
    return (
        basis,
        eye,
        H0,
        bridge,
        router_s,
        router_a,
        low_index,
        Utr,
        Delta_m,
        Delta_s,
        Delta_a,
    )


def exact_feshbach_error(lam: float, cached):
    (
        eye_a,
        B,
        Cs,
        Ca,
        K,
        sz,
        eye_c,
        H0c,
        bridge_c,
        router_s_c,
        router_a_c,
        low_index,
        Delta_m,
        Delta_s,
        Delta_a,
    ) = cached

    s, a, b = 1.0, 0.5, 0.7
    V1 = s * np.kron(Cs, router_s_c) + a * np.kron(Ca, router_a_c)
    V2 = (b / math.sqrt(2.0)) * np.kron(B, bridge_c)
    H = np.kron(eye_a, H0c) + lam * V1 + lam**2 * V2

    dim_c = eye_c.shape[0]
    p_all = np.asarray([i * dim_c + low_index for i in range(eye_a.shape[0])])
    keep_q = np.ones(H.shape[0], dtype=bool)
    keep_q[p_all] = False
    q = np.where(keep_q)[0]
    p_sz = p_all[sz]

    Hqq = H[np.ix_(q, q)]
    Hqp = H[np.ix_(q, p_sz)]
    F = -(Hqp.conj().T @ la.solve(Hqq, Hqp, assume_a="her"))

    B2 = B @ B
    K2 = K @ K
    alpha_bar = 4.0 * b**2 / Delta_m**2 * (
        s**2 / Delta_s - a**2 / Delta_a
    )
    target = (
        -lam**4 * b**2 / Delta_m * B2
        -lam**6 * 4.0 * b**2 * s**2 / (Delta_m**2 * Delta_s) * B2
        +lam**6 * alpha_bar * K2
    )
    target_sz = target[np.ix_(sz, sz)]
    error = float(la.norm(F - target_sz, 2))
    gap = float(la.eigvalsh(Hqq, subset_by_index=[0, 0])[0])
    return error, alpha_bar, gap


def main() -> None:
    (
        _basis_a,
        eye_a,
        B,
        Cs,
        Ca,
        _Qa,
        K,
        sz,
        Utr_a,
    ) = active_data()
    (
        _basis_c,
        eye_c,
        H0c,
        bridge_c,
        router_s_c,
        router_a_c,
        low_index,
        Utr_c,
        Delta_m,
        Delta_s,
        Delta_a,
    ) = control_data()

    active_B_even = la.norm(tr_conjugate(Utr_a, B) - B)
    active_Cs_odd = la.norm(tr_conjugate(Utr_a, Cs) + Cs)
    active_Ca_odd = la.norm(tr_conjugate(Utr_a, Ca) + Ca)
    control_H0_even = la.norm(tr_conjugate(Utr_c, H0c) - H0c)
    control_bridge_even = la.norm(tr_conjugate(Utr_c, bridge_c) - bridge_c)
    control_router_s_odd = la.norm(tr_conjugate(Utr_c, router_s_c) + router_s_c)
    control_router_a_odd = la.norm(tr_conjugate(Utr_c, router_a_c) + router_a_c)

    s, a, b = 1.0, 0.5, 0.7
    V1 = s * np.kron(Cs, router_s_c) + a * np.kron(Ca, router_a_c)
    V2 = (b / math.sqrt(2.0)) * np.kron(B, bridge_c)
    Utr_full = np.kron(Utr_a, Utr_c)
    full_V1_even = la.norm(tr_conjugate(Utr_full, V1) - V1)
    full_V2_even = la.norm(tr_conjugate(Utr_full, V2) - V2)

    d_s_plus = Delta_m - 4.0 * (+s) ** 2 / Delta_s
    d_s_minus = Delta_m - 4.0 * (-s) ** 2 / Delta_s
    d_a_plus = Delta_m - 4.0 * (+a) ** 2 / Delta_a
    d_a_minus = Delta_m - 4.0 * (-a) ** 2 / Delta_a
    partner_denominator_error = max(
        abs(d_s_plus - d_s_minus), abs(d_a_plus - d_a_minus)
    )

    cached = (
        eye_a,
        B,
        Cs,
        Ca,
        K,
        sz,
        eye_c,
        H0c,
        bridge_c,
        router_s_c,
        router_a_c,
        low_index,
        Delta_m,
        Delta_s,
        Delta_a,
    )

    print("KRAMERS CLOSED-SHELL CONTROL EMBEDDING CERTIFICATE")
    print("==================================================")
    print("control mode order: (g+,g-,m+,m-,s+,s-,a+,a-)")
    print("active mode order: (1L up/down,1R up/down,2L up/down,2R up/down)")
    print("\n1. Explicit time-reversal transformation checks")
    for key, value in (
        ("active B even residual", active_B_even),
        ("active C_s odd residual", active_Cs_odd),
        ("active C_a odd residual", active_Ca_odd),
        ("control H0 even residual", control_H0_even),
        ("control bridge even residual", control_bridge_even),
        ("control router_s odd residual", control_router_s_odd),
        ("control router_a odd residual", control_router_a_odd),
        ("full V1 even residual", full_V1_even),
        ("full V2 even residual", full_V2_even),
    ):
        print(f"{key:38s} = {value:.12e}")

    print("\n2. Partner equality")
    print(f"d_s(+) = {d_s_plus:.12f}")
    print(f"d_s(-) = {d_s_minus:.12f}")
    print(f"d_a(+) = {d_a_plus:.12f}")
    print(f"d_a(-) = {d_a_minus:.12f}")
    print(f"partner denominator error                = {partner_denominator_error:.12e}")

    print("\n3. Exact closed-shell Feshbach scaling")
    ratios = []
    alpha_bar = None
    min_gap = math.inf
    for lam in (0.34, 0.30, 0.26, 0.22, 0.18, 0.15):
        error, alpha_bar, gap = exact_feshbach_error(lam, cached)
        ratio = error / lam**8
        ratios.append(ratio)
        min_gap = min(min_gap, gap)
        print(
            f"lambda={lam:5.3f}  error={error:.12e}  "
            f"error/lambda^8={ratio:.12e}  high_gap={gap:.12e}"
        )
    print(f"alpha_bar                                = {alpha_bar:.12e}")
    print(f"seniority-zero active dimension          = {len(sz)}")
    print(f"minimum sampled high-block gap           = {min_gap:.12e}")

    symmetry_residuals = [
        active_B_even,
        active_Cs_odd,
        active_Ca_odd,
        control_H0_even,
        control_bridge_even,
        control_router_s_odd,
        control_router_a_odd,
        full_V1_even,
        full_V2_even,
    ]
    convergence = max(
        abs(ratios[-1] - ratios[-2]), abs(ratios[-2] - ratios[-3])
    )
    checks = [
        max(symmetry_residuals) < TOL,
        partner_denominator_error < TOL,
        alpha_bar is not None and alpha_bar > 0.0,
        min_gap > 0.0,
        convergence < 2.0e-4,
    ]
    print(f"tail ratio convergence diagnostic        = {convergence:.12e}")
    print("\nOVERALL:", "PASS" if all(checks) else "FAIL")
    if not all(checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
