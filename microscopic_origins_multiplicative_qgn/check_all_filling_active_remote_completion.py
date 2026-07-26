#!/usr/bin/env python3
"""All-filling electron-only active--remote completion certificate.

This script checks four exact statements on the complete local Fock space of
four spinful active orbitals (1L,1R,2L,2R).

1. The ordinary crossed hopping B maps the complete seniority-zero space into
   one fixed seniority-two source space R=ran(B P0), for every local charge.
2. The block-current quadratures J1,J2 resolve simultaneous endpoint parity on
   that source space:
       J1^2 = J2^2 = I,   J1 J2 = W  on R.
   Equivalently C_s=J1+J2 and C_a=J1-J2 obey
       C_s^2 = 4 Q_s,     C_a^2 = 4 Q_a  on R.
3. A four-level one-electron control chain with ordinary two-electron exchange
   vertices B, C_s, C_a has an exact zero-energy Schur complement whose
   restriction to every seniority-zero filling is
       -t_B^2/d_s B^2 + t_B^2(1/d_s-1/d_a) K^2.
   Thus the reduced source is the full crossed B, uniformly over all fillings.
4. The charge polynomial accompanying a conventional molecular-shell B^2
   channel projects onto the product-AGP composition basis as an explicit
   diagonal selection potential.  The finite-size combinatorial formula is
   checked by exact enumeration.

Dependencies: numpy, scipy.
"""
from __future__ import annotations

import itertools
import math
from typing import Sequence

import numpy as np
import scipy.linalg as la

TOL = 2.0e-11


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


def mode(site: int, spin: int) -> int:
    return 2 * site + spin


def active_operators():
    basis = list(range(1 << 8))
    dim = len(basis)
    eye = np.eye(dim, dtype=complex)

    def cc(dst: int, src: int) -> np.ndarray:
        return cdagger_c(basis, dst, src)

    # Crossed hopping B: (1L,2R) and (1R,2L).
    B = np.zeros((dim, dim), dtype=complex)
    for left, right in ((0, 3), (1, 2)):
        for spin in (0, 1):
            B += cc(mode(left, spin), mode(right, spin))
            B += cc(mode(right, spin), mode(left, spin))

    # Block currents J_a = i(c_L^dag c_R-c_R^dag c_L).
    currents = []
    for left, right in ((0, 1), (2, 3)):
        J = np.zeros_like(B)
        for spin in (0, 1):
            J += 1j * cc(mode(left, spin), mode(right, spin))
            J -= 1j * cc(mode(right, spin), mode(left, spin))
        currents.append(J)
    J1, J2 = currents
    Cs, Ca = J1 + J2, J1 - J2

    # Individual and simultaneous endpoint swaps.
    swaps = []
    for left, right in ((0, 1), (2, 3)):
        perm = list(range(8))
        for spin in (0, 1):
            a, b = mode(left, spin), mode(right, spin)
            perm[a], perm[b] = perm[b], perm[a]
        swaps.append(induced_permutation(basis, perm))
    W1, W2 = swaps
    W = W1 @ W2
    Qs, Qa = 0.5 * (eye + W), 0.5 * (eye - W)

    # Seniority-zero and exactly-two-unpaired projectors.
    p0_diag = []
    p2_diag = []
    number_diag = []
    for state in basis:
        seniority = 0
        for site in range(4):
            up = (state >> mode(site, 0)) & 1
            dn = (state >> mode(site, 1)) & 1
            seniority += (up - dn) ** 2
        p0_diag.append(1.0 if seniority == 0 else 0.0)
        p2_diag.append(1.0 if seniority == 2 else 0.0)
        number_diag.append(float(state.bit_count()))
    P0 = np.diag(p0_diag)
    P2 = np.diag(p2_diag)
    N = np.diag(number_diag)

    return basis, eye, B, J1, J2, Cs, Ca, W1, W2, W, Qs, Qa, P0, P2, N


def source_range_projector(B: np.ndarray, P0: np.ndarray) -> np.ndarray:
    U, singular, _ = la.svd(B @ P0, full_matrices=False)
    rank = int(np.sum(singular > 1.0e-11))
    R = U[:, :rank]
    return R @ R.conj().T


def all_filling_local_test() -> dict[str, float]:
    (
        _basis,
        eye,
        B,
        J1,
        J2,
        Cs,
        Ca,
        W1,
        W2,
        W,
        Qs,
        Qa,
        P0,
        P2,
        _N,
    ) = active_operators()

    PR = source_range_projector(B, P0)
    K = Qa @ B

    errors = {
        "seniority_two_source_error": float(la.norm((eye - P2) @ B @ P0)),
        "B_swap_commutator": float(la.norm(B @ W - W @ B)),
        "low_J1_square_error": float(la.norm(P0 @ J1 @ J1 @ P0 - 4.0 * P0 @ (eye - W1) / 2.0 @ P0)),
        "low_J2_square_error": float(la.norm(P0 @ J2 @ J2 @ P0 - 4.0 * P0 @ (eye - W2) / 2.0 @ P0)),
        "low_cross_current_error": float(la.norm(P0 @ J1 @ J2 @ P0)),
        "range_J1_square_error": float(la.norm((J1 @ J1 - eye) @ B @ P0)),
        "range_J2_square_error": float(la.norm((J2 @ J2 - eye) @ B @ P0)),
        "range_cross_current_error": float(la.norm((J1 @ J2 - W) @ B @ P0)),
        "range_Cs_square_error": float(la.norm((Cs @ Cs - 4.0 * Qs) @ B @ P0)),
        "range_Ca_square_error": float(la.norm((Ca @ Ca - 4.0 * Qa) @ B @ P0)),
        "range_invariance_Cs2": float(la.norm((eye - PR) @ Cs @ Cs @ PR)),
        "range_invariance_Ca2": float(la.norm((eye - PR) @ Ca @ Ca @ PR)),
        "Cs_preserves_seniority_two_source": float(la.norm((eye - P2) @ Cs @ PR)),
        "Ca_preserves_seniority_two_source": float(la.norm((eye - P2) @ Ca @ PR)),
        "K_definition_error": float(la.norm(K - B @ Qa)),
    }

    # Exact nested active--remote/control Schur complement.
    Delta_m = 10.0
    Delta_s_remote = 8.0
    Delta_a_remote = 8.0
    t_s = 1.0
    t_a = 0.5
    t_B = 0.7

    D = (
        Delta_m * eye
        - (t_s**2 / Delta_s_remote) * (Cs @ Cs)
        - (t_a**2 / Delta_a_remote) * (Ca @ Ca)
    )
    eig_D = la.eigvalsh(D)
    d_s = Delta_m - 4.0 * t_s**2 / Delta_s_remote
    d_a = Delta_m - 4.0 * t_a**2 / Delta_a_remote
    alpha = t_B**2 * (1.0 / d_s - 1.0 / d_a)

    exact = -t_B**2 * P0 @ B @ la.inv(D) @ B @ P0
    formula = P0 @ (
        -(t_B**2 / d_s) * (B @ B)
        + alpha * (K @ K)
    ) @ P0
    errors["control_high_block_min_eigenvalue"] = float(eig_D[0])
    errors["exact_all_filling_schur_error"] = float(la.norm(exact - formula))
    errors["effective_even_denominator"] = d_s
    errors["effective_odd_denominator"] = d_a
    errors["positive_K2_coefficient"] = alpha

    # Direct B^2 compensation and its controlled feedback into the high block.
    beta = t_B**2 / d_s
    E = beta * (B @ B)
    D_beta = D + E
    exact_compensated = P0 @ (E - t_B**2 * B @ la.inv(D_beta) @ B) @ P0
    target_compensated = P0 @ (alpha * K @ K) @ P0
    feedback_error = float(la.norm(exact_compensated - target_compensated))
    d_min = float(eig_D[0])
    e_norm = float(la.norm(E, 2))
    b_norm = float(la.norm(B, 2))
    feedback_bound = (
        t_B**2 * b_norm**2 * e_norm / (d_min * (d_min - e_norm))
        if e_norm < d_min
        else math.inf
    )
    errors["B_operator_norm"] = b_norm
    errors["B2_feedback_error"] = feedback_error
    errors["B2_feedback_bound"] = feedback_bound

    # Scaling check: the uncompensated high-block feedback is fourth order in t_B
    # when beta=t_B^2/d_s.
    scaling = []
    for trial_t in (0.70, 0.50, 0.35, 0.25, 0.175):
        trial_beta = trial_t**2 / d_s
        trial_alpha = trial_t**2 * (1.0 / d_s - 1.0 / d_a)
        trial_E = trial_beta * (B @ B)
        trial_exact = P0 @ (
            trial_E - trial_t**2 * B @ la.inv(D + trial_E) @ B
        ) @ P0
        trial_target = P0 @ (trial_alpha * K @ K) @ P0
        trial_error = float(la.norm(trial_exact - trial_target))
        scaling.append((trial_t, trial_error, trial_error / trial_t**4))
    errors["B2_feedback_scaling"] = scaling

    return errors


def product_dicke_vector(volume: int, n1: int, n2: int):
    basis = []
    amplitudes = []
    # Hard-core pair bitstrings: first V bits block 1, next V block 2.
    for occ1 in itertools.combinations(range(volume), n1):
        mask1 = sum(1 << x for x in occ1)
        for occ2 in itertools.combinations(range(volume), n2):
            mask2 = sum(1 << (volume + x) for x in occ2)
            basis.append(mask1 | mask2)
            amplitudes.append(1.0)
    vec = np.asarray(amplitudes, dtype=float)
    vec /= la.norm(vec)
    return basis, vec


def charge_energy_for_state(volume: int, state: int, A: float, U: float) -> float:
    total = 0.0
    for x in range(volume):
        y = (x + 1) % volume
        pairs = 0
        for block in (0, 1):
            pairs += (state >> (block * volume + x)) & 1
            pairs += (state >> (block * volume + y)) & 1
        electron_charge = 2.0 * pairs
        total += A * electron_charge**2 - 0.5 * U * electron_charge
    return total


def charge_formula(volume: int, n1: int, n2: int, A: float, U: float) -> float:
    n = n1 + n2
    quadratic = (
        8.0 * n
        + 8.0 * (n1 * (n1 - 1) + n2 * (n2 - 1)) / (volume - 1)
        + 32.0 * n1 * n2 / volume
    )
    return A * quadratic - 2.0 * U * n


def charge_projection_test() -> dict[str, float]:
    U = 2.0
    V_m = 1.25
    A = 0.25 * (U + V_m)
    beta = 0.25 * (U - V_m)
    max_error = 0.0
    max_rearranged_error = 0.0
    examples = []

    for volume in range(2, 7):
        for n1 in range(volume + 1):
            for n2 in range(volume + 1):
                basis, vec = product_dicke_vector(volume, n1, n2)
                diagonal = np.asarray(
                    [charge_energy_for_state(volume, state, A, U) for state in basis]
                )
                exact = float(np.sum(np.abs(vec) ** 2 * diagonal))
                formula = charge_formula(volume, n1, n2, A, U)
                max_error = max(max_error, abs(exact - formula))

                n = n1 + n2
                constant = A * (
                    8.0 * n + 8.0 * (n * n - n) / (volume - 1)
                ) - 2.0 * U * n
                selector = 16.0 * A * (volume - 2) * n1 * n2 / (
                    volume * (volume - 1)
                )
                max_rearranged_error = max(
                    max_rearranged_error, abs(formula - constant - selector)
                )
                if volume == 5 and (n1, n2) in ((4, 0), (2, 2), (1, 3)):
                    examples.append((n1, n2, exact, selector))

    return {
        "charge_projection_max_error": max_error,
        "charge_rearrangement_max_error": max_rearranged_error,
        "charge_channel_A": A,
        "orbital_channel_beta": beta,
        "example_count": float(len(examples)),
        "example_selector_min": min(v[3] for v in examples),
        "example_selector_max": max(v[3] for v in examples),
    }


def main() -> None:
    local = all_filling_local_test()
    charge = charge_projection_test()

    print("ALL-FILLING ACTIVE--REMOTE / CHARGE-CHANNEL CERTIFICATE")
    print("=======================================================")
    print("\n1. Exact local source and parity-routing identities")
    for key in (
        "seniority_two_source_error",
        "B_swap_commutator",
        "low_J1_square_error",
        "low_J2_square_error",
        "low_cross_current_error",
        "range_J1_square_error",
        "range_J2_square_error",
        "range_cross_current_error",
        "range_Cs_square_error",
        "range_Ca_square_error",
        "range_invariance_Cs2",
        "range_invariance_Ca2",
        "Cs_preserves_seniority_two_source",
        "Ca_preserves_seniority_two_source",
        "K_definition_error",
    ):
        print(f"{key}: {local[key]:.12e}")

    print("\n2. Exact nested control-chain Schur complement")
    for key in (
        "control_high_block_min_eigenvalue",
        "effective_even_denominator",
        "effective_odd_denominator",
        "positive_K2_coefficient",
        "exact_all_filling_schur_error",
    ):
        print(f"{key}: {local[key]:.12e}")

    print("\n3. Direct B^2 compensation feedback")
    for key in (
        "B_operator_norm",
        "B2_feedback_error",
        "B2_feedback_bound",
    ):
        print(f"{key}: {local[key]:.12e}")
    print("t_B     feedback_error      error/t_B^4")
    for trial_t, trial_error, ratio in local["B2_feedback_scaling"]:
        print(f"{trial_t:5.3f}   {trial_error:.12e}   {ratio:.12e}")

    print("\n4. Product-AGP projection of Coulomb charge terms")
    for key, value in charge.items():
        print(f"{key}: {value:.12e}")

    identity_keys = [
        key for key in local
        if key.endswith("error") and key != "B2_feedback_error"
    ]
    checks = [local[key] < TOL for key in identity_keys]
    checks += [
        local["control_high_block_min_eigenvalue"] > 0.0,
        local["effective_odd_denominator"] > local["effective_even_denominator"],
        local["positive_K2_coefficient"] > 0.0,
        local["B2_feedback_error"] <= local["B2_feedback_bound"] + 1.0e-12,
        charge["charge_projection_max_error"] < TOL,
        charge["charge_rearrangement_max_error"] < TOL,
    ]
    print("\nOVERALL:", "PASS" if all(checks) else "FAIL")
    if not all(checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
