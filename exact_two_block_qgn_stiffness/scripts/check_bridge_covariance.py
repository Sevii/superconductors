#!/usr/bin/env python3
"""Audit bridge twist covariance and the conditional two-block Jacobi operator.

This script separates two logically distinct calculations:

1. A fixed multi-site bridge kernel B(k) with only the active projector twisted.
   This reproduces the nonzero derivative used in the superseded verifier.
2. A covariantly Peierls-twisted bridge kernel B(k+A), twisted together with the
   active projector.  The projected zero row then vanishes for every A and its
   derivative cancels to numerical precision.

It also checks the displacement-balance criterion for four-fermion vertices,
the exact transfer-dark zero mode of the bridge-only Jacobi matrix, and the
large-volume middle-filling stiffness formula quoted in the revised draft.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np

Array = np.ndarray


def haar_filter(z: complex) -> Array:
    return 0.5 * np.array(
        [[1.0 + z, 1.0 - z], [1.0 - z, 1.0 + z]], dtype=complex
    )


D_PHASE = np.diag([1.0, 1.0j])


def frame(kx: float, ky: float, block: int) -> Array:
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


def projected_bridge_matrix(
    kx: float,
    ky: float,
    ay: float,
    quadrature: str,
    *,
    covariant_kernel: bool,
) -> Array:
    """Project either a fixed or a covariantly twisted bridge kernel."""
    u1_a = active(kx, ky + ay, 1)
    u2_a = active(kx, ky + ay, 2)

    kernel_shift = ay if covariant_kernel else 0.0
    r1_b = bright(kx, ky + kernel_shift, 1)
    u2_b = active(kx, ky + kernel_shift, 2)
    m12 = np.outer(r1_b, np.conjugate(u2_b))

    coefficient_12 = np.vdot(u1_a, m12 @ u2_a)
    coefficient_21 = np.conjugate(coefficient_12)
    if quadrature == "X":
        return np.array(
            [[0.0, coefficient_12], [coefficient_21, 0.0]], dtype=complex
        )
    if quadrature == "Y":
        return np.array(
            [[0.0, -1.0j * coefficient_12], [1.0j * coefficient_21, 0.0]],
            dtype=complex,
        )
    raise ValueError("quadrature must be X or Y")


def bridge_covariance_audit(step: float = 1.0e-6, grid: int = 19) -> dict[str, float]:
    expected_fixed_x = 0.5 * np.array(
        [[0.0, 1.0j], [-1.0j, 0.0]], dtype=complex
    )
    expected_fixed_y = 0.5 * np.array(
        [[0.0, 1.0], [1.0, 0.0]], dtype=complex
    )
    fixed_error_x = 0.0
    fixed_error_y = 0.0
    covariant_derivative_max = 0.0
    covariant_value_max = 0.0

    for ix in range(grid):
        kx = 2.0 * np.pi * ix / grid
        for iy in range(grid):
            ky = 2.0 * np.pi * iy / grid
            for quadrature, expected in (
                ("X", expected_fixed_x),
                ("Y", expected_fixed_y),
            ):
                fixed = (
                    projected_bridge_matrix(
                        kx, ky, step, quadrature, covariant_kernel=False
                    )
                    - projected_bridge_matrix(
                        kx, ky, -step, quadrature, covariant_kernel=False
                    )
                ) / (2.0 * step)
                error = float(np.max(np.abs(fixed - expected)))
                if quadrature == "X":
                    fixed_error_x = max(fixed_error_x, error)
                else:
                    fixed_error_y = max(fixed_error_y, error)

                plus = projected_bridge_matrix(
                    kx, ky, step, quadrature, covariant_kernel=True
                )
                minus = projected_bridge_matrix(
                    kx, ky, -step, quadrature, covariant_kernel=True
                )
                zero = projected_bridge_matrix(
                    kx, ky, 0.0, quadrature, covariant_kernel=True
                )
                derivative = (plus - minus) / (2.0 * step)
                covariant_derivative_max = max(
                    covariant_derivative_max, float(np.max(np.abs(derivative)))
                )
                covariant_value_max = max(
                    covariant_value_max,
                    float(np.max(np.abs(plus))),
                    float(np.max(np.abs(minus))),
                    float(np.max(np.abs(zero))),
                )

    return {
        "fixed_kernel_derivative_error_X": fixed_error_x,
        "fixed_kernel_derivative_error_Y": fixed_error_y,
        "covariant_projected_value_max": covariant_value_max,
        "covariant_projected_derivative_max": covariant_derivative_max,
    }


def net_displacement(
    created: Iterable[tuple[float, float]],
    annihilated: Iterable[tuple[float, float]],
) -> Array:
    return np.sum(np.asarray(list(created), dtype=float), axis=0) - np.sum(
        np.asarray(list(annihilated), dtype=float), axis=0
    )


def peierls_phase(delta_r: Array, vector_potential: tuple[float, float]) -> complex:
    return np.exp(-1.0j * float(np.dot(delta_r, np.asarray(vector_potential))))


def displacement_balance_audit() -> dict[str, object]:
    x = (0.0, 0.0)
    xp = (1.0, 0.0)
    pair_delta = net_displacement((x, xp), (xp, x))
    bridge_delta = net_displacement((xp,), (x,))
    test_a = (0.37, -0.21)
    return {
        "position_balanced_pair_hop_delta_R": pair_delta.tolist(),
        "position_balanced_pair_hop_phase": {
            "real": float(peierls_phase(pair_delta, test_a).real),
            "imag": float(peierls_phase(pair_delta, test_a).imag),
        },
        "one_body_bridge_delta_R": bridge_delta.tolist(),
        "one_body_bridge_phase": {
            "real": float(peierls_phase(bridge_delta, test_a).real),
            "imag": float(peierls_phase(bridge_delta, test_a).imag),
        },
    }


def rho(capacity: int, pair_number: int) -> float:
    if capacity <= 1:
        return 0.0
    return pair_number * (capacity - pair_number) / (capacity - 1)


def middle_jacobi(
    capacity: int, interaction_u: float, gx2: float, gy2: float
) -> Array:
    """Conditional y-direction Jacobi matrix at n=V from the draft."""
    volume = capacity
    n_pairs = capacity
    j_total = gx2 + gy2
    p_ph = gy2 - gx2
    matrix = np.zeros((volume + 1, volume + 1), dtype=float)
    for r in range(volume + 1):
        matrix[r, r] = (
            2.0 * interaction_u * rho(volume, r)
            + 2.5 * interaction_u * rho(volume, n_pairs - r)
            + j_total
            / (2.0 * volume)
            * ((n_pairs - r) * (volume - r) + r * (volume - n_pairs + r))
        )
        if r < volume:
            root = math.sqrt(
                (n_pairs - r)
                * (volume - r)
                * (r + 1)
                * (volume - n_pairs + r + 1)
            )
            matrix[r, r + 1] = p_ph * root / (2.0 * volume)
            matrix[r + 1, r] = matrix[r, r + 1]
    return matrix


def transfer_dark_zero_mode(capacity: int, gx2: float, gy2: float) -> Array:
    p_ph = gy2 - gx2
    if min(abs(gx2), abs(gy2)) > 1.0e-14:
        raise ValueError("an exact transfer-dark vector requires gx2=0 or gy2=0")
    if abs(p_ph) < 1.0e-14:
        raise ValueError("both quadratures cannot vanish in this check")
    sign = -1.0 if p_ph > 0 else 1.0
    vector = np.array(
        [(sign**r) * math.comb(capacity, r) for r in range(capacity + 1)],
        dtype=float,
    )
    return vector / np.linalg.norm(vector)


def asymptotic_middle_stiffness(interaction_u: float, gx2: float, gy2: float) -> float:
    """Limit of lambda_min(Q_y^(V))/(4V) for the conditional Jacobi family."""
    j_total = gx2 + gy2
    g_min2 = min(gx2, gy2)
    return 0.125 * min(j_total, g_min2 + 2.25 * interaction_u)


def jacobi_audit(
    maximum_capacity: int = 120,
    interaction_u: float = 1.0,
    gx: float = 0.8,
    gy: float = 0.3,
) -> dict[str, object]:
    gx2 = gx * gx
    gy2 = gy * gy
    zero_mode_rows = []
    for dark in ("X", "Y"):
        zx2, zy2 = (0.0, 1.0) if dark == "X" else (1.0, 0.0)
        q = middle_jacobi(maximum_capacity, 0.0, zx2, zy2)
        vector = transfer_dark_zero_mode(maximum_capacity, zx2, zy2)
        zero_mode_rows.append(
            {
                "vanishing_quadrature": dark,
                "capacity": maximum_capacity,
                "residual_norm": float(np.linalg.norm(q @ vector)),
                "smallest_eigenvalue": float(np.linalg.eigvalsh(q)[0]),
            }
        )

    capacities = sorted(
        set(v for v in (8, 16, 32, 64, maximum_capacity) if v <= maximum_capacity)
    )
    predicted = asymptotic_middle_stiffness(interaction_u, gx2, gy2)
    convergence = []
    for capacity in capacities:
        q = middle_jacobi(capacity, interaction_u, gx2, gy2)
        stiffness = float(np.linalg.eigvalsh(q)[0] / (4.0 * capacity))
        convergence.append(
            {
                "capacity": capacity,
                "stiffness": stiffness,
                "difference_from_limit": stiffness - predicted,
            }
        )

    endpoint_condition = interaction_u >= 4.0 * max(gx2, gy2) / 9.0
    return {
        "parameters": {
            "U": interaction_u,
            "gX": gx,
            "gY": gy,
            "gX2": gx2,
            "gY2": gy2,
            "J": gx2 + gy2,
            "P_ph": gy2 - gx2,
        },
        "transfer_dark_zero_modes": zero_mode_rows,
        "thermodynamic_limit": predicted,
        "endpoint_localization_condition": endpoint_condition,
        "endpoint_J_over_8": (gx2 + gy2) / 8.0,
        "finite_capacity_convergence": convergence,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=int, default=19)
    parser.add_argument("--step", type=float, default=1.0e-6)
    parser.add_argument("--max-capacity", type=int, default=120)
    parser.add_argument("--u", type=float, default=1.0)
    parser.add_argument("--gx", type=float, default=0.8)
    parser.add_argument("--gy", type=float, default=0.3)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    results = {
        "bridge_covariance": bridge_covariance_audit(args.step, args.grid),
        "displacement_balance": displacement_balance_audit(),
        "conditional_jacobi": jacobi_audit(
            args.max_capacity, args.u, args.gx, args.gy
        ),
    }

    text = json.dumps(results, indent=2)
    print(text)
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text + "\n", encoding="utf-8")

    tolerance = 5.0e-8
    cov = results["bridge_covariance"]
    if cov["fixed_kernel_derivative_error_X"] > tolerance:
        raise SystemExit("fixed-kernel X derivative regression failed")
    if cov["fixed_kernel_derivative_error_Y"] > tolerance:
        raise SystemExit("fixed-kernel Y derivative regression failed")
    if cov["covariant_projected_derivative_max"] > tolerance:
        raise SystemExit("covariant bridge source failed to cancel")
    for row in results["conditional_jacobi"]["transfer_dark_zero_modes"]:
        if row["residual_norm"] > tolerance:
            raise SystemExit("transfer-dark exact zero mode check failed")


if __name__ == "__main__":
    main()
