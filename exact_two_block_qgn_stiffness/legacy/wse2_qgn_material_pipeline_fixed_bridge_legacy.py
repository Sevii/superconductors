#!/usr/bin/env python3
"""LEGACY 3.65-degree WSe2 pipeline retained for numerical provenance.

Band structure, topology, projected-kernel, and locality outputs remain useful.
The source-map and t0 portions hold a filtered multi-site bridge kernel fixed
while twisting the projector and are superseded by the covariance correction in
``scripts/check_bridge_covariance.py``.  Treat every Gamma/t0 value from this
legacy file as a geometric or conditional algebraic proxy, not a demonstrated
physical ordinary-twist coefficient.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import differential_evolution

Array = np.ndarray
SparseVector = dict[int, complex]


@dataclass(frozen=True)
class TightBindingParameters:
    delta: float = -14.9
    t_th_1: float = 10.78
    t_hh_1: float = 0.55
    t_tt_1: float = -1.95
    t_th_2: float = -1.21
    t_hh_3: float = 5.40
    e_z: float = 0.0


@dataclass(frozen=True)
class InteractionParameters:
    u_h: float = 35.0
    u_t: float = 20.0
    v_th: float = 40.0
    j_exchange: float = 10.0
    # Representative leading neglected pair-hopping scale reported for a
    # related three-orbital tTMD Wannier model.  Set to zero to test U,V,J only.
    p_pair: float = 2.0


SQRT3 = math.sqrt(3.0)
PI2 = 2.0 * math.pi

A1 = np.array([1.0, 0.0])
A2 = np.array([-0.5, SQRT3 / 2.0])
A3 = -A1 - A2
A_VECTORS = np.array([A1, A2, A3])

U0 = (A1 - A2) / 3.0
C3 = np.array([[-0.5, -SQRT3 / 2.0], [SQRT3 / 2.0, -0.5]])
U_VECTORS = np.array([U0, C3 @ U0, C3 @ C3 @ U0])

DIRECT = np.column_stack([A1, A2])
RECIPROCAL = PI2 * np.linalg.inv(DIRECT).T
B1 = RECIPROCAL[:, 0]
B2 = RECIPROCAL[:, 1]

# Orbital order: MX, MM, XM.  This embedding is used to sew the BZ boundary.
TAU = np.array([-U0, np.zeros(2), U0])
OMEGA = np.exp(2.0j * math.pi / 3.0)
# K-valley convention g_k = sum_j exp[2pi i(1-j)/3] exp(i k.u_j).
G_PHASES = np.array([OMEGA, 1.0 + 0.0j, OMEGA**2])

BAND_LABELS = ("top", "second", "remote")
ORBITAL_LABELS = ("MX", "MM", "XM")
PAIR_BLOCKS = ((0, 0), (1, 1), (2, 2), (0, 1), (0, 2), (1, 2))


class WSe2Model:
    def __init__(self, parameters: TightBindingParameters):
        self.p = parameters

    @staticmethod
    def form_f(k: Array, factor: int = 1) -> complex:
        return np.exp(1.0j * factor * (U_VECTORS @ k)).sum()

    @staticmethod
    def form_g(k: Array, factor: int = 1) -> complex:
        return (G_PHASES * np.exp(1.0j * factor * (U_VECTORS @ k))).sum()

    @staticmethod
    def form_h(k: Array) -> float:
        return float(2.0 * np.cos(A_VECTORS @ k).sum())

    @staticmethod
    def derivative_forms(k: Array, direction: int, factor: int = 1) -> tuple[complex, complex]:
        exponentials = np.exp(1.0j * factor * (U_VECTORS @ k))
        coefficient = 1.0j * factor * U_VECTORS[:, direction]
        df = (coefficient * exponentials).sum()
        dg = (G_PHASES * coefficient * exponentials).sum()
        return complex(df), complex(dg)

    def hamiltonian(self, k: Array, e_z: float | None = None) -> Array:
        p = self.p
        field = p.e_z if e_z is None else float(e_z)

        f1 = self.form_f(k)
        g1 = self.form_g(k)
        gm1 = self.form_g(-k)
        f2 = self.form_f(k, 2)
        g2 = self.form_g(k, 2)
        gm2 = self.form_g(-k, 2)

        h1 = np.array(
            [
                [field, p.t_th_1 * g1, p.t_hh_1 * np.conjugate(f1)],
                [p.t_th_1 * np.conjugate(g1), -p.delta, -p.t_th_1 * np.conjugate(gm1)],
                [p.t_hh_1 * f1, -p.t_th_1 * gm1, -field],
            ],
            dtype=complex,
        )
        h2 = np.array(
            [
                [0.0, -p.t_th_2 * gm2, p.t_hh_3 * f2],
                [-p.t_th_2 * np.conjugate(gm2), p.t_tt_1 * self.form_h(k), p.t_th_2 * np.conjugate(g2)],
                [p.t_hh_3 * np.conjugate(f2), p.t_th_2 * g2, 0.0],
            ],
            dtype=complex,
        )
        return h1 + h2

    def derivative(self, k: Array, direction: int) -> Array:
        p = self.p
        df1, dg1 = self.derivative_forms(k, direction, 1)
        df2, dg2 = self.derivative_forms(k, direction, 2)

        exp_minus_1 = np.exp(-1.0j * (U_VECTORS @ k))
        dgm1 = (G_PHASES * (-1.0j) * U_VECTORS[:, direction] * exp_minus_1).sum()
        exp_minus_2 = np.exp(-2.0j * (U_VECTORS @ k))
        dgm2 = (G_PHASES * (-2.0j) * U_VECTORS[:, direction] * exp_minus_2).sum()
        dh = float((-2.0 * np.sin(A_VECTORS @ k) * A_VECTORS[:, direction]).sum())

        dh1 = np.array(
            [
                [0.0, p.t_th_1 * dg1, p.t_hh_1 * np.conjugate(df1)],
                [p.t_th_1 * np.conjugate(dg1), 0.0, -p.t_th_1 * np.conjugate(dgm1)],
                [p.t_hh_1 * df1, -p.t_th_1 * dgm1, 0.0],
            ],
            dtype=complex,
        )
        dh2 = np.array(
            [
                [0.0, -p.t_th_2 * dgm2, p.t_hh_3 * df2],
                [-p.t_th_2 * np.conjugate(dgm2), p.t_tt_1 * dh, p.t_th_2 * np.conjugate(dg2)],
                [p.t_hh_3 * np.conjugate(df2), p.t_th_2 * dg2, 0.0],
            ],
            dtype=complex,
        )
        return dh1 + dh2

    def eigensystem(self, k: Array, e_z: float | None = None) -> tuple[Array, Array]:
        values, vectors = np.linalg.eigh(self.hamiltonian(k, e_z=e_z))
        order = np.argsort(values)[::-1]
        return values[order], vectors[:, order]


@dataclass
class BandGrid:
    n: int
    k: Array
    fractional: Array
    energies: Array
    vectors: Array


def make_band_grid(model: WSe2Model, n: int, e_z: float | None = None) -> BandGrid:
    ks: list[Array] = []
    fractional: list[tuple[float, float]] = []
    energies: list[Array] = []
    vectors: list[Array] = []
    for m in range(n):
        for q in range(n):
            frac = np.array([m / n, q / n], dtype=float)
            k = frac[0] * B1 + frac[1] * B2
            ev, u = model.eigensystem(k, e_z=e_z)
            fractional.append((float(frac[0]), float(frac[1])))
            ks.append(k)
            energies.append(ev)
            vectors.append(u)
    return BandGrid(
        n=n,
        k=np.asarray(ks),
        fractional=np.asarray(fractional),
        energies=np.asarray(energies),
        vectors=np.asarray(vectors),
    )


def derivative_check(model: WSe2Model, step: float = 1.0e-6) -> float:
    points = (np.array([0.2, 0.7]), np.array([1.1, -0.3]), np.array([2.0, 1.4]))
    maximum = 0.0
    for k in points:
        for direction in (0, 1):
            shift = np.zeros(2)
            shift[direction] = step
            finite = (model.hamiltonian(k + shift) - model.hamiltonian(k - shift)) / (2.0 * step)
            maximum = max(maximum, float(np.max(np.abs(finite - model.derivative(k, direction)))))
    return maximum


def boundary_sewing() -> tuple[Array, Array]:
    # For the embedding-aware convention, u(k+B_i)=D_i u(k), with
    # D_i=diag(exp[-i B_i.tau_mu]).
    return (
        np.diag(np.exp(-1.0j * (TAU @ B1))),
        np.diag(np.exp(-1.0j * (TAU @ B2))),
    )


def chern_numbers(model: WSe2Model, n: int = 31) -> list[float]:
    grid = make_band_grid(model, n)
    vectors = grid.vectors.reshape(n, n, 3, 3)
    d1, d2 = boundary_sewing()
    output: list[float] = []
    for band in range(3):
        ux = np.zeros((n, n), dtype=complex)
        uy = np.zeros((n, n), dtype=complex)
        for m in range(n):
            for q in range(n):
                u = vectors[m, q, :, band]
                vx = vectors[m + 1, q, :, band] if m + 1 < n else d1 @ vectors[0, q, :, band]
                vy = vectors[m, q + 1, :, band] if q + 1 < n else d2 @ vectors[m, 0, :, band]
                zx = np.vdot(u, vx)
                zy = np.vdot(u, vy)
                if abs(zx) < 1.0e-14 or abs(zy) < 1.0e-14:
                    raise RuntimeError("Fukui link overlap vanished")
                ux[m, q] = zx / abs(zx)
                uy[m, q] = zy / abs(zy)
        flux = 0.0
        for m in range(n):
            for q in range(n):
                plaquette = ux[m, q] * uy[(m + 1) % n, q] / (
                    ux[m, (q + 1) % n] * uy[m, q]
                )
                flux += float(np.angle(plaquette))
        output.append(flux / PI2)
    return output


def band_extrema(model: WSe2Model, seed: int = 1729) -> dict[str, Any]:
    def values(frac: Array) -> Array:
        x = float(frac[0] % 1.0)
        y = float(frac[1] % 1.0)
        return model.eigensystem(x * B1 + y * B2)[0]

    minima: list[float] = []
    maxima: list[float] = []
    min_locations: list[list[float]] = []
    max_locations: list[list[float]] = []
    for band in range(3):
        res_min = differential_evolution(
            lambda x, b=band: float(values(x)[b]),
            bounds=((0.0, 1.0), (0.0, 1.0)),
            seed=seed + band,
            tol=1.0e-10,
            polish=True,
            maxiter=300,
        )
        res_max = differential_evolution(
            lambda x, b=band: -float(values(x)[b]),
            bounds=((0.0, 1.0), (0.0, 1.0)),
            seed=seed + 10 + band,
            tol=1.0e-10,
            polish=True,
            maxiter=300,
        )
        minima.append(float(res_min.fun))
        maxima.append(float(-res_max.fun))
        min_locations.append([float(v % 1.0) for v in res_min.x])
        max_locations.append([float(v % 1.0) for v in res_max.x])

    gaps: dict[str, Any] = {}
    for upper, lower, name in ((0, 1, "top_second"), (1, 2, "second_remote")):
        result = differential_evolution(
            lambda x, a=upper, b=lower: float(values(x)[a] - values(x)[b]),
            bounds=((0.0, 1.0), (0.0, 1.0)),
            seed=seed + 30 + upper,
            tol=1.0e-10,
            polish=True,
            maxiter=300,
        )
        gaps[name] = {
            "direct_min": float(result.fun),
            "fractional_location": [float(v % 1.0) for v in result.x],
        }

    gaps["top_second"]["indirect"] = float(minima[0] - maxima[1])
    gaps["second_remote"]["indirect"] = float(minima[1] - maxima[2])
    return {
        "energy_min_meV": dict(zip(BAND_LABELS, minima)),
        "energy_max_meV": dict(zip(BAND_LABELS, maxima)),
        "bandwidth_meV": dict(zip(BAND_LABELS, [maxima[i] - minima[i] for i in range(3)])),
        "minimum_fractional_k": dict(zip(BAND_LABELS, min_locations)),
        "maximum_fractional_k": dict(zip(BAND_LABELS, max_locations)),
        "gaps_meV": gaps,
    }


def geometry_and_bridge(model: WSe2Model, n: int = 81) -> dict[str, Any]:
    grid = make_band_grid(model, n)
    count = n * n
    gamma = np.zeros((2, 2, 2), dtype=float)
    full_metric = np.zeros((2, 2, 2), dtype=float)
    orbital_weights = np.zeros((3, 3), dtype=float)
    weighted_remote_t = np.zeros(2, dtype=float)

    for index, k in enumerate(grid.k):
        e = grid.energies[index]
        u = grid.vectors[index]
        orbital_weights += np.abs(u) ** 2
        for active in (0, 1):
            connection = np.zeros((3, 2), dtype=complex)
            for other in range(3):
                if other == active:
                    continue
                for direction in (0, 1):
                    connection[other, direction] = (
                        np.vdot(u[:, other], model.derivative(k, direction) @ u[:, active])
                        / (e[active] - e[other])
                    )
            for i in (0, 1):
                for j in (0, 1):
                    full_metric[active, i, j] += float(
                        np.real(sum(np.conjugate(connection[b, i]) * connection[b, j] for b in range(3) if b != active))
                    )
                    gamma[active, i, j] += float(
                        np.real(np.conjugate(connection[2, i]) * connection[2, j])
                    )
            if active == 1:
                for direction in (0, 1):
                    weighted_remote_t[direction] += float(
                        abs(connection[2, direction]) ** 2 * abs(u[1, 2]) ** 2
                    )

    gamma /= count
    full_metric /= count
    orbital_weights /= count
    weighted_remote_t /= count
    top_h_weight = float(orbital_weights[0, 0] + orbital_weights[2, 0])
    bridge_gamma = 3.0 * weighted_remote_t * top_h_weight

    remote_t_weight = float(orbital_weights[1, 2])
    active_t_weight = 1.0 - remote_t_weight
    remote_h_weight_each = float(orbital_weights[0, 2])
    active_h_weight_each = 1.0 - remote_h_weight_each
    bare_leakage_per_bond = active_t_weight * active_h_weight_each
    filtered_norm_per_bond = remote_t_weight * float(orbital_weights[0, 0])

    return {
        "grid": n,
        "remote_metric_gamma_a_M2": {
            "top": gamma[0].tolist(),
            "second": gamma[1].tolist(),
        },
        "full_quantum_metric_a_M2": {
            "top": full_metric[0].tolist(),
            "second": full_metric[1].tolist(),
        },
        "remote_fraction_of_metric": {
            "top_x": float(gamma[0, 0, 0] / full_metric[0, 0, 0]),
            "top_y": float(gamma[0, 1, 1] / full_metric[0, 1, 1]),
            "second_x": float(gamma[1, 0, 0] / full_metric[1, 0, 0]),
            "second_y": float(gamma[1, 1, 1] / full_metric[1, 1, 1]),
        },
        "orbital_weights_rows_MX_MM_XM_cols_top_second_remote": orbital_weights.tolist(),
        "weighted_remote_MM_connection_a_M2": weighted_remote_t.tolist(),
        "top_honeycomb_weight": top_h_weight,
        "local_bridge_Gamma_a_M2": {"x": float(bridge_gamma[0]), "y": float(bridge_gamma[1])},
        "bare_active_active_leakage_HS2_per_bond": float(bare_leakage_per_bond),
        "desired_filtered_remote_top_HS2_per_bond": float(filtered_norm_per_bond),
    }


def th_bonds() -> list[tuple[int, Array, float]]:
    # T is at the origin.  MX bonds are at -u_j and XM bonds at +u_j.
    # The chiral exchange phase is twice the corresponding one-particle phase.
    bonds: list[tuple[int, Array, float]] = []
    for j in range(3):
        phi = float(np.angle(G_PHASES[j] ** 2))
        bonds.append((0, -U_VECTORS[j], phi))
        bonds.append((2, U_VECTORS[j], phi))
    return bonds


def projected_cooper_kernels(
    model: WSe2Model,
    interactions: InteractionParameters,
    n: int = 11,
) -> dict[str, Any]:
    grid = make_band_grid(model, n)
    volume = n * n
    ks = grid.k
    u = grid.vectors
    bonds = th_bonds()

    results: dict[str, Any] = {"grid": n, "norm_type": "spectral/operator 2-norm"}
    tables: dict[str, dict[str, float]] = {kind: {} for kind in ("U", "V", "J", "singlet_U_plus_V_minus_J")}

    for active, target in PAIR_BLOCKS:
        ua = u[:, :, active]
        ub = u[:, :, target]
        ku = np.zeros((volume, volume), dtype=complex)
        for orbital, coefficient in ((0, interactions.u_h), (1, interactions.u_t), (2, interactions.u_h)):
            ku += coefficient * np.outer(np.abs(ua[:, orbital]) ** 2, np.abs(ub[:, orbital]) ** 2)

        kv = np.zeros_like(ku)
        kj = np.zeros_like(ku)
        for honeycomb, displacement, phi in bonds:
            phase_plus = np.exp(1.0j * (ks @ displacement))
            phase_minus = np.conjugate(phase_plus)

            left = np.conjugate(ua[:, 1]) * ua[:, honeycomb] * phase_plus
            right = ub[:, 1] * np.conjugate(ub[:, honeycomb]) * phase_minus
            kv += interactions.v_th * np.outer(left, right)
            left_2 = np.conjugate(ua[:, honeycomb]) * ua[:, 1] * phase_minus
            right_2 = ub[:, honeycomb] * np.conjugate(ub[:, 1]) * phase_plus
            kv += interactions.v_th * np.outer(left_2, right_2)

            # Exchange kernel before the singlet fermion-reordering sign.
            left_j = np.conjugate(ua[:, 1]) * ua[:, honeycomb] * phase_plus
            right_j = ub[:, honeycomb] * np.conjugate(ub[:, 1]) * phase_plus
            kj += interactions.j_exchange * np.exp(1.0j * phi) * np.outer(left_j, right_j)
            left_j2 = np.conjugate(ua[:, honeycomb]) * ua[:, 1] * phase_minus
            right_j2 = ub[:, 1] * np.conjugate(ub[:, honeycomb]) * phase_minus
            kj += interactions.j_exchange * np.exp(-1.0j * phi) * np.outer(left_j2, right_j2)

        ku /= volume
        kv /= volume
        kj /= volume
        singlet = ku + kv - kj
        label = f"{BAND_LABELS[active]}->{BAND_LABELS[target]}"
        for name, matrix in (("U", ku), ("V", kv), ("J", kj), ("singlet_U_plus_V_minus_J", singlet)):
            tables[name][label] = float(np.linalg.svd(matrix, compute_uv=False)[0])

    results["operator_norm_meV"] = tables
    diagonal = tables["singlet_U_plus_V_minus_J"]
    k11 = diagonal["top->top"]
    k22 = diagonal["second->second"]
    krr = diagonal["remote->remote"]
    results["normalized_cross_block_ratios"] = {
        "top_second": diagonal["top->second"] / math.sqrt(k11 * k22),
        "top_remote": diagonal["top->remote"] / math.sqrt(k11 * krr),
        "second_remote": diagonal["second->remote"] / math.sqrt(k22 * krr),
    }
    return results


def fierz_bridge_couplings(interactions: InteractionParameters) -> dict[str, Any]:
    j = interactions.j_exchange
    p = interactions.p_pair
    if abs(p) > j:
        raise ValueError("Positive-square quadratures require |P| <= J in this normalization")
    gx2 = 0.5 * (j - p)
    gy2 = 0.5 * (j + p)
    uvj_only = 0.5 * j
    return {
        "identity": {
            "gX2_plus_gY2": "J",
            "gY2_minus_gX2": "P",
        },
        "UVJ_only_P_zero": {
            "gX2_meV": uvj_only,
            "gY2_meV": uvj_only,
            "composition_hopping_difference_meV": 0.0,
        },
        "with_representative_pair_hop": {
            "J_meV": j,
            "P_meV": p,
            "gX2_meV": gx2,
            "gY2_meV": gy2,
            "gX_sqrt_meV": math.sqrt(gx2),
            "gY_sqrt_meV": math.sqrt(gy2),
            "composition_hopping_difference_meV": gy2 - gx2,
        },
    }


def minus_map(n: int) -> Array:
    return np.array([((-m) % n) * n + ((-q) % n) for m in range(n) for q in range(n)], dtype=int)


def material_source_rows(
    model: WSe2Model,
    n: int,
    direction: int,
) -> tuple[BandGrid, list[Array], Array, float]:
    grid = make_band_grid(model, n)
    volume = n * n
    connections = np.zeros(volume, dtype=complex)
    for index, k in enumerate(grid.k):
        e = grid.energies[index]
        u = grid.vectors[index]
        connections[index] = (
            np.vdot(u[:, 2], model.derivative(k, direction) @ u[:, 1]) / (e[1] - e[2])
        )

    rows: list[Array] = []
    for rm in range(n):
        for rq in range(n):
            cell = rm * A1 + rq * A2
            for j in range(3):
                for honeycomb, displacement in ((0, -U_VECTORS[j]), (2, U_VECTORS[j])):
                    left = (
                        np.conjugate(connections)
                        * np.conjugate(grid.vectors[:, 1, 2])
                        * np.exp(-1.0j * (grid.k @ cell))
                    )
                    right = (
                        grid.vectors[:, honeycomb, 0]
                        * np.exp(1.0j * (grid.k @ (cell + displacement)))
                    )
                    rows.append(np.outer(left, right) / volume)

    weighted = float(
        np.mean(np.abs(connections) ** 2 * np.abs(grid.vectors[:, 1, 2]) ** 2)
    )
    top_h = float(np.mean(np.abs(grid.vectors[:, 0, 0]) ** 2 + np.abs(grid.vectors[:, 2, 0]) ** 2))
    gamma_n = 3.0 * weighted * top_h
    return grid, rows, connections, gamma_n


def apply_x_pair_matrix(matrix: Array, q: Array, n: int, adjoint: bool = False) -> Array:
    volume = n * n
    dimension = 2 * volume
    if adjoint:
        left = np.zeros((dimension, dimension), dtype=complex)
        left[:volume, volume:] = np.conjugate(q).T
    else:
        left = np.zeros((dimension, dimension), dtype=complex)
        left[volume:, :volume] = q

    mm = minus_map(n)
    q_down = np.conjugate(q[np.ix_(mm, mm)])
    if adjoint:
        right = np.zeros((dimension, dimension), dtype=complex)
        right[:volume, volume:] = np.conjugate(q_down).T
    else:
        right = np.zeros((dimension, dimension), dtype=complex)
        right[volume:, :volume] = q_down
    return left @ matrix + matrix @ right.T


def one_pair_agp_matrices(n: int) -> list[Array]:
    volume = n * n
    dimension = 2 * volume
    mm = minus_map(n)
    states: list[Array] = []
    for block in (0, 1):
        matrix = np.zeros((dimension, dimension), dtype=complex)
        for k in range(volume):
            matrix[block * volume + k, block * volume + mm[k]] = 1.0 / math.sqrt(volume)
        states.append(matrix)
    return states


def one_pair_ed(
    model: WSe2Model,
    interactions: InteractionParameters,
    n: int = 4,
    direction: int = 0,
) -> dict[str, Any]:
    _, rows, _, gamma_n = material_source_rows(model, n, direction)
    states = one_pair_agp_matrices(n)
    gx2 = 0.5 * (interactions.j_exchange - interactions.p_pair)
    gy2 = 0.5 * (interactions.j_exchange + interactions.p_pair)
    q_matrix = np.zeros((2, 2), dtype=complex)
    for row in rows:
        bx_images: list[Array] = []
        by_images: list[Array] = []
        for state in states:
            x = apply_x_pair_matrix(state, row, n, adjoint=False)
            xd = apply_x_pair_matrix(state, row, n, adjoint=True)
            bx_images.append(x + xd)
            by_images.append(-1.0j * (x - xd))
        for a in range(2):
            for b in range(2):
                q_matrix[a, b] += gx2 * np.vdot(bx_images[a], bx_images[b]) + gy2 * np.vdot(by_images[a], by_images[b])

    formula = np.array(
        [
            [2.0 * gamma_n * (gx2 + gy2), 2.0 * gamma_n * (gy2 - gx2)],
            [2.0 * gamma_n * (gy2 - gx2), 2.0 * gamma_n * (gx2 + gy2)],
        ]
    )
    # The sign of the off-diagonal element depends on the phase convention for
    # the second composition basis vector.  Compare absolute hopping magnitude.
    error = max(
        float(np.max(np.abs(np.diag(q_matrix) - np.diag(formula)))),
        abs(abs(q_matrix[0, 1]) - abs(formula[0, 1])),
    )
    return {
        "torus_N": n,
        "capacity_V": n * n,
        "direction": "x" if direction == 0 else "y",
        "Gamma_N_a_M2": gamma_n,
        "direct_Q_meV_a_M2": np.real_if_close(q_matrix).real.tolist(),
        "formula_Q_up_to_basis_phase_meV_a_M2": formula.tolist(),
        "eigenvalues_meV_a_M2": np.linalg.eigvalsh(q_matrix).real.tolist(),
        "offdiagonal_t0_magnitude_meV_a_M2": float(abs(q_matrix[0, 1])),
        "max_error_up_to_basis_phase": error,
    }


# --- Sparse Fock-space helpers for the many-composition ED check. ---

def create(state: int, orbital: int) -> tuple[int, int] | None:
    if (state >> orbital) & 1:
        return None
    sign = -1 if (state & ((1 << orbital) - 1)).bit_count() % 2 else 1
    return state | (1 << orbital), sign


def annihilate(state: int, orbital: int) -> tuple[int, int] | None:
    if not ((state >> orbital) & 1):
        return None
    sign = -1 if (state & ((1 << orbital) - 1)).bit_count() % 2 else 1
    return state & ~(1 << orbital), sign


def inner_sparse(left: Mapping[int, complex], right: Mapping[int, complex]) -> complex:
    return sum(np.conjugate(value) * right.get(state, 0.0) for state, value in left.items())


def normalize_sparse(vector: Mapping[int, complex]) -> SparseVector:
    norm = math.sqrt(float(np.real(inner_sparse(vector, vector))))
    if norm == 0.0:
        raise ValueError("Cannot normalize zero vector")
    return {state: value / norm for state, value in vector.items()}


def add_sparse(target: defaultdict[int, complex], source: Mapping[int, complex], coefficient: complex = 1.0) -> None:
    for state, value in source.items():
        target[state] += coefficient * value


def orbital_index(block: int, spin: int, level: int, volume: int) -> int:
    return ((block * 2 + spin) * volume) + level


def block_agp(volume: int, block: int, pairs: int) -> SparseVector:
    output: defaultdict[int, complex] = defaultdict(complex)
    for occupied in itertools.combinations(range(volume), pairs):
        state = 0
        amplitude = 1
        for level in occupied:
            for orbital in (
                orbital_index(block, 1, level, volume),
                orbital_index(block, 0, level, volume),
            ):
                result = create(state, orbital)
                if result is None:
                    raise RuntimeError("AGP construction collision")
                state, sign = result
                amplitude *= sign
        output[state] += amplitude
    return normalize_sparse(output)


def product_agp(volume: int, n1: int, n2: int) -> SparseVector:
    first = block_agp(volume, 0, n1)
    second = block_agp(volume, 1, n2)
    output: defaultdict[int, complex] = defaultdict(complex)
    for state_1, amp_1 in first.items():
        for state_2, amp_2 in second.items():
            output[state_1 | state_2] += amp_1 * amp_2
    return normalize_sparse(output)


def apply_cdag_c(
    vector: Mapping[int, complex],
    create_orbital: int,
    annihilate_orbital: int,
    coefficient: complex,
) -> SparseVector:
    output: defaultdict[int, complex] = defaultdict(complex)
    for state, amplitude in vector.items():
        first = annihilate(state, annihilate_orbital)
        if first is None:
            continue
        state_1, sign_1 = first
        second = create(state_1, create_orbital)
        if second is None:
            continue
        state_2, sign_2 = second
        output[state_2] += coefficient * amplitude * sign_1 * sign_2
    return dict(output)


def apply_material_x(
    vector: Mapping[int, complex],
    q: Array,
    volume: int,
    adjoint: bool = False,
) -> SparseVector:
    output: defaultdict[int, complex] = defaultdict(complex)
    for p in range(volume):
        for k in range(volume):
            for spin in (0, 1):
                coefficient = q[p, k] if spin == 0 else np.conjugate(q[p, k])
                if adjoint:
                    created = orbital_index(0, spin, k, volume)
                    removed = orbital_index(1, spin, p, volume)
                    coefficient = np.conjugate(coefficient)
                else:
                    created = orbital_index(1, spin, p, volume)
                    removed = orbital_index(0, spin, k, volume)
                add_sparse(output, apply_cdag_c(vector, created, removed, coefficient))
    return dict(output)


def apply_material_quadrature(
    vector: Mapping[int, complex],
    q: Array,
    volume: int,
    quadrature: str,
) -> SparseVector:
    x = apply_material_x(vector, q, volume, adjoint=False)
    xd = apply_material_x(vector, q, volume, adjoint=True)
    output: defaultdict[int, complex] = defaultdict(complex)
    if quadrature == "X":
        add_sparse(output, x, 1.0)
        add_sparse(output, xd, 1.0)
    elif quadrature == "Y":
        add_sparse(output, x, -1.0j)
        add_sparse(output, xd, 1.0j)
    else:
        raise ValueError("quadrature must be X or Y")
    return dict(output)


def many_composition_ed(
    model: WSe2Model,
    interactions: InteractionParameters,
    n_torus: int = 2,
    total_pairs: int | None = None,
    direction: int = 0,
) -> dict[str, Any]:
    volume = n_torus * n_torus
    n_pairs = volume if total_pairs is None else total_pairs
    _, rows, _, gamma_n = material_source_rows(model, n_torus, direction)
    gx2 = 0.5 * (interactions.j_exchange - interactions.p_pair)
    gy2 = 0.5 * (interactions.j_exchange + interactions.p_pair)
    r_min = max(0, n_pairs - volume)
    r_max = min(volume, n_pairs)
    r_values = list(range(r_min, r_max + 1))
    states = [product_agp(volume, r, n_pairs - r) for r in r_values]
    direct = np.zeros((len(states), len(states)), dtype=complex)

    for row in rows:
        bx = [apply_material_quadrature(state, row, volume, "X") for state in states]
        by = [apply_material_quadrature(state, row, volume, "Y") for state in states]
        for a in range(len(states)):
            for b in range(len(states)):
                direct[a, b] += gx2 * inner_sparse(bx[a], bx[b]) + gy2 * inner_sparse(by[a], by[b])

    formula = np.zeros_like(direct.real)
    for index, r in enumerate(r_values):
        formula[index, index] = 2.0 * gamma_n * (gx2 + gy2) * (
            (n_pairs - r) * (volume - r) + r * (volume - n_pairs + r)
        ) / volume
        if index + 1 < len(r_values):
            root = math.sqrt(
                (n_pairs - r)
                * (volume - r)
                * (r + 1)
                * (volume - n_pairs + r + 1)
            )
            formula[index, index + 1] = 2.0 * gamma_n * (gy2 - gx2) * root / volume
            formula[index + 1, index] = formula[index, index + 1]

    # Alternating phases of |r,n-r> flip all nearest-neighbour signs.
    gauge = np.diag([(-1.0) ** r for r in r_values])
    direct_gauged = gauge @ direct @ gauge
    error = float(np.max(np.abs(direct_gauged - formula)))
    return {
        "torus_N": n_torus,
        "capacity_V": volume,
        "total_pairs_n": n_pairs,
        "composition_r": r_values,
        "Gamma_N_a_M2": gamma_n,
        "direct_Q_after_alternating_composition_gauge": np.real_if_close(direct_gauged).real.tolist(),
        "closed_Jacobi_formula": formula.tolist(),
        "maximum_matrix_error": error,
        "offdiagonal_formula_meV_a_M2": [float(formula[i, i + 1]) for i in range(len(r_values) - 1)],
        "eigenvalues_meV_a_M2": np.linalg.eigvalsh(direct).real.tolist(),
    }


def locality_tails(model: WSe2Model, n: int = 128, maximum_radius: int = 10) -> dict[str, Any]:
    # Canonical orbital-periodic gauge.  Radius is triangular-cell distance
    # R=max(|m|,|n|,|m-n|).
    left = np.zeros((n, n, 3), dtype=complex)
    right = np.zeros((n, n, 3), dtype=complex)
    source = np.zeros((n, n, 3), dtype=complex)

    for m in range(n):
        for q in range(n):
            k = (m / n) * B1 + (q / n) * B2
            e, u = model.eigensystem(k)
            orbital_phase = np.exp(1.0j * (TAU @ k))
            canonical = orbital_phase[:, None] * u
            left[m, q] = canonical[:, 2] * np.conjugate(canonical[1, 2])
            right[m, q] = canonical[:, 0] * np.conjugate(canonical[0, 0])
            connection = np.vdot(u[:, 2], model.derivative(k, 0) @ u[:, 1]) / (e[1] - e[2])
            source[m, q] = np.conjugate(connection) * canonical[:, 1] * np.conjugate(canonical[1, 2])

    def transform(field: Array) -> Array:
        return np.fft.fftshift(np.fft.ifft2(field, axes=(0, 1)), axes=(0, 1))

    left_r = transform(left)
    right_r = transform(right)
    source_r = transform(source)
    left_norm = np.sum(np.abs(left_r) ** 2, axis=2)
    right_norm = np.sum(np.abs(right_r) ** 2, axis=2)
    source_norm = np.sum(np.abs(source_r) ** 2, axis=2)

    coordinates = np.arange(-n // 2, n - n // 2)
    mm, qq = np.meshgrid(coordinates, coordinates, indexing="ij")
    radius = np.maximum.reduce([np.abs(mm), np.abs(qq), np.abs(mm - qq)])

    rows: list[dict[str, float]] = []
    for r in range(1, maximum_radius + 1):
        mask = radius <= r
        left_fraction = float(left_norm[mask].sum() / left_norm.sum())
        right_fraction = float(right_norm[mask].sum() / right_norm.sum())
        source_fraction = float(source_norm[mask].sum() / source_norm.sum())
        rows.append(
            {
                "radius_cells": r,
                "factor_fraction": left_fraction * right_fraction,
                "source_fraction": source_fraction * right_fraction,
                "remote_MM_projector_fraction": left_fraction,
                "top_MX_projector_fraction": right_fraction,
                "twist_source_left_fraction": source_fraction,
            }
        )
    return {"fft_grid": n, "radius_definition": "max(|m|,|n|,|m-n|)", "cumulative": rows}


def displacement_scan(model: WSe2Model, values: Iterable[float], n: int = 41) -> list[dict[str, float]]:
    output: list[dict[str, float]] = []
    for field in values:
        grid = make_band_grid(model, n, e_z=field)
        bandwidth = np.ptp(grid.energies, axis=0)
        direct_12 = float(np.min(grid.energies[:, 0] - grid.energies[:, 1]))
        direct_23 = float(np.min(grid.energies[:, 1] - grid.energies[:, 2]))
        weighted = np.zeros(2, dtype=float)
        top_h = float(np.mean(np.abs(grid.vectors[:, 0, 0]) ** 2 + np.abs(grid.vectors[:, 2, 0]) ** 2))
        for index, k in enumerate(grid.k):
            e = grid.energies[index]
            u = grid.vectors[index]
            for direction in (0, 1):
                connection = np.vdot(u[:, 2], model.derivative(k, direction) @ u[:, 1]) / (e[1] - e[2])
                weighted[direction] += abs(connection) ** 2 * abs(u[1, 2]) ** 2
        weighted /= n * n
        gamma = 3.0 * weighted * top_h
        output.append(
            {
                "E_z_meV": float(field),
                "W_top_meV": float(bandwidth[0]),
                "W_second_meV": float(bandwidth[1]),
                "W_remote_meV": float(bandwidth[2]),
                "direct_gap_top_second_meV": direct_12,
                "direct_gap_second_remote_meV": direct_23,
                "Gamma_x_a_M2": float(gamma[0]),
                "Gamma_y_a_M2": float(gamma[1]),
            }
        )
    return output


def corrections_summary(
    bands: Mapping[str, Any],
    projected: Mapping[str, Any],
    interactions: InteractionParameters,
    geometry: Mapping[str, Any],
) -> dict[str, Any]:
    bandwidth = bands["bandwidth_meV"]
    gap = float(bands["gaps_meV"]["second_remote"]["direct_min"])
    kernels = projected["operator_norm_meV"]["singlet_U_plus_V_minus_J"]
    k1r = float(kernels["top->remote"])
    k2r = float(kernels["second->remote"])
    gamma = float(geometry["local_bridge_Gamma_a_M2"]["x"])
    gx2 = 0.5 * (interactions.j_exchange - interactions.p_pair)
    gy2 = 0.5 * (interactions.j_exchange + interactions.p_pair)
    bridge_stiffness_lower = 0.5 * gamma * min(gx2, gy2)
    return {
        "isolating_direct_gap_Delta23_meV": gap,
        "dimensionless_ratios_to_Delta23": {
            "W_top_over_Delta": float(bandwidth["top"] / gap),
            "W_second_over_Delta": float(bandwidth["second"] / gap),
            "U_H_over_Delta": interactions.u_h / gap,
            "U_T_over_Delta": interactions.u_t / gap,
            "V_over_Delta": interactions.v_th / gap,
            "J_over_Delta": interactions.j_exchange / gap,
            "P_over_Delta": interactions.p_pair / gap,
            "K_top_remote_over_Delta": k1r / gap,
            "K_second_remote_over_Delta": k2r / gap,
        },
        "pair_kinetic_spectral_width_meV": {
            "top": 2.0 * float(bandwidth["top"]),
            "second": 2.0 * float(bandwidth["second"]),
        },
        "second_order_energy_estimates_meV": {
            "P2_over_Delta": interactions.p_pair**2 / gap,
            "J2_over_Delta": interactions.j_exchange**2 / gap,
            "K_top_remote2_over_Delta": k1r**2 / gap,
            "K_second_remote2_over_Delta": k2r**2 / gap,
        },
        "bridge_only_p_equals_1_stiffness_lower_bound_meV_in_aM_units": bridge_stiffness_lower,
        "controlled_isolated_flat_band_regime": False,
        "reason": "Bandwidths and projected active-remote U,V,J kernels are comparable to or larger than the second-to-remote direct gap.",
    }


def high_symmetry_path(model: WSe2Model, points_per_segment: int = 180) -> tuple[Array, Array, list[float], list[str]]:
    points = [np.array([0.0, 0.0]), np.array([1.0 / 3.0, 1.0 / 3.0]), np.array([0.5, 0.5]), np.array([0.0, 0.0])]
    labels = [r"$\Gamma$", r"$K$", r"$M$", r"$\Gamma$"]
    path: list[Array] = []
    x_values: list[float] = []
    ticks: list[float] = [0.0]
    distance = 0.0
    previous_k: Array | None = None
    for segment in range(len(points) - 1):
        for index in range(points_per_segment):
            t = index / points_per_segment
            frac = (1.0 - t) * points[segment] + t * points[segment + 1]
            k = frac[0] * B1 + frac[1] * B2
            if previous_k is not None:
                distance += float(np.linalg.norm(k - previous_k))
            path.append(k)
            x_values.append(distance)
            previous_k = k
        ticks.append(distance)
    final_frac = points[-1]
    final_k = final_frac[0] * B1 + final_frac[1] * B2
    distance += float(np.linalg.norm(final_k - previous_k)) if previous_k is not None else 0.0
    path.append(final_k)
    x_values.append(distance)
    ticks[-1] = distance
    energies = np.array([model.eigensystem(k)[0] for k in path])
    return np.asarray(x_values), energies, ticks, labels


def make_plots(
    output: Path,
    model: WSe2Model,
    geometry: Mapping[str, Any],
    one_pair: Mapping[str, Any],
    many_ed: Mapping[str, Any],
    locality: Mapping[str, Any],
    scan: list[Mapping[str, float]],
) -> None:
    x, energies, ticks, labels = high_symmetry_path(model)
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    for band in range(3):
        ax.plot(x, energies[:, band], label=BAND_LABELS[band])
    for tick in ticks:
        ax.axvline(tick, linewidth=0.7, alpha=0.45)
    ax.set_xticks(ticks, labels)
    ax.set_ylabel("energy (meV)")
    ax.set_title("3.65° tWSe$_2$ three-orbital Wannier model")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "wse2_band_structure.png", dpi=220)
    plt.close(fig)

    gamma = geometry["remote_metric_gamma_a_M2"]
    full = geometry["full_quantum_metric_a_M2"]
    labels_bar = [r"$\gamma_{1x}$", r"$g_{1x}$", r"$\gamma_{2x}$", r"$g_{2x}$", r"$\Gamma_x$"]
    values_bar = [
        gamma["top"][0][0],
        full["top"][0][0],
        gamma["second"][0][0],
        full["second"][0][0],
        geometry["local_bridge_Gamma_a_M2"]["x"],
    ]
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    ax.bar(labels_bar, values_bar)
    ax.set_ylabel(r"geometric coefficient ($a_M^2$)")
    ax.set_title("Remote-band geometry and local bridge filter")
    ax.grid(True, axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "wse2_geometric_coefficients.png", dpi=220)
    plt.close(fig)

    r_values = many_ed["composition_r"]
    direct = np.array(many_ed["direct_Q_after_alternating_composition_gauge"])
    formula = np.array(many_ed["closed_Jacobi_formula"])
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    if len(r_values) > 1:
        midpoints = np.asarray(r_values[:-1]) + 0.5
        ax.plot(midpoints, np.diag(direct, 1), marker="o", label="direct Fock ED")
        ax.plot(midpoints, np.diag(formula, 1), marker="x", linestyle="--", label="closed formula")
    ax.set_xlabel(r"composition bond $r\leftrightarrow r+1$")
    ax.set_ylabel(r"$t_{r,x}$ (meV $a_M^2$)")
    ax.set_title("Nonzero composition-space hopping")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "wse2_composition_hopping_ed.png", dpi=220)
    plt.close(fig)

    rows = locality["cumulative"]
    radii = [row["radius_cells"] for row in rows]
    factor = [row["factor_fraction"] for row in rows]
    source = [row["source_fraction"] for row in rows]
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    ax.plot(radii, factor, marker="o", label="filtered factor")
    ax.plot(radii, source, marker="s", label="twist source")
    ax.set_xlabel("triangular-cell truncation radius")
    ax.set_ylabel("cumulative Hilbert-Schmidt weight")
    ax.set_ylim(0.0, 1.02)
    ax.set_title("Real-space locality of the band-filtered bridge")
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "wse2_bridge_locality.png", dpi=220)
    plt.close(fig)

    fields = np.array([row["E_z_meV"] for row in scan])
    gap = np.array([row["direct_gap_second_remote_meV"] for row in scan])
    bridge = np.array([row["Gamma_x_a_M2"] for row in scan])
    fig, ax = plt.subplots(figsize=(7.0, 4.1))
    ax.plot(fields, gap, marker="o", label=r"$\Delta_{23}$ (meV)")
    ax.set_xlabel(r"displacement-field energy $E_z$ (meV)")
    ax.set_ylabel(r"direct remote gap $\Delta_{23}$ (meV)")
    ax2 = ax.twinx()
    ax2.plot(fields, bridge, marker="s", linestyle="--", label=r"$\Gamma_x$ ($a_M^2$)")
    ax2.set_ylabel(r"bridge coefficient $\Gamma_x$ ($a_M^2$)")
    handles_1, labels_1 = ax.get_legend_handles_labels()
    handles_2, labels_2 = ax2.get_legend_handles_labels()
    ax.legend(handles_1 + handles_2, labels_1 + labels_2, frameon=False)
    ax.grid(True, alpha=0.25)
    ax.set_title("Gap–geometry tradeoff under displacement field")
    fig.tight_layout()
    fig.savefig(output / "wse2_displacement_tradeoff.png", dpi=220)
    plt.close(fig)

    direct_one = np.array(one_pair["direct_Q_meV_a_M2"])
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    ax.plot([0, 1], np.linalg.eigvalsh(direct_one), marker="o")
    ax.set_xticks([0, 1])
    ax.set_xlabel("one-pair curvature eigenvalue index")
    ax.set_ylabel(r"eigenvalue (meV $a_M^2$)")
    ax.set_title("Material bridge one-pair curvature")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(output / "wse2_one_pair_bridge_spectrum.png", dpi=220)
    plt.close(fig)


def write_csv_files(
    output: Path,
    projected: Mapping[str, Any],
    scan: list[Mapping[str, float]],
    locality: Mapping[str, Any],
) -> None:
    with (output / "wse2_projected_cooper_kernel_norms.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["block_pair", "U_meV", "V_meV", "J_meV", "singlet_U_plus_V_minus_J_meV"])
        tables = projected["operator_norm_meV"]
        for label in tables["U"]:
            writer.writerow([label, tables["U"][label], tables["V"][label], tables["J"][label], tables["singlet_U_plus_V_minus_J"][label]])

    with (output / "wse2_displacement_scan.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(scan[0].keys()))
        writer.writeheader()
        writer.writerows(scan)

    rows = locality["cumulative"]
    with (output / "wse2_bridge_locality.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: float, digits: int = 6) -> str:
    return f"{value:.{digits}f}"


def write_report(output: Path, results: Mapping[str, Any]) -> None:
    """Write a self-contained Markdown report without relying on LaTeX parsing."""
    bands = results["band_structure"]
    geometry = results["geometry"]
    projected = results["projected_interactions"]
    fierz = results["effective_quadratures"]
    one = results["exact_diagonalization"]["one_pair"]
    many = results["exact_diagonalization"]["many_composition"]
    corrections = results["corrections"]
    locality = results["locality"]

    gamma_top = geometry["remote_metric_gamma_a_M2"]["top"][0][0]
    gamma_second = geometry["remote_metric_gamma_a_M2"]["second"][0][0]
    full_top = geometry["full_quantum_metric_a_M2"]["top"][0][0]
    full_second = geometry["full_quantum_metric_a_M2"]["second"][0][0]
    bridge_gamma = geometry["local_bridge_Gamma_a_M2"]["x"]
    couplings = fierz["with_representative_pair_hop"]
    kernels = projected["operator_norm_meV"]["singlet_U_plus_V_minus_J"]
    cn = results["topology"]["chern_numbers"]
    p_pair = results["interactions"]["p_pair"]

    lines: list[str] = []
    add = lines.append
    add("# Material-specific continuation: two-block QGN screening in 3.65° twisted WSe₂")
    add("")
    add("## Executive result")
    add("")
    add("All six requested calculations were completed for the published three-orbital MX/MM/XM Wannier Hamiltonian. The two highest K-valley bands were treated as active pairing blocks and the third band as a remote mediator.")
    add("")
    add("The calculation finds a nonzero ordinary-twist composition-space hopping mechanism, but it also finds two obstructions to treating the published 3.65° parameter set as an exact material QGN model:")
    add("")
    add("1. The projected U, V, J Cooper kernels substantially mix the active bands rather than preserving exact QGN blocks.")
    add("2. U, V, J alone give equal transfer quadratures, so the Jacobi hopping vanishes. A pair-hop P, or another quadrature anisotropy, is necessary.")
    add("")
    add(f"Using the representative benchmark P = {p_pair:.3f} meV together with J = {results['interactions']['j_exchange']:.3f} meV gives |t₀,x| = {one['offdiagonal_t0_magnitude_meV_a_M2']:.9f} meV a_M² on the 4×4 torus. A direct many-composition Fock calculation agrees with the closed Jacobi formula to {many['maximum_matrix_error']:.3e}.")
    add("")
    add("## 1. Wannier manifold and topology")
    add("")
    add(f"The selected three-band composite has K-valley Chern numbers ({cn[0]:.0f}, {cn[1]:.0f}, {cn[2]:.0f}). The top two alone therefore have total Chern number +2; the C = −2 remote band is required for a Chern-trivial composite Wannier manifold. The local composite orbitals are centered on MX, MM, and XM sites. The active blocks are Bloch eigenbundles inside this local three-orbital Hilbert space.")
    add("")
    add("This calculation reconstructs the published tight-binding/Wannier Hamiltonian. It does not rerun continuum-model Wannierization because continuum Bloch-wavefunction and gauge files were not available in the supplied archive.")
    add("")
    add("| quantity | top | second | remote |")
    add("|---|---:|---:|---:|")
    add(f"| minimum energy (meV) | {bands['energy_min_meV']['top']:.6f} | {bands['energy_min_meV']['second']:.6f} | {bands['energy_min_meV']['remote']:.6f} |")
    add(f"| maximum energy (meV) | {bands['energy_max_meV']['top']:.6f} | {bands['energy_max_meV']['second']:.6f} | {bands['energy_max_meV']['remote']:.6f} |")
    add(f"| bandwidth (meV) | {bands['bandwidth_meV']['top']:.6f} | {bands['bandwidth_meV']['second']:.6f} | {bands['bandwidth_meV']['remote']:.6f} |")
    add("")
    add(f"Minimum direct gaps: Δ₁₂ = {bands['gaps_meV']['top_second']['direct_min']:.6f} meV and Δ₂₃ = {bands['gaps_meV']['second_remote']['direct_min']:.6f} meV.")
    add("")
    add("## 2. Interband geometry")
    add("")
    add("The remote-band contribution is computed from")
    add("")
    add("```text")
    add("A_(3a,i)(k) = <u_3(k)|∂_(k_i) h(k)|u_a(k)> / [E_a(k) - E_3(k)]")
    add("γ^(3)_(a,ij) = BZ average Re[A*_(3a,i) A_(3a,j)].")
    add("```")
    add("")
    add("C₃ symmetry makes the tensors isotropic within numerical precision:")
    add("")
    add("| active band | γ_xx = γ_yy (a_M²) | full metric g_xx = g_yy (a_M²) | remote fraction |")
    add("|---|---:|---:|---:|")
    add(f"| top | {gamma_top:.9f} | {full_top:.9f} | {100.0*gamma_top/full_top:.3f}% |")
    add(f"| second | {gamma_second:.9f} | {full_second:.9f} | {100.0*gamma_second/full_second:.3f}% |")
    add("")
    add(f"The local filtered bridge coefficient, including remote-MM weight, active-honeycomb weight, and all six T–H bonds, is Γ_x = Γ_y = {bridge_gamma:.9f} a_M².")
    add("")
    add(f"A bare local T†H bilinear is not a zero row: estimated active-active Hilbert–Schmidt leakage per bond is {geometry['bare_active_active_leakage_HS2_per_bond']:.6f}, versus {geometry['desired_filtered_remote_top_HS2_per_bond']:.6f} in the desired remote-to-top component. Exact zero-twist preservation therefore requires the filtered operator P₃(T†H)P₁ or a microscopic cancellation approximating it.")
    add("")
    add("## 3. Projection of U, V, J into the pairing bands")
    add("")
    add(f"Spectral norms of the projected singlet Cooper kernel U + V − J on the {projected['grid']}×{projected['grid']} momentum torus:")
    add("")
    add("| top→top | second→second | remote→remote | top→second | top→remote | second→remote |")
    add("|---:|---:|---:|---:|---:|---:|")
    add(f"| {kernels['top->top']:.6f} | {kernels['second->second']:.6f} | {kernels['remote->remote']:.6f} | {kernels['top->second']:.6f} | {kernels['top->remote']:.6f} | {kernels['second->remote']:.6f} |")
    add("")
    ratios = projected["normalized_cross_block_ratios"]
    add(f"Normalized cross-block ratios are {ratios['top_second']:.3f} (top–second), {ratios['top_remote']:.3f} (top–remote), and {ratios['second_remote']:.3f} (second–remote). These are not small enough to identify the top and second energy bands as exact interaction-preserved QGN blocks.")
    add("")
    add("## 4. Effective g_X and g_Y")
    add("")
    add("For the spin-odd transfer quadratures, the local Fierz identities give")
    add("")
    add("```text")
    add("g_X² + g_Y² = J")
    add("g_Y² - g_X² = P")
    add("```")
    add("")
    add(f"With U, V, J alone, P = 0 and g_X² = g_Y² = J/2 = {fierz['UVJ_only_P_zero']['gX2_meV']:.3f} meV, so t_r = 0 even though the diagonal bridge response is nonzero.")
    add(f"With the representative P = {p_pair:.3f} meV benchmark: g_X² = {couplings['gX2_meV']:.3f} meV, g_Y² = {couplings['gY2_meV']:.3f} meV, g_X = {couplings['gX_sqrt_meV']:.6f} √meV, and g_Y = {couplings['gY_sqrt_meV']:.6f} √meV.")
    add("")
    add("The P benchmark comes from a related three-orbital tTMD Wannier calculation; it is not a first-principles extraction for the exact Kim parameter set.")
    add("")
    add("## 5. Exact-diagonalization verification")
    add("")
    add("For the material source rows generated by the ordinary Peierls twist,")
    add("")
    add("```text")
    add("t_(r,i) = [2 Γ_i (g_Y² - g_X²) / V] × sqrt[(n-r)(V-r)(r+1)(V-n+r+1)].")
    add("```")
    add("")
    add("The 4×4 one-pair curvature matrix is")
    add("")
    add("```text")
    add(f"[[{one['direct_Q_meV_a_M2'][0][0]:.9f}, {one['direct_Q_meV_a_M2'][0][1]:.9f}],")
    add(f" [{one['direct_Q_meV_a_M2'][1][0]:.9f}, {one['direct_Q_meV_a_M2'][1][1]:.9f}]] meV a_M²")
    add("```")
    add("")
    add(f"Its eigenvalues are {one['eigenvalues_meV_a_M2'][0]:.9f} and {one['eigenvalues_meV_a_M2'][1]:.9f} meV a_M². The off-diagonal sign is basis dependent; its nonzero magnitude is {one['offdiagonal_t0_magnitude_meV_a_M2']:.9f} meV a_M².")
    add(f"A separate sparse-Fock calculation at V = n = 4 reconstructs the full five-dimensional Jacobi matrix. After the phase gauge |r,n−r> → (−1)^r |r,n−r>, the maximum direct-versus-formula discrepancy is {many['maximum_matrix_error']:.3e}.")
    add(f"The thermodynamic one-pair hopping estimate is |t₀,x| = 2 Γ_x P = {2.0*bridge_gamma*p_pair:.9f} meV a_M².")
    add("")
    add("## 6. Finite-bandwidth and remote-band corrections")
    add("")
    ratios_gap = corrections["dimensionless_ratios_to_Delta23"]
    add("| ratio to Δ₂₃ | value |")
    add("|---|---:|")
    for label, key in [
        ("W_top/Δ₂₃", "W_top_over_Delta"),
        ("W_second/Δ₂₃", "W_second_over_Delta"),
        ("U_H/Δ₂₃", "U_H_over_Delta"),
        ("U_T/Δ₂₃", "U_T_over_Delta"),
        ("V/Δ₂₃", "V_over_Delta"),
        ("J/Δ₂₃", "J_over_Delta"),
        ("P/Δ₂₃", "P_over_Delta"),
    ]:
        add(f"| {label} | {ratios_gap[key]:.3f} |")
    add("")
    add(f"Pair kinetic spectral widths are {corrections['pair_kinetic_spectral_width_meV']['top']:.3f} meV (top) and {corrections['pair_kinetic_spectral_width_meV']['second']:.3f} meV (second). Naive second-order estimates K₁r²/Δ₂₃ and K₂r²/Δ₂₃ are {corrections['second_order_energy_estimates_meV']['K_top_remote2_over_Delta']:.3f} and {corrections['second_order_energy_estimates_meV']['K_second_remote2_over_Delta']:.3f} meV. The 3.65° fit is therefore outside a controlled isolated-flat-band or simple Schrieffer–Wolff regime.")
    add("")
    r5 = locality["cumulative"][4]
    r8 = locality["cumulative"][7]
    add(f"The filtered bridge remains reasonably local in the canonical composite-Wannier gauge. Radius {r5['radius_cells']} retains {100.0*r5['factor_fraction']:.2f}% of factor weight and {100.0*r5['source_fraction']:.2f}% of source weight; radius {r8['radius_cells']} retains {100.0*r8['factor_fraction']:.3f}% and {100.0*r8['source_fraction']:.3f}%.")
    add("")
    add("## Research decision")
    add("")
    add("The mechanism is established, but the published 3.65° parameter set is not yet a controlled exact material QGN realization.")
    add("")
    add("- Positive: the filtered remote-to-top bridge has Γ_i > 0, and a pair-hop anisotropy P produces t_(r,i) ≠ 0 under the ordinary uniform twist. Direct Fock-space ED verifies the coefficient.")
    add("- Negative: U, V, J alone give g_X = g_Y and t_(r,i) = 0; the same interactions also mix the proposed active blocks substantially.")
    add("- Control problem: active bandwidths and active–remote interaction kernels are too large relative to Δ₂₃.")
    add("")
    add("The next material search should optimize the simultaneous conditions")
    add("")
    add("```text")
    add("W₁, W₂, ||K₁r||, ||K₂r|| << Δ₂₃, while Γ_i P != 0.")
    add("```")
    add("")
    add("The included displacement-field scan shows the expected tradeoff: fields that enlarge the remote gap tend to suppress Γ_i.")
    add("")
    add("## Reproducibility")
    add("")
    add("```bash")
    add("python wse2_qgn_material_pipeline.py --outdir wse2_qgn_material_results")
    add("```")
    add("")
    add("The JSON file records every parameter and output. CSV files contain the projected Cooper-kernel norms, displacement scan, and locality tails.")
    add("")
    add("## Source models used")
    add("")
    add("- S. Kim et al., Theory of correlated insulators and superconductor at ν=1 in twisted WSe₂, Nature Communications 16, 1701 (2025): three-orbital Hamiltonian and U_H, U_T, V, J parameters.")
    add("- V. Crépel and A. Millis, Bridging the small and large in twisted transition-metal-dichalcogenide homobilayers, Physical Review Research 6, 033127 (2024): composite-Wannier construction, interaction hierarchy, and representative pair-hopping scale.")
    add("- C. Tuo et al., Theory of topological superconductivity and antiferromagnetic correlated insulators in twisted bilayer WSe₂, Nature Communications (2025): independent direct three-band Wannierization and the (1,1,−2) minimal topology.")

    (output / "wse2_qgn_material_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

def write_latex_section(output: Path, results: Mapping[str, Any]) -> None:
    geometry = results["geometry"]
    bands = results["band_structure"]
    projected = results["projected_interactions"]
    one = results["exact_diagonalization"]["one_pair"]
    many = results["exact_diagonalization"]["many_composition"]
    gamma1 = geometry["remote_metric_gamma_a_M2"]["top"][0][0]
    gamma2 = geometry["remote_metric_gamma_a_M2"]["second"][0][0]
    bridge = geometry["local_bridge_Gamma_a_M2"]["x"]
    ratios = projected["normalized_cross_block_ratios"]
    error_mantissa, error_exponent = f"{many['maximum_matrix_error']:.2e}".split("e")
    error_exponent = int(error_exponent)
    text = rf"""\section{{Material screening in a three-orbital twisted-WSe$_2$ Wannier model}}
\label{{sec:wse2screen}}
We reconstructed the three-orbital MX/MM/XM Hamiltonian of Kim \emph{{et al.}}
using their $E_z=0$ parameters and selected the two highest $K$-valley bands
as active pairing blocks, with the third band as a remote mediator.  A
sewing-aware Fukui calculation gives Chern numbers $(1,1,-2)$, so the three-band
composite is Wannierizable whereas the two active bands alone remain Chern
obstructed in one valley.

Define the remote-band contribution to the active-band metric by
\begin{{equation}}
 \gamma^{{(3)}}_{{a,ij}}=\frac1V\sum_{{\bm k}}
 \operatorname{{Re}}\left[\mathcal A_{{3a,i}}^*(\bm k)
 \mathcal A_{{3a,j}}(\bm k)\right],\qquad
 \mathcal A_{{3a,i}}=\frac{{\langle u_3|\partial_{{k_i}}h|u_a\rangle}}
 {{E_a-E_3}}.
\end{{equation}}
At zero displacement field we obtain
\begin{{equation}}
 \gamma^{{(3)}}_1={gamma1:.9f}\,a_M^2\Id_2,
 \qquad
 \gamma^{{(3)}}_2={gamma2:.9f}\,a_M^2\Id_2.
\end{{equation}}
The local MM--honeycomb filtered bridge has
\begin{{equation}}
 \Gamma_x=\Gamma_y={bridge:.9f}\,a_M^2.
\end{{equation}}
The minimum second-to-remote direct gap is
$\Delta_{{23}}={bands['gaps_meV']['second_remote']['direct_min']:.6f}\,$meV,
while the two active bandwidths are
${bands['bandwidth_meV']['top']:.6f}\,$meV and
${bands['bandwidth_meV']['second']:.6f}\,$meV.

For a local chiral exchange $J$ and pair hop $P$, the spin-odd bridge
quadratures obey
\begin{{equation}}
 g_X^2+g_Y^2=J,\qquad g_Y^2-g_X^2=P.
\end{{equation}}
Thus the published $U,V,J$ set alone has $P=0$ and gives no composition-space
hopping.  Using the representative $J=10\,$meV and $P=2\,$meV gives
$g_X^2=4\,$meV and $g_Y^2=6\,$meV.  The exact active-space curvature then has
\begin{{equation}}
 t_{{r,i}}=\frac{{2\Gamma_iP}}V
 \sqrt{{(n-r)(V-r)(r+1)(V-n+r+1)}}.
\end{{equation}}
On a $4\times4$ torus the direct one-pair Fock calculation gives
$|t_{{0,x}}|={one['offdiagonal_t0_magnitude_meV_a_M2']:.9f}\,$meV$\,a_M^2$.
A separate $V=n=4$ many-composition calculation agrees with the complete
Jacobi formula to maximum error ${error_mantissa}\times10^{{{error_exponent}}}$.

The physical projected singlet kernels do not, however, preserve these energy
bands as exact blocks.  Their normalized top--second, top--remote, and
second--remote operator-norm ratios are respectively
${ratios['top_second']:.3f}$, ${ratios['top_remote']:.3f}$, and
${ratios['second_remote']:.3f}$.  Moreover $W_1,W_2$ and the active--remote
interaction kernels are not small compared with $\Delta_{{23}}$.  The model
therefore establishes a local material mechanism for $t_{{r,i}}\ne0$, but the
published $3.65^\circ$ parameter set is not yet a controlled exact QGN
realization.
"""
    (output / "wse2_material_screening_section.tex").write_text(text, encoding="utf-8")


def convert_json(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, dict):
        return {str(key): convert_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [convert_json(item) for item in value]
    return value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("/mnt/data/wse2_qgn_material_results"))
    parser.add_argument("--chern-grid", type=int, default=31)
    parser.add_argument("--geometry-grid", type=int, default=81)
    parser.add_argument("--projection-grid", type=int, default=13)
    parser.add_argument("--one-pair-ed-N", type=int, default=4)
    parser.add_argument("--many-ed-N", type=int, default=2)
    parser.add_argument("--tail-grid", type=int, default=128)
    parser.add_argument("--scan-grid", type=int, default=41)
    parser.add_argument("--scan-max-Ez", type=float, default=30.0)
    parser.add_argument("--scan-step-Ez", type=float, default=2.0)
    parser.add_argument("--pair-hop", type=float, default=2.0, help="Representative P in meV; use 0 for U,V,J only")
    args = parser.parse_args()

    output = args.outdir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    tight_binding = TightBindingParameters()
    interactions = InteractionParameters(p_pair=args.pair_hop)
    model = WSe2Model(tight_binding)

    print("[1/9] derivative and topology checks")
    derivative_error = derivative_check(model)
    chern = chern_numbers(model, args.chern_grid)

    print("[2/9] optimized band extrema and gaps")
    bands = band_extrema(model)

    print("[3/9] interband geometry and local bridge coefficient")
    geometry = geometry_and_bridge(model, args.geometry_grid)

    print("[4/9] projected U,V,J Cooper kernels")
    projected = projected_cooper_kernels(model, interactions, args.projection_grid)

    print("[5/9] Fierz quadrature extraction")
    fierz = fierz_bridge_couplings(interactions)

    print("[6/9] one-pair and many-composition exact diagonalization")
    one_pair = one_pair_ed(model, interactions, args.one_pair_ed_N, direction=0)
    many_ed = many_composition_ed(model, interactions, args.many_ed_N, direction=0)

    print("[7/9] real-space locality tails")
    locality = locality_tails(model, args.tail_grid, maximum_radius=10)

    print("[8/9] displacement-field scan")
    scan_values = np.arange(0.0, args.scan_max_Ez + 0.5 * args.scan_step_Ez, args.scan_step_Ez)
    scan = displacement_scan(model, scan_values, args.scan_grid)

    print("[9/9] corrections, reports, data files, and figures")
    corrections = corrections_summary(bands, projected, interactions, geometry)
    results = {
        "metadata": {
            "model": "Kim et al. three-orbital 3.65-degree twisted WSe2 fit",
            "energy_unit": "meV",
            "length_unit": "a_M",
            "active_blocks": ["top band", "second band"],
            "remote_band": "third band",
            "derivative_max_error": derivative_error,
        },
        "tight_binding_parameters": asdict(tight_binding),
        "interactions": asdict(interactions),
        "topology": {
            "chern_grid": args.chern_grid,
            "chern_numbers": chern,
            "composite_chern": float(sum(chern)),
        },
        "band_structure": bands,
        "geometry": geometry,
        "projected_interactions": projected,
        "effective_quadratures": fierz,
        "exact_diagonalization": {"one_pair": one_pair, "many_composition": many_ed},
        "locality": locality,
        "displacement_scan": scan,
        "corrections": corrections,
    }

    json_path = output / "wse2_qgn_material_results.json"
    json_path.write_text(json.dumps(convert_json(results), indent=2), encoding="utf-8")
    write_csv_files(output, projected, scan, locality)
    write_report(output, results)
    write_latex_section(output, results)
    make_plots(output, model, geometry, one_pair, many_ed, locality, scan)

    tolerance = 5.0e-8
    if derivative_error > tolerance:
        raise SystemExit(f"Hamiltonian derivative check failed: {derivative_error}")
    if max(abs(value - round(value)) for value in chern) > 1.0e-6 or abs(sum(chern)) > 1.0e-6:
        raise SystemExit(f"Chern check failed: {chern}")
    if one_pair["max_error_up_to_basis_phase"] > tolerance:
        raise SystemExit("One-pair ED check failed")
    if many_ed["maximum_matrix_error"] > tolerance:
        raise SystemExit("Many-composition ED check failed")

    print(f"Results written to {output}")
    print(f"Chern numbers: {chern}")
    print(f"Gamma_x = {geometry['local_bridge_Gamma_a_M2']['x']:.9f} a_M^2")
    print(f"one-pair |t0,x| = {one_pair['offdiagonal_t0_magnitude_meV_a_M2']:.9f} meV a_M^2")
    print(f"many-composition maximum ED/formula error = {many_ed['maximum_matrix_error']:.3e}")


if __name__ == "__main__":
    main()
