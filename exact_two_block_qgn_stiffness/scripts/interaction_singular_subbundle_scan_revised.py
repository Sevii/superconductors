#!/usr/bin/env python3
"""Scan twisted-bilayer WSe2 continuum parameters using interaction singular subbundles.

The model is the conventional lowest-order two-layer continuum Hamiltonian for
the +K valence valley.  Instead of declaring the two highest energy bands to be
the two QGN blocks, this program:

1. restricts local density and interlayer-transfer interaction channels to the
   top-three-energy-band composite manifold;
2. diagonalizes the positive channel Gram operator S(k);
3. treats its isolated right-singular eigenlines as candidate pairing blocks;
4. selects two active singular subbundles by low-energy capture, singular
   isolation, active/remote expectation gap, interaction leakage, and a
   finite-bandwidth control ratio;
5. evaluates the active-remote geometric bridge and the resulting one-pair
   composition-hopping coefficient t0 = 2 P Gamma.

The archived scan used:
  exploratory: 2,800 coarse points, 120 shell-3 refinements;
  conservative: 1,800 coarse points, all valid shell-3 refinements;
  seed: 20260725.

Energies are in meV and lengths in Angstrom unless explicitly normalized by
the moire lattice constant a_M.

Examples
--------
Quick verification of the three archived candidates:
    python interaction_singular_subbundle_scan.py --mode verify

Rerun the full random scans:
    python interaction_singular_subbundle_scan.py --mode full \
        --output interaction_singular_scan_rerun

A reduced exploratory scan:
    python interaction_singular_subbundle_scan.py --mode full \
        --exploratory-points 300 --conservative-points 200 --refine-count 30
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.linalg as la

HBAR2_2ME = 3809.98212  # meV Angstrom^2


@dataclass
class LOParams:
    theta_deg: float = 3.65
    V: float = 9.0
    psi_deg: float = 128.0
    w: float = 18.0
    mstar: float = 0.45
    vz: float = 0.0
    a0: float = 3.317
    shell: int = 2


class LOContinuum:
    """Lowest-order continuum model in a hexagonal plane-wave basis."""

    def __init__(self, p: LOParams):
        self.p = p
        theta = np.deg2rad(p.theta_deg)
        k_mono = 4.0 * np.pi / (3.0 * p.a0)
        q = 2.0 * k_mono * np.sin(theta / 2.0)
        self.qmag = q
        self.q1 = q * np.array([0.0, -1.0])
        self.q2 = q * np.array([np.sqrt(3.0) / 2.0, 0.5])
        self.q3 = q * np.array([-np.sqrt(3.0) / 2.0, 0.5])
        self.b1 = self.q2 - self.q1
        self.b2 = self.q3 - self.q1

        indices: list[tuple[int, int]] = []
        s = p.shell
        for m in range(-s, s + 1):
            for n in range(-s, s + 1):
                if max(abs(m), abs(n), abs(m + n)) <= s:
                    indices.append((m, n))
        self.inds = indices
        self.index = {mn: i for i, mn in enumerate(indices)}
        self.N = len(indices)
        self.G = np.array([m * self.b1 + n * self.b2 for m, n in indices])
        self.dim = 2 * self.N
        self.alpha = HBAR2_2ME / p.mstar
        self.pot_q = [(1, 0), (-1, 1), (0, -1)]

    def hamiltonian(self, k: np.ndarray) -> np.ndarray:
        p = self.p
        h = np.zeros((self.dim, self.dim), dtype=complex)

        p_top = k[None, :] + self.G
        p_bottom = k[None, :] + self.G + self.q1
        ids = np.arange(self.N)
        h[ids, ids] = -self.alpha * np.sum(p_top * p_top, axis=1) + p.vz / 2.0
        h[self.N + ids, self.N + ids] = (
            -self.alpha * np.sum(p_bottom * p_bottom, axis=1) - p.vz / 2.0
        )

        psi = np.deg2rad(p.psi_deg)
        for layer, sign in ((0, +1), (1, -1)):
            offset = layer * self.N
            for i, (m, n) in enumerate(self.inds):
                for dm, dn in self.pot_q:
                    j = self.index.get((m - dm, n - dn))
                    if j is not None:
                        h[offset + i, offset + j] += p.V * np.exp(1.0j * sign * psi)
                    j2 = self.index.get((m + dm, n + dn))
                    if j2 is not None:
                        h[offset + i, offset + j2] += p.V * np.exp(-1.0j * sign * psi)

        for i, (m, n) in enumerate(self.inds):
            for dm, dn in ((0, 0), (1, 0), (0, 1)):
                j = self.index.get((m + dm, n + dn))
                if j is not None:
                    h[i, self.N + j] += p.w
                    h[self.N + j, i] += p.w
        return h

    # Backward-compatible name used in the archived notebook.
    H = hamiltonian

    def eig(self, k: np.ndarray, n: int = 6) -> tuple[np.ndarray, np.ndarray]:
        if n <= 0 or n > self.dim:
            raise ValueError(f"requested {n} bands for a {self.dim}-dimensional basis")
        values, vectors = la.eigh(
            self.hamiltonian(k), subset_by_index=[self.dim - n, self.dim - 1]
        )
        return values[::-1], vectors[:, ::-1]


def kmesh(model: LOContinuum, nk: int) -> np.ndarray:
    if nk < 2:
        raise ValueError("nk must be at least 2")
    points = []
    for i in range(nk):
        for j in range(nk):
            u = (i + 0.5) / nk - 0.5
            v = (j + 0.5) / nk - 0.5
            points.append(u * model.b1 + v * model.b2)
    return np.asarray(points)


def channel_operators(
    model: LOContinuum,
    eta: float = 0.3,
    density_asym: float = 0.0,
    path_anis: Sequence[float] = (1.0, 1.0, 1.0),
) -> list[tuple[float, np.ndarray, str]]:
    """Local channel operators entering the positive interaction Gram map."""
    if len(path_anis) != 3 or sum(path_anis) <= 0:
        raise ValueError("path_anis must contain three nonnegative weights")
    n, dim = model.N, model.dim
    p_top = np.zeros((dim, dim), dtype=complex)
    p_bottom = np.zeros((dim, dim), dtype=complex)
    p_top[:n, :n] = np.eye(n)
    p_bottom[n:, n:] = np.eye(n)

    operators: list[tuple[float, np.ndarray, str]] = [
        (1.0 + density_asym, p_top, "rho_t"),
        (1.0 - density_asym, p_bottom, "rho_b"),
    ]

    norm = float(sum(path_anis))
    for j, ((dm, dn), path_weight) in enumerate(
        zip(((0, 0), (1, 0), (0, 1)), path_anis)
    ):
        transfer = np.zeros((dim, dim), dtype=complex)
        for i, (m, n_idx) in enumerate(model.inds):
            jj = model.index.get((m + dm, n_idx + dn))
            if jj is not None:
                transfer[i, model.N + jj] = 1.0
        bx = (transfer + transfer.conj().T) / np.sqrt(2.0)
        by = -1.0j * (transfer - transfer.conj().T) / np.sqrt(2.0)
        weight = eta * float(path_weight) / norm
        operators.extend(((weight, bx, f"X{j}"), (weight, by, f"Y{j}")))
    return operators


def block_at(
    model: LOContinuum,
    k: np.ndarray,
    eta: float = 0.3,
    density_asym: float = 0.0,
    nb: int = 5,
) -> dict[str, Any]:
    """Interaction singular subbundles inside the top-three energy manifold."""
    values, energy_vectors = model.eig(k, n=nb)
    values3 = values[:3]
    u3 = energy_vectors[:, :3]

    gram = np.zeros((3, 3), dtype=complex)
    frames: list[tuple[float, np.ndarray, str, np.ndarray]] = []
    for weight, operator, name in channel_operators(model, eta, density_asym):
        frame = u3.conj().T @ operator @ u3
        gram += weight * (frame.conj().T @ frame)
        frames.append((weight, frame, name, operator))

    singular_values, q = np.linalg.eigh(gram)
    order = np.argsort(singular_values)[::-1]
    singular_values = singular_values[order]
    q = q[:, order]
    singular_vectors = u3 @ q
    h_block = q.conj().T @ np.diag(values3) @ q

    weights = np.abs(q) ** 2
    energy_expectation = np.sum(weights * values3[:, None], axis=0)
    energy_second = np.sum(weights * values3[:, None] ** 2, axis=0)
    energy_variance = np.maximum(0.0, energy_second - energy_expectation**2)

    leakage = np.zeros((3, 3))
    leakage_energy = np.zeros((3, 3))
    for weight, frame, _, _ in frames:
        transformed = q.conj().T @ frame @ q
        leakage += weight * np.abs(transformed) ** 2
        leakage_energy += weight * np.abs(frame) ** 2

    return {
        "vals": values,
        "U": energy_vectors,
        "s": singular_values,
        "q": q,
        "Q": singular_vectors,
        "Hb": h_block,
        "eps": energy_expectation,
        "var": energy_variance,
        "L": leakage,
        "Le": leakage_energy,
        "Fs": frames,
    }


def coarse_metrics(
    params: Mapping[str, float],
    nk: int = 5,
    eta: float = 0.3,
    density_asym: float = 0.0,
    shell: int = 2,
) -> dict[str, Any]:
    """Evaluate all three possible two-active/one-remote singular partitions."""
    p = LOParams(
        theta_deg=float(params["theta"]),
        V=float(params["V"]),
        psi_deg=float(params["psi"]),
        w=float(params["w"]),
        mstar=float(params["mstar"]),
        vz=float(params["vz"]),
        shell=shell,
    )
    model = LOContinuum(p)
    data = [block_at(model, k, eta, density_asym, nb=5) for k in kmesh(model, nk)]

    singular = np.array([d["s"] for d in data])
    eps = np.array([d["eps"] for d in data])
    variance = np.array([d["var"] for d in data])
    h_block = np.array([d["Hb"] for d in data])
    values = np.array([d["vals"] for d in data])
    leakage = np.array([d["L"] for d in data])
    leakage_energy = np.array([d["Le"] for d in data])
    q = np.array([d["q"] for d in data])

    relative_gaps = np.stack(
        (
            (singular[:, 0] - singular[:, 1]) / np.maximum(singular[:, 0], 1e-12),
            (singular[:, 1] - singular[:, 2]) / np.maximum(singular[:, 1], 1e-12),
        ),
        axis=1,
    )
    isolation = float(np.min(relative_gaps))

    pair_metrics: list[dict[str, Any]] = []
    for pair in combinations(range(3), 2):
        remote = next(i for i in range(3) if i not in pair)
        a, b = pair
        capture = float(
            np.mean(
                (
                    np.sum(np.abs(q[:, :2, a]) ** 2, axis=1)
                    + np.sum(np.abs(q[:, :2, b]) ** 2, axis=1)
                )
                / 2.0
            )
        )
        widths = np.ptp(eps[:, [a, b]], axis=0)
        remote_gap = float(
            np.min(np.minimum(eps[:, a] - eps[:, remote], eps[:, b] - eps[:, remote]))
        )
        mix_ab = float(np.max(np.abs(h_block[:, a, b])))
        mix_ar = float(
            np.max(
                np.maximum(
                    np.abs(h_block[:, a, remote]),
                    np.abs(h_block[:, b, remote]),
                )
            )
        )
        variance_rms = np.sqrt(np.mean(variance[:, [a, b]], axis=0))

        def normalized_leakage(x: int, y: int) -> np.ndarray:
            return leakage[:, x, y] / np.sqrt(
                np.maximum(singular[:, x] * singular[:, y], 1e-15)
            )

        leak_ab = float(np.mean(normalized_leakage(a, b)))
        leak_ar = float(
            max(
                np.mean(normalized_leakage(a, remote)),
                np.mean(normalized_leakage(b, remote)),
            )
        )
        denominator = 1.0 + widths.mean() + 2.0 * mix_ab + 2.0 * mix_ar + 5.0 * leak_ar
        score = (
            max(isolation, 0.0)
            * capture**2
            * max(remote_gap, 0.0)
            * max(0.0, 1.0 - leak_ar)
            / denominator
        )
        control_ratio = (
            max(widths.max(), 2.0 * mix_ab, 2.0 * mix_ar, variance_rms.max())
            / max(remote_gap, 1e-6)
            if remote_gap > 0
            else 1e6
        )
        pair_metrics.append(
            {
                "active_pair": pair,
                "remote_index": remote,
                "capture_top2": capture,
                "block_width1": float(widths[0]),
                "block_width2": float(widths[1]),
                "remote_gap": remote_gap,
                "energy_mix12_max": mix_ab,
                "energy_mix_ar_max": mix_ar,
                "energy_var1_rms": float(variance_rms[0]),
                "energy_var2_rms": float(variance_rms[1]),
                "interaction_leak12": leak_ab,
                "interaction_leak_ar": leak_ar,
                "control_ratio": float(control_ratio),
                "score": float(score),
            }
        )

    best = max(pair_metrics, key=lambda item: item["score"])
    if best["score"] == 0:
        best = max(
            pair_metrics,
            key=lambda item: (
                item["remote_gap"],
                item["capture_top2"],
                -item["interaction_leak_ar"],
            ),
        )

    diagonal_energy = np.diagonal(leakage_energy, axis1=1, axis2=2)
    energy_leaks: dict[str, float] = {}
    for a, b in combinations(range(3), 2):
        energy_leaks[f"{a}{b}"] = float(
            np.mean(
                leakage_energy[:, a, b]
                / np.sqrt(np.maximum(diagonal_energy[:, a] * diagonal_energy[:, b], 1e-15))
            )
        )

    energy_widths = np.ptp(values[:, :4], axis=0)
    output: dict[str, Any] = {
        "valid": True,
        "theta": p.theta_deg,
        "V": p.V,
        "psi": p.psi_deg,
        "w": p.w,
        "mstar": p.mstar,
        "vz": p.vz,
        "eta": eta,
        "density_asym": density_asym,
        "singular_gap12_min": float(np.min(relative_gaps[:, 0])),
        "singular_gap23_min": float(np.min(relative_gaps[:, 1])),
        "singular_gap_min": isolation,
        "energy_gap12": float(np.min(values[:, 0] - values[:, 1])),
        "energy_gap23": float(np.min(values[:, 1] - values[:, 2])),
        "energy_gap34": float(np.min(values[:, 2] - values[:, 3])),
        "energy_bw_top": float(energy_widths[0]),
        "energy_bw_second": float(energy_widths[1]),
        "energy_bw_third": float(energy_widths[2]),
        "s_mean": np.mean(singular, axis=0),
        "s_min": np.min(singular, axis=0),
        "s_max": np.max(singular, axis=0),
        "overlap_matrix": np.mean(np.abs(q) ** 2, axis=0),
        "pair_metrics": pair_metrics,
        "energybasis_leaks": energy_leaks,
    }
    output.update(best)
    return output


def transfer_operators(model: LOContinuum) -> list[np.ndarray]:
    transfers: list[np.ndarray] = []
    for dm, dn in ((0, 0), (1, 0), (0, 1)):
        operator = np.zeros((model.dim, model.dim), dtype=complex)
        for i, (m, n) in enumerate(model.inds):
            j = model.index.get((m + dm, n + dn))
            if j is not None:
                operator[i, model.N + j] = 1.0
        transfers.append(operator)
    return transfers


def fine_geometry(
    params: Mapping[str, float],
    eta: float,
    density_asym: float,
    active_pair: Sequence[int] = (0, 1),
    shell: int = 3,
    nk: int = 9,
    delta_frac: float = 2e-3,
    p_pair: float = 2.0,
) -> dict[str, Any]:
    """Gauge-invariant finite-difference geometry of the singular subbundles.

    The archived implementation phase-aligned singular vectors and then
    differentiated them.  This revision finite-differences rank-one projectors,
    which is gauge invariant and remains stable when the singular gap is small.

    ``Gamma_proxy`` is a screening diagnostic, not by itself a physical twist
    source.  For directed mechanism (source,target) it is

      < Tr[P_r (d_i P_source)^2] * max_l(
          |<q_r|T_l|q_target>|^2 + |<q_r|T_l^dag|q_target>|^2 ) >_k,

    summed over (a,b) and (b,a), and divided by a_M^2.
    """
    p = LOParams(
        theta_deg=float(params["theta"]),
        V=float(params["V"]),
        psi_deg=float(params["psi"]),
        w=float(params["w"]),
        mstar=float(params["mstar"]),
        vz=float(params["vz"]),
        shell=shell,
    )
    model = LOContinuum(p)
    a, b = int(active_pair[0]), int(active_pair[1])
    remote = next(i for i in range(3) if i not in (a, b))
    transfers = transfer_operators(model)
    area = abs(float(model.b1[0] * model.b2[1] - model.b1[1] * model.b2[0]))
    delta = delta_frac * np.linalg.norm(model.b1)

    gamma = np.zeros((2, 2))
    metric = np.zeros((3, 2))
    berry = np.zeros(3)
    local_gamma = np.zeros((2, 2))
    path_weights = np.zeros((2, 3))
    singular_gaps: list[tuple[float, float]] = []

    for k in kmesh(model, nk):
        center = block_at(model, k, eta, density_asym, nb=5)
        q0 = center["Q"]
        p0 = [np.outer(q0[:, j], np.conjugate(q0[:, j])) for j in range(3)]
        singular = center["s"]
        singular_gaps.append(
            (
                float((singular[0] - singular[1]) / singular[0]),
                float((singular[1] - singular[2]) / singular[1]),
            )
        )

        dprojectors: list[list[np.ndarray]] = []
        for direction in range(2):
            shift = np.zeros(2)
            shift[direction] = delta
            plus = block_at(model, k + shift, eta, density_asym, nb=5)
            minus = block_at(model, k - shift, eta, density_asym, nb=5)
            dps: list[np.ndarray] = []
            for j in range(3):
                pp = np.outer(plus["Q"][:, j], np.conjugate(plus["Q"][:, j]))
                pm = np.outer(minus["Q"][:, j], np.conjugate(minus["Q"][:, j]))
                dp = (pp - pm) / (2.0 * delta)
                dps.append(dp)
                metric[j, direction] += 0.5 * float(np.trace(dp @ dp).real)
            dprojectors.append(dps)

        dx, dy = dprojectors
        for j in range(3):
            commutator = dx[j] @ dy[j] - dy[j] @ dx[j]
            berry[j] += float(np.real(-1.0j * np.trace(p0[j] @ commutator)))

        for mechanism, (source, target) in enumerate(((a, b), (b, a))):
            weights = []
            for transfer in transfers:
                z1 = np.vdot(q0[:, remote], transfer @ q0[:, target])
                z2 = np.vdot(q0[:, remote], transfer.conj().T @ q0[:, target])
                weights.append(float(abs(z1) ** 2 + abs(z2) ** 2))
            max_path_weight = max(weights)

            for direction in range(2):
                dp = dprojectors[direction][source]
                geometric = float(np.trace(p0[remote] @ dp @ dp).real)
                # Roundoff can produce tiny negative values although the exact
                # projector expression is nonnegative.
                geometric = max(geometric, 0.0)
                gamma[mechanism, direction] += geometric
                local_gamma[mechanism, direction] += geometric * max_path_weight

            for path_index, weight in enumerate(weights):
                path_weights[mechanism, path_index] += weight

    npoints = nk * nk
    gamma /= npoints
    metric /= npoints
    local_gamma /= npoints
    path_weights /= npoints
    berry = np.asarray(berry) * area / npoints / (2.0 * np.pi)

    a_m = p.a0 / (2.0 * np.sin(np.deg2rad(p.theta_deg) / 2.0))
    gamma_dimless = gamma / a_m**2
    metric_dimless = metric / a_m**2
    local_dimless = local_gamma / a_m**2

    gamma_sum = np.sum(local_dimless, axis=0)
    gamma_max = np.max(local_dimless, axis=0)
    ideal_sum = np.sum(gamma_dimless, axis=0)
    return {
        "geometry_method": "central finite differences of rank-one projectors",
        "Gamma_definition": (
            "sum over directed active mechanisms of BZ-average "
            "Tr[P_r(d_i P_source)^2] times max over local transfer paths of "
            "(|<q_r|T_l|q_target>|^2+|<q_r|T_l^dag|q_target>|^2), divided by a_M^2"
        ),
        "gamma_active_remote_aM2": gamma_dimless,
        "metric_aM2": metric_dimless,
        "local_gamma_aM2": local_dimless,
        "Gamma_proxy_sum_aM2": gamma_sum,
        "Gamma_proxy_max_aM2": gamma_max,
        # Backward-compatible aliases.  They remain geometric screening proxies.
        "Gamma_sum_aM2": gamma_sum,
        "Gamma_max_aM2": gamma_max,
        "ideal_gamma_sum_aM2": ideal_sum,
        "path_weights": path_weights,
        "chern_approx": berry,
        "conditional_t0_proxy_sum_meV_aM2": 2.0 * p_pair * gamma_sum,
        "conditional_t0_proxy_max_meV_aM2": 2.0 * p_pair * gamma_max,
        "t0_sum_meV_aM2": 2.0 * p_pair * gamma_sum,
        "t0_max_meV_aM2": 2.0 * p_pair * gamma_max,
        "singular_gap_min": float(np.min(singular_gaps)),
        "aM_A": a_m,
        "delta_Ainv": delta,
        "delta_frac": delta_frac,
    }

def augmented_score(result: Mapping[str, Any]) -> float:
    gap34 = max(float(result["energy_gap34"]), 0.0)
    control = float(result["control_ratio"])
    return (
        float(result["score"])
        * gap34
        / (1.0 + gap34)
        / (1.0 + max(control - 1.0, 0.0))
    )


def is_valid(result: Mapping[str, Any]) -> bool:
    return (
        float(result["remote_gap"]) > 0
        and float(result["singular_gap_min"]) > 0.02
        and float(result["capture_top2"]) > 0.9
    )


def is_controlled(result: Mapping[str, Any]) -> bool:
    return (
        is_valid(result)
        and float(result["control_ratio"]) < 1.0
        and float(result["energy_gap34"]) > 1.0
        and float(result["singular_gap_min"]) > 0.05
        and float(result["capture_top2"]) > 0.95
    )


def parameter_dict(result: Mapping[str, Any]) -> dict[str, float]:
    return {
        key: float(result[key])
        for key in ("theta", "V", "psi", "w", "mstar", "vz")
    }


def angle_scan(
    template: Mapping[str, float],
    eta: float,
    density_asym: float,
    thetas: Iterable[float],
    shell: int = 3,
    nk: int = 7,
    geometry: bool = True,
) -> list[dict[str, Any]]:
    output = []
    for theta in thetas:
        params = dict(template)
        params["theta"] = float(theta)
        result = coarse_metrics(params, nk=nk, eta=eta, density_asym=density_asym, shell=shell)
        if geometry:
            geom = fine_geometry(
                params,
                eta,
                density_asym,
                result["active_pair"],
                shell=2,
                nk=5,
            )
            result["Gamma_min"] = float(np.min(geom["Gamma_sum_aM2"]))
            result["t0_P2_min"] = float(np.min(geom["t0_sum_meV_aM2"]))
        output.append(result)
    return output


def jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, tuple):
        return [jsonable(item) for item in value]
    if isinstance(value, list):
        return [jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (str, bool, int, float)) or value is None:
        return value
    return repr(value)


def flat_row(result: Mapping[str, Any]) -> dict[str, Any]:
    row: dict[str, Any] = {}
    for key, value in result.items():
        if isinstance(value, (str, bool, int, float, np.integer, np.floating)):
            row[key] = jsonable(value)
        elif isinstance(value, tuple):
            row[key] = "-".join(str(item) for item in value)
        else:
            row[key] = json.dumps(jsonable(value), separators=(",", ":"))
    return row


ARCHIVED_CANDIDATES: dict[str, dict[str, float]] = {
    "balanced": {
        "theta": 2.0711737480196906,
        "V": 14.015497868832437,
        "psi": 147.75480281472278,
        "w": 8.978910378844677,
        "mstar": 0.6095970039206522,
        "vz": -5.726330509642874,
        "eta": 0.29118751954084376,
        "density_asym": 0.27456065782160355,
    },
    "max_hopping": {
        "theta": 2.2855807969513804,
        "V": 13.20085339450338,
        "psi": 149.30808764913382,
        "w": 17.650847080610927,
        "mstar": 0.6984796553180537,
        "vz": -0.7837588646244242,
        "eta": 0.3879556759227445,
        "density_asym": 0.10927877619375825,
    },
    "most_controlled": {
        "theta": 2.012049592666495,
        "V": 13.993781518700583,
        "psi": 130.07772318208907,
        "w": 9.839545951526041,
        "mstar": 0.6798350180504138,
        "vz": 2.2939117094026393,
        "eta": 0.7829188168564329,
        "density_asym": -0.15310960598951856,
    },
}


def verify_archived(output: Path) -> dict[str, Any]:
    output.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {}
    rows = []
    for name, archived in ARCHIVED_CANDIDATES.items():
        params = {key: archived[key] for key in ("theta", "V", "psi", "w", "mstar", "vz")}
        eta = archived["eta"]
        density_asym = archived["density_asym"]

        metric3 = coarse_metrics(params, nk=9, eta=eta, density_asym=density_asym, shell=3)
        geometry3 = fine_geometry(
            params,
            eta,
            density_asym,
            metric3["active_pair"],
            shell=3,
            nk=7,
        )
        payload[name] = {
            "parameters": archived,
            "shell3_metrics": jsonable(metric3),
            "shell3_geometry": jsonable(geometry3),
        }
        rows.append(
            {
                "candidate": name,
                **archived,
                "active_pair": "-".join(map(str, metric3["active_pair"])),
                "remote_gap": metric3["remote_gap"],
                "singular_gap_min": metric3["singular_gap_min"],
                "capture_top2": metric3["capture_top2"],
                "control_ratio": metric3["control_ratio"],
                "energy_gap34": metric3["energy_gap34"],
                "Gamma_x_aM2": geometry3["Gamma_sum_aM2"][0],
                "Gamma_y_aM2": geometry3["Gamma_sum_aM2"][1],
                "t0_x_P2_meV_aM2": geometry3["t0_sum_meV_aM2"][0],
                "t0_y_P2_meV_aM2": geometry3["t0_sum_meV_aM2"][1],
            }
        )
    (output / "verified_archived_candidates.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )
    pd.DataFrame(rows).to_csv(output / "verified_archived_candidates.csv", index=False)
    return payload


def draw_parameter(rng: np.random.Generator, conservative: bool) -> tuple[dict[str, float], float, float]:
    if conservative:
        params = {
            "theta": rng.uniform(2.0, 5.0),
            "V": rng.uniform(6.75, 11.25),
            "psi": rng.uniform(115.0, 141.0),
            "w": rng.uniform(13.5, 22.5),
            "mstar": rng.uniform(0.36, 0.56),
            "vz": rng.uniform(-10.0, 10.0),
        }
        eta = rng.uniform(0.1, 0.5)
        density_asym = rng.uniform(-0.2, 0.2)
    else:
        params = {
            "theta": rng.uniform(2.0, 3.3),
            "V": rng.uniform(8.0, 16.0),
            "psi": rng.uniform(108.0, 150.0),
            "w": rng.uniform(8.0, 23.0),
            "mstar": rng.uniform(0.40, 0.80),
            "vz": rng.uniform(-12.0, 12.0),
        }
        eta = rng.uniform(0.0, 1.0)
        density_asym = rng.uniform(-0.4, 0.4)
    return params, float(eta), float(density_asym)


def run_random_scan(
    rng: np.random.Generator,
    count: int,
    conservative: bool,
    progress_every: int = 100,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for index in range(count):
        params, eta, density_asym = draw_parameter(rng, conservative=conservative)
        result = coarse_metrics(
            params,
            nk=5,
            eta=eta,
            density_asym=density_asym,
            shell=2,
        )
        result["aug_score"] = augmented_score(result)
        rows.append(result)
        if progress_every and (index + 1) % progress_every == 0:
            print(f"  completed {index + 1}/{count}", file=sys.stderr)
    return rows


def refine_scan(
    rows: Sequence[Mapping[str, Any]],
    refine_count: int | None,
) -> list[dict[str, Any]]:
    candidates = [row for row in rows if is_valid(row)]
    candidates = sorted(candidates, key=lambda item: float(item["aug_score"]), reverse=True)
    if refine_count is not None:
        candidates = candidates[:refine_count]

    output: list[dict[str, Any]] = []
    for index, row in enumerate(candidates):
        result = coarse_metrics(
            parameter_dict(row),
            nk=9,
            eta=float(row["eta"]),
            density_asym=float(row["density_asym"]),
            shell=3,
        )
        result["aug_score"] = augmented_score(result)
        output.append(result)
        if (index + 1) % 10 == 0:
            print(f"  refined {index + 1}/{len(candidates)}", file=sys.stderr)
    return output


def add_geometry(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for index, row in enumerate(rows):
        if not is_valid(row):
            continue
        geometry = fine_geometry(
            parameter_dict(row),
            float(row["eta"]),
            float(row["density_asym"]),
            row["active_pair"],
            shell=2,
            nk=7,
        )
        merged = dict(row)
        merged.update(
            {
                "Gamma_x": float(geometry["Gamma_sum_aM2"][0]),
                "Gamma_y": float(geometry["Gamma_sum_aM2"][1]),
                "Gamma_min": float(np.min(geometry["Gamma_sum_aM2"])),
                "ideal_gamma_min": float(np.min(geometry["ideal_gamma_sum_aM2"])),
                "t0_P2_min": float(np.min(geometry["t0_sum_meV_aM2"])),
                "chern1": float(geometry["chern_approx"][row["active_pair"][0]]),
                "chern2": float(geometry["chern_approx"][row["active_pair"][1]]),
                "chern_remote": float(
                    geometry["chern_approx"][
                        next(i for i in range(3) if i not in row["active_pair"])
                    ]
                ),
            }
        )
        output.append(merged)
        if (index + 1) % 10 == 0:
            print(f"  geometry {index + 1}/{len(rows)}", file=sys.stderr)
    return output


def select_representatives(rows: Sequence[Mapping[str, Any]]) -> list[tuple[str, Mapping[str, Any]]]:
    controlled = [row for row in rows if is_controlled(row)]
    if not controlled:
        raise RuntimeError("scan produced no controlled candidate")
    balanced = max(
        controlled,
        key=lambda row: float(row["Gamma_min"]) / (1.0 + float(row["control_ratio"])),
    )
    max_hopping = max(controlled, key=lambda row: float(row["Gamma_min"]))
    most_controlled = min(controlled, key=lambda row: float(row["control_ratio"]))
    return [
        ("balanced", balanced),
        ("max_hopping", max_hopping),
        ("most_controlled", most_controlled),
    ]


def refine_representatives(
    representatives: Sequence[tuple[str, Mapping[str, Any]]],
) -> list[tuple[str, dict[str, Any], dict[str, Any], dict[str, Any], dict[str, Any]]]:
    output = []
    for name, row in representatives:
        params = parameter_dict(row)
        eta = float(row["eta"])
        density_asym = float(row["density_asym"])
        metric3 = coarse_metrics(params, nk=11, eta=eta, density_asym=density_asym, shell=3)
        metric4 = coarse_metrics(params, nk=9, eta=eta, density_asym=density_asym, shell=4)
        geometry3 = fine_geometry(
            params, eta, density_asym, metric3["active_pair"], shell=3, nk=11
        )
        geometry4 = fine_geometry(
            params, eta, density_asym, metric4["active_pair"], shell=4, nk=7
        )
        output.append((name, metric3, metric4, geometry3, geometry4))
    return output


def write_rows(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    pd.DataFrame([flat_row(row) for row in rows]).to_csv(path, index=False)


def make_basic_plots(
    output: Path,
    exploratory_geometry: Sequence[Mapping[str, Any]],
    scan_baseline: Sequence[Mapping[str, Any]],
    scan_balanced: Sequence[Mapping[str, Any]],
    refined: Sequence[tuple[str, Mapping[str, Any], Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]],
) -> None:
    figures = output / "figures"
    figures.mkdir(exist_ok=True)

    theta_base = np.array([row["theta"] for row in scan_baseline])
    control_base = np.array(
        [row["control_ratio"] if row["remote_gap"] > 0 else np.nan for row in scan_baseline]
    )
    theta_bal = np.array([row["theta"] for row in scan_balanced])
    control_bal = np.array(
        [row["control_ratio"] if row["remote_gap"] > 0 else np.nan for row in scan_balanced]
    )

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(theta_base, control_base, label="standard parameter set")
    ax.plot(theta_bal, control_bal, label="balanced candidate parameters")
    ax.plot(theta_base, np.ones(len(theta_base)), label="control threshold")
    ax.set_xlabel("twist angle theta (degrees)")
    ax.set_ylabel("control ratio")
    ax.set_ylim(0, 1.5)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figures / "twist_angle_control.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    valid = [
        row
        for row in exploratory_geometry
        if row["remote_gap"] > 0
        and row["singular_gap_min"] > 0.02
        and row["capture_top2"] > 0.9
        and row["energy_gap34"] > 1
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.scatter(
        [row["control_ratio"] for row in valid],
        [row["t0_P2_min"] for row in valid],
        label="refined exploratory candidates",
    )
    for name, _, metric4, _, geometry4 in refined:
        ax.scatter(
            [metric4["control_ratio"]],
            [float(np.min(geometry4["t0_sum_meV_aM2"]))],
            label=name.replace("_", " "),
        )
    ax.set_xlabel("control ratio")
    ax.set_ylabel("one-pair hopping proxy for P=2 meV")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(figures / "pareto_control_hopping.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def run_full(args: argparse.Namespace) -> None:
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    # The archived exploratory scan followed an earlier 1,000-point pilot scan.
    # Advancing 8,000 uniforms reproduces that archived RNG position.
    if args.reproduce_archive_rng:
        rng.random(1000 * 8)

    start = time.time()
    print("Exploratory coarse scan", file=sys.stderr)
    exploratory = run_random_scan(
        rng, args.exploratory_points, conservative=False, progress_every=args.progress_every
    )
    write_rows(output / "exploratory_coarse_scan.csv", exploratory)

    print("Exploratory shell-3 refinement", file=sys.stderr)
    exploratory_refined = refine_scan(exploratory, args.refine_count)
    exploratory_refined_valid = [row for row in exploratory_refined if is_valid(row)]
    write_rows(output / "exploratory_refined_scan.csv", exploratory_refined_valid)

    print("Exploratory geometry", file=sys.stderr)
    exploratory_geometry = add_geometry(exploratory_refined_valid)
    write_rows(output / "exploratory_geometry_scan.csv", exploratory_geometry)

    representatives = select_representatives(exploratory_geometry)
    refined = refine_representatives(representatives)

    print("Conservative coarse scan", file=sys.stderr)
    conservative = run_random_scan(
        rng, args.conservative_points, conservative=True, progress_every=args.progress_every
    )
    write_rows(output / "conservative_coarse_scan.csv", conservative)

    print("Conservative shell-3 refinement", file=sys.stderr)
    conservative_refined = refine_scan(conservative, None)
    conservative_refined_valid = [row for row in conservative_refined if is_valid(row)]
    write_rows(output / "conservative_refined_scan.csv", conservative_refined_valid)

    balanced = dict(representatives[0][1])
    balanced_template = parameter_dict(balanced)
    baseline_template = {
        "theta": 3.65,
        "V": 9.0,
        "psi": 128.0,
        "w": 18.0,
        "mstar": 0.45,
        "vz": 0.0,
    }
    thetas = np.linspace(2.0, 5.0, 31)
    print("Twist-angle line scans", file=sys.stderr)
    scan_baseline = angle_scan(
        baseline_template, eta=10.0 / 35.0, density_asym=0.0, thetas=thetas
    )
    scan_balanced = angle_scan(
        balanced_template,
        eta=float(balanced["eta"]),
        density_asym=float(balanced["density_asym"]),
        thetas=thetas,
    )
    write_rows(output / "angle_scan_baseline.csv", scan_baseline)
    write_rows(output / "angle_scan_balanced.csv", scan_balanced)

    payload: dict[str, Any] = {
        "configuration": {
            "seed": args.seed,
            "exploratory_points": args.exploratory_points,
            "conservative_points": args.conservative_points,
            "refine_count": args.refine_count,
            "pair_hop_benchmark_meV": 2.0,
        },
        "counts": {
            "exploratory_valid": sum(is_valid(row) for row in exploratory),
            "exploratory_refined_valid": len(exploratory_refined_valid),
            "exploratory_controlled": sum(is_controlled(row) for row in exploratory_geometry),
            "conservative_valid": sum(is_valid(row) for row in conservative),
            "conservative_refined_valid": len(conservative_refined_valid),
        },
        "refined_candidates": {},
    }
    for name, metric3, metric4, geometry3, geometry4 in refined:
        payload["refined_candidates"][name] = {
            "shell3_metrics": jsonable(metric3),
            "shell4_metrics": jsonable(metric4),
            "shell3_geometry": jsonable(geometry3),
            "shell4_geometry": jsonable(geometry4),
        }
    payload["elapsed_seconds"] = time.time() - start
    (output / "scan_results.json").write_text(
        json.dumps(jsonable(payload), indent=2), encoding="utf-8"
    )
    make_basic_plots(output, exploratory_geometry, scan_baseline, scan_balanced, refined)
    print(f"Completed in {payload['elapsed_seconds']:.1f} s", file=sys.stderr)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode",
        choices=("verify", "full"),
        default="verify",
        help="verify archived candidates or rerun random scans",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("interaction_singular_scan_rerun"),
    )
    parser.add_argument("--seed", type=int, default=20260725)
    parser.add_argument("--exploratory-points", type=int, default=2800)
    parser.add_argument("--conservative-points", type=int, default=1800)
    parser.add_argument("--refine-count", type=int, default=120)
    parser.add_argument("--progress-every", type=int, default=100)
    parser.add_argument(
        "--reproduce-archive-rng",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="advance the RNG past the archived 1,000-point pilot scan",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.mode == "verify":
        verify_archived(args.output)
    else:
        run_full(args)


if __name__ == "__main__":
    main()
