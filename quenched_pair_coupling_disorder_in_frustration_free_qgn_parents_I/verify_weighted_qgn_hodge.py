#!/usr/bin/env python3
"""Numerical certificates for weighted swap-parent QGN robustness.

The script checks:
  1. exact finite-graph weighted Hodge curvature against direct many-body least squares;
  2. exact weighted dynamical filter against the full target-space response;
  3. gauge invariance and deterministic ellipticity bounds;
  4. the exact effective-resistance formula on a disordered ring;
  5. finite-volume random-torus effective conductivity statistics;
  6. soft-source leakage for macroscopically modulated bounded conductances.

Only the hard-core-pair (seniority-zero) sector is needed: complete site-swap
rows act there as transpositions of occupations, and the Dicke/AGP state is the
uniform vector on n-subsets.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import platform
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import scipy
import scipy.linalg as la


TOL = 5.0e-9


@dataclass(frozen=True)
class Graph:
    n_vertices: int
    edges: tuple[tuple[int, int], ...]  # oriented tail -> head
    displacements: np.ndarray | None = None  # E x d covering-lattice displacement

    @property
    def incidence(self) -> np.ndarray:
        """Return B: edge fields -> site divergences, shape N x E."""
        b = np.zeros((self.n_vertices, len(self.edges)), dtype=float)
        for e, (tail, head) in enumerate(self.edges):
            b[tail, e] = -1.0
            b[head, e] = +1.0
        return b


def ring_graph(n: int) -> Graph:
    edges = tuple((x, (x + 1) % n) for x in range(n))
    r = np.ones((n, 1), dtype=float)
    return Graph(n, edges, r)


def torus_graph(shape: Sequence[int]) -> Graph:
    shape = tuple(int(v) for v in shape)
    d = len(shape)
    n = math.prod(shape)

    def idx(coord: Sequence[int]) -> int:
        out = 0
        stride = 1
        for c, ell in zip(reversed(coord), reversed(shape)):
            out += c * stride
            stride *= ell
        return out

    edges: list[tuple[int, int]] = []
    displacements: list[np.ndarray] = []
    for coord in itertools.product(*(range(ell) for ell in shape)):
        x = idx(coord)
        for mu, ell in enumerate(shape):
            nbr = list(coord)
            nbr[mu] = (nbr[mu] + 1) % ell
            edges.append((x, idx(nbr)))
            disp = np.zeros(d)
            disp[mu] = 1.0
            displacements.append(disp)
    return Graph(n, tuple(edges), np.asarray(displacements, dtype=float))


def path_graph(n: int) -> Graph:
    return Graph(n, tuple((x, x + 1) for x in range(n - 1)), None)


def pinv_psd(a: np.ndarray, rtol: float = 1e-12) -> np.ndarray:
    vals, vecs = la.eigh((a + a.T.conj()) / 2)
    cutoff = rtol * max(1.0, float(np.max(np.abs(vals))))
    inv = np.zeros_like(vals)
    inv[vals > cutoff] = 1.0 / vals[vals > cutoff]
    return (vecs * inv) @ vecs.T.conj()


def weighted_hodge_energy(
    graph: Graph, conductances: np.ndarray, a: np.ndarray
) -> tuple[float, np.ndarray, np.ndarray]:
    """Return min_phi sum_e j_e(a_e + grad phi_e)^2 and minimizer."""
    b = graph.incidence
    j = np.asarray(conductances, dtype=float)
    a = np.asarray(a, dtype=float)
    lap = (b * j[None, :]) @ b.T
    g = b @ (j * a)
    phi = -pinv_psd(lap) @ g
    residual = a + b.T @ phi
    energy = float(np.dot(j * residual, residual))
    return energy, phi, residual


def hodge_matrix(graph: Graph, conductances: np.ndarray) -> np.ndarray:
    b = graph.incidence
    j = np.asarray(conductances, dtype=float)
    J = np.diag(j)
    lap = b @ J @ b.T
    return J - J @ b.T @ pinv_psd(lap) @ b @ J


def gamma_factor(n_vertices: int, n_pairs: int) -> float:
    if n_vertices <= 1:
        return 0.0
    return 2.0 * n_pairs * (n_vertices - n_pairs) / (
        n_vertices * (n_vertices - 1)
    )


def subset_basis(n_vertices: int, n_pairs: int) -> tuple[tuple[int, ...], ...]:
    return tuple(itertools.combinations(range(n_vertices), n_pairs))


def swap_matrix(
    basis: Sequence[tuple[int, ...]], tail: int, head: int
) -> np.ndarray:
    index = {state: k for k, state in enumerate(basis)}
    m = len(basis)
    w = np.zeros((m, m), dtype=float)
    for col, state in enumerate(basis):
        occ = set(state)
        tail_in = tail in occ
        head_in = head in occ
        if tail_in != head_in:
            if tail_in:
                occ.remove(tail)
                occ.add(head)
            else:
                occ.remove(head)
                occ.add(tail)
        row = index[tuple(sorted(occ))]
        w[row, col] = 1.0
    return w


def occupation_vector(
    basis: Sequence[tuple[int, ...]], vertex: int
) -> np.ndarray:
    return np.asarray([1.0 if vertex in s else 0.0 for s in basis])


def many_body_rows_and_source(
    graph: Graph, conductances: np.ndarray, a: np.ndarray, n_pairs: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[tuple[int, ...], ...]]:
    """Build stacked D and S|Omega> in the hard-core-pair sector.

    The physical electronic Peierls phase is used. A pair at an edge head has
    charge two, so d/dt Pi_e(t a)|Omega> = i a_e(n_head-n_tail)|Omega>.
    """
    basis = subset_basis(graph.n_vertices, n_pairs)
    dim = len(basis)
    omega = np.ones(dim, dtype=complex) / math.sqrt(dim)
    blocks: list[np.ndarray] = []
    source_blocks: list[np.ndarray] = []
    for (tail, head), je, ae in zip(graph.edges, conductances, a):
        w = swap_matrix(basis, tail, head)
        pi = (np.eye(dim) - w) / 2.0
        blocks.append(math.sqrt(float(je)) * pi)
        n_head = occupation_vector(basis, head)
        n_tail = occupation_vector(basis, tail)
        z = (n_head - n_tail) * omega
        source_blocks.append(1j * math.sqrt(float(je)) * float(ae) * z)
    D = np.vstack(blocks).astype(complex)
    source = np.concatenate(source_blocks)
    return D, source, omega, basis


def direct_static_curvature(D: np.ndarray, source: np.ndarray) -> float:
    # min_x ||source + D x||^2; sign of x is immaterial.
    x, *_ = la.lstsq(D, -source, cond=1e-12, lapack_driver="gelsd")
    residual = source + D @ x
    return float(np.vdot(residual, residual).real)


def direct_dynamic_kernel(
    D: np.ndarray, source: np.ndarray, zeta: float
) -> float:
    target_lap = 0.5 * D @ D.conj().T
    vals, vecs = la.eigh((target_lap + target_lap.conj().T) / 2)
    coeff = vecs.conj().T @ source
    multiplier = zeta**2 / (vals**2 + zeta**2)
    return float(np.sum(multiplier * np.abs(coeff) ** 2).real)


def graph_dynamic_kernel(
    graph: Graph,
    conductances: np.ndarray,
    a: np.ndarray,
    n_pairs: int,
    zeta: float,
) -> float:
    b = graph.incidence
    sqrt_j = np.sqrt(conductances)
    m = 0.25 * (sqrt_j[:, None] * (b.T @ b) * sqrt_j[None, :])
    vals, vecs = la.eigh((m + m.T) / 2)
    s = sqrt_j * a
    coeff = vecs.T @ s
    multiplier = zeta**2 / (vals**2 + zeta**2)
    return gamma_factor(graph.n_vertices, n_pairs) * float(
        np.sum(multiplier * coeff**2)
    )


def effective_tensor(graph: Graph, conductances: np.ndarray) -> np.ndarray:
    if graph.displacements is None:
        raise ValueError("Graph has no covering-lattice displacement matrix")
    p = hodge_matrix(graph, conductances)
    return graph.displacements.T @ p @ graph.displacements


def soft_source_mass(
    graph: Graph,
    conductances: np.ndarray,
    a: np.ndarray,
    cutoff: float,
) -> tuple[float, float, float]:
    """Normalized source mass and H^-1 mass below a target-laplacian cutoff."""
    b = graph.incidence
    sqrt_j = np.sqrt(conductances)
    m = 0.25 * (sqrt_j[:, None] * (b.T @ b) * sqrt_j[None, :])
    vals, vecs = la.eigh((m + m.T) / 2)
    s = sqrt_j * a
    coeff2 = np.abs(vecs.T @ s) ** 2
    positive = vals > 1e-12
    low = positive & (vals <= cutoff)
    mass = float(np.sum(coeff2[low]) / graph.n_vertices)
    hminus1 = float(np.sum(coeff2[low] / vals[low]) / graph.n_vertices)
    total_positive = float(np.sum(coeff2[positive]) / graph.n_vertices)
    return mass, hminus1, total_positive


def assert_close(name: str, x: float, y: float, tol: float = TOL) -> None:
    err = abs(x - y)
    scale = max(1.0, abs(x), abs(y))
    if err > tol * scale:
        raise AssertionError(f"{name}: {x} != {y}; err={err}")
    print(f"PASS {name}: value={x:.12g}, abs_err={err:.3e}")


def run_exact_checks(rng: np.random.Generator) -> dict[str, float]:
    results: dict[str, float] = {}
    cases = [
        (ring_graph(7), 3),
        (torus_graph((2, 3)), 2),
        (path_graph(6), 2),
    ]
    for case_id, (graph, n_pairs) in enumerate(cases):
        e = len(graph.edges)
        j = rng.uniform(0.35, 2.4, size=e)
        a = rng.normal(size=e)
        D, source, _, _ = many_body_rows_and_source(graph, j, a, n_pairs)
        direct = direct_static_curvature(D, source)
        network, _, _ = weighted_hodge_energy(graph, j, a)
        predicted = gamma_factor(graph.n_vertices, n_pairs) * network
        assert_close(f"static_many_body_vs_hodge_case_{case_id}", direct, predicted)
        results[f"static_case_{case_id}"] = direct

        for zeta in (0.03, 0.21, 1.1):
            full = direct_dynamic_kernel(D, source, zeta)
            graph_value = graph_dynamic_kernel(graph, j, a, n_pairs, zeta)
            assert_close(
                f"dynamic_full_target_vs_graph_case_{case_id}_zeta_{zeta}",
                full,
                graph_value,
                tol=2e-8,
            )
        results[f"dynamic_case_{case_id}"] = full

        # Gauge invariance of the static cell problem.
        phi = rng.normal(size=graph.n_vertices)
        shifted, _, _ = weighted_hodge_energy(graph, j, a + graph.incidence.T @ phi)
        assert_close(f"gauge_invariance_case_{case_id}", shifted, network)

        # Deterministic ellipticity comparison with unit weights.
        clean, _, _ = weighted_hodge_energy(graph, np.ones(e), a)
        lower = float(np.min(j) * clean)
        upper = float(np.max(j) * clean)
        if not (direct / gamma_factor(graph.n_vertices, n_pairs) >= lower - 2e-8):
            raise AssertionError("ellipticity lower bound failed")
        if not (direct / gamma_factor(graph.n_vertices, n_pairs) <= upper + 2e-8):
            raise AssertionError("ellipticity upper bound failed")
        print(
            f"PASS ellipticity_case_{case_id}: "
            f"{lower:.8g} <= {network:.8g} <= {upper:.8g}"
        )

    # Exact ring resistance formula.
    graph = ring_graph(31)
    j = rng.uniform(0.4, 2.1, size=31)
    A = 0.37
    energy, _, _ = weighted_hodge_energy(graph, j, A * np.ones(31))
    resistance_formula = (31 * A) ** 2 / np.sum(1.0 / j)
    assert_close("ring_effective_resistance", energy, resistance_formula)
    results["ring_energy"] = energy
    return results


def run_random_torus_scan(
    base_seed: int, out_csv: Path, sizes: Iterable[int], samples: int
) -> list[dict[str, float]]:
    """Run a reproducible per-size Monte Carlo scan.

    Each size receives an independent SeedSequence derived from the recorded
    base seed.  This prevents changes in the exact-check draw count or loop
    ordering from changing the published torus rows.
    """
    rows: list[dict[str, float]] = []
    for L in sizes:
        rng = np.random.default_rng(
            np.random.SeedSequence([int(base_seed), 200, int(L), int(samples)])
        )
        graph = torus_graph((L, L))
        vals_x: list[float] = []
        vals_y: list[float] = []
        vals_xy: list[float] = []
        for _ in range(samples):
            # Uniformly elliptic i.i.d. conductances.
            j = rng.uniform(0.5, 2.0, size=len(graph.edges))
            k = effective_tensor(graph, j) / graph.n_vertices
            eig = la.eigvalsh((k + k.T) / 2)
            if eig[0] < 0.5 - 2e-8 or eig[-1] > 2.0 + 2e-8:
                raise AssertionError("effective tensor violates ellipticity bounds")
            vals_x.append(float(k[0, 0]))
            vals_y.append(float(k[1, 1]))
            vals_xy.append(float(k[0, 1]))
        row = {
            "L": float(L),
            "N": float(graph.n_vertices),
            "samples": float(samples),
            "mean_Kxx": float(np.mean(vals_x)),
            "std_Kxx": float(np.std(vals_x, ddof=1)),
            "mean_Kyy": float(np.mean(vals_y)),
            "std_Kyy": float(np.std(vals_y, ddof=1)),
            "mean_Kxy": float(np.mean(vals_xy)),
            "std_Kxy": float(np.std(vals_xy, ddof=1)),
            "N_var_Kxx": float(graph.n_vertices * np.var(vals_x, ddof=1)),
        }
        rows.append(row)
        print(
            "PASS random_torus "
            f"L={L}: mean Kxx={row['mean_Kxx']:.8f}, "
            f"std={row['std_Kxx']:.4g}, N*var={row['N_var_Kxx']:.4g}"
        )
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return rows


def run_macroscopic_modulation_scan(out_csv: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    eps = 0.35
    for N in (32, 64, 128, 256, 512):
        graph = ring_graph(N)
        x = np.arange(N)
        j = 1.0 + eps * np.cos(2.0 * np.pi * x / N)
        a = np.ones(N)
        cutoff = 8.0 * math.pi**2 / (N**2)
        mass, hminus1, total_positive = soft_source_mass(
            graph, j, a, cutoff=cutoff
        )
        static, _, _ = weighted_hodge_energy(graph, j, a)
        dynamic_zeta = graph_dynamic_kernel(
            graph, j, a, n_pairs=max(1, N // 3), zeta=1.0 / N
        )
        static_curvature = gamma_factor(N, max(1, N // 3)) * static
        row = {
            "N": float(N),
            "epsilon": eps,
            "cutoff": cutoff,
            "soft_source_mass_per_site": mass,
            "soft_Hminus1_mass_per_site": hminus1,
            "total_positive_source_mass_per_site": total_positive,
            "static_network_energy_per_site": static / N,
            "dynamic_minus_static_per_site_at_zeta_1_over_N": (
                dynamic_zeta - static_curvature
            )
            / N,
        }
        rows.append(row)
        print(
            "PASS macroscopic_modulation "
            f"N={N}: soft_mass/N={mass:.8g}, "
            f"positive_mass/N={total_positive:.8g}, "
            f"(K-C)/N={row['dynamic_minus_static_per_site_at_zeta_1_over_N']:.8g}"
        )
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # The positive soft source mass should not vanish for this L-dependent,
    # macroscopically modulated deterministic sequence.
    if rows[-1]["soft_source_mass_per_site"] < 1e-3:
        raise AssertionError("macroscopic modulation did not produce soft leakage")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=Path("."))
    parser.add_argument("--seed", type=int, default=20260727)
    parser.add_argument("--samples", type=int, default=24)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    exact_rng = np.random.default_rng(
        np.random.SeedSequence([int(args.seed), 100])
    )

    exact = run_exact_checks(exact_rng)
    torus = run_random_torus_scan(
        args.seed,
        args.outdir / "random_torus_effective_tensor.csv",
        sizes=(4, 6, 8, 12),
        samples=args.samples,
    )
    modulation = run_macroscopic_modulation_scan(
        args.outdir / "macroscopic_modulation_soft_source.csv"
    )

    summary = {
        "seed": args.seed,
        "rng": {
            "generator": "numpy.random.default_rng(PCG64)",
            "exact_seed_sequence": [args.seed, 100],
            "torus_seed_sequence_template": [args.seed, 200, "L", args.samples],
        },
        "runtime": {
            "python": sys.version.split()[0],
            "python_implementation": platform.python_implementation(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
        },
        "samples_per_size": args.samples,
        "tolerance": TOL,
        "exact_values": exact,
        "random_torus": torus,
        "macroscopic_modulation": modulation,
        "overall": "PASS",
    }
    (args.outdir / "weighted_qgn_hodge_certificate.json").write_text(
        json.dumps(summary, indent=2)
    )
    print("OVERALL: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
