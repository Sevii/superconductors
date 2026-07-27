#!/usr/bin/env python3
"""Independent full-Fock certificate for the exact swap-QGN Meissner paper.

Unlike the compressed hard-core-pair verifier, this script builds canonical
spinful-fermion Fock spaces and the complete second-quantized site-swap
operators.  It checks the fermionic layer of the construction on a path,
triangle, and four-site cycle:

  * W_e and W_e[A] are Hermitian involutions;
  * H_G[0] is positive and its even-sector kernel is the paired Dicke state
    with the canonical fermionic signs, while odd sectors have no zero modes;
  * U_phi H[A] U_phi^dagger = H[A+B^*phi] on the full Fock space;
  * seniority-zero compression gives the charge-two Peierls pair generator;
  * commutator current/stress derivatives agree with finite differences;
  * the full-Fock Kubo kernel reproduces the graph-Hodge static formula and
    the imaginary-frequency graph filter, including exact transverse
    frequency independence.

The calculation is a finite-dimensional certificate, not a substitute for the
analytic proofs in the paper.
"""
from __future__ import annotations

import itertools
import math
import platform
import sys
from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import scipy
import scipy.linalg

TOL = 8.0e-10


@dataclass(frozen=True)
class Graph:
    name: str
    nsites: int
    edges: tuple[tuple[int, int], ...]
    incidence: np.ndarray


def make_graph(name: str, nsites: int, edges: Sequence[tuple[int, int]]) -> Graph:
    b = np.zeros((nsites, len(edges)), dtype=float)
    for e, (tail, head) in enumerate(edges):
        b[tail, e] = -1.0
        b[head, e] = +1.0
    return Graph(name, nsites, tuple(edges), b)


def mode(site: int, spin: int) -> int:
    return 2 * site + spin


def inversion_parity(values: Sequence[int]) -> int:
    inversions = 0
    for i in range(len(values)):
        for j in range(i + 1, len(values)):
            inversions += values[i] > values[j]
    return -1 if inversions % 2 else 1


def site_swap_matrix(nsites: int, x: int, y: int) -> np.ndarray:
    """Second quantization of the one-particle permutation x <-> y."""
    nmodes = 2 * nsites
    dim = 1 << nmodes
    perm = list(range(nmodes))
    for spin in (0, 1):
        perm[mode(x, spin)], perm[mode(y, spin)] = perm[mode(y, spin)], perm[mode(x, spin)]
    out = np.zeros((dim, dim), dtype=complex)
    for state in range(dim):
        occupied = [i for i in range(nmodes) if (state >> i) & 1]
        images = [perm[i] for i in occupied]
        target = sum(1 << i for i in images)
        out[target, state] = inversion_parity(images)
    return out


def number_diagonal(nsites: int, site: int) -> np.ndarray:
    dim = 1 << (2 * nsites)
    return np.array(
        [((state >> mode(site, 0)) & 1) + ((state >> mode(site, 1)) & 1) for state in range(dim)],
        dtype=float,
    )


def onsite_seniority_penalty(nsites: int) -> np.ndarray:
    dim = 1 << (2 * nsites)
    diag = np.zeros(dim, dtype=float)
    for state in range(dim):
        for x in range(nsites):
            up = (state >> mode(x, 0)) & 1
            down = (state >> mode(x, 1)) & 1
            diag[state] += (up - down) ** 2
    return np.diag(diag)


def twisted_swap(w: np.ndarray, n_head: np.ndarray, field: float) -> np.ndarray:
    phase = np.exp(-1.0j * field * n_head)
    return phase[:, None] * w * phase.conj()[None, :]


def full_hamiltonian(graph: Graph, field: np.ndarray, U: float, j: float) -> np.ndarray:
    field = np.asarray(field, dtype=float)
    if field.shape != (len(graph.edges),):
        raise ValueError("wrong field shape")
    dim = 1 << (2 * graph.nsites)
    h = 0.5 * U * onsite_seniority_penalty(graph.nsites)
    ident = np.eye(dim, dtype=complex)
    for e, (tail, head) in enumerate(graph.edges):
        del tail
        w = site_swap_matrix(graph.nsites, graph.edges[e][0], head)
        wa = twisted_swap(w, number_diagonal(graph.nsites, head), field[e])
        pi = 0.5 * (ident - wa)
        h = h + 0.5 * j * (pi @ pi)
    return h


def current_stress(graph: Graph, profile: np.ndarray, j: float) -> tuple[np.ndarray, np.ndarray]:
    """Analytic derivatives of the complete full-Fock Peierls family."""
    dim = 1 << (2 * graph.nsites)
    current = np.zeros((dim, dim), dtype=complex)
    stress = np.zeros((dim, dim), dtype=complex)
    for e, (tail, head) in enumerate(graph.edges):
        w = site_swap_matrix(graph.nsites, tail, head)
        n = np.diag(number_diagonal(graph.nsites, head))
        comm = n @ w - w @ n
        double_comm = n @ comm - comm @ n
        current += 0.25j * j * profile[e] * comm
        stress += 0.25 * j * profile[e] ** 2 * double_comm
    return current, stress


def sector_indices(nsites: int, particles: int) -> np.ndarray:
    return np.array(
        [state for state in range(1 << (2 * nsites)) if state.bit_count() == particles],
        dtype=int,
    )


def restrict(a: np.ndarray, indices: np.ndarray) -> np.ndarray:
    return a[np.ix_(indices, indices)]


def pair_state_index(nsites: int, subset: Iterable[int]) -> int:
    state = 0
    for x in subset:
        state |= 1 << mode(x, 0)
        state |= 1 << mode(x, 1)
    return state


def dicke_full(nsites: int, npairs: int) -> np.ndarray:
    dim = 1 << (2 * nsites)
    out = np.zeros(dim, dtype=complex)
    norm = math.sqrt(math.comb(nsites, npairs))
    for subset in itertools.combinations(range(nsites), npairs):
        out[pair_state_index(nsites, subset)] = 1.0 / norm
    return out


def seniority_zero_isometry(nsites: int, npairs: int) -> tuple[np.ndarray, tuple[int, ...], dict[int, int]]:
    subsets = tuple(sum(1 << x for x in occ) for occ in itertools.combinations(range(nsites), npairs))
    lookup = {bits: i for i, bits in enumerate(subsets)}
    full_indices = sector_indices(nsites, 2 * npairs)
    full_lookup = {state: i for i, state in enumerate(full_indices)}
    p = np.zeros((len(full_indices), len(subsets)), dtype=complex)
    for col, bits in enumerate(subsets):
        sites = [x for x in range(nsites) if (bits >> x) & 1]
        p[full_lookup[pair_state_index(nsites, sites)], col] = 1.0
    return p, subsets, lookup


def pair_expected(graph: Graph, npairs: int, field: np.ndarray, j: float) -> np.ndarray:
    _p, basis, lookup = seniority_zero_isometry(graph.nsites, npairs)
    dim = len(basis)
    h = np.zeros((dim, dim), dtype=complex)
    for e, (x, y) in enumerate(graph.edges):
        for col, bits in enumerate(basis):
            h[col, col] += j / 4.0
            ox, oy = (bits >> x) & 1, (bits >> y) & 1
            target = bits if ox == oy else bits ^ (1 << x) ^ (1 << y)
            row = lookup[target]
            if ox == oy:
                phase = 1.0
            elif ox == 1:
                phase = np.exp(-2.0j * field[e])
            else:
                phase = np.exp(+2.0j * field[e])
            h[row, col] -= (j / 4.0) * phase
    return h


def hodge_projector(graph: Graph) -> np.ndarray:
    b = graph.incidence
    return np.eye(len(graph.edges)) - b.T @ np.linalg.pinv(b @ b.T, rcond=1.0e-13) @ b


def spectral_kubo(
    h: np.ndarray,
    current: np.ndarray,
    stress: np.ndarray,
    omega: np.ndarray,
    zeta: float | None,
) -> float:
    evals, evecs = scipy.linalg.eigh(h)
    direct = float(np.real(np.vdot(omega, stress @ omega)))
    amps = evecs.conj().T @ (current @ omega)
    param = 0.0
    for energy, amp in zip(evals, amps):
        if energy < 1.0e-10:
            continue
        if zeta is None:
            param += 2.0 * abs(amp) ** 2 / energy
        else:
            param += 2.0 * energy * abs(amp) ** 2 / (energy * energy + zeta * zeta)
    return direct - param


def exact_graph_filter(graph: Graph, npairs: int, j: float, profile: np.ndarray, zeta: float) -> float:
    lap1 = graph.incidence.T @ graph.incidence
    vals, vecs = scipy.linalg.eigh(lap1)
    filt = (vecs * (zeta * zeta / ((j * vals / 4.0) ** 2 + zeta * zeta))) @ vecs.T
    gamma = 2.0 * npairs * (graph.nsites - npairs) / (graph.nsites * (graph.nsites - 1))
    return j * gamma * float(profile @ filt @ profile)


def max_abs(a: np.ndarray) -> float:
    return float(np.max(np.abs(a))) if a.size else 0.0


def check(name: str, error: float, tol: float = TOL) -> None:
    status = "PASS" if error <= tol else "FAIL"
    print(f"{name:<78s} {status}  error={error:.3e}  tol={tol:.1e}")
    if error > tol:
        raise AssertionError(name)


def report(name: str, value: float) -> None:
    print(f"{name:<78s} VALUE {value:.12g}")


def graph_algebra_checks(graph: Graph, rng: np.random.Generator, U: float, j: float) -> None:
    dim = 1 << (2 * graph.nsites)
    ident = np.eye(dim, dtype=complex)
    field = rng.normal(scale=0.21, size=len(graph.edges))
    max_w = 0.0
    max_wa = 0.0
    for e, (tail, head) in enumerate(graph.edges):
        w = site_swap_matrix(graph.nsites, tail, head)
        wa = twisted_swap(w, number_diagonal(graph.nsites, head), field[e])
        max_w = max(max_w, max_abs(w - w.conj().T), max_abs(w @ w - ident))
        max_wa = max(max_wa, max_abs(wa - wa.conj().T), max_abs(wa @ wa - ident))
    check(f"{graph.name}: complete swaps are Hermitian involutions", max_w)
    check(f"{graph.name}: twisted swaps are Hermitian involutions", max_wa)

    h0 = full_hamiltonian(graph, np.zeros(len(graph.edges)), U, j)
    check(f"{graph.name}: full Hamiltonian Hermiticity", max_abs(h0 - h0.conj().T))
    min_global = float(np.min(scipy.linalg.eigvalsh(h0)))
    check(f"{graph.name}: full Hamiltonian positivity", max(0.0, -min_global))

    for particles in range(2 * graph.nsites + 1):
        idx = sector_indices(graph.nsites, particles)
        hs = restrict(h0, idx)
        evals, evecs = scipy.linalg.eigh(hs)
        if particles % 2:
            check(f"{graph.name}: odd sector N={particles} has no zero mode", max(0.0, 1.0e-8 - float(evals[0])), 1.1e-8)
        else:
            npairs = particles // 2
            expected_full = dicke_full(graph.nsites, npairs)
            expected = expected_full[idx]
            check(f"{graph.name}: Dicke state N={particles} is exactly null", float(np.linalg.norm(hs @ expected)))
            if len(evals) > 1:
                check(f"{graph.name}: even-sector kernel N={particles} is one-dimensional", abs(float(evals[0])) + max(0.0, 1.0e-8 - float(evals[1])), 1.1e-8)
            overlap = abs(np.vdot(expected, evecs[:, 0]))
            check(f"{graph.name}: fermionic Dicke signs N={particles}", abs(1.0 - overlap), 2.0e-10)

    phi = rng.normal(scale=0.18, size=graph.nsites)
    h_a = full_hamiltonian(graph, field, U, j)
    h_ag = full_hamiltonian(graph, field + graph.incidence.T @ phi, U, j)
    phases = np.empty(dim, dtype=complex)
    for state in range(dim):
        charge_phase = 0.0
        for x in range(graph.nsites):
            charge_phase += phi[x] * (((state >> mode(x, 0)) & 1) + ((state >> mode(x, 1)) & 1))
        phases[state] = np.exp(-1.0j * charge_phase)
    u = np.diag(phases)
    check(f"{graph.name}: full-Fock gauge covariance", max_abs(u @ h_a @ u.conj().T - h_ag), 2.0e-10)

    for npairs in range(graph.nsites + 1):
        idx = sector_indices(graph.nsites, 2 * npairs)
        p_sz, _basis, _lookup = seniority_zero_isometry(graph.nsites, npairs)
        compressed = p_sz.conj().T @ restrict(h_a, idx) @ p_sz
        expected = pair_expected(graph, npairs, field, j)
        check(f"{graph.name}: charge-two paired compression n={npairs}", max_abs(compressed - expected), 2.0e-10)


def response_checks(graph: Graph, rng: np.random.Generator, U: float, j: float) -> None:
    npairs = 2
    particles = 2 * npairs
    idx = sector_indices(graph.nsites, particles)
    h0_full = full_hamiltonian(graph, np.zeros(len(graph.edges)), U, j)
    h0 = restrict(h0_full, idx)
    omega = dicke_full(graph.nsites, npairs)[idx]
    profile = rng.normal(size=len(graph.edges))
    current_full, stress_full = current_stress(graph, profile, j)
    current = restrict(current_full, idx)
    stress = restrict(stress_full, idx)

    eps = 5.0e-4
    h_plus = restrict(full_hamiltonian(graph, eps * profile, U, j), idx)
    h_minus = restrict(full_hamiltonian(graph, -eps * profile, U, j), idx)
    h_plus2 = restrict(full_hamiltonian(graph, 2.0 * eps * profile, U, j), idx)
    h_minus2 = restrict(full_hamiltonian(graph, -2.0 * eps * profile, U, j), idx)
    j_fd = (-h_plus2 + 8.0 * h_plus - 8.0 * h_minus + h_minus2) / (12.0 * eps)
    t_fd = (-h_plus2 + 16.0 * h_plus - 30.0 * h0 + 16.0 * h_minus - h_minus2) / (12.0 * eps * eps)
    check("C4: commutator current equals five-point finite difference", max_abs(current - j_fd), 2.0e-9)
    check("C4: commutator stress equals five-point finite difference", max_abs(stress - t_fd), 2.0e-8)

    gamma = 2.0 * npairs * (graph.nsites - npairs) / (graph.nsites * (graph.nsites - 1))
    static = spectral_kubo(h0, current, stress, omega, None)
    expected_static = j * gamma * float(profile @ hodge_projector(graph) @ profile)
    check("C4: full-Fock static Kubo equals Hodge formula", abs(static - expected_static), 3.0e-10)

    zeta = 0.61
    dynamic = spectral_kubo(h0, current, stress, omega, zeta)
    expected_dynamic = exact_graph_filter(graph, npairs, j, profile, zeta)
    check("C4: full-Fock imaginary-frequency graph filter", abs(dynamic - expected_dynamic), 3.0e-10)

    transverse = np.ones(len(graph.edges), dtype=float)
    check("C4: oriented cycle profile is divergence free", float(np.linalg.norm(graph.incidence @ transverse)))
    jt_full, tt_full = current_stress(graph, transverse, j)
    jt, tt = restrict(jt_full, idx), restrict(tt_full, idx)
    check("C4: transverse current annihilates full-Fock Dicke state", float(np.linalg.norm(jt @ omega)))
    exact = j * gamma * float(transverse @ transverse)
    check("C4: transverse static response is diamagnetic", abs(spectral_kubo(h0, jt, tt, omega, None) - exact))
    check("C4: transverse response is frequency independent", abs(spectral_kubo(h0, jt, tt, omega, zeta) - exact))
    report("C4 transverse full-Fock response", exact)


def main() -> int:
    np.set_printoptions(precision=9, suppress=True)
    rng = np.random.default_rng(2026072701)
    U, j = 2.3, 1.4
    graphs = (
        make_graph("P3", 3, ((0, 1), (1, 2))),
        make_graph("C3", 3, ((0, 1), (1, 2), (2, 0))),
        make_graph("C4", 4, ((0, 1), (1, 2), (2, 3), (3, 0))),
    )

    print("INDEPENDENT FULL-FOCK SWAP-QGN CHECK")
    print("python", platform.python_version())
    print("numpy", np.__version__)
    print("scipy", scipy.__version__)
    print("-" * 116)
    for graph in graphs:
        graph_algebra_checks(graph, rng, U, j)
    response_checks(graphs[-1], rng, U, j)
    print("-" * 116)
    print("OVERALL: PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"OVERALL: FAIL ({type(exc).__name__}: {exc})", file=sys.stderr)
        raise
