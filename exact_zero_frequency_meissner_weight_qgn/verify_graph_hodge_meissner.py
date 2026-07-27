#!/usr/bin/env python3
"""Finite-volume certificate for the exact graph-Hodge Meissner theorem.

The verifier works in a fixed-n seniority-zero hard-core-pair sector and checks:
  * the unique Dicke/AGP zero state and exact pair ODLRO;
  * H|f> = (j/4)|Delta f> and the site-function inner product;
  * complete bondwise gauge covariance;
  * the current and stress source identities;
  * the static Kubo/Hodge formula for arbitrary link profiles;
  * the full imaginary-frequency edge-Laplacian filter;
  * pure-gauge cancellation and transverse frequency independence;
  * the L^{-2} fixed-number gap, including a multipair rectangular torus;
  * the selected two-block scalar floor;
  * the operator-valued backbone floor with nontrivial residual rows.
"""
from __future__ import annotations

import itertools
import math
import platform
import sys
from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np
import scipy
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg

TOL = 5.0e-11


@dataclass(frozen=True)
class Torus:
    L: int
    d: int
    sites: Tuple[Tuple[int, ...], ...]
    index: Dict[Tuple[int, ...], int]
    edges: Tuple[Tuple[int, int, int], ...]  # (tail, head, direction)
    incidence: np.ndarray  # site x edge, - at tail and + at head


def make_torus(L: int, d: int) -> Torus:
    if L < 3:
        raise ValueError("Use L >= 3 so positive-direction bonds are unique.")
    sites = tuple(itertools.product(range(L), repeat=d))
    index = {x: i for i, x in enumerate(sites)}
    edges: List[Tuple[int, int, int]] = []
    for x in sites:
        tail = index[x]
        for mu in range(d):
            y = list(x)
            y[mu] = (y[mu] + 1) % L
            head = index[tuple(y)]
            edges.append((tail, head, mu))
    b = np.zeros((len(sites), len(edges)), dtype=float)
    for e, (tail, head, _mu) in enumerate(edges):
        b[tail, e] = -1.0
        b[head, e] = +1.0
    return Torus(L, d, sites, index, tuple(edges), b)




def make_rect_torus(shape: Tuple[int, ...]) -> Torus:
    """Periodic rectangular torus; all side lengths must be at least three."""
    if any(length < 3 for length in shape):
        raise ValueError("Use side lengths >= 3 so positive-direction bonds are unique.")
    sites = tuple(itertools.product(*[range(length) for length in shape]))
    index = {x: i for i, x in enumerate(sites)}
    edges: List[Tuple[int, int, int]] = []
    for x in sites:
        tail = index[x]
        for mu, length in enumerate(shape):
            y = list(x)
            y[mu] = (y[mu] + 1) % length
            edges.append((tail, index[tuple(y)], mu))
    b = np.zeros((len(sites), len(edges)), dtype=float)
    for e, (tail, head, _mu) in enumerate(edges):
        b[tail, e] = -1.0
        b[head, e] = +1.0
    return Torus(max(shape), len(shape), sites, index, tuple(edges), b)


def fixed_n_basis(M: int, n: int) -> Tuple[Tuple[int, ...], Dict[int, int]]:
    bits: List[int] = []
    for occ in itertools.combinations(range(M), n):
        state = 0
        for x in occ:
            state |= 1 << x
        bits.append(state)
    return tuple(bits), {state: i for i, state in enumerate(bits)}


def occupied(state: int, x: int) -> int:
    return (state >> x) & 1


def swapped_bit(state: int, x: int, y: int) -> int:
    if occupied(state, x) == occupied(state, y):
        return state
    return state ^ ((1 << x) | (1 << y))


def pair_hamiltonian(
    torus: Torus,
    basis: Sequence[int],
    lookup: Dict[int, int],
    j: float,
    link_field: np.ndarray | None = None,
) -> np.ndarray:
    """H[A]=(j/4) sum_e (1-W_e[A]) in the paired sector."""
    ecount = len(torus.edges)
    field = np.zeros(ecount) if link_field is None else np.asarray(link_field, dtype=float)
    if field.shape != (ecount,):
        raise ValueError("wrong link-field shape")
    dim = len(basis)
    h = np.zeros((dim, dim), dtype=complex)
    for e, (x, y, _mu) in enumerate(torus.edges):
        for col, state in enumerate(basis):
            h[col, col] += j / 4.0
            ox, oy = occupied(state, x), occupied(state, y)
            row = lookup[swapped_bit(state, x, y)]
            if ox == oy:
                phase = 1.0
            elif ox == 1 and oy == 0:
                phase = np.exp(-2.0j * field[e])
            else:
                phase = np.exp(+2.0j * field[e])
            h[row, col] -= (j / 4.0) * phase
    return h




def pair_hamiltonian_sparse(
    torus: Torus,
    basis: Sequence[int],
    lookup: Dict[int, int],
    j: float,
) -> scipy.sparse.csr_matrix:
    """Sparse zero-field paired Hamiltonian for multipair gap checks."""
    dim = len(basis)
    rows: List[int] = []
    cols: List[int] = []
    data: List[complex] = []
    diagonal = np.zeros(dim, dtype=float)
    for x, y, _mu in torus.edges:
        for col, state in enumerate(basis):
            diagonal[col] += j / 4.0
            row = lookup[swapped_bit(state, x, y)]
            rows.append(row)
            cols.append(col)
            data.append(-j / 4.0)
    rows.extend(range(dim))
    cols.extend(range(dim))
    data.extend(diagonal.astype(complex))
    return scipy.sparse.coo_matrix((data, (rows, cols)), shape=(dim, dim)).tocsr()


def matrix_function_psd(a: np.ndarray, func) -> np.ndarray:
    vals, vecs = scipy.linalg.eigh((a + a.conj().T) / 2.0)
    return (vecs * func(vals)) @ vecs.conj().T


def verify_backbone_operator_floor(rng: np.random.Generator) -> Tuple[float, float, float]:
    """Nontrivial direct-sum target test of Theorem 10.1.

    The complete target Laplacian has nonzero swap/residual off-diagonal blocks,
    so the test exercises the kernel-eigenvector step rather than a block-diagonal
    special case.
    """
    ground_dim, excited_dim = 2, 3
    state_dim = ground_dim + excited_dim
    m_sw, m_r = 4, 3

    d_sw = np.zeros((m_sw, state_dim), dtype=complex)
    d_sw[:2, ground_dim:] = rng.normal(size=(2, excited_dim)) + 1j * rng.normal(size=(2, excited_dim))
    d_r = np.zeros((m_r, state_dim), dtype=complex)
    d_r[:, ground_dim:] = rng.normal(size=(m_r, excited_dim)) + 1j * rng.normal(size=(m_r, excited_dim))

    s_sw = np.zeros((m_sw, ground_dim), dtype=complex)
    s_sw[2:, :] = np.array([[1.2, 0.1j], [-0.3j, 0.8]], dtype=complex)
    s_r = rng.normal(size=(m_r, ground_dim)) + 1j * rng.normal(size=(m_r, ground_dim))

    d = np.vstack([d_sw, d_r])
    source = np.vstack([s_sw, s_r])
    lap = 0.5 * d @ d.conj().T
    floor = s_sw.conj().T @ s_sw
    mixing_norm = float(np.linalg.norm(d_sw @ d_r.conj().T))
    kernel_residual = float(np.linalg.norm(d.conj().T @ np.vstack([s_sw, np.zeros_like(s_r)])))

    worst_dynamic = math.inf
    for zeta in (0.03, 0.2, 1.1, 7.0):
        filt = matrix_function_psd(lap, lambda x: zeta * zeta / (x * x + zeta * zeta))
        full = source.conj().T @ filt @ source
        worst_dynamic = min(worst_dynamic, float(np.min(np.linalg.eigvalsh((full - floor + (full - floor).conj().T) / 2.0))))

    pker = matrix_function_psd(lap, lambda x: (x < 1.0e-10).astype(float))
    static = source.conj().T @ pker @ source
    worst_static = float(np.min(np.linalg.eigvalsh((static - floor + (static - floor).conj().T) / 2.0)))
    return mixing_norm, kernel_residual, min(worst_dynamic, worst_static)


def current_and_stress(
    torus: Torus,
    basis: Sequence[int],
    lookup: Dict[int, int],
    j: float,
    profile: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    """First and second derivatives of H[t a] at t=0."""
    a = np.asarray(profile, dtype=float)
    if a.shape != (len(torus.edges),):
        raise ValueError("wrong profile shape")
    dim = len(basis)
    current = np.zeros((dim, dim), dtype=complex)
    stress = np.zeros((dim, dim), dtype=complex)
    for e, (x, y, _mu) in enumerate(torus.edges):
        for col, state in enumerate(basis):
            ox, oy = occupied(state, x), occupied(state, y)
            if ox == oy:
                continue
            row = lookup[swapped_bit(state, x, y)]
            if ox == 1 and oy == 0:
                current[row, col] += +0.5j * j * a[e]
            else:
                current[row, col] += -0.5j * j * a[e]
            stress[row, col] += j * a[e] ** 2
    return current, stress


def omega_state(dim: int) -> np.ndarray:
    return np.ones(dim, dtype=complex) / math.sqrt(dim)


def alpha(M: int, n: int) -> float:
    return n * (M - n) / (M * (M - 1))


def gamma(M: int, n: int) -> float:
    return 2.0 * alpha(M, n)


def expectation(v: np.ndarray, a: np.ndarray, w: np.ndarray | None = None) -> complex:
    if w is None:
        w = v
    return np.vdot(v, a @ w)


def spectral_kubo(
    h: np.ndarray,
    current: np.ndarray,
    stress: np.ndarray,
    omega: np.ndarray,
    zeta: float | None,
) -> float:
    evals, evecs = scipy.linalg.eigh(h)
    direct = float(np.real(expectation(omega, stress)))
    amps = evecs.conj().T @ (current @ omega)
    param = 0.0
    for energy, amp in zip(evals, amps):
        if energy < 1.0e-10:
            continue
        weight = abs(amp) ** 2
        if zeta is None:
            param += 2.0 * weight / energy
        else:
            param += 2.0 * energy * weight / (energy * energy + zeta * zeta)
    return direct - param


def hodge_projector(torus: Torus) -> np.ndarray:
    b = torus.incidence
    lap0 = b @ b.T
    return np.eye(b.shape[1]) - b.T @ np.linalg.pinv(lap0, rcond=1.0e-13) @ b


def exact_imaginary_kernel(
    torus: Torus, j: float, M: int, n: int, a: np.ndarray, zeta: float
) -> float:
    b = torus.incidence
    lap1 = b.T @ b
    vals, vecs = scipy.linalg.eigh(lap1)
    multipliers = zeta * zeta / (zeta * zeta + (j * vals / 4.0) ** 2)
    filt = (vecs * multipliers) @ vecs.T
    return j * gamma(M, n) * float(a @ filt @ a)


def site_vector(basis: Sequence[int], f: np.ndarray) -> np.ndarray:
    out = np.zeros(len(basis), dtype=complex)
    norm = math.sqrt(len(basis))
    for i, state in enumerate(basis):
        out[i] = sum(f[x] for x in range(len(f)) if occupied(state, x)) / norm
    return out


def pair_density_matrix(
    basis: Sequence[int], lookup: Dict[int, int], omega: np.ndarray, M: int
) -> np.ndarray:
    g = np.zeros((M, M), dtype=complex)
    for y in range(M):
        for col, state in enumerate(basis):
            if not occupied(state, y):
                continue
            g[y, y] += np.conj(omega[col]) * omega[col]
            for x in range(M):
                if x == y or occupied(state, x):
                    continue
                target = state ^ (1 << y) ^ (1 << x)
                row = lookup[target]
                g[x, y] += np.conj(omega[row]) * omega[col]
    return g


def max_abs(a: np.ndarray) -> float:
    return float(np.max(np.abs(a))) if a.size else 0.0


def check(name: str, error: float, tol: float = TOL) -> None:
    status = "PASS" if error <= tol else "FAIL"
    print(f"{name:<72s} {status}  error={error:.3e}  tol={tol:.1e}")
    if error > tol:
        raise AssertionError(name)


def report(name: str, value: float) -> None:
    print(f"{name:<72s} VALUE {value:.12g}")


def main() -> int:
    np.set_printoptions(precision=8, suppress=True)
    rng = np.random.default_rng(2026072603)

    L, d, n, j = 3, 2, 4, 1.7
    torus = make_torus(L, d)
    M = len(torus.sites)
    basis, lookup = fixed_n_basis(M, n)
    dim = len(basis)
    omega = omega_state(dim)
    h = pair_hamiltonian(torus, basis, lookup, j)

    print("EXACT GRAPH-HODGE MEISSNER VERIFIER")
    print("python", platform.python_version())
    print("numpy", np.__version__)
    print("scipy", scipy.__version__)
    print(f"torus={L}^{d}, M={M}, n={n}, dim={dim}, j={j}")
    print("-" * 108)

    check("H(0) Hermiticity", max_abs(h - h.conj().T))
    check("uniform Dicke/AGP state is an exact zero state", float(np.linalg.norm(h @ omega)))
    evals = scipy.linalg.eigvalsh(h)
    check("unique zero eigenvalue", abs(evals[0]) + max(0.0, 1.0e-8 - evals[1]), 1.1e-8)
    report("first positive eigenvalue", float(evals[1]))

    g = pair_density_matrix(basis, lookup, omega, M)
    g_expected = np.full((M, M), alpha(M, n), dtype=float)
    np.fill_diagonal(g_expected, n / M)
    check("pair-density matrix equals exact combinatorial formula", max_abs(g - g_expected))
    lam = scipy.linalg.eigvalsh(g)
    lambda_max_exact = n * (M - n + 1) / M
    check("largest pair-density eigenvalue", abs(float(lam[-1]) - lambda_max_exact))
    report("ODLRO off-diagonal amplitude alpha(M,n)", alpha(M, n))

    f = rng.normal(size=M)
    f -= np.mean(f)
    gfun = rng.normal(size=M)
    gfun -= np.mean(gfun)
    vf = site_vector(basis, f)
    vg = site_vector(basis, gfun)
    check(
        "site-function inner product alpha <f,g>",
        abs(np.vdot(vf, vg) - alpha(M, n) * np.dot(f, gfun)),
    )
    lap0 = torus.incidence @ torus.incidence.T
    check(
        "H|f> = (j/4)|Delta f>",
        float(np.linalg.norm(h @ vf - (j / 4.0) * site_vector(basis, lap0 @ f))),
    )

    # Complete gauge covariance for a random finite field.
    field = rng.normal(scale=0.17, size=len(torus.edges))
    phi = rng.normal(scale=0.23, size=M)
    dphi = torus.incidence.T @ phi
    h_a = pair_hamiltonian(torus, basis, lookup, j, field)
    h_ag = pair_hamiltonian(torus, basis, lookup, j, field + dphi)
    phases = np.empty(dim, dtype=complex)
    for i, state in enumerate(basis):
        phases[i] = np.exp(-2.0j * sum(phi[x] for x in range(M) if occupied(state, x)))
    u = np.diag(phases)
    check(
        "bondwise gauge covariance U H[A] U^dagger = H[A+B^*phi]",
        max_abs(u @ h_a @ u.conj().T - h_ag),
    )

    # Arbitrary profile: source, stress, static Hodge, and dynamical filter.
    a = rng.normal(size=len(torus.edges))
    current, stress = current_and_stress(torus, basis, lookup, j, a)
    check("current Hermiticity", max_abs(current - current.conj().T))
    check("stress Hermiticity", max_abs(stress - stress.conj().T))
    div_a = torus.incidence @ a
    source_expected = 0.5j * j * site_vector(basis, div_a)
    check("J[a]|Omega> = i j/2 |B a>", float(np.linalg.norm(current @ omega - source_expected)))
    stress_expected = j * gamma(M, n) * float(a @ a)
    check("<T[a]> = j gamma ||a||^2", abs(float(np.real(expectation(omega, stress))) - stress_expected))

    pcyc = hodge_projector(torus)
    kubo_static = spectral_kubo(h, current, stress, omega, zeta=None)
    hodge_static = j * gamma(M, n) * float(a @ pcyc @ a)
    check("static Kubo curvature equals exact Hodge projection", abs(kubo_static - hodge_static), 2.0e-10)

    zeta = 0.73
    kubo_zeta = spectral_kubo(h, current, stress, omega, zeta=zeta)
    exact_zeta = exact_imaginary_kernel(torus, j, M, n, a, zeta)
    check("imaginary-frequency Kubo kernel equals graph filter", abs(kubo_zeta - exact_zeta), 2.0e-10)

    phi2 = rng.normal(size=M)
    phi2 -= np.mean(phi2)
    a_long = torus.incidence.T @ phi2
    j_long, t_long = current_and_stress(torus, basis, lookup, j, a_long)
    k_long = spectral_kubo(h, j_long, t_long, omega, zeta=None)
    check("pure-gradient static response vanishes", abs(k_long), 2.0e-10)

    # y-polarized profile varying in x: an exactly co-closed lattice Fourier mode.
    a_trans = np.zeros(len(torus.edges), dtype=float)
    q = 2.0 * math.pi / L
    for e, (tail, _head, mu) in enumerate(torus.edges):
        xcoord = torus.sites[tail][0]
        if mu == 1:
            a_trans[e] = math.cos(q * xcoord)
    check("chosen finite-q profile is exactly divergence free", float(np.linalg.norm(torus.incidence @ a_trans)))
    j_trans, t_trans = current_and_stress(torus, basis, lookup, j, a_trans)
    check("transverse current annihilates the ground state", float(np.linalg.norm(j_trans @ omega)))
    k_t0 = spectral_kubo(h, j_trans, t_trans, omega, zeta=None)
    k_tz = spectral_kubo(h, j_trans, t_trans, omega, zeta=zeta)
    exact_t = j * gamma(M, n) * float(a_trans @ a_trans)
    check("transverse static kernel equals diamagnetic value", abs(k_t0 - exact_t))
    check("transverse finite-frequency kernel is frequency independent", abs(k_tz - exact_t))
    report("transverse curvature for certificate profile", exact_t)

    qhat2_min = 4.0 * math.sin(math.pi / L) ** 2
    gap_exact = (j / 4.0) * qhat2_min
    check("n-pair gap equals one-pair random-walk gap on 3x3 torus", abs(float(evals[1]) - gap_exact), 2.0e-10)
    report("gap formula j sin^2(pi/L)", gap_exact)

    for L2 in (3, 4, 5, 6, 7):
        t2 = make_torus(L2, 2)
        b2, l2 = fixed_n_basis(len(t2.sites), 1)
        h2 = pair_hamiltonian(t2, b2, l2, j)
        e2 = scipy.linalg.eigvalsh(h2)
        target = j * math.sin(math.pi / L2) ** 2
        check(f"one-pair gap at L={L2}", abs(float(e2[1]) - target), 2.0e-10)

    # A genuine Aldous-content check: 4x3 torus, six hard-core pairs.
    rect = make_rect_torus((4, 3))
    b_rect, l_rect = fixed_n_basis(len(rect.sites), 6)
    h_rect = pair_hamiltonian_sparse(rect, b_rect, l_rect, j)
    evals_rect = np.sort(scipy.sparse.linalg.eigsh(h_rect, k=3, which="SA", tol=1.0e-11, return_eigenvectors=False))
    target_rect = (j / 4.0) * min(4.0 * math.sin(math.pi / 4.0) ** 2, 4.0 * math.sin(math.pi / 3.0) ** 2)
    check("multipair half-filled gap on 4x3 torus", abs(float(evals_rect[1]) - target_rect), 2.0e-9)
    report("4x3 half-filled sector dimension", float(len(b_rect)))

    # Selected two-block floor at fixed total pair number.
    M2 = 12
    j1, j2 = 0.9, 1.6
    jmin = min(j1, j2)
    for ntot in (2, 7, 12, 16, 22):
        lo, hi = max(0, ntot - M2), min(M2, ntot)
        values = []
        for n1 in range(lo, hi + 1):
            n2 = ntot - n1
            values.append(2.0 * (j1 * alpha(M2, n1) + j2 * alpha(M2, n2)))
        actual = min(values)
        if ntot <= M2:
            endpoint = 2.0 * jmin * ntot * (M2 - ntot) / (M2 * (M2 - 1))
        else:
            endpoint = 2.0 * jmin * (ntot - M2) * (2 * M2 - ntot) / (M2 * (M2 - 1))
        check(f"two-block selected-composition floor at total n={ntot}", max(0.0, endpoint - actual))

    mixing_norm, kernel_residual, floor_minimum = verify_backbone_operator_floor(rng)
    report("nonzero swap/residual target-Laplacian mixing norm", mixing_norm)
    check("swap source lies in complete target kernel", kernel_residual)
    check("dynamic and static operator backbone floors", max(0.0, -floor_minimum), 2.0e-10)

    print("-" * 108)
    print("OVERALL: PASS")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"OVERALL: FAIL ({type(exc).__name__}: {exc})", file=sys.stderr)
        raise
