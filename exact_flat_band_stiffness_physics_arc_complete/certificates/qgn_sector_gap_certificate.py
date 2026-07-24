#!/usr/bin/env python3
"""Certificates for the QGN sector-gap reduction and one-pair Fourier formula.

The script checks four claims used in the companion note:

1. For a projected PSD-Hubbard Hamiltonian built from rank-one frame matrices,
   the finite-size sector gaps observed in connected examples equal the one-pair
   gap (numerical evidence for the rank-one Aldous property).
2. The exact one-pair Hamiltonian of a translation-invariant UPC band is reduced
   to a small orbital Gram matrix in each momentum-transfer block.
3. The Gram formula agrees with direct dense commutator diagonalization on small
   Gao-Han-Khalaf Model-I/Model-II tori.
4. Generic Hermitian square factors do NOT have filling-independent gaps, showing
   that rank-one projected-density structure is essential.

The equality in item 1 is evidence, not a proof. The representation-theoretic
reduction in the note identifies the remaining statement precisely.
"""
from __future__ import annotations

import csv
import itertools
import json
import math
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy.linalg import eigvalsh

Array = NDArray[np.complex128]
ROOT = Path(__file__).resolve().parent


def masks(L: int, n: int) -> tuple[int, ...]:
    return tuple(sum(1 << i for i in c) for c in itertools.combinations(range(L), n))


def dgamma(a: Array, n: int) -> Array:
    """Second quantization dGamma_n(a) on Lambda^n C^L."""
    L = a.shape[0]
    basis = masks(L, n)
    index = {m: i for i, m in enumerate(basis)}
    out = np.zeros((len(basis), len(basis)), dtype=np.complex128)
    for q in range(L):
        for p in range(L):
            z = a[p, q]
            if abs(z) < 1e-15:
                continue
            for col, mask in enumerate(basis):
                if not ((mask >> q) & 1):
                    continue
                sign = -1 if (mask & ((1 << q) - 1)).bit_count() % 2 else 1
                removed = mask ^ (1 << q)
                if (removed >> p) & 1:
                    continue
                if (removed & ((1 << p) - 1)).bit_count() % 2:
                    sign *= -1
                out[index[removed | (1 << p)], col] += sign * z
    return out


def commutator_hamiltonian(generators: Sequence[Array], n: int, U: float = 1.0) -> Array:
    """H_n=(U/2) sum_x ad_{dGamma_n(a_x)}^2 on End(Lambda^n C^L)."""
    L = generators[0].shape[0]
    d = math.comb(L, n)
    ident = np.eye(d, dtype=np.complex128)
    H = np.zeros((d * d, d * d), dtype=np.complex128)
    for a in generators:
        A = dgamma(a, n)
        # Column-major vectorization: vec(AT-TA)=(I kron A-A^T kron I)vec(T).
        S = np.kron(ident, A) - np.kron(A.T, ident)
        H += (U / 2) * (S.conj().T @ S)
    return (H + H.conj().T) / 2


def gap_and_nullity(H: Array, tol: float = 2e-9) -> tuple[float, int]:
    vals = eigvalsh(H)
    nullity = int(np.sum(np.abs(vals) < tol))
    positive = vals[vals >= tol]
    return (float(positive[0]) if len(positive) else float("nan"), nullity)


def random_parseval_frame(L: int, m: int, rng: np.random.Generator) -> list[Array]:
    z = rng.normal(size=(m, L)) + 1j * rng.normal(size=(m, L))
    q, _ = np.linalg.qr(z)
    return [np.conjugate(q[x, :]) for x in range(m)]


def rank_one_generators(frame: Sequence[Array]) -> list[Array]:
    return [np.outer(v, np.conjugate(v)) for v in frame]


def grid_2d(Nx: int, Ny: int) -> Array:
    return np.array(
        [(2 * np.pi * ix / Nx, 2 * np.pi * iy / Ny)
         for iy in range(Ny) for ix in range(Nx)],
        dtype=float,
    )


def model_I_d(k: Array) -> NDArray[np.float64]:
    x, y = k[:, 0], k[:, 1]
    t1, t2 = 1.0, 1 / np.sqrt(2)
    z = -t1 * (
        np.exp(1j * np.pi / 4)
        + np.exp(1j * (x + y + np.pi / 4))
        + np.exp(1j * (x - np.pi / 4))
        + np.exp(1j * (y - np.pi / 4))
    )
    return np.column_stack((z.real, z.imag, -2 * t2 * (np.cos(y) - np.cos(x))))


def model_II_d(k: Array, xi: float) -> NDArray[np.float64]:
    x, y = k[:, 0], k[:, 1]
    alpha = xi * (np.cos(x) + np.cos(y))
    return np.column_stack((-np.sin(alpha), -np.cos(alpha), np.zeros_like(alpha)))


def lower_spinors(d: NDArray[np.float64]) -> Array:
    """Deterministic eigenvectors of d.sigma with the lower eigenvalue."""
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.array([[1, 0], [0, -1]], complex)
    out = np.empty((len(d), 2), dtype=np.complex128)
    for i, (x, y, z) in enumerate(d):
        norm = math.sqrt(float(x*x + y*y + z*z))
        if norm < 1e-12:
            raise ValueError("sampled band touching")
        vals, vecs = np.linalg.eigh((x * sx + y * sy + z * sz) / norm)
        u = vecs[:, np.argmin(vals)]
        pivot = int(np.argmax(np.abs(u)))
        u *= np.exp(-1j * np.angle(u[pivot]))
        out[i] = u
    return out


def frame_from_spinors(k: Array, spinor: Array, Nx: int, Ny: int) -> list[Array]:
    V = Nx * Ny
    frame: list[Array] = []
    for Ry in range(Ny):
        for Rx in range(Nx):
            phase = np.exp(-1j * (k[:, 0] * Rx + k[:, 1] * Ry)) / np.sqrt(V)
            for alpha in range(spinor.shape[1]):
                frame.append(phase * np.conjugate(spinor[:, alpha]))
    return frame


def one_pair_gram_gap(spinor: Array, Nx: int, Ny: int, U: float = 1.0) -> tuple[float, tuple[int, int], float]:
    """Exact UPC one-pair gap from orbital Gram matrices.

    Returns (gap, minimizing momentum-transfer label, epsilon).  The formula is
      H_Q = U [epsilon I - (1/V) sum_alpha |a_alpha,Q><a_alpha,Q|],
      a_alpha,Q(k)=conj(u_alpha(k+Q)) u_alpha(k).
    """
    V, norb = spinor.shape
    weights = np.mean(np.abs(spinor) ** 2, axis=0)
    if np.max(weights) - np.min(weights) > 2e-10:
        raise ValueError(f"UPC failed: orbital weights {weights}")
    epsilon = float(np.mean(weights))
    grid = spinor.reshape(Ny, Nx, norb)
    best = float("inf")
    best_q = (0, 0)
    for qy in range(Ny):
        for qx in range(Nx):
            shifted = np.roll(np.roll(grid, -qy, axis=0), -qx, axis=1)
            a = np.conjugate(shifted) * grid
            avecs = a.reshape(V, norb)
            G = (avecs.conj().T @ avecs) / V
            lambdas = eigvalsh((G + G.conj().T) / 2)
            energies = U * (epsilon - lambdas)
            for e in energies:
                if e > 2e-10 and e < best:
                    best = float(e)
                    best_q = (qx, qy)
            # The orthogonal complement of the orbital span has energy U*epsilon.
            if V > norb and U * epsilon < best:
                best = U * epsilon
                best_q = (qx, qy)
    return best, best_q, epsilon


def dense_one_pair_gap(spinor: Array, Nx: int, Ny: int, U: float = 1.0) -> float:
    k = grid_2d(Nx, Ny)
    frame = frame_from_spinors(k, spinor, Nx, Ny)
    H = commutator_hamiltonian(rank_one_generators(frame), 1, U=U)
    gap, nullity = gap_and_nullity(H)
    if nullity != 1:
        raise AssertionError(f"unexpected one-pair nullity {nullity}")
    return gap


def run_rank_one_tests() -> list[dict]:
    rng = np.random.default_rng(20260723)
    rows: list[dict] = []
    cases = [(4, 7, 120), (5, 8, 40)]
    worst = 0.0
    retained = 0
    for L, m, count in cases:
        for sample in range(count):
            gens = rank_one_generators(random_parseval_frame(L, m, rng))
            gaps = []
            nullities = []
            for n in range(1, L):
                gap, nullity = gap_and_nullity(commutator_hamiltonian(gens, n))
                gaps.append(gap)
                nullities.append(nullity)
            if min(gaps) < 1e-7:
                continue
            retained += 1
            deviation = max(abs(g / gaps[0] - 1.0) for g in gaps)
            worst = max(worst, deviation)
            rows.append({
                "L": L, "frame_vectors": m, "sample": sample,
                "gap_1": gaps[0], "min_gap": min(gaps), "max_gap": max(gaps),
                "max_relative_deviation": deviation,
                "nullities": "+".join(map(str, nullities)),
            })
    if retained < 100:
        raise AssertionError(f"too few retained random cases: {retained}")
    if worst > 3e-7:
        raise AssertionError(f"rank-one gap mismatch above numerical tolerance: {worst}")
    print(f"rank-one tests: retained={retained}, worst relative deviation={worst:.3e}")
    return rows


def run_ghk_scaling() -> tuple[list[dict], list[dict]]:
    scaling: list[dict] = []
    validations: list[dict] = []
    model_specs: list[tuple[str, Callable[[Array], NDArray[np.float64]]]] = [
        ("GHK-Model-I", model_I_d),
        ("GHK-Model-II-xi0.5", lambda k: model_II_d(k, 0.5)),
        ("GHK-Model-II-xi1", lambda k: model_II_d(k, 1.0)),
        ("GHK-Model-II-xi1.5", lambda k: model_II_d(k, 1.5)),
    ]
    for name, model in model_specs:
        for N in (4, 6, 8, 10, 12, 16, 20, 24, 32, 40):
            k = grid_2d(N, N)
            spinor = lower_spinors(model(k))
            gap, qmin, epsilon = one_pair_gram_gap(spinor, N, N)
            scaling.append({
                "model": name, "N": N, "V": N*N, "gap": gap,
                "N2_gap": N*N*gap, "qmin_x": qmin[0], "qmin_y": qmin[1],
                "epsilon": epsilon,
            })
        for N in (2, 3, 4):
            k = grid_2d(N, N)
            spinor = lower_spinors(model(k))
            gram_gap, qmin, _ = one_pair_gram_gap(spinor, N, N)
            dense_gap = dense_one_pair_gap(spinor, N, N)
            err = abs(gram_gap - dense_gap)
            validations.append({
                "model": name, "N": N, "gram_gap": gram_gap,
                "dense_gap": dense_gap, "absolute_error": err,
                "qmin_x": qmin[0], "qmin_y": qmin[1],
            })
            if err > 3e-9:
                raise AssertionError((name, N, gram_gap, dense_gap, err))
    # Model-II asymptotic coefficient: N^2 Delta -> pi^2 xi^2 / 4.
    for xi in (0.5, 1.0, 1.5):
        name = f"GHK-Model-II-xi{xi:g}"
        last = [r for r in scaling if r["model"] == name and r["N"] == 40][0]
        target = math.pi**2 * xi**2 / 4
        rel = abs(last["N2_gap"] - target) / target
        if rel > 0.012:
            raise AssertionError((name, last["N2_gap"], target, rel))
    print("GHK Gram/dense validations: PASS")
    return scaling, validations


def run_generic_hermitian_counterexample() -> dict:
    """Find a deterministic non-rank-one Hermitian-factor gap mismatch."""
    rng = np.random.default_rng(20260723)
    best: dict | None = None
    for trial in range(600):
        L = 4
        gens: list[Array] = []
        for _ in range(3):
            z = rng.normal(size=(L, L)) + 1j * rng.normal(size=(L, L))
            a = (z + z.conj().T) / 2
            a /= np.linalg.norm(a)
            gens.append(a)
        gap1, null1 = gap_and_nullity(commutator_hamiltonian(gens, 1))
        gap2, null2 = gap_and_nullity(commutator_hamiltonian(gens, 2))
        if null1 != 1 or null2 != 1 or gap1 < 1e-8:
            continue
        ratio = gap2 / gap1
        if best is None or ratio < best["ratio_gap2_gap1"]:
            best = {
                "trial": trial, "gap_1": gap1, "gap_2": gap2,
                "ratio_gap2_gap1": ratio,
                "generators_real_imag": [
                    {"real": a.real.tolist(), "imag": a.imag.tolist()} for a in gens
                ],
            }
    if best is None or best["ratio_gap2_gap1"] >= 0.8:
        raise AssertionError(f"failed to find generic Hermitian mismatch: {best}")
    print(
        "generic Hermitian control: gap2/gap1=",
        f"{best['ratio_gap2_gap1']:.6f}",
    )
    return best


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    rank_rows = run_rank_one_tests()
    scaling, validations = run_ghk_scaling()
    generic = run_generic_hermitian_counterexample()

    write_csv(ROOT / "qgn_sector_gap_rank_one_tests.csv", rank_rows)
    write_csv(ROOT / "qgn_sector_gap_scaling.csv", scaling)
    write_csv(ROOT / "qgn_sector_gap_gram_validation.csv", validations)
    (ROOT / "qgn_sector_gap_generic_hermitian_control.json").write_text(
        json.dumps(generic, indent=2)
    )
    summary = {
        "rank_one_cases": len(rank_rows),
        "rank_one_worst_relative_deviation": max(r["max_relative_deviation"] for r in rank_rows),
        "gram_validation_worst_absolute_error": max(r["absolute_error"] for r in validations),
        "generic_hermitian_gap_ratio": generic["ratio_gap2_gap1"],
        "model_II_asymptotic_targets": {
            str(xi): math.pi**2 * xi**2 / 4 for xi in (0.5, 1.0, 1.5)
        },
    }
    (ROOT / "qgn_sector_gap_certificate_summary.json").write_text(
        json.dumps(summary, indent=2)
    )
    print(json.dumps(summary, indent=2))
    print("PASS")


if __name__ == "__main__":
    main()
