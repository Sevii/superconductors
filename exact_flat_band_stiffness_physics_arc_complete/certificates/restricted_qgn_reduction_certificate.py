#!/usr/bin/env python3
"""Finite-dimensional certificate for the restricted QGN reduction theorem.

The script checks, in a genuinely multiband basis, the identity

    C_n = rho_n C_1 + (n^2-rho_n) Gamma,
    rho_n = n(L-n)/(L-1),

for a positive-square, Hermitian-QGN Hamiltonian.  It uses:

* a random full-rank skew-unitary pairing matrix J on 2L one-particle modes;
* the complete Hermitian one-body commutant of the associated pseudospin SU(2);
* a frustration-free Hamiltonian H_0 = (1/2) sum_lambda S_lambda^2;
* two independent generic number-conserving first twist derivatives B_{lambda,i};
* the exact least-squares curvature formula, plus a direct eigenvalue
  finite-difference cross-check.

The random seed is fixed, and all headline claims are asserted.  This is a
finite-dimensional algebra certificate, not a proof of the lattice locality
corollary.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import sympy as sp
from scipy.linalg import eigh, null_space

Array = np.ndarray


@dataclass(frozen=True)
class SectorResult:
    L: int
    n: int
    sector_dimension: int
    ground_energy: float
    gap: float
    agp_overlap: float
    rho_n: float
    tensor_error: float
    finite_difference_error: float


def fermion_operators(number_modes: int) -> tuple[list[Array], list[Array]]:
    """Dense creation/annihilation matrices in occupation-number ordering."""
    dim = 1 << number_modes
    creators: list[Array] = []
    annihilators: list[Array] = []
    for p in range(number_modes):
        cd = np.zeros((dim, dim), dtype=np.complex128)
        lower_mask = (1 << p) - 1
        for mask in range(dim):
            if (mask >> p) & 1:
                continue
            parity = (mask & lower_mask).bit_count() & 1
            cd[mask | (1 << p), mask] = -1.0 if parity else 1.0
        creators.append(cd)
        annihilators.append(cd.conj().T)
    return creators, annihilators


def haar_unitary(n: int, rng: np.random.Generator) -> Array:
    z = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    q, r = np.linalg.qr(z)
    diag = np.diag(r)
    phases = np.where(np.abs(diag) > 0, diag / np.abs(diag), 1.0)
    return q @ np.diag(np.conjugate(phases))


def canonical_pairing_matrix(L: int) -> Array:
    J = np.zeros((2 * L, 2 * L), dtype=np.complex128)
    for a in range(L):
        J[2 * a, 2 * a + 1] = 1.0
        J[2 * a + 1, 2 * a] = -1.0
    return J


def hermitian_matrix_basis(n: int) -> list[Array]:
    """Hilbert-Schmidt orthonormal basis of Hermitian n-by-n matrices."""
    out: list[Array] = []
    for i in range(n):
        a = np.zeros((n, n), dtype=np.complex128)
        a[i, i] = 1.0
        out.append(a)
    for i in range(n):
        for j in range(i + 1, n):
            a = np.zeros((n, n), dtype=np.complex128)
            a[i, j] = a[j, i] = 1.0 / np.sqrt(2.0)
            out.append(a)
            a = np.zeros((n, n), dtype=np.complex128)
            a[i, j] = -1j / np.sqrt(2.0)
            a[j, i] = 1j / np.sqrt(2.0)
            out.append(a)
    return out


def pseudospin_scalar_basis(J: Array, tolerance: float = 1e-11) -> list[Array]:
    """Basis of Hermitian s satisfying sJ + J s^T = 0.

    For a skew-unitary J on 2L modes this is the compact symplectic Lie
    algebra in Hermitian convention, of real dimension L(2L+1).
    """
    vectors: list[Array] = []
    matrices: list[Array] = []
    for h in hermitian_matrix_basis(J.shape[0]):
        # Fixed-point projection of h -> -J h^T J^dagger.
        projected = 0.5 * (h - J @ h.T @ J.conj().T)
        projected = 0.5 * (projected + projected.conj().T)
        v = projected.reshape(-1).copy()
        for q in vectors:
            v -= np.vdot(q, v) * q
        norm = np.linalg.norm(v)
        if norm > tolerance:
            v /= norm
            vectors.append(v)
            matrices.append(v.reshape(J.shape))
    return matrices


def second_quantize(one_body: Array, creators: list[Array], annihilators: list[Array]) -> Array:
    out = np.zeros_like(creators[0])
    for p in range(one_body.shape[0]):
        for q in range(one_body.shape[1]):
            coefficient = one_body[p, q]
            if abs(coefficient) > 1e-14:
                out += coefficient * (creators[p] @ annihilators[q])
    return 0.5 * (out + out.conj().T)


def pair_creation(J: Array, creators: list[Array]) -> Array:
    eta_plus = np.zeros_like(creators[0])
    for p in range(J.shape[0]):
        for q in range(J.shape[1]):
            coefficient = J[p, q]
            if abs(coefficient) > 1e-14:
                eta_plus += 0.5 * coefficient * (creators[p] @ creators[q])
    return eta_plus


def sector_indices(number_modes: int, particle_number: int) -> Array:
    return np.asarray(
        [mask for mask in range(1 << number_modes) if mask.bit_count() == particle_number],
        dtype=np.int64,
    )


def restrict(operator: Array, indices: Array) -> Array:
    return operator[np.ix_(indices, indices)]


def normalized_agp(eta_plus: Array, vacuum: Array, n: int, indices: Array) -> Array:
    state = vacuum.copy()
    for _ in range(n):
        state = eta_plus @ state
    state = state[indices]
    norm = np.linalg.norm(state)
    if norm == 0:
        raise RuntimeError("AGP state vanished")
    return state / norm


def orthogonal_basis(psi: Array) -> Array:
    return null_space(np.conjugate(psi)[None, :])


def least_squares_hessian(
    factors: list[Array],
    derivatives_by_direction: list[list[Array]],
    psi: Array,
) -> Array:
    """Curvature tensor from the frustration-free least-squares formula."""
    D = np.vstack(factors)
    Q = orthogonal_basis(psi)
    DQ = D @ Q
    residuals: list[Array] = []
    for derivative_factors in derivatives_by_direction:
        source = np.concatenate([B @ psi for B in derivative_factors])
        coefficients, *_ = np.linalg.lstsq(DQ, -source, rcond=1e-12)
        residuals.append(source + DQ @ coefficients)
    dimension = len(residuals)
    tensor = np.empty((dimension, dimension), dtype=float)
    for i in range(dimension):
        for j in range(dimension):
            tensor[i, j] = float(np.vdot(residuals[i], residuals[j]).real)
    return tensor


def direct_curvature(
    factors: list[Array],
    derivative_factors: list[Array],
    step: float,
) -> float:
    """Five-point curvature of the exact lowest eigenvalue of H(t)."""

    def energy(t: float) -> float:
        H = sum(
            (S + t * B).conj().T @ (S + t * B)
            for S, B in zip(factors, derivative_factors)
        ) / 2.0
        return float(eigh(H, eigvals_only=True, subset_by_index=[0, 0])[0].real)

    values = [energy(s * step) for s in (-2, -1, 0, 1, 2)]
    em2, em1, e0, ep1, ep2 = values
    return float((-ep2 + 16 * ep1 - 30 * e0 + 16 * em1 - em2) / (12 * step**2))


def symbolic_agp_ratio() -> tuple[sp.Expr, sp.Symbol, sp.Symbol]:
    n, L = sp.symbols("n L", integer=True, positive=True)
    norm_agp = sp.factorial(n) * sp.factorial(L) / sp.factorial(L - n)
    norm_lowest_weight = sp.factorial(n - 1) * sp.factorial(L - 2) / sp.factorial(L - n - 1)
    ratio = sp.simplify(n**2 * norm_lowest_weight * L / norm_agp)
    return sp.factor(ratio), n, L


def build_generic_perturbations(
    number_directions: int,
    number_factors: int,
    number_modes: int,
    creators: list[Array],
    annihilators: list[Array],
    rng: np.random.Generator,
    scale: float,
) -> tuple[list[list[Array]], list[list[Array]]]:
    one_body_by_direction: list[list[Array]] = []
    many_body_by_direction: list[list[Array]] = []
    for _ in range(number_directions):
        one_body_direction: list[Array] = []
        many_body_direction: list[Array] = []
        for _ in range(number_factors):
            z = rng.normal(size=(number_modes, number_modes)) + 1j * rng.normal(
                size=(number_modes, number_modes)
            )
            b = 0.5 * (z + z.conj().T)
            b *= scale / np.linalg.norm(b)
            one_body_direction.append(b)
            many_body_direction.append(second_quantize(b, creators, annihilators))
        one_body_by_direction.append(one_body_direction)
        many_body_by_direction.append(many_body_direction)
    return one_body_by_direction, many_body_by_direction


def run_certificate(L: int, seed: int, step: float) -> tuple[dict, list[SectorResult]]:
    if L < 3:
        raise ValueError("L must be at least 3")
    rng = np.random.default_rng(seed)
    number_modes = 2 * L
    full_dimension = 1 << number_modes

    J0 = canonical_pairing_matrix(L)
    W = haar_unitary(number_modes, rng)
    J = W @ J0 @ W.T
    skew_error = float(np.linalg.norm(J + J.T))
    unitary_error = float(np.linalg.norm(J.conj().T @ J - np.eye(number_modes)))

    scalar_one_body = pseudospin_scalar_basis(J)
    expected_scalar_dimension = L * (2 * L + 1)
    scalar_condition_error = max(
        float(np.linalg.norm(s @ J + J @ s.T)) for s in scalar_one_body
    )

    creators, annihilators = fermion_operators(number_modes)
    eta_plus = pair_creation(J, creators)
    scalar_factors_full = [second_quantize(s, creators, annihilators) for s in scalar_one_body]
    one_body_B, many_body_B = build_generic_perturbations(
        number_directions=2,
        number_factors=len(scalar_factors_full),
        number_modes=number_modes,
        creators=creators,
        annihilators=annihilators,
        rng=rng,
        scale=0.02,
    )

    vacuum = np.zeros(full_dimension, dtype=np.complex128)
    vacuum[0] = 1.0

    agp_states: dict[int, Array] = {}
    sectors: dict[int, Array] = {}
    factor_sectors: dict[int, list[Array]] = {}
    B_sectors: dict[int, list[list[Array]]] = {}
    spectra: dict[int, tuple[float, float, float]] = {}

    for n in range(1, L):
        indices = sector_indices(number_modes, 2 * n)
        psi = normalized_agp(eta_plus, vacuum, n, indices)
        factors = [restrict(S, indices) for S in scalar_factors_full]
        H0 = sum(S.conj().T @ S for S in factors) / 2.0
        eigenvalues, eigenvectors = eigh(H0)
        gap = float(eigenvalues[1] - eigenvalues[0])
        overlap = float(abs(np.vdot(psi, eigenvectors[:, 0])))
        agp_states[n] = psi
        sectors[n] = indices
        factor_sectors[n] = factors
        B_sectors[n] = [
            [restrict(B, indices) for B in direction] for direction in many_body_B
        ]
        spectra[n] = (float(eigenvalues[0].real), gap, overlap)

    # One-pair longitudinal tensor Gamma, computed independently from expectation values.
    psi_one = agp_states[1]
    beta = np.zeros((2, len(scalar_factors_full)), dtype=np.complex128)
    for i in range(2):
        for lam, B in enumerate(B_sectors[1][i]):
            beta[i, lam] = np.vdot(psi_one, B @ psi_one)
    gamma = np.real(beta @ beta.conj().T)

    curvature_tensors: dict[int, Array] = {}
    for n in range(1, L):
        curvature_tensors[n] = least_squares_hessian(
            factor_sectors[n], B_sectors[n], agp_states[n]
        )

    test_direction = np.asarray([1.0, -0.37])
    results: list[SectorResult] = []
    maximum_tensor_error = 0.0
    maximum_fd_error = 0.0
    for n in range(1, L):
        rho = n * (L - n) / (L - 1)
        prediction = rho * curvature_tensors[1] + (n * n - rho) * gamma
        tensor_error = float(np.linalg.norm(curvature_tensors[n] - prediction))
        maximum_tensor_error = max(maximum_tensor_error, tensor_error)

        derivative_factors = [
            sum(test_direction[i] * B_sectors[n][i][lam] for i in range(2))
            for lam in range(len(factor_sectors[n]))
        ]
        fd = direct_curvature(factor_sectors[n], derivative_factors, step)
        exact_directional = float(test_direction @ curvature_tensors[n] @ test_direction)
        fd_error = abs(fd - exact_directional)
        maximum_fd_error = max(maximum_fd_error, fd_error)

        e0, gap, overlap = spectra[n]
        results.append(
            SectorResult(
                L=L,
                n=n,
                sector_dimension=len(sectors[n]),
                ground_energy=e0,
                gap=gap,
                agp_overlap=overlap,
                rho_n=float(rho),
                tensor_error=tensor_error,
                finite_difference_error=float(fd_error),
            )
        )

    symbolic_ratio, symbolic_n, symbolic_L = symbolic_agp_ratio()
    expected_ratio = symbolic_n * (symbolic_L - symbolic_n) / (symbolic_L - 1)

    summary = {
        "L": L,
        "number_modes": number_modes,
        "seed": seed,
        "pairing_skew_error": skew_error,
        "pairing_unitarity_error": unitary_error,
        "scalar_basis_dimension": len(scalar_one_body),
        "expected_scalar_basis_dimension": expected_scalar_dimension,
        "scalar_condition_error": scalar_condition_error,
        "symbolic_agp_ratio": str(symbolic_ratio),
        "gamma_tensor": gamma.tolist(),
        "one_pair_curvature_tensor": curvature_tensors[1].tolist(),
        "maximum_tensor_identity_error": maximum_tensor_error,
        "maximum_finite_difference_error": maximum_fd_error,
    }

    # Headline assertions.
    assert skew_error < 1e-12
    assert unitary_error < 1e-12
    assert len(scalar_one_body) == expected_scalar_dimension
    assert scalar_condition_error < 1e-11
    assert sp.simplify(symbolic_ratio - expected_ratio) == 0
    assert all(result.ground_energy < 1e-10 for result in results)
    assert all(result.gap > 1.0 for result in results)
    assert all(abs(result.agp_overlap - 1.0) < 1e-10 for result in results)
    assert maximum_tensor_error < 1e-11
    assert maximum_fd_error < 1e-6

    return summary, results


def write_outputs(output_directory: Path, summary: dict, results: Iterable[SectorResult]) -> None:
    output_directory.mkdir(parents=True, exist_ok=True)
    (output_directory / "certificate_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    rows = [asdict(result) for result in results]
    with (output_directory / "certificate_sectors.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--L", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260723)
    parser.add_argument("--step", type=float, default=5e-4)
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "restricted_qgn_certificate_output",
    )
    args = parser.parse_args()

    summary, results = run_certificate(args.L, args.seed, args.step)
    write_outputs(args.out, summary, results)

    print("Restricted QGN reduction certificate")
    print(json.dumps(summary, indent=2))
    for result in results:
        print(asdict(result))
    print("PASS")
    print(f"Wrote outputs to {args.out}")


if __name__ == "__main__":
    main()
