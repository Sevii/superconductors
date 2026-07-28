#!/usr/bin/env python3
"""Finite-matrix verifier for the weighted row-space Hodge identities.

Checks:
  * range/cokernel criterion;
  * explicit minimizing target field;
  * exact Schur-complement cost;
  * weighted target-space projection;
  * Pythagorean increase under a target-cokernel perturbation;
  * approximate-witness inequality;
  * nonuniqueness of the operator-norm minimizer on a rank-two target domain.

The script uses deterministic random seeds and exits nonzero on failure.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Tuple

import numpy as np


def hermitian(a: np.ndarray) -> np.ndarray:
    return (a + a.conj().T) / 2


def psd_sqrt_and_inv(a: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    vals, vecs = np.linalg.eigh(hermitian(a))
    if np.min(vals) <= 0:
        raise ValueError("matrix is not strictly positive")
    root = (vecs * np.sqrt(vals)) @ vecs.conj().T
    inv_root = (vecs * (1 / np.sqrt(vals))) @ vecs.conj().T
    return root, inv_root


def projector_from_basis(v: np.ndarray) -> np.ndarray:
    if v.size == 0:
        return np.zeros((v.shape[0], v.shape[0]), dtype=complex)
    q, _ = np.linalg.qr(v)
    return q @ q.conj().T


def nullspace(a: np.ndarray, tol: float = 1e-11) -> np.ndarray:
    _, s, vh = np.linalg.svd(a, full_matrices=True)
    rank = int(np.sum(s > tol * max(1.0, s[0] if s.size else 1.0)))
    return vh[rank:].conj().T


def opnorm(a: np.ndarray) -> float:
    return float(np.linalg.norm(a, 2))


@dataclass
class CheckResult:
    name: str
    error: float
    tolerance: float

    @property
    def passed(self) -> bool:
        return self.error <= self.tolerance


def run(seed: int = 20260727, tol: float = 5e-9) -> list[CheckResult]:
    rng = np.random.default_rng(seed)

    # D: state (n) -> target (m), with a prescribed two-dimensional kernel.
    n, m, nullity, srcdim = 9, 7, 2, 3
    raw = rng.normal(size=(m, n - nullity)) + 1j * rng.normal(size=(m, n - nullity))
    # Make full column rank in the active state block.
    u, _, vh = np.linalg.svd(raw, full_matrices=False)
    active = u[:, : n - nullity] @ np.diag(np.linspace(0.8, 2.0, n - nullity)) @ vh
    D = np.concatenate([active, np.zeros((m, nullity), dtype=complex)], axis=1)

    x = rng.normal(size=(m, m)) + 1j * rng.normal(size=(m, m))
    G = x.conj().T @ x + 0.7 * np.eye(m)
    G_root, G_inv_root = psd_sqrt_and_inv(G)
    G_inv = np.linalg.inv(G)

    M = hermitian(D.conj().T @ G @ D)
    H = M / 2
    M_plus = np.linalg.pinv(M, rcond=1e-12)
    H_plus = np.linalg.pinv(H, rcond=1e-12)

    Zbasis = nullspace(D)
    Z = projector_from_basis(Zbasis)
    Q = np.eye(n) - Z

    # Source columns live in the complete zero manifold domain, while their
    # images are forced into Ran(D^dagger).
    domain = Zbasis[:, :srcdim] if Zbasis.shape[1] >= srcdim else Zbasis
    srcdim = domain.shape[1]
    Xi = rng.normal(size=(m, srcdim)) + 1j * rng.normal(size=(m, srcdim))
    K = D.conj().T @ Xi

    eta_star = G @ D @ M_plus @ K
    B = G_root @ D
    Pi = B @ M_plus @ B.conj().T

    # An arbitrary feasible witness is eta_star plus target-cokernel columns.
    ker_Bdag = nullspace(B.conj().T)
    if ker_Bdag.shape[1]:
        coeff = rng.normal(size=(ker_Bdag.shape[1], srcdim)) + 1j * rng.normal(
            size=(ker_Bdag.shape[1], srcdim)
        )
        zeta = G_root @ (ker_Bdag @ coeff)  # D^dagger zeta = 0
    else:
        zeta = np.zeros_like(eta_star)
    eta0 = eta_star + zeta

    cost_star = hermitian(eta_star.conj().T @ G_inv @ eta_star)
    schur = hermitian(0.5 * K.conj().T @ H_plus @ K)
    cost0 = hermitian(eta0.conj().T @ G_inv @ eta0)
    pyth_rhs = hermitian(cost_star + zeta.conj().T @ G_inv @ zeta)

    z0 = G_inv_root @ eta0
    zstar = G_inv_root @ eta_star

    results = [
        CheckResult("range factorization", opnorm(D.conj().T @ eta_star - K), tol),
        CheckResult("zero-space orthogonality", opnorm(Z @ K), tol),
        CheckResult("Schur cost identity", opnorm(cost_star - schur), 20 * tol),
        CheckResult("target projection is Hermitian", opnorm(Pi - Pi.conj().T), 20 * tol),
        CheckResult("target projection is idempotent", opnorm(Pi @ Pi - Pi), 50 * tol),
        CheckResult("projected witness", opnorm(zstar - Pi @ z0), 50 * tol),
        CheckResult("Pythagorean identity", opnorm(cost0 - pyth_rhs), 50 * tol),
    ]

    # Approximate-witness inequality. Remove a compatible piece from eta_star.
    delta = 0.15 * (rng.normal(size=eta_star.shape) + 1j * rng.normal(size=eta_star.shape))
    eta_app = eta_star - delta
    c = Q @ (K - D.conj().T @ eta_app)
    theta = 0.7
    lhs = hermitian(0.5 * K.conj().T @ H_plus @ K)
    rhs = hermitian(
        (1 + theta) * eta_app.conj().T @ G_inv @ eta_app
        + 0.5 * (1 + 1 / theta) * c.conj().T @ H_plus @ c
    )
    min_eig = float(np.min(np.linalg.eigvalsh(hermitian(rhs - lhs))))
    results.append(CheckResult("approximate-witness Loewner bound", max(0.0, -min_eig), 100 * tol))

    # A nonprojected witness can be an operator-norm minimizer when rank(P)>=2.
    # Here the Loewner-minimal cost is diag(4,1), while a positive target-cokernel
    # increment diag(0,1) fits inside the spectral slack and leaves the norm 4.
    Dn = np.array(
        [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]],
        dtype=complex,
    )
    Kn = np.array(
        [[2.0, 0.0], [0.0, 1.0], [0.0, 0.0], [0.0, 0.0]],
        dtype=complex,
    )
    eta_n = np.array(
        [[2.0, 0.0], [0.0, 1.0], [0.0, 0.0]],
        dtype=complex,
    )
    kappa_n = np.array(
        [[0.0, 0.0], [0.0, 0.0], [0.0, 1.0]],
        dtype=complex,
    )
    Cn = hermitian(eta_n.conj().T @ eta_n)
    Ck = hermitian(kappa_n.conj().T @ kappa_n)
    lam = opnorm(Cn)
    slack = lam * np.eye(2) - Cn
    results.extend(
        [
            CheckResult("norm-example factorization", opnorm(Dn.conj().T @ eta_n - Kn), tol),
            CheckResult("norm-example target cokernel", opnorm(Dn.conj().T @ kappa_n), tol),
            CheckResult(
                "nonprojected norm equality",
                abs(opnorm(Cn + Ck) - opnorm(Cn)),
                tol,
            ),
            CheckResult(
                "spectral-slack criterion",
                max(0.0, -float(np.min(np.linalg.eigvalsh(hermitian(slack - Ck))))),
                tol,
            ),
            CheckResult(
                "nonprojected perturbation is nonzero",
                max(0.0, 0.5 - opnorm(kappa_n)),
                tol,
            ),
        ]
    )

    # Strict selected subspace: PJ P = 0 can coexist with mixing into Z-P.
    if Zbasis.shape[1] >= 2:
        pvec = Zbasis[:, [0]]
        rvec = Zbasis[:, [1]]
        Psel = pvec @ pvec.conj().T
        # Hermitian J with J p = r and P J P = 0.
        Jmix = rvec @ pvec.conj().T + pvec @ rvec.conj().T
        selected_source = Jmix @ pvec
        pjpp = pvec.conj().T @ Jmix @ pvec
        obstruction = Z @ selected_source
        results.extend(
            [
                CheckResult("selected PJP is zero", opnorm(pjpp), tol),
                # This should be nonzero. Record error as distance below a threshold.
                CheckResult(
                    "selected-branch cokernel detected",
                    max(0.0, 0.5 - opnorm(obstruction)),
                    tol,
                ),
            ]
        )

    return results


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20260727)
    parser.add_argument("--tol", type=float, default=5e-9)
    args = parser.parse_args()

    results = run(args.seed, args.tol)
    failed = [r for r in results if not r.passed]
    print("Weighted row-space Hodge finite-matrix verification")
    print(f"seed={args.seed} tolerance={args.tol:.2e}")
    for r in results:
        status = "PASS" if r.passed else "FAIL"
        print(f"{status:4s}  {r.name:40s} error={r.error:.3e} limit={r.tolerance:.3e}")
    if failed:
        print(f"\n{len(failed)} check(s) failed.")
        return 1
    print("\nAll checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
