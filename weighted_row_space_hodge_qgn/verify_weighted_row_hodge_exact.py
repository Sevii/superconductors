#!/usr/bin/env python3
"""Exact rational check of the weighted row-space Hodge identity."""
from __future__ import annotations

import sympy as sp


def assert_zero(name: str, matrix: sp.Matrix) -> None:
    matrix = sp.simplify(matrix)
    if matrix != sp.zeros(*matrix.shape):
        raise AssertionError(f"{name} failed:\n{matrix}")
    print(f"PASS  {name}")


def main() -> int:
    # State dimension 3, target dimension 2, ker(D)=span(e_3).
    D = sp.Matrix([[1, 0, 0], [0, 2, 0]])
    G = sp.Matrix([[2, 1], [1, 3]])
    Ginv = G.inv()
    M = D.T * G * D
    H = M / 2
    Mplus = M.pinv()
    Hplus = H.pinv()

    e3 = sp.Matrix([0, 0, 1])
    K = sp.Matrix([1, 2, 0])
    J = K * e3.T + e3 * K.T
    P = e3 * e3.T

    eta = G * D * Mplus * K
    factor_residual = D.T * eta - K
    cost = (eta.T * Ginv * eta)[0]
    schur = (sp.Rational(1, 2) * K.T * Hplus * K)[0]

    assert_zero("D^T eta_* = K", factor_residual)
    if sp.simplify(cost - schur) != 0:
        raise AssertionError(f"cost mismatch: {cost} versus {schur}")
    print("PASS  eta_*^T G^{-1} eta_* = (1/2) K^T H^+ K")
    if sp.simplify((P * J * P)[2, 2]) != 0:
        raise AssertionError("PJP safety failed")
    print("PASS  full-zero-space PJP condition")

    print("\nExact data")
    print("D =")
    sp.pprint(D)
    print("G =")
    sp.pprint(G)
    print("eta_* =")
    sp.pprint(eta)
    print(f"minimum cost = {sp.simplify(cost)}")

    # Exact nonuniqueness of the operator-norm minimizer for rank(P)=2.
    Dn = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 0]])
    Kn = sp.Matrix([[2, 0], [0, 1], [0, 0], [0, 0]])
    eta_n = sp.Matrix([[2, 0], [0, 1], [0, 0]])
    kappa_n = sp.Matrix([[0, 0], [0, 0], [0, 1]])
    assert_zero("norm example D^T eta = K", Dn.T * eta_n - Kn)
    assert_zero("norm example D^T kappa = 0", Dn.T * kappa_n)
    Cn = eta_n.T * eta_n
    Ck = kappa_n.T * kappa_n
    if Cn != sp.diag(4, 1) or Ck != sp.diag(0, 1):
        raise AssertionError("unexpected norm-example costs")
    if max((Cn + Ck).eigenvals()) != max(Cn.eigenvals()):
        raise AssertionError("nonprojected witness did not preserve operator norm")
    slack = 4 * sp.eye(2) - Cn - Ck
    if any(v < 0 for v in slack.eigenvals()):
        raise AssertionError("spectral-slack criterion failed")
    print("PASS  nonprojected rank-two witness has the same operator norm")

    # Strict selected branch counterexample: ker(D2)=span(e3,e4).
    D2 = sp.Matrix([[1, 0, 0, 0], [0, 1, 0, 0]])
    e3_4 = sp.Matrix([0, 0, 1, 0])
    e4_4 = sp.Matrix([0, 0, 0, 1])
    Psel = e3_4 * e3_4.T
    Jmix = e4_4 * e3_4.T + e3_4 * e4_4.T
    Z2 = sp.diag(0, 0, 1, 1)
    assert_zero("selected P J P = 0", Psel * Jmix * Psel)
    obstruction = Z2 * Jmix * e3_4
    if obstruction == sp.zeros(4, 1):
        raise AssertionError("selected-branch obstruction was not detected")
    print("PASS  selected branch has nonzero Z J P obstruction")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
