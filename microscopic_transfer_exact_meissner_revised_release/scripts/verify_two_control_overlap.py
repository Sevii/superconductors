#!/usr/bin/env python3
"""Verify the first connected two-control Feshbach overlap coefficient."""
from __future__ import annotations
import itertools
import math
import sys
import numpy as np
import scipy
import scipy.linalg as la

RNG = np.random.default_rng(17290726)


def herm(a: np.ndarray) -> np.ndarray:
    return (a + a.conj().T) / 2.0


def opnorm(a: np.ndarray) -> float:
    return float(la.norm(a, 2))


def feshbach(A: np.ndarray, B: np.ndarray, t: float, Delta: float) -> np.ndarray:
    n = A.shape[0]
    I, Z = np.eye(n), np.zeros_like(A)
    hq = np.block(
        [[Delta * I, Z, t * B], [Z, Delta * I, t * A], [t * B, t * A, 2.0 * Delta * I]]
    )
    c = np.hstack([t * A, t * B, Z])
    return herm(-c @ la.solve(hq, c.conj().T, assume_a="her"))


def fermion_overlap_ground_compression(L: int = 3) -> tuple[np.ndarray, float]:
    n_modes = 4 * L

    def mode(a: int, x: int, spin: int) -> int:
        return ((a * L + x) * 2 + spin)

    basis = [sum(1 << i for i in comb) for comb in itertools.combinations(range(n_modes), 2)]
    lookup = {state: k for k, state in enumerate(basis)}
    dim = len(basis)

    def annihilate(state: int, i: int):
        if not ((state >> i) & 1):
            return None
        sign = -1 if (state & ((1 << i) - 1)).bit_count() % 2 else 1
        return state ^ (1 << i), sign

    def create(state: int, i: int):
        if (state >> i) & 1:
            return None
        sign = -1 if (state & ((1 << i) - 1)).bit_count() % 2 else 1
        return state | (1 << i), sign

    def cdagc(i: int, j: int) -> np.ndarray:
        out = np.zeros((dim, dim), dtype=complex)
        for col, state in enumerate(basis):
            step1 = annihilate(state, j)
            if step1 is None:
                continue
            state1, sign1 = step1
            step2 = create(state1, i)
            if step2 is None:
                continue
            state2, sign2 = step2
            out[lookup[state2], col] += sign1 * sign2
        return out

    def bridge(x: int) -> np.ndarray:
        y = (x + 1) % L
        out = np.zeros((dim, dim), dtype=complex)
        for spin in (0, 1):
            for i, j in [
                (mode(0, y, spin), mode(1, x, spin)),
                (mode(0, x, spin), mode(1, y, spin)),
            ]:
                out += cdagc(i, j) + cdagc(j, i)
        return out

    b0, b1 = bridge(0), bridge(1)
    anticom = b0 @ b1 + b1 @ b0
    r = herm(anticom @ anticom)
    agps = []
    for a in (0, 1):
        v = np.zeros(dim, dtype=complex)
        for x in range(L):
            state = (1 << mode(a, x, 0)) | (1 << mode(a, x, 1))
            v[lookup[state]] = 1.0 / math.sqrt(L)
        agps.append(v)
    z = np.column_stack(agps)
    return herm(z.conj().T @ r @ z), opnorm(r @ z)


def main() -> int:
    n, Delta = 5, 2.3
    A = herm(RNG.normal(size=(n, n)) + 1j * RNG.normal(size=(n, n)))
    B = herm(RNG.normal(size=(n, n)) + 1j * RNG.normal(size=(n, n)))
    A, B = A / opnorm(A), B / opnorm(B)
    anti = A @ B + B @ A
    c2 = -(A @ A + B @ B) / Delta
    c4 = -(anti @ anti) / (2.0 * Delta**3)

    print("TWO-CONTROL CONNECTED OVERLAP CERTIFICATE")
    print("numpy", np.__version__)
    print("scipy", scipy.__version__)
    print(f"Delta={Delta}")
    print("Formula: F=-t^2(A^2+B^2)/Delta - t^4{A,B}^2/(2Delta^3)+O(t^6)")
    print()
    ok, prev = True, None
    for t in [0.20, 0.14, 0.10, 0.07, 0.05, 0.035]:
        F = feshbach(A, B, t, Delta)
        extracted = (F - t**2 * c2) / t**4
        err = opnorm(extracted - c4)
        ratio = err / prev if prev is not None else float("nan")
        print(f"t={t:6.3f} ||c4(t)-c4||={err:.6e} err/t^2={err/t**2:.6e} ratio={ratio:.4f}")
        if prev is not None:
            ok &= err < prev
        prev = err
    ok &= prev is not None and prev < 2e-5
    print(f"||c2||={opnorm(c2):.10f}")
    print(f"||c4||={opnorm(c4):.10f}")

    print()
    print("One-pair product-AGP compression of the bridge-overlap operator")
    compression, source_norm = fermion_overlap_ground_compression(L=3)
    print("P_Z {B_0,B_1}^2 P_Z =")
    for row in compression:
        print("  [" + ", ".join(f"{z.real:.10f}" for z in row) + "]")
    print(f"||{{B_0,B_1}}^2 P_Z||={source_norm:.10f}")
    expected = (16.0 / 3.0) * np.ones((2, 2))
    ok &= opnorm(compression - expected) < 5e-10 and source_norm > 1.0
    print("Conclusion: the isolated bridge-overlap coefficient is not in the local two-sided parent-row ideal.")

    print()
    print("OVERALL:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
