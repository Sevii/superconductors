#!/usr/bin/env python3
"""Exact finite-size counterexample to the UPC-only filling law.

The model is the two-orbital, M=4 construction described in the accompanying
answer.  SymPy evaluates the finite-dimensional Kohn curvature exactly.

Spin-up flat-band vector:
    u_up(k) = ( cos(pi/4 + c sin 4k),
                exp[i b sin k] sin(pi/4 + c sin 4k) )^T
with c=1/8 and b=pi/2.  Spin down is u_down(k)=u_up(-k)^*.
The parent Hamiltonian h_sigma(k)=I-|u_sigma(k)><u_sigma(k)| has exactly-flat
bands at energies 0 and 1.

At the M=4 momenta k_m=m*pi/2, sin(4 k_m)=0, so the finite-torus UPC holds.
After a twist A, sin(4(k_m+A))=sin(4A), and the projected Hubbard Hamiltonian
factorizes as
    H_c(A) = f(A) H_0(A),
    f(A) = cos^2(2 c sin 4A).
For c=1/8, f''(0)=-2.  Since E_n(0)=-n/2, this adds +n to every energy
curvature.  The base c=0 curvature has the desired hard-core form, but the
added linear-in-n term does not.
"""

from __future__ import annotations

import itertools
from functools import lru_cache

import sympy as sp

M = 4
I = sp.I
PI = sp.pi
B = PI / 2
C = sp.Rational(1, 8)


@lru_cache(maxsize=None)
def masks(n: int) -> tuple[int, ...]:
    return tuple(
        sum(1 << i for i in combination)
        for combination in itertools.combinations(range(M), n)
    )


@lru_cache(maxsize=None)
def transition_table(n: int) -> dict[tuple[int, int], tuple[tuple[int, int], ...]]:
    """Matrix elements of d_p^dagger d_q in an n-particle spin sector."""
    basis = masks(n)
    index = {mask: i for i, mask in enumerate(basis)}
    table: dict[tuple[int, int], tuple[tuple[int, int], ...]] = {}

    for p in range(M):
        for q in range(M):
            entries: list[tuple[int, int]] = []
            for mask in basis:
                if not ((mask >> q) & 1):
                    entries.append((-1, 0))
                    continue

                sign = -1 if (mask & ((1 << q) - 1)).bit_count() % 2 else 1
                after_annihilation = mask ^ (1 << q)
                if (after_annihilation >> p) & 1:
                    entries.append((-1, 0))
                    continue

                if (after_annihilation & ((1 << p) - 1)).bit_count() % 2:
                    sign *= -1
                final_mask = after_annihilation | (1 << p)
                entries.append((index[final_mask], sign))
            table[(p, q)] = tuple(entries)

    return table


# For the base model c=0, the second orbital phase is
# z_m(A)=exp[i b sin(k_m+A)].  At b=pi/2 and k_m=m*pi/2, the values and
# derivatives below are exact.
Z = (sp.Integer(1), I, sp.Integer(1), -I)
Z_PRIME = (I * B, sp.Integer(0), -I * B, sp.Integer(0))
Z_SECOND = (-B**2, B, -B**2, B)
Z_SERIES = tuple((Z[m], Z_PRIME[m], Z_SECOND[m] / 2) for m in range(M))


def conjugate_series(series: tuple[sp.Expr, sp.Expr, sp.Expr]) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    return tuple(sp.conjugate(value) for value in series)


def multiply_series(
    left: tuple[sp.Expr, sp.Expr, sp.Expr],
    right: tuple[sp.Expr, sp.Expr, sp.Expr],
) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Multiply Taylor series through A^2; entries are ordinary coefficients."""
    return tuple(
        sp.expand(sum(left[j] * right[order - j] for j in range(order + 1)))
        for order in range(3)
    )


def phase_product_series(p: int, q: int, r: int, s: int) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    result: tuple[sp.Expr, sp.Expr, sp.Expr] = (sp.Integer(1), sp.Integer(0), sp.Integer(0))
    factors = (
        conjugate_series(Z_SERIES[p]),
        Z_SERIES[q],
        conjugate_series(Z_SERIES[r]),
        Z_SERIES[s],
    )
    for factor in factors:
        result = multiply_series(result, factor)
    return tuple(sp.simplify(value) for value in result)


def base_hamiltonian_derivatives(n: int) -> tuple[sp.Matrix, sp.Matrix, sp.Matrix]:
    """Return exact H_0(0), dH_0/dA, and d^2H_0/dA^2 for c=0."""
    basis = masks(n)
    table = transition_table(n)
    dim_spin = len(basis)
    dimension = dim_spin**2
    coefficients = [sp.zeros(dimension, dimension) for _ in range(3)]
    occupied = [[q for q in range(M) if (mask >> q) & 1] for mask in basis]
    cache: dict[tuple[int, int, int, int], tuple[sp.Expr, sp.Expr, sp.Expr]] = {}

    for up_index in range(dim_spin):
        for down_index in range(dim_spin):
            column = up_index * dim_spin + down_index
            for q in occupied[up_index]:
                for s in occupied[down_index]:
                    for p in range(M):
                        r = (q + s - p) % M
                        up_dest, up_sign = table[(p, q)][up_index]
                        down_dest, down_sign = table[(r, s)][down_index]
                        if up_dest < 0 or down_dest < 0:
                            continue

                        key = (p, q, r, s)
                        if key not in cache:
                            phase = phase_product_series(*key)
                            # Each orbital has weight 1/sqrt(2), so the onsite
                            # form factor is (1 + phase_product)/4.  The Fourier
                            # interaction contributes the additional factor 1/M.
                            cache[key] = (
                                -(1 + phase[0]) / 16,
                                -phase[1] / 16,
                                -phase[2] / 16,
                            )

                        row = up_dest * dim_spin + down_dest
                        sign = up_sign * down_sign
                        for order in range(3):
                            coefficients[order][row, column] += cache[key][order] * sign

    h0 = sp.simplify((coefficients[0] + coefficients[0].H) / 2)
    h1 = sp.simplify((coefficients[1] + coefficients[1].H) / 2)
    # coefficients[2] multiplies A^2, so the second derivative is 2*coefficients[2].
    h2 = sp.simplify(coefficients[2] + coefficients[2].H)
    return h0, h1, h2


def exact_kohn_curvature(n: int) -> tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, sp.Expr]:
    """Exact base-model energy, diamagnetic, paramagnetic, and total curvature."""
    h0, h1, h2 = base_hamiltonian_derivatives(n)
    dimension = h0.rows
    energy = -sp.Rational(n, 2)

    ground_vectors = (h0 - energy * sp.eye(dimension)).nullspace()
    if len(ground_vectors) != 1:
        raise RuntimeError(f"Expected a unique ground state for n={n}; got nullity {len(ground_vectors)}")

    ground = ground_vectors[0]
    ground /= sp.sqrt((ground.H * ground)[0])
    source = sp.simplify(h1 * ground)
    if sp.simplify((ground.H * source)[0]) != 0:
        raise RuntimeError("The first derivative of the ground-state energy did not vanish")

    # Solve (H-E)x = H'|g> with <g|x>=0.  Then the paramagnetic contribution is
    # -2 <g|H' x>.
    shifted = h0 - energy * sp.eye(dimension)
    augmented = shifted.row_join(ground)
    augmented = augmented.col_join(ground.H.row_join(sp.zeros(1, 1)))
    rhs = source.col_join(sp.zeros(1, 1))
    solution = next(iter(sp.linsolve((augmented, rhs))))
    response_vector = sp.Matrix(solution[:dimension])

    diamagnetic = sp.simplify((ground.H * h2 * ground)[0])
    paramagnetic = sp.simplify(-2 * (source.H * response_vector)[0])
    total = sp.simplify(diamagnetic + paramagnetic)

    exact_levels = sorted(h0.eigenvals(), key=lambda value: float(sp.N(value)))
    gap = sp.simplify(exact_levels[1] - exact_levels[0])
    return energy, diamagnetic, paramagnetic, total, gap


def main() -> None:
    print("M=4, two orbitals, U=1, b=pi/2, c=1/8")
    print("At k_m=m*pi/2: sin(4k_m)=0, so orbital weights are exactly (1/2, 1/2).")
    print("The parent h(k)=I-P(k) has exactly-flat eigenvalues (0,1).")
    print("Alias factor: f(A)=cos^2(2c sin 4A), with f''(0)=-2.")
    print()

    full_curvatures: dict[int, sp.Expr] = {}
    for n in (1, 2, 3):
        energy, dia, para, base_total, gap = exact_kohn_curvature(n)
        full_total = sp.simplify(base_total - 2 * energy)  # f''(0) E_n(0) = (-2)(-n/2)=n
        expected_base = sp.simplify(PI**2 * n * (M - n) / 24)
        assert sp.simplify(base_total - expected_base) == 0
        full_curvatures[n] = full_total
        print(
            f"n={n}: E(0)={energy}, gap={gap}, "
            f"dia={dia}, para={para}, base E''={base_total}, full E''={full_total}"
        )

    actual_e2 = full_curvatures[2]
    predicted_e2 = sp.simplify(sp.Rational(4, 3) * full_curvatures[1])
    raw_discrepancy = sp.simplify(actual_e2 - predicted_e2)

    # The proposal's convention is K_n=E_n''/(4V), with V=M=4 here.
    k1 = sp.simplify(full_curvatures[1] / 16)
    k2 = sp.simplify(full_curvatures[2] / 16)
    one_pair_coefficient = sp.simplify(k1 / (sp.Rational(1, 4) * sp.Rational(3, 4)))
    predicted_k2 = sp.simplify(sp.Rational(1, 2) ** 2 * one_pair_coefficient)
    k_discrepancy = sp.simplify(k2 - predicted_k2)

    print()
    print("Filling-law test at n=2:")
    print("  predicted raw E_2'' = (4/3) E_1'' =", predicted_e2)
    print("  actual raw E_2''                    =", actual_e2)
    print("  exact raw discrepancy               =", raw_discrepancy)
    print()
    print("In the proposal's K=E''/(4V) convention, V=4:")
    print("  K_1                   =", k1)
    print("  one-pair coefficient =", one_pair_coefficient)
    print("  predicted K_2         =", predicted_k2)
    print("  actual K_2            =", k2)
    print("  exact K discrepancy   =", k_discrepancy)


if __name__ == "__main__":
    main()
