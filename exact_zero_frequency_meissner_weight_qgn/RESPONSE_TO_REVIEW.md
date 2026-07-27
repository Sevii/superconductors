# Response to the reviewer

The reviewer found the central mathematics correct and recommended minor revision. This document records the disposition of every substantive comment.

## 1. Literature positioning

**Implemented.** The introduction now positions the exact response theorem relative to:

- Peotta–Törmä and the flat-band quantum-geometric superfluid-weight literature;
- the Törmä–Peotta–Bernevig review;
- Tada–Koma’s order-of-limits/no-go analysis;
- Sewell and Nieh–Su–Zhao on ODLRO and the Meissner effect under additional hypotheses; and
- Yang’s eta-pairing construction.

The text now states explicitly that the present result evaluates the finite-size current kernel directly rather than inferring transport from ODLRO. The DOI for Gao–Han–Khalaf has also been added.

## 2. Formal finite-q definition

**Implemented.** Section 2 now defines `J_L[a]`, `T_L[a]`, `C_L[a]`, and `K_L(i zeta; a)` for an arbitrary real bond profile before any model specialization. It then defines the finite-`q` transverse objects and fixes the convention

```text
K_{L,T}(q,0;v) := C_L[a_{L,q,v}].
```

Thus `zeta=0` always means the static finite-volume curvature. The paper also defines admissible `q -> 0` sequences on the volume-dependent reciprocal grids and distinguishes the thermodynamic-first finite-`q` frequency order.

## 3. Space-time summability proof

**Expanded.** Proposition 6.3 now includes:

- the connected subtraction through `Q_L`;
- the exact time/spectral representation;
- the uniform low-energy tail estimate
  `1/E <= e integral_{1/delta}^infinity exp(-tE) dt` for `E <= delta`;
- uniform `H^{-1}` equi-integrability from the envelope; and
- Fourier-space dominated convergence for the static paramagnetic term, followed by the assumed diamagnetic continuity.

The proposition remains a proved sufficient criterion rather than a proof sketch.

## 4. Multiblock cross term

**Corrected.** The proof now emphasizes that `s_0=(S_sw,0)` lies in `ker L`, so `f_zeta(L)s_0=s_0`. Only after using self-adjointness and this eigenvector identity does direct-sum orthogonality eliminate the cross term. No block-diagonal assumption on `L=DD^†/2` is made.

## 5. Clean-gas temperature bookkeeping

**Corrected.** The proof first uses a smooth finite-beta Fermi occupation, takes `q -> 0` by dominated convergence and integration by parts, and then sends `beta -> infinity` at fixed density. The paper now states that the numerical verifier checks this smooth regularization and that the proposition is its zero-temperature limit.

## 6. Verification gaps

**Closed and strengthened.** A new `independent_fock_check.py` constructs the spinful-fermion Fock space on a path, triangle, and four-site cycle. It verifies complete fermionic site swaps, twisted involutions, even/odd kernel statements with fermionic signs, full-Fock gauge covariance, charge-two seniority-zero compression, analytic current/stress derivatives against five-point finite differences, and full-Fock Kubo/Hodge response.

The graph verifier now also includes:

- a half-filled six-pair check on a `4 x 3` torus, testing the genuine multipair content of Aldous’ gap theorem; and
- an operator-level backbone-floor test with nonzero swap/residual target-Laplacian mixing.

All three certificates end with `OVERALL: PASS`.

## Minor comments

Also implemented:

- explicit statement that `E_L=0` and `Hhat_L=H_L(0)` for positive-square parents;
- the parenthetical term “uniform H^{-1} equi-integrability”;
- the exact spin-1/2 ferromagnet representation;
- the observation that `U>0` only selects the paired kernel and drops out of response formulas;
- a prominent graph-filter remark;
- explicit notice that the exact transverse current annihilates the ground state, making the response purely diamagnetic and non-dispersive at finite size; and
- updated reproducibility and claim-status documents.
