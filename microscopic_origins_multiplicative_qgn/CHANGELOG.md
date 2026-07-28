# Changelog

## Leakage-integrated revision — July 28, 2026

### Exact periodic leakage theorem

The paper now defines
\[
R_{{\rm ch},L}^{(n)}=(1-P_{\mathcal Z_n})H_{\rm ch}P_{\mathcal Z_n}
\]
and gives its exact action on each composition state,
\[
R_{{\rm ch},L}^{(n)}|c\rangle
=A_c\bigl(D_L-\langle D_L\rangle_c\bigr)|c\rangle.
\]
Normalized centered images provide the full singular-value decomposition. The revision states the exact Gram operator, general norm, rank, closed variance formula, complete periodic zero set, and the benchmark
\[
\|R_{4,4}^{(2)}|2,0\rangle\|=\frac{8\sqrt2}{3}A_4.
\]

### Open-boundary theorem

A new appendix corollary treats the standard open path with `L-1` bonds. It derives the endpoint decomposition
\[
H_{{\rm ch},L}^{\rm op}=(8A_c-2U_c)n+8A_cZ-V_cB,
\]
the open compressed selector, the covariance formula for the singular values, the exact norm, and the complete physical-shell zero set. The paper now explicitly distinguishes the periodic theorem-level lattice from the open truncation used by finite-system DMRG.

### Isolation gate

A new remark introduces
\[
\Xi_L^{(n)}(\lambda)
=\frac{\lambda^4\|R_{{\rm ch},L}^{(n)}\|}
       {\lambda^6\gamma_L^{(n)}}.
\]
It makes explicit that nonzero order-four leakage is not suppressed relative to an order-six isolation gap by taking \(\lambda\to0\). An autonomous composition Hamiltonian therefore needs order-four protection or exact cancellation.

### Corrected claim discipline

- Replaced every implication that the unscreened charge operator itself is “diagonal in composition space” by the precise compression statement.
- Recast the generic unscreened family as a charge-selection problem rather than an automatically isolated selector-plus-mixing model.
- Corrected the one-pair discussion: both the periodic selector and periodic leakage vanish at one pair, so the order-six one-pair curvature formula is unchanged.
- Distinguished three valid uses of the unscreened result: projected selector algebra, the true charge-selected spectrum, or a modified model with order-four isolation.
- Linked the off-manifold block to the “state source” vocabulary of the companion microscopic-transfer analysis without importing any theorem from that paper.

### Reviewer-polish changes

- Replaced “order-one off-manifold singular values” by the exact scale \(A_c\sigma_L(c)\), of the same shell order as the selector shift.
- Added two main-text references to the open-path corollary.
- Clarified that exact enumeration checks centered moments, while the operator SVD follows analytically from occupation-basis diagonality and block-charge orthogonality.
- Extended the default exact scan through \(L=9\).
- Clarified that tabulated decimal open-path singular values already include \(A_4\).

### Verification and editorial changes

- The exact leakage verifier now reports 15/15 PASS over 760 sectors and 699,040 configurations.
- Updated the abstract, claim hierarchy, main-results list, numerical-certificates section, discussion, conclusion, and bibliography.
- Renamed the canonical paper files without reviewer-round or draft suffixes.
- Made `lmodern` optional for portable compilation.
