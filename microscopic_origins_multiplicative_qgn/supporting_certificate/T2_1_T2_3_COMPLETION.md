# Completion record for T2.1 and T2.3

## T2.1 — complete

The periodic proposition now contains:

1. the exact composition-resolved action
   \[
   R|c\rangle=(\mathcal C-\langle\mathcal C\rangle_c)|c\rangle;
   \]
2. normalized left singular vectors obtained from the centered charge observable;
3. the exact SVD
   \[
   R=\sum_c\tau_c|\ell_c\rangle\langle c|;
   \]
4. the diagonal Gram operator
   \[
   R^\dagger R=\sum_c\tau_c^2|c\rangle\langle c|;
   \]
5. the general norm and rank formulas
   \[
   \|R\|=\max_c\tau_c,
   \qquad
   \operatorname{rank}R=\#\{c:\tau_c>0\}.
   \]

For the periodic order-four shell, \(\tau_c=A_4\sigma_L(c)\), with the closed rational formula for \(\sigma_L^2(c)\) stated in `paper_insert.tex`.

## T2.3 — complete after scope split

- **Periodic cycle:** complete across all fillings, including the exact projected selector, off-manifold singular values, general norm, and full zero set.
- **Standard open path:** complete across all fillings for general physical shell parameters \(A_4>0,V_4>0\), including the endpoint correction, compressed selector, exact singular values, norm, and zero set.

The result applies to the standard numerical truncation with \(L\) active sites and \(L-1\) bond terms. A different boundary implementation—for example explicit edge counterterms chosen to restore degree two—would define a different operator and would require a separate short certificate.

## Verification status

The exact verifier checks both boundary conditions for every composition and every hard-core pair configuration through \(L=9\):

- 760 boundary/composition sectors;
- 699,040 configuration evaluations;
- 15 exact PASS checks under the default shell parameters;
- additional PASS runs at nondefault values of \(V_4/A_4\).

No floating-point threshold is used in any theorem-level PASS/FAIL decision.
