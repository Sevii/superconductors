# Deferred Mathematical Note: Rank-One Reflection Gap Ordering

**Status:** outline only; not part of the physics proof.

## Proposed title

**Aldous-Type Gap Ordering for Rank-One Householder Reflection Laplacians**

## Problem

For a connected weighted Parseval frame of rank-one projectors \(q_x\) on \(V\simeq\mathbb C^L\), define Householder reflections \(r_x=I-2q_x\) and representation Laplacians

\[
\mathcal L_k=\frac{U}{4}\sum_x w_x^2[I-\rho_k(r_x)]
\]

on the irreducible components \(\mathcal W_k\subset\operatorname{End}(\Lambda^nV)\), with \(\mathcal W_1\) the adjoint sector. Determine whether

\[
\lambda_{\min}(\mathcal L_k)\ge \lambda_{\min}(\mathcal L_1)
\qquad (k\ge2)
\]

for every connected frame.

If true, it implies \(\Delta_{L,n}=\Delta_{L,1}\) for every nonsingular filling.

## What is already proved

- Exact decomposition:
  \[
  \Delta_{L,n}=\min_{1\le k\le\min(n,L-n)}\gamma_k.
  \]
- Therefore:
  \[
  \Delta_{L,n}\le\Delta_{L,1}.
  \]
- Extensive small-\(L\) numerics found equality for rank-one projected-density models.
- Equality fails for generic Hermitian one-body square factors, so rank one is essential.
- None of this is needed for the stiffness theorem.

## Related reference

G. Alon and D. Puder, “Aldous-type Spectral Gaps in Unitary Groups,” arXiv:2603.00353 (2026), formulates unitary-group analogues and proves several cases for hypergraph-induced measures. The present Householder-reflection Laplacian is related but has not been shown to fall under their proved families.

## Candidate proof routes

1. Specialized octopus inequality for rank-one reflections.
2. A completely positive intertwiner from \(\mathcal W_k\) to \(\mathcal W_1\).
3. Frame-operator comparison exploiting \(\sum_xw_xq_x=I\).
4. Highest-weight/Casimir estimates adapted to Householder generators.
5. Reduction to a hypergraph-induced unitary measure, if an exact identification exists.

## Publication boundary

This note should be developed independently. The physics paper uses only the exact finite-size curvature reduction, H4 connectivity theorem, and the proved inequality \(\Delta_{L,n}\le\Delta_{L,1}\).
