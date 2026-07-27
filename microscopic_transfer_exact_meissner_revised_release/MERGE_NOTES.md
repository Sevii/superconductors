# Merge notes

This release replaces the two prior manuscripts with one 20-page microscopic-transfer paper.

## Structural merge

- **Part I** retains the lasting transfer framework from the Schrieffer-Wolff note: target-metric stability, the gap-independent current \(H_{-1}\) estimate, the anisotropic counterexample, the point-group isotropy theorem, and the finite-cluster two-sided row-ideal criterion.
- **Part II** replaces the prior note's isolated two-control overlap section with the complete degree-eight resolvent enumeration, closed source formula, two-pair all-filling obstruction, and explicit \(\mathcal C_8\) counterterm.
- **Part III** contains the exact dark branch, the positive-semidefinite finite-coupling completion, and the microscopic Meissner floor.

## Reconciliation of the two order-eight completion claims

The earlier Schrieffer-Wolff note gave an abstract completion: remove the unsafe cluster projection and replace the surviving target metric by a point-group average. That statement predated the complete microscopic source formula.

The merged paper uses the later and stronger result:

\[
\mathcal C_8=
\frac{16b^2s^4}{\bar\Delta_m^3\bar\Delta_s^2}\sum_eB_e^2
-\frac{b^4}{\bar\Delta_m^3}\mathcal M.
\]

This explicit finite-range counterterm cancels the full degree-eight source, not only the isolated bridge-overlap cluster. The old generic completion is therefore not retained as a separate theorem. Its useful locality and Peierls-covariance argument is incorporated as Proposition 10.2.

The merged claim is also stricter: \(\mathcal C_8\) proves the two-sided row-ideal/state-source condition through weighted degree eight. Full truncated Kohn-Abelian-Meissner equality additionally requires hydrodynamic closure of the current-created soft sector; point-group symmetry then supplies the scalar zero-momentum symbol. The paper does not assume that closure automatically follows from row-relative form.

## Claim-boundary changes

- Removed the obsolete possibility that unenumerated degree-eight clusters might cancel the isolated bridge overlap. The complete single-flavor source and Kramers correction are now included, and the two-pair witness proves the no-go.
- Added the analytic scope qualification: the nonzero degree-eight coefficient rules out all-composition invariance for every sufficiently small coupling, but does not formally exclude an isolated accidental finite-coupling zero involving higher orders.
- Preserved the exact dark branch as an eigenbranch of the original screened family, while keeping the ground-state claim restricted to the modified strictly overcompensated positive-semidefinite family.
- Kept the transport result as a strict finite-coupling floor. Equality of the total response remains conditional on the Part I regularity audit for the extra gadget rows.

## Expert-review revision

The accept-with-minor-revisions report was incorporated after the structural merge. The revision closes the missing commutator-square step in the two-pair witness proof, formalizes hydrodynamic closure by an explicit momentum-fiber intertwiner, adds an exact scalar-convolution example and a concrete audit protocol for operator-valued remainders, attaches all prior-work citations, defines the quasi-local interaction norm and response-density notation, removes the apparent translation-invariance restriction from the dark-kernel uniqueness proof, and clarifies the termwise gauge argument for Peierls positivity. The reviewer-safe headline and the later degree-eight claim boundary are unchanged.
