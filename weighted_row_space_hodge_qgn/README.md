# Weighted row-space Hodge theorem package

This reviewer-integrated package develops O1 of the variational QGN functional: replacing the residual-current many-body pseudoinverse by a proof-carrying target-row Hodge problem.

## Main finite-volume result

For

`H_s = D^† G D / 2`

and a current source satisfying the full range condition `P_ker(D) J_R P = 0`, the weighted Hodge representative is

`eta_* P = G D (D^† G D)^+ J_R P`.

It obeys the operator identity

`P eta^† G^{-1} eta P = (1/2) P J_R H_s^+ J_R P + P kappa^† G^{-1} kappa P`

for every feasible witness `eta = eta_* + kappa`. Hence the Schur-complement current weight is the exact Loewner minimum and the operator-norm minimization has no duality gap.

The revised note corrects an earlier overclaim: the projected witness is the unique **Loewner-order** minimizer, but for `rank(P) >= 2` a nonprojected witness can have the same operator norm. The exact norm-equality condition is

`C_kappa <= ||C_*|| P - C_*`,

where `C_* = P eta_*^† G^{-1} eta_* P` and `C_kappa = P kappa^† G^{-1} kappa P`.

## Locality and decision procedures

The package now distinguishes three regimes.

1. A buffered cluster theorem gives an exponentially local witness when `p_X P = P`, `p_X j_X p_X = 0`, and finite-cluster pseudoinverses grow polynomially. The decay bound explicitly includes the parent-row support buffer.
2. Finite-range translation-invariant factorization is a Laurent-module membership problem. Smith normal form decides the one-dimensional case; Gröbner-basis saturation and syzygy computation give a terminating exact decision procedure in every dimension over an effective coefficient field.
3. At genuine rank drops, a channel-resolved theorem gives the precise vanishing-order condition for bounded Hodge cost. The note does **not** claim a general exponentially local, non-polynomial soft-mode factorization theorem.

## Files

- `weighted_row_space_hodge_qgn.pdf` - reviewer-integrated research note.
- `weighted_row_space_hodge_qgn.tex` - editable LaTeX source.
- `REVIEW_RESPONSE.md` - point-by-point response to the reviewer comments.
- `CLAIM_STATUS.md` - exact claims, conditional claims, and theorem boundary.
- `O1_IMPLEMENTATION_SPEC.md` - cluster, Laurent-module, and cost-certificate workflow.
- `verify_weighted_row_hodge.py` - deterministic finite-matrix checker.
- `verify_weighted_row_hodge_exact.py` - exact rational checks, including norm-minimizer nonuniqueness.
- `VERIFICATION_OUTPUT.txt` and `EXACT_VERIFICATION_OUTPUT.txt` - verifier outputs.
- `o1_certificate_schema.json` and `example_o1_certificate.json` - machine-readable witness schema and illustrative instance.
- `BUILD_INFO.txt`, `PDF_PREFLIGHT.txt`, `MANIFEST.txt`, and `SHA256SUMS.txt` - build, inspection, inventory, and integrity records.

## Recommended next use

Run the complete weighted-degree-eight microscopic current audit. For every connected cluster `X`, verify the load-bearing frustration-free condition `p_X P = P` and test `p_X j_X p_X`. Safe terms yield `eta_X=(D_X^+)^†j_X`; their buffered local cost can be inserted into the existing QGN stiffness floor. Translation-invariant orbit sums should also be passed through the Laurent-module decision when exact Fourier data are available.

## Verification

```bash
python verify_weighted_row_hodge.py
python verify_weighted_row_hodge_exact.py
```

The scripts check range factorization, the Schur-complement identity, target-space projection, the operator Pythagorean decomposition, the approximate-witness inequality, selected-branch cokernel mixing, and the explicit nonprojected operator-norm minimizer.
