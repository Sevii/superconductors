# Microscopic Origins of Multiplicative QGN Interactions

## Leakage-integrated revision — July 28, 2026

This release contains the revised paper

**Microscopic Origins of Multiplicative Null Interactions in Quantum-Geometric-Nesting Superconductors**

and the exact periodic/open-boundary unscreened leakage certificate supporting the new theorem statements.

## Main files

- `paper/microscopic_origins_multiplicative_qgn.pdf` — compiled 31-page paper.
- `paper/microscopic_origins_multiplicative_qgn.tex` — canonical LaTeX source.
- `paper/microscopic_origins_multiplicative_qgn_from_draft4.patch` — unified diff against the uploaded fourth draft.
- `CHANGELOG.md` — scientific and editorial changes.
- `REVIEW_RESPONSE.md` — point-by-point response to the leakage-certificate review.
- `BUILD_REPORT.md` — compilation, source-audit, exact-verifier, and render checks.
- `supporting_certificate/` — derivation, exact verifier, certificate output, CSV tables, focused Phase-3 matrices, and reproducibility metadata.

## Central correction

The periodic order-four charge operator has an exactly diagonal **compression** to the product-AGP composition space, but it does not generally preserve that space. The revised paper gives its exact composition-resolved action, singular-value decomposition, Gram operator, norm, rank, and complete periodic zero set. A separate open-path corollary gives the endpoint-corrected selector and leakage formulas required by open-chain DMRG calculations.

The generic unscreened family is therefore presented as an order-four charge-selection problem. Its compression contains the selector, but a closed selector-plus-mixing composition Hamiltonian requires projection, exact cancellation of the centered charge action, or an additional order-four-or-stronger isolation mechanism.

The new isolation diagnostic is

\[
\Xi_L^{(n)}(\lambda)
=\frac{\lambda^4\|R_{{\rm ch},L}^{(n)}\|}
       {\lambda^6\gamma_L^{(n)}}
=\frac{A_c\max\sigma_L}{\lambda^2\gamma_L^{(n)}}.
\]

## Reproduce the leakage certificate

```bash
cd supporting_certificate
bash run_certificate.sh
```

The default exact scan covers both boundaries and all compositions for `L=2,...,9`: 760 boundary/composition sectors and 699,040 hard-core configurations. All 15 theorem-level checks pass using integer and `fractions.Fraction` arithmetic.

## Build the paper

```bash
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error \
  microscopic_origins_multiplicative_qgn.tex
```

The source conditionally loads Latin Modern, so it remains portable to TeX installations where `lmodern` is unavailable.
