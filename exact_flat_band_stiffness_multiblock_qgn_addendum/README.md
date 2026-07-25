# Multi-block QGN addendum — review-integrated draft

This package contains the review-integrated draft of **Multi-Block Reduction and Irreducibility Laws for QGN Superconductors**, a companion to *Exact Flat-Band Stiffness from Pair Mobility in a QGN Model*.

## Primary files

- `multiblock_qgn_addendum_reviewed.pdf` — typeset 25-page draft.
- `multiblock_qgn_addendum_reviewed.tex` — authoritative LaTeX source.
- `multiblock_qgn_addendum_reviewed.md` — repository-oriented Markdown conversion.
- `multiblock_qgn_addendum_reviewed_CHANGELOG.md` — mapping from the external review to the applied changes.

## Verification files

- `verify_multiblock_reduction.py`
- `multiblock_reduction_verification_reviewed.txt`
- `verify_multiblock_connectivity_reviewed.py`
- `multiblock_connectivity_verification_reviewed.txt`

The connectivity script is independent of the manuscript formulas at the exact-diagonalization level. It constructs the one-body generators, closes their associative and Lie algebras, builds the many-body positive-square Hamiltonians, resolves every spin/block charge sector, and compares the observed nullities with the complete kernel theorem. The two new adversarial models test the odd-parity clause and a synchronization class embedded in a larger center-resonant system.

## Supporting files

- `multiblock_connectivity_open_problem.md` — historical problem note with a resolved-status header.
- `multiblock_connectivity_solution.md` — detailed solution note underlying Section 4.
- `multiblock_reduction_law.md` — detailed source-reduction note underlying Sections 5–9.
- `claude_review_multiblock_v2_2026-07-24.md` — the review applied in this revision.

## Build

From the package directory:

```bash
latexmk -pdf -pdflatex='pdflatex %O %S' \
  -interaction=nonstopmode -halt-on-error \
  multiblock_qgn_addendum_reviewed.tex
```

The retained build was clean: 25 pages, no LaTeX warnings, no undefined references, and no overfull or underfull boxes.

## Re-run the certificates

```bash
python verify_multiblock_reduction.py
python verify_multiblock_connectivity_reviewed.py
```

Expected final lines:

```text
PASS
PASS: all six multi-block connectivity and kernel certificates reproduced.
```

The retained numerical outputs are included for exact comparison.
