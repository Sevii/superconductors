# Build and verification report

## Paper build

- Engine: pdfLaTeX through `latexmk`.
- Output: 31 pages, US Letter.
- Undefined references: none.
- Undefined citations: none.
- Duplicate labels: none.
- Overfull boxes: none.
- Fatal LaTeX errors: none.
- Remaining layout messages: one harmless underfull-box notice.

## Source audit

- 139 labels; all unique.
- 91 reference occurrences covering 56 distinct targets; every target resolves.
- 25 citation occurrences covering 14 distinct keys; every cited key resolves.
- No occurrence of the superseded phrases “is diagonal in composition space,” “order-one off-manifold singular values,” or an unqualified “controlled avoided-crossing problem” for the generic unscreened family.
- The open-boundary corollary is referenced twice from the main text.
- The order-four isolation gate is referenced from the discussion.

## Exact leakage verifier

- Result: OVERALL PASS.
- Exact checks: 15/15.
- Volumes: \(L=2,\ldots,9\).
- Boundaries: periodic cycle and standard open path.
- Boundary/composition sectors: 760.
- Hard-core configuration evaluations: 699,040.
- Arithmetic: integers and `fractions.Fraction` for theorem-level decisions.

## PDF preflight and visual audit

- PDF opens successfully and is not encrypted or scan-only.
- 31 pages at 612 x 792 points.
- Fonts are embedded.
- The final PDF was rendered at 160 dpi and all pages were inspected in contact sheets, with direct inspection of the title/abstract, periodic leakage theorem, isolation gate, numerical-certificate section, open-path corollary, discussion, conclusion, and bibliography.
- No clipped text, overlapping equations, broken glyphs, black boxes, or margin overflows were observed.
