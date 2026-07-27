# Build and preflight record

## LaTeX build

- Engine: `pdfTeX` through `pdflatex`.
- Passes: two successful passes with `-interaction=nonstopmode -halt-on-error`.
- Output: `exact_zero_frequency_meissner_weight_qgn.pdf`.
- Length: 19 pages, US Letter.
- PDF metadata: title, author, and subject embedded.
- Final-log scan: no LaTeX errors, undefined references, overfull boxes, underfull boxes, or package warnings.
- Source audit: 108 unique labels; no duplicate labels or missing internal references.

## Numerical certificates

- `verify_dynamical_response.py`: `OVERALL: PASS`.
- `verify_graph_hodge_meissner.py`: `OVERALL: PASS`.
- `independent_fock_check.py`: `OVERALL: PASS`.

The reviewer-requested verification additions are present:

1. a full-Fock spinful-fermion implementation on `P3`, `C3`, and `C4`;
2. a genuine multipair Aldous-gap check in the half-filled six-pair sector of a `4 x 3` torus (dimension 924); and
3. an operator-level multiblock backbone-floor check with nonzero swap/residual target-Laplacian mixing.

The scripts support reproducibility but do not replace review of the analytic proofs.

## PDF preflight

- PDF opens successfully and has 19 pages.
- No encryption, JavaScript, XFA, or scanned-page dependency.
- Text is extractable.
- All reported fonts are embedded and subsetted.
- Extracted-text scan found no unresolved `??`, broken control sequences, or missing-reference markers.

## Render and visual inspection

The original 17-page draft and the reviewer-revised 19-page PDF were both rendered. The revised PDF was rendered again after the final metadata build at 160 dpi. All 19 final pages were visually inspected for clipping, overlap, broken glyphs, malformed equations, table overflow, and reference-layout defects. No visible problems were found. Intermediate render and diff files are omitted from the release archive.
