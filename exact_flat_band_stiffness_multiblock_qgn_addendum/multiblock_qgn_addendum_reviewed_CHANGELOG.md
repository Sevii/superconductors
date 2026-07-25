# Review-integrated addendum changelog

Date: July 24, 2026

This draft applies the recommendations in `claude_review_multiblock_v2_2026-07-24.md` to the v2 addendum.

## Manuscript changes

1. Promoted the central conceptual correction to the abstract and introduction: second quantization `dΓ` is a Lie-algebra representation, not an associative-algebra homomorphism, so associative fullness does not by itself generate relative block-charge constraints.
2. Closed the implicit active-label point in the joint Lie-algebra proof. A synchronized identification cannot send a nonzero centered rank-one projector to zero, centered spectra force matching traces, and the implementing unitary therefore intertwines the full labeled projectors.
3. Renamed Theorem 4.4 as the exact criterion for fixed-composition uniqueness (MB3), defined “the MB3 criterion,” and removed ambiguous later uses of “under MB3.”
4. Added an explicit paragraph on one-dimensional blocks: they carry no semisimple factor, are harmless for fixed-composition uniqueness, and are controlled entirely by the center test.
5. Added the odd-capacity resonance discussion. Odd `F(δ)` resonances leave every even pairing sector intact but create exact fermionic determinant zero modes in odd particle-number sectors.
6. Expanded the reproducibility section from four even-sector examples to six sector-resolved examples, including:
   - `L=(1,2)`, `Z=span(2,1)`, testing an odd-capacity resonance;
   - `L=(1,2,2)` with synchronized `{B,C}` and simultaneous even/odd center resonances.
7. Added direct Lie-span membership tests for the available total center and the missing relative center, plus an explicit factor-residual check of the two determinant dark states.
8. Marked `multiblock_connectivity_open_problem.md` as resolved while recording which half of the original conjecture survived and which half was replaced.

## Certificate status

- Source-reduction audit: maximum identity error `3.553e-14`; minimal interference spectrum `{0,4}`.
- Connectivity/kernel audit: six models, exact agreement with Theorem 4.3 in every charge sector; total-center residual `1.592e-15`, missing-relative-center residual `2.000`, determinant-state residual `1.110e-16`.
- LaTeX: clean 25-page build with no warnings, undefined references, overfull boxes, or underfull boxes.
- PDF: rendered with PDFium and visually inspected page by page; no clipping, overlap, or broken glyphs found.
