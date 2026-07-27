# Microscopic transfer of exact Meissner weight

This archive contains the merged microscopic-transfer manuscript and its reproducibility materials.

## Paper

`paper/microscopic_transfer_exact_meissner.pdf` is the 20-page expert-revised merged paper. Its narrative is:

1. infrared-stable target-metric and row-ideal transfer theory;
2. complete weighted-degree-eight audit and no-go for the unmodified all-composition family;
3. explicit finite-order counterterm plus an exact dark-branch positive completion with a finite-coupling Meissner floor.

The reviewer-safe claim boundary is in `CLAIM_STATUS.md`. `REVIEW_RESPONSE.md` records the response to the expert report, while `MERGE_NOTES.md` records how the two prior manuscripts and their differing order-eight claims were reconciled.

## Reproduce the certificates and PDF

From the release root:

```bash
chmod +x RUN_ALL.sh
./RUN_ALL.sh
```

The scripts require Python 3, NumPy, and SciPy. The paper requires a standard LaTeX installation with `pdflatex` or `latexmk`.

## Directory map

- `paper/`: final LaTeX source, PDF, PDF metadata, and preflight report.
- `scripts/`: eight deterministic scientific verification programs plus one release-hygiene checker.
- `certificates/`: retained fresh outputs from every verification program.
- `data/`: machine-readable all-filling, tuning, and Kramers audit data.
- `context/`: the two superseded source manuscripts and their original claim-status files for provenance.
- `CLAIM_STATUS.md`: exact established/not-established boundary.
- `FINDINGS.md`: consolidated scientific decision and next research tasks.
- `MERGE_NOTES.md`: editorial and theorem-level reconciliation performed in the merge.
- `REVIEW_RESPONSE.md`: itemized incorporation of the expert reviewer comments.
- `review/`: the expert report included for provenance.

All retained scientific and release-hygiene certificate outputs end in `PASS`, `OVERALL: PASS`, or `AUDIT IDENTITIES: PASS`.
