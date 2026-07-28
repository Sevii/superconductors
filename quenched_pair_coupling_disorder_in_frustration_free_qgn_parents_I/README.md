# Reviewer package: Quenched Pair-Coupling Disorder in Frustration-Free QGN Parents I

**Subtitle:** Exact weighted Hodge response, random-conductance homogenization, and annealed dynamical bounds

This directory is a self-contained reviewer package for Paper I.

## Manuscript

- `quenched_pair_coupling_disorder_in_frustration_free_qgn_parents_I.pdf` — submission-ready PDF
- `quenched_pair_coupling_disorder_in_frustration_free_qgn_parents_I.tex` — complete standalone LaTeX source with inline bibliography

Build the manuscript from this directory with:

```bash
pdflatex -interaction=nonstopmode -halt-on-error quenched_pair_coupling_disorder_in_frustration_free_qgn_parents_I.tex
pdflatex -interaction=nonstopmode -halt-on-error quenched_pair_coupling_disorder_in_frustration_free_qgn_parents_I.tex
pdflatex -interaction=nonstopmode -halt-on-error quenched_pair_coupling_disorder_in_frustration_free_qgn_parents_I.tex
```

No files outside this directory are required to compile the paper.

## Paper-specific support

- `CLAIM_STATUS.md`
- `STOCHASTIC_HOMOGENIZATION_SOURCE_MAP.md`
- `VERIFICATION_OUTPUT.txt`
- `FINITE_TORUS_VERIFICATION_OUTPUT.txt`
- Reproducibility scripts, CSV data, figures, and JSON certificates used for the numerical and algebraic audits

To rerun the principal checks:

```bash
python verify_weighted_qgn_hodge.py --outdir . --seed 20260727 --samples 24
python verify_finite_torus_source_tail.py --outdir .
```

The random-torus Monte Carlo diagnostic uses reproducible independent seed sequences. Runtime versions are recorded in `weighted_qgn_hodge_certificate.json`.

## Shared reviewer record

The package also includes the review response, original review, changelog, paper-split map, and source/build audit records. The audit records retain the paths and filenames from the combined working archive for provenance.

`MANIFEST.txt` lists every delivered file and `CHECKSUMS.sha256` provides a SHA-256 integrity check.
