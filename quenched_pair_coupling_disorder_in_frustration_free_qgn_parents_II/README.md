# Reviewer package: Quenched Pair-Coupling Disorder in Frustration-Free QGN Parents II

**Subtitle:** Quenched dynamical-source fluctuations, sharp cutoffs, and disorder-noise geometry

This directory is a self-contained reviewer package for Paper II.

## Manuscript

- `quenched_pair_coupling_disorder_in_frustration_free_qgn_parents_II.pdf` — submission-ready PDF
- `quenched_pair_coupling_disorder_in_frustration_free_qgn_parents_II.tex` — complete standalone LaTeX source with inline bibliography

Build the manuscript from this directory with:

```bash
pdflatex -interaction=nonstopmode -halt-on-error quenched_pair_coupling_disorder_in_frustration_free_qgn_parents_II.tex
pdflatex -interaction=nonstopmode -halt-on-error quenched_pair_coupling_disorder_in_frustration_free_qgn_parents_II.tex
pdflatex -interaction=nonstopmode -halt-on-error quenched_pair_coupling_disorder_in_frustration_free_qgn_parents_II.tex
```

No files outside this directory are required to compile the paper.

## Paper-specific support

- `CLAIM_STATUS.md`
- `QUENCHED_FLUCTUATION_SOURCE_MAP.md`
- Human-readable verification outputs
- Reproducibility scripts, CSV data, figures, and JSON certificates used for the fluctuation, cutoff, coarea, and boundary-layer audits

To rerun representative exact checks:

```bash
python verify_abel_tail_identities.py --output-dir .
python verify_maximal_shell_coarea.py --quick
python verify_disorder_noise_boundary_layer.py --quick --outdir .
```

This paper is deliberately a partial-results and theorem-boundary paper. It does not claim the still-open mesoscopic nonlinear sharp-cutoff variance theorem.

## Shared reviewer record

The package also includes the review response, original review, changelog, paper-split map, and source/build audit records. The audit records retain the paths and filenames from the combined working archive for provenance.

`MANIFEST.txt` lists every delivered file and `CHECKSUMS.sha256` provides a SHA-256 integrity check.
