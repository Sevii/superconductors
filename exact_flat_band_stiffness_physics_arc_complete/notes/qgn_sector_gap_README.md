# QGN Sector-Gap Reduction Package

**Date:** July 23, 2026  
**Author:** Nicholas Sledgianowski  
**Computational assistance:** ChatGPT Sol 5.6

## Purpose

This package follows the H4 connectivity theorem by investigating the observed filling-independence of the first positive fixed-number sector gap.

The main exact result is the representation-theoretic reduction

\[
\Delta_{L,n}=\min_{1\le k\le \min(n,L-n)}\gamma_k,
\qquad
\Delta_{L,1}=\gamma_1,
\]

where \(\gamma_k\) is the lowest eigenvalue of a common rank-one reflection Laplacian in the irreducible \(U(L)\) sector \(\mathcal W_k\). Therefore the empirical identity \(\Delta_{L,n}=\Delta_{L,1}\) is equivalent to the still-unproved ordering \(\gamma_k\ge\gamma_1\) for all \(k\ge2\), called the **rank-one Aldous property** in the note.

The package also derives an exact one-pair orbital-Gram formula for translation-invariant UPC models and shows that the complete finite-size sector gap cannot have a positive thermodynamic lower bound. For every smooth fixed local projector, the one-pair gap has an inverse-square upper bound in the smallest linear dimension. In GHK Model II, the asymptotic formula is

\[
\Delta_{N,1}\sim \frac{U\pi^2\xi^2}{4N^2}.
\]

## Scientific status

- **Proved:** reflection-Laplacian identity; multiplicity-free sector reduction; \(\Delta_{L,n}\le\Delta_{L,1}\); one-pair Gram formula; inverse-square upper bound; Model-II asymptotic.
- **Strongly supported, not proved:** rank-one Aldous ordering \(\gamma_k\ge\gamma_1\), equivalently filling-independent sector gap.
- **Ruled out:** a size-independent lower bound on the complete fixed-number sector gap for smooth translation-invariant thermodynamic sequences.
- **Unaffected:** finite-size H4 simplicity and analyticity, and the exact curvature-reduction theorem at each finite size.

## Files

### Main note

- `qgn_sector_gap_reduction_note.md` — editable manuscript source.
- `qgn_sector_gap_reduction_note.pdf` — publication-style PDF.
- `qgn_sector_gap_reduction_note.docx` — editable Word version.

### Main theorem amendment

- `restricted_qgn_reduction_theorem_h4_gap_amended.md` — theorem manuscript with the gap correction integrated.
- `restricted_qgn_reduction_theorem_h4_gap_amended.pdf` — amended PDF.
- `restricted_qgn_reduction_theorem_h4_gap_amended_fixed_layout.docx` — fixed-layout Word view copy; pages are embedded to preserve the verified mathematical layout.

### Reproducibility

- `qgn_sector_gap_certificate.py` — deterministic certificate and numerical tests.
- `qgn_sector_gap_rank_one_tests.csv` — 160 random connected Parseval-frame tests.
- `qgn_sector_gap_gram_validation.csv` — direct dense/Gram comparisons for GHK Models I and II.
- `qgn_sector_gap_scaling.csv` — finite-size one-pair gap data.
- `qgn_sector_gap_generic_hermitian_control.json` — connected non-rank-one control showing gap equality can fail.
- `qgn_sector_gap_certificate_summary.json` — headline test summary.
- `qgn_sector_gap_scaled.png` — scaled-gap figure.
- `qgn_sector_gap_loglog.png` — log-log gap-closing figure.
- `TEST_RESULTS.txt` — archived output from a clean rerun.
- `SHA256SUMS.txt` — hashes of package files.

## Reproduction

From a Python environment with NumPy and SciPy installed:

```bash
python qgn_sector_gap_certificate.py
```

Expected final lines:

```text
rank-one tests: retained=160, worst relative deviation=1.019e-13
GHK Gram/dense validations: PASS
generic Hermitian control: gap2/gap1= 0.432639
PASS
```

The script rewrites the CSV/JSON outputs in its working output directory. It uses deterministic random seeds.

## Interpretation

The expert suggestion that the all-filling sector-gap question should reduce to a one-pair calculation was directionally correct, but the one-pair calculation reveals a closing phase/pseudospin-wave branch rather than a uniform thermodynamic gap. The next mathematical target is the rank-one Aldous ordering. The next physics targets are pair-mass nonvanishing and a gapless order-of-limits argument relating canonical flat-connection curvature to physical zero-frequency superfluid stiffness.
