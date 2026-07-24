# Independent rerun of all certificates against the audited paper (July 23, 2026)

All nine certificate scripts in `certificates/` were rerun from scratch in a clean
Python 3.11.15 environment (numpy 2.4.4, scipy 1.17.1, sympy 1.14.0) and their
outputs were compared line by line against the archived outputs in
`certificates/data/` and against the numerical claims in the audited paper
(`paper/exact_flat_band_stiffness_physics_arc.md`, revised July 23, 2026).
Raw rerun logs are in `certificates/data/rerun_2026-07-23/`.

## Paper-version check

Applying `revision_history/exact_flat_band_stiffness_physics_arc_corrections.diff`
to the pre-audit manuscript reproduces the audited manuscript **byte for byte**,
and the changes correspond one-to-one to the nine items in
`revision_history/exact_flat_band_stiffness_revision_notes.md`. The corrected PDF
contains the corrected text (spot-checked on the T3 lifted-mass statement and the
four-cell certificate parameters). The audited version is therefore the latest and
is shipped as the canonical `paper/` copy.

## Certificate results versus paper claims

| Certificate | Paper / archive claim | Rerun result | Status |
|---|---|---|---|
| Four-cell winding (`flat_band_counterexample_exact.py`) | raw defect 2/3, normalized 1/24 at g=1 (exact SymPy) | raw 2/3, normalized 1/24, exact | MATCH |
| Resonance controls (`flat_band_resonance_phase12.py`) | M=8 harmonic-4 control obeys filling law; resonant defects match analytic values | control defect 9.9e-10 (FD floor); resonant M=4 defects 2/3, 2 and M=8 defect 8/7 reproduce analytic predictions | MATCH |
| Reduction theorem (`restricted_qgn_reduction_certificate.py`) | max tensor-identity error 1.25e-18; max finite-difference error 5.17e-08 | 1.2500e-18; 4.51e-08 (FD floor) | MATCH |
| H4 connectivity (`qgn_h4_connectivity_certificate.py`) | full Lie dimensions 16 and 25; exact connected/disconnected nullity counts; all 55 GHK graphs connected; max Parseval-frame error 8.10e-16 | dimensions 16/25; all nullity counts exact; 55/55 connected; frame error 8.100e-16; output CSVs identical to archived copies | MATCH (exact) |
| Model-II closure (`model_II_pair_dispersion_certificate.py`) | direct lifted pair matrix 4.26e-16; twist identity 5.88e-16; finite-grid mass 3.56e-09; many-body curvature 2.50e-09 | 4.259e-16; 5.881e-16; 3.558e-09; curvature 5.10e-09 (FD floor) | MATCH |
| Connectivity counterexample (`qgn_connectivity_counterexample.py`) | half-filled even-torus zero mode exact; positive one-pair inverse mass | zero-mode residuals 0.0 at all tested twists; m_pair^{-1}=1; CSV agrees with archive to the FD floor (<3e-09) | MATCH |
| Sector-gap evidence (`qgn_sector_gap_certificate.py`) | 160 rank-one cases, worst relative deviation 1.02e-13; non-rank-one control ratio 0.432639 | 160 cases, 1.019e-13; 0.432639 | MATCH (exact) |
| GHK models ED (`ghk_models_ed_search.py`) | independent projected-Hubbard ED, stiffness ratios 1 | all Model-I/II ratio defects at or below ~2.5e-10 | PASS |
| Longitudinal defect (`qgn_longitudinal_defect_certificate.py`) | flagged in revision notes as not rerunnable (missing `qgn_search_v3` module and input JSON) | **now fixed**: dependencies restored from `legacy_packages/qgn_connected_irreducible_search_package.zip`; rerun passes with worst residual 1.41e-08 across M=5,6,7, all fillings | PASS (restored) |

Values labeled "FD floor" are finite-difference noise floors; they vary at the
same order of magnitude between machines/BLAS builds and do not affect any
verified identity. All machine-precision identities (1e-13 to 1e-18) and all
exact symbolic results reproduce exactly.

## Changes made to make this folder self-contained

1. `certificates/qgn_search_v3.py` and `certificates/search_connected_candidates_m6.json`
   were extracted from `legacy_packages/qgn_connected_irreducible_search_package.zip`
   so `qgn_longitudinal_defect_certificate.py` runs standalone (it previously
   failed on a missing import; see the revision notes).
2. Four hardcoded `/mnt/data/...` paths (in `qgn_longitudinal_defect_certificate.py`,
   `ghk_models_ed_search.py`, `qgn_connectivity_counterexample.py`,
   `restricted_qgn_reduction_certificate.py`) were replaced with paths relative to
   the script's own directory. These affect only input lookup and default output
   locations, not any computation; all scripts were recompiled and rerun after
   patching.
3. The audited manuscript replaces the pre-audit copy as the canonical `paper/`
   files; the pre-audit manuscript, the fixed-layout DOCX view of it, the audit
   corrections diff, and the revision notes moved to `revision_history/`.
4. `SHA256SUMS.txt` was regenerated for the new layout.

## How to rerun everything

```bash
pip install -r requirements.txt
cd certificates
for s in *.py; do python3 "$s"; done
```

Every script prints PASS (the four-cell and resonance scripts print their exact
defect tables). Runtime is a few minutes total on a laptop; the reduction and
sector-gap certificates dominate.
