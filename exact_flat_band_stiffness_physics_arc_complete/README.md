# Exact Flat-Band Stiffness Physics Arc - Complete Package (audited revision)

This archive accompanies the paper:

> Nicholas Sledgianowski, **Exact Flat-Band Stiffness from Pair Mobility in a QGN Model: Finite-size obstructions, a many-body reduction theorem, and a proof of the Gao-Han-Khalaf formula for Model II** (revised July 23, 2026; audited corrections applied).

## Primary deliverables

- `paper/exact_flat_band_stiffness_physics_arc.pdf` - authoritative fixed-layout paper (audited revision).
- `paper/exact_flat_band_stiffness_physics_arc.md` - editable canonical manuscript source (audited revision).

The audited revision incorporates the external audit's mathematical corrections:
the lifted finite-volume pair family replacing the invalid momentum-block
identification, the explicit common Fock-space and Hermitian positive-kernel
hypotheses, the expanded spanning-tree connectivity proof, Lemma 2's validity
range, corrected locality constants, the displayed trace proof of the exact
Model-II cancellation Gamma = 0, precise reproducible counterexample statements,
and a reproducibility section synchronized with this archive. See
`revision_history/exact_flat_band_stiffness_revision_notes.md` for the itemized
list and `revision_history/exact_flat_band_stiffness_physics_arc_corrections.diff`
for the exact diff. The pre-audit manuscript, its PDF, and its fixed-layout DOCX
view are retained in `revision_history/` for provenance; they are superseded.

## Verification

`VERIFICATION_RERUN_2026-07-23.md` records a full independent rerun of every
certificate in this archive against the audited paper's claims: all pass, all
machine-precision identities and exact symbolic values reproduce exactly, and
the previously non-rerunnable longitudinal-defect certificate has been restored
(its missing helper module now ships in `certificates/`). Raw rerun logs are in
`certificates/data/rerun_2026-07-23/`.

## July 23 mathematical-review revision

The earlier same-day revision added the missing all-size Model-II connectivity
proof for `0 < xi < pi/2`, removed non-load-bearing boxed sector-gap assertions,
corrected the H2 derivation wording, separated pair-orbital and torus-side
notation, regularized equation numbering, and added the requested presentation
clarifications. The unused all-size Model-I connectivity assertion has been
withdrawn; its finite graph audit remains reproducibility evidence only.

Review records and the itemized responses are in:

- `reviews/review_flat_band_stiffness_math_2026-07-23.md`
- `audit/RESPONSE_TO_MATH_REVIEW.md` and `.pdf`
- `audit/math_review_revision.diff`
- `audit/FINAL_REVISION_QA.md`

## Scope and earlier review records

- `audit/prior_art_and_scope_audit.pdf` and `.md` - targeted literature and theorem-boundary audit, including the non-Hermitian QGN boundary and Alon-Puder verification.
- `audit/RESPONSE_TO_FINAL_CONSOLIDATION_REVIEW.md` - itemized response to the earlier consolidation review.
- `audit/rank_one_aldous_note_outline.md` - deferred, non-load-bearing mathematics-note outline.
- `reviews/` - independent review records supplied during the project.

## Reproducibility

The `certificates/` directory is now fully self-contained: it holds all nine
standalone scripts, the helper module `qgn_search_v3.py` with its input records,
the archived outputs (`certificates/data/`), and the July 23 rerun logs:

- exact finite-size winding obstruction;
- resonant/nonresonant aliasing controls;
- fixed-local reducibility obstruction;
- connected-model search and longitudinal-defect certificate;
- full multiband reduction certificate;
- H4 connectivity and GHK graph audit;
- sector-gap reduction controls;
- Model-II pair dispersion, `Q=2A`, positive mass, and finite-size stiffness checks;
- independent projected-Hubbard ED for GHK Models I and II.

The `legacy_packages/` directory retains the archived stage packages as
provenance. The consolidated paper supersedes their separate manuscripts but not
their numerical records.

## Key rerun commands

From an environment with Python 3.11+ and `numpy`, `scipy`, and `sympy`
(`pip install -r requirements.txt`):

```bash
python certificates/model_II_pair_dispersion_certificate.py
python certificates/restricted_qgn_reduction_certificate.py
python certificates/qgn_h4_connectivity_certificate.py
python certificates/qgn_longitudinal_defect_certificate.py
```

All scripts write any outputs next to themselves (no absolute paths). Every
script is standalone apart from standard scientific Python dependencies.

## Verification status (archived certificate values)

The final Model-II certificate reports:

- direct lifted pair block versus closed formula: maximum error `4.259e-16`;
- electronic twist versus `Q=2A`: maximum error `5.881e-16`;
- finite-grid mass Hessian: maximum error `3.558e-09`;
- many-body finite-size curvature: maximum error `2.496e-09` (rerun: `5.10e-09`, finite-difference floor).

The general reduction certificate reports a maximum tensor-identity error of
`1.249e-18`. The H4 certificate verifies the connected and disconnected
zero-mode counts, attains the full Lie dimensions 16 and 25, and audits 55
finite GHK projected-frame graphs (all connected, maximum Parseval-frame error
`8.10e-16`). The exact four-cell winding certificate gives raw defect `2/3` and
normalized defect `1/24`. The sector-gap certificate retains 160 rank-one cases
to worst relative deviation `1.02e-13` with non-rank-one control ratio
`0.432639`. The restored longitudinal-defect certificate passes with worst
residual `1.41e-08`.

## Integrity

`SHA256SUMS.txt` contains hashes for every archived file except itself.
