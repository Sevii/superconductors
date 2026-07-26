# Two-block QGN reviewer-corrected package

Primary deliverables:

- `two_block_qgn_stiffness.pdf` -- compiled reviewer-corrected manuscript.
- `two_block_qgn_stiffness.tex` -- LaTeX source.
- `REVIEWER_RESPONSE.md` -- map from the review comments to manuscript and code changes.

Corrected verification:

- `scripts/check_bridge_covariance.py`
- `scripts/interaction_singular_subbundle_scan_revised.py`
- `scripts/check_projector_geometry_convergence.py`
- `scripts/wse2_interaction_subbundle_scan.py`
- `legacy/verify_local_material_qgn_legacy.py` -- retained only to document the superseded fixed-kernel convention.

Machine-readable outputs:

- `data/bridge_covariance_audit.json`
- `data/verified_archived_candidates.json`
- `data/verified_archived_candidates.csv`
- `data/projector_geometry_convergence.json`
- `data/projector_geometry_comparison.csv`

The revised manuscript compiles without undefined references, overfull boxes, or PDF preflight warnings.
