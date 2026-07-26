MICROSCOPIC ORIGINS OF MULTIPLICATIVE QGN INTERACTIONS - DRAFT PACKAGE
=================================================

Primary manuscript
------------------
microscopic_origins_multiplicative_qgn_draft.tex
microscopic_origins_multiplicative_qgn_draft.pdf

Review documentation
--------------------
review_microscopic_qgn_revised.md       reviewer report addressed here
RESPONSE_TO_REVIEWER.md                 point-by-point response
REVISION_NOTES.md                       concise change log

Verification programs
---------------------
check_multiplet_uv_completion.py
check_electron_only_parity_channels.py
check_all_filling_active_remote_completion.py
check_global_lattice_downfolding.py
check_kramers_closed_shell_embedding.py
run_all_certificates.sh

Raw certificate outputs
-----------------------
multiplet_uv_completion_certificate.txt
electron_only_parity_channels_certificate.txt
all_filling_active_remote_completion_certificate.txt
global_lattice_downfolding_certificate.txt
kramers_closed_shell_embedding_certificate.txt

Stable summary
--------------
PLATFORM_STABLE_CERTIFICATE_SUMMARY.txt

The raw floating-point residuals may vary in their final digits with BLAS/LAPACK
implementations. The stable summary reports PASS/FAIL thresholds and rounded
physical coefficients suitable for cross-platform comparison.

Reproduction
------------
Python dependencies:
  python 3.11+
  numpy
  scipy

Run all certificates:
  bash run_all_certificates.sh

Build the paper:
  latexmk -pdf -interaction=nonstopmode -halt-on-error \
    microscopic_origins_multiplicative_qgn_draft.tex

Important scope notes
---------------------
- The all-filling local router and local Schur identities are checked on the
  complete four-orbital active Fock space.
- The global L=2 and L=3 lattice certificate is restricted to the one-pair
  active sector. It does not numerically certify the all-filling or
  volume-uniform theorem.
- The Kramers closed-shell certificate performs an exact one-bond Feshbach
  calculation in the one-pair active sector and checks the predicted
  lambda^8 residual scaling. The time-reversal transformation and bridge-count
  selection rule are proved analytically in the manuscript.
- The screened exact parent is tuned. The paper gives finite-lambda mismatch
  tolerances and treats the unscreened selector-plus-mixing family as the
  generic physical alternative.

Integrity
---------
MANIFEST.txt lists the distributed files.
SHA256SUMS records their SHA-256 hashes.
ENVIRONMENT.txt records the build environment and verification status.
