# Exact Zero-Frequency Meissner Weight in a Gapless Frustration-Free QGN Parent

This reviewer-revised package contains the standalone paper devoted to the zero-frequency Meissner/superfluid-weight problem. It is separate from the broader dynamical-response manuscript and from the microscopic Schrieffer–Wolff transfer investigation.

## Central result

For a gauge-covariant swap-QGN/AGP parent on a periodic lattice of dimension `d >= 2`, the complete finite-volume transverse response is solved exactly by graph Hodge decomposition. Every divergence-free lattice field lies in the target kernel, so the static and imaginary-frequency transverse kernels agree at every finite size, allowed wave vector, and positive frequency regulator.

At fixed pair density `n/N -> rho` with `0 < rho < 1`, the electronic-link convention gives

```text
kappa_K = D_A = D_M = 2 j rho (1-rho) > 0,
```

and the pair-phase convention gives

```text
D_s = (j/2) rho (1-rho).
```

The theorem does not use a uniform many-body gap; the fixed-number gap closes as `L^{-2}`.

## Reviewer-integrated changes

The July 27 revision:

1. defines the general profile-dependent Kohn and Abelian kernels before specializing to the graph model;
2. defines the finite-`q` transverse kernel, fixes the meaning of `zeta=0`, and specifies admissible `q -> 0` sequences on finite reciprocal grids;
3. expands the space-time summability argument into a complete spectral/time-integral proof;
4. makes the kernel-eigenvector step explicit in the multiblock source-floor theorem;
5. aligns the clean-gas statement, proof, and verifier through a smooth finite-temperature regularization followed by `beta -> infinity`;
6. positions the result relative to flat-band quantum geometry, no-go/order-of-limits work, ODLRO-to-Meissner theorems, and eta pairing;
7. adds the isotropic spin-1/2 ferromagnet representation and explains that the exactly transverse response is purely diamagnetic at finite size;
8. adds an independent full-Fock fermionic verifier;
9. strengthens the compressed verifier with a half-filled 4x3 multipair gap check and an operator-valued multiblock-floor test with nonzero target-block mixing.

A detailed point-by-point response is in `RESPONSE_TO_REVIEW.md`.

## Files

- `exact_zero_frequency_meissner_weight_qgn.pdf` — reviewer-revised 19-page paper.
- `exact_zero_frequency_meissner_weight_qgn.tex` — complete LaTeX source.
- `verify_dynamical_response.py` — abstract response identities, sharp temporal bounds, and counterexamples.
- `verify_graph_hodge_meissner.py` — seniority-zero graph-Hodge model, multipair gap, and operator-floor checks.
- `independent_fock_check.py` — independent spinful-fermion full-Fock construction and response checks.
- `DYNAMICAL_RESPONSE_CERTIFICATE.txt`, `GRAPH_HODGE_MEISSNER_CERTIFICATE.txt`, and `FULL_FOCK_CERTIFICATE.txt` — retained passing outputs; `CERTIFICATE_SUMMARY.txt` gives a compact audit.
- `REVIEW_exact_zero_frequency_meissner_weight_qgn.md` — reviewer report supplied for this revision.
- `RESPONSE_TO_REVIEW.md` and `REVISION_NOTES.md` — comment-by-comment disposition and concise change log.
- `CLAIM_STATUS.md` — theorem scope, qualifications, and nonclaims.
- `BUILD_AND_PREFLIGHT.md`, `PDF_PREFLIGHT.txt`, `PDF_FONTS.txt`, and `build_pass2.log` — build and PDF quality records.
- `ENVIRONMENT.txt` — local tool and library versions.
- `MANIFEST.txt` and `SHA256SUMS.txt` — inventory and integrity hashes.

## Reproduction

From a Python environment with NumPy and SciPy installed:

```bash
python verify_dynamical_response.py
python verify_graph_hodge_meissner.py
python independent_fock_check.py
```

All three retained runs terminate with `OVERALL: PASS`.

To rebuild the paper with a standard TeX Live installation:

```bash
pdflatex -interaction=nonstopmode -halt-on-error exact_zero_frequency_meissner_weight_qgn.tex
pdflatex -interaction=nonstopmode -halt-on-error exact_zero_frequency_meissner_weight_qgn.tex
```

## Scope boundary

The paper proves a zero-temperature **matter** Meissner/superfluid weight for the exact local parent before coupling to a dynamical Maxwell field. It does not claim a finite-temperature phase theorem, a complete real-axis conductivity decomposition, or unconditional transfer to the unmodified finite-coupling microscopic Schrieffer–Wolff realization. That transfer problem requires separate uniform response control in the current-source `H^{-1}` topology.
