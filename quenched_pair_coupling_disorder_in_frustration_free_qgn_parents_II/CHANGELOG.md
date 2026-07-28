# Changelog: monolithic draft to two-paper reviewer revision

## Structural split

- Moved the exact weighted parent, Hodge formula, deterministic bounds, dynamical graph filter, deterministic counterexample, homogenization, almost-sure source tightness, annealed finite-torus rates, and multiblock swap backbone to Paper I.
- Moved the Abel variance theorem, Born concentration, sharp-cutoff closure routes, maximal-shell analysis, one-bond Krein geometry, participation coarea, acoustic obstruction, and disorder-noise target to Paper II.
- Replaced the approximately 65-page monolith with two 21-page papers, each with a focused abstract, claim hierarchy, theorem boundary, bibliography, and reproducibility section.

## Self-containedness and analytic sourcing

- Defined the imaginary-frequency kernel in Paper I and derived the positive-square dynamical residual formula in an appendix.
- Added a precise definition of Kohn-Abelian equality and references for the terminology.
- Defined the torus heat kernel, parabolic weight, horizontal generator, periodic vertical derivative, and period-uniform spectral-gap condition.
- Restated the exact GNO long-version inputs by theorem/equation number and displayed the weighted-Holder exponent cancellation.
- Fixed the finite-volume probability spaces: space-periodization for almost-sure limits and an independent product law on torus coefficient blocks for quantitative estimates.

## Claim-boundary and notation fixes

- Added the explicit caveat that the bounded-effective-mode sharp variance follows only from `Var T <= E T^2`.
- Renamed the generic heat exponent to `sigma`, the heat edge field to `w_e(t)`, and failure probability to `vartheta`.
- Added monotonicity and mollification details to the relative Wegner proof.
- Added the necessary condition `m_e(r)<0` for a positive one-bond crossing.
- Removed package-history wording and replaced it with self-contained statements.
- Cited or removed every bibliography entry.

## Reproducibility

- Replaced a shared Monte Carlo random stream with deterministic per-size `SeedSequence` instances.
- Fixed the default random-torus sample count at 24.
- Recorded Python, NumPy, and SciPy versions in the machine-readable certificate.
- Recompiled both PDFs without undefined references, citation warnings, or overfull boxes.
