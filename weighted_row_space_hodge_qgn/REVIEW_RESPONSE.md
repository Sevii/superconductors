# Response to reviewer comments

## 1. Corollary 2.4 norm-equality overclaim

Corrected. The paper no longer states that Hodge projection is necessary for equality of operator norms. It now distinguishes:

- the unique Loewner-order minimizer `eta_*`;
- potentially nonunique operator-norm minimizers when `rank(P) >= 2`.

The exact condition is

`C_kappa <= ||C_*|| P - C_*`.

A concrete rank-two counterexample is included in the text and in both verification scripts.

## 2. Gapless branch scope and multichannel proof

The abstract and claim hierarchy now state explicitly that the general exponentially local soft-mode branch is not solved. The analytic-strip theorem is labeled as a uniformly reduced-gapped route and followed by a scope remark explaining why it does not cover genuine real-axis rank drops.

The former single-channel scaling proof has been replaced by a channel-resolved proposition using mutually orthogonal spectral projectors of `M(k)`. The result gives an exact block decomposition of the cost, an upper estimate from channelwise source orders, and a converse necessary bound. Approximate channel projectors are explicitly excluded unless off-block couplings are controlled separately.

## 3. Higher-dimensional Laurent-module decision

Added a theorem reducing Laurent-module membership to the saturated polynomial module `N:h^infinity`. A module Gröbner-basis elimination and syzygy computation now gives a terminating existence/nonexistence decision in every dimension over an effective coefficient field. The text is explicit that worst-case complexity may be poor and that floating-point data do not support exact nonmembership without reconstruction or enclosure.

## 4. Cluster theorem hypotheses and support bookkeeping

Revised as follows:

- `p_X P = P` is now a numbered, boxed hypothesis of the cluster theorem.
- The parent-row support radius `R_D` and enlarged support `X^{+R_D}` are defined explicitly.
- The locality theorem includes the buffer factor `N_(R_D) exp(2 mu' R_D)`.
- The cluster theorem points forward to the locality theorem for infinite-volume convergence.
- The interaction-norm definition now explains that only summability and one triangle inequality are used; standard F-function norms are cited as alternatives.

## 5. Minor corrections and references

- Corrected the mangled dagger placement in the proof of the Hodge identity.
- Added the note after the `b` definition that `Z J_R P = 0` makes it exactly the O1 target expression.
- Cross-referenced the normalization appendix from the stiffness-floor equations.
- Added references for Moore-Penrose/SVD conventions, Smith form, Gröbner bases and syzygies, Chebyshev inverse approximation, banded inverse decay, and Lieb-Robinson/F-function interaction norms.
- Added pointers to the companion QGN functional, Peierls escape, microscopic transfer, and zero-frequency Meissner manuscripts.
