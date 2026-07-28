# Response to `REVIEW_exact_unscreened_leakage_v2.md`

## Overall assessment

The reviewer independently reran the shipped verifier, reimplemented every closed formula from scratch, verified all zero sets and tables, and confirmed that the theorem insertion was mathematically correct. This revision preserves that verified core and applies every requested polish item.

## Comment 1 — scale of the off-manifold singular values

**Reviewer concern:** “order-one off-manifold singular values” could be read literally even though the singular values are \(A_c\sigma_L(c)\).

**Revision:** The paper now states exactly that the singular values are \(A_c\sigma_L(c)\), nonvanishing at the same shell order as the selector shift before restoring the common microscopic factor \(\lambda^4\Delta\).

## Comment 2 — open corollary not referenced from the main text

**Revision:** Added two explicit main-text references to the open-path corollary: one immediately after the periodic proposition and one in the numerical-certificates discussion.

## Comment 3 — certificate self-description

**Revision:** Updated the verifier header, support README, platform-stable summary, and paper. They now state that exact enumeration verifies the centered first and second moments. The operator SVD then follows analytically from occupation-basis diagonality and orthogonality of distinct block-charge sectors; no floating-point Fock-space SVD is claimed.

## Comment 4 — extend the default enumeration range

**Revision:** Extended the default scan from \(L=2,\ldots,8\) to \(L=2,\ldots,9\). The exact verifier now checks 760 boundary/composition sectors and 699,040 hard-core configurations, with all 15 checks passing.

## Comment 5 — whether decimal open norms include \(A_4\)

**Revision:** `supporting_certificate/RESULTS.md` now says explicitly that the displayed decimal \(\tau^{\rm op}\) values at \(A_4=0.098\) already include \(A_4\); only the external microscopic factor \(\lambda^4\Delta\) remains.

## Optional connection to microscopic Meissner transfer

**Revision:** Added a concise terminology link. The off-manifold block
\[
(1-P_{\mathcal Z_n})H_{\rm ch}P_{\mathcal Z_n}
\]
is identified as the leading unscreened state source in the two-sided row-ideal language of the companion transfer paper. The leakage theorem remains self-contained.

## Additional strengthening

The revision adds an explicit order-four isolation gate \(\Xi_L^{(n)}(\lambda)\). This converts the qualitative warning about order-four leakage versus an order-six parent gap into a named, auditable diagnostic before unscreened full-model spectra are interpreted.
