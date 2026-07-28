# Review: exact unscreened leakage certificate v2 and manuscript update

**Package:** `exact_unscreened_leakage_certificate_v2` (update to *Microscopic Origins of Multiplicative QGN*)
**Review date:** July 28, 2026
**Verification performed:** fresh rerun of the shipped verifier; independent from-scratch re-implementation of every closed formula; patch application and full LaTeX compile of the revised manuscript.

---

## 1. Verdict

The mathematics is correct, the certificate is genuine, and the manuscript edits do exactly what the patch notes promise. I found no errors in any formula, any table entry, or any claimed zero set. The wording changes consistently and correctly narrow the old "diagonal in composition space" and "selector-plus-mixing" language to compression statements, and the new open-boundary corollary is both correct and genuinely necessary for the DMRG plan. This is ready to integrate, modulo the small comments in Section 5.

## 2. What I verified and how

**Fresh rerun.** `check_unscreened_agp_leakage.py` reruns cleanly in a clean environment and reports 15/15 PASS. The regenerated certificate and every regenerated CSV are byte-identical to the shipped files.

**Independent re-implementation.** I wrote a from-scratch brute-force checker sharing no code with the package: it enumerates all subset pairs for L = 2,...,8, both boundaries, all (n₁,n₂), computes exact Fraction moments of the diagonal charge operator, and compares against the closed formulas as transcribed directly from `DERIVATION.md`. Every check passed with zero failures, specifically:

- the periodic mean ⟨D_L⟩ (both forms) and the closed variance σ_L², including the equality Var(D_L) = σ_L² in every sector;
- the periodic zero set (one block empty/full, other in {0,1,L−1,L}; all sectors at L=2) — exact match, no extra or missing zeros;
- the open configuration-wise identities ΣN_x = 4n−2B and ΣN_x² = 8n−4B+8Z on every one of the 87,376 open configurations;
- the moment formulas b_L, h_L, e_L, w_L and the identifications Γ₀ = Var(Z), Γ₁ = Cov(Z,B), Γ₂ = Var(B), plus the sign structure Γ₁ ≤ 0 everywhere;
- the full open variance (τ^op)² = 64A₄²Γ₀ − 16A₄V₄Γ₁ + V₄²Γ₂ against the directly computed Var(A₄D − U₄S/2) at the project default V₄/A₄ = 1;
- the open compressed selector μ^op = C_L^op(n) + 16A₄(L−2)n₁n₂/L² and the (L−1)/L ratio claim;
- the open zero set {0,L}×{0,L} (and the disappearance of the periodic L=4, n=3 accidental zeros);
- all twelve Phase-3 Gram-diagonal tables in `RESULTS.md` (periodic and open), rank counts, and the decimal norm table at A₄ = 0.098 (all entries agree to the printed precision);
- the L=4 (2,0) benchmark ‖R|2,0⟩‖ = (8√2/3)A₄.

The claimed audit counts also check out arithmetically: 2·Σ_{L=2}^{8}(L+1)² = 560 sectors and 2·Σ 4^L = 174,752 configurations.

**Structural soundness of the SVD claim.** Because C₄ is diagonal in the hard-core occupation basis, everything reduces to classical moments of uniform fixed-size subsets, and the operator-level SVD statements follow rigorously from two facts the derivation states correctly: the centered image is orthogonal to |c⟩ within its composition sector (Z_n intersects each sector in exactly the one product-AGP line), and images from distinct compositions are orthogonal by block-charge conservation. The proofs of the mixed-covariance steps (E[C|X₁] constant on the cycle; E[C|X₁] = (2n₂/L)(2n₁−B₁) on the path) are correct.

**Manuscript patch.** The unified diff applies cleanly to `microscopic_origins_multiplicative_qgn_draft.tex` and produces output byte-identical to the shipped `..._leakage_updated.tex`. The updated manuscript compiles to a 30-page PDF with zero undefined references, zero multiply-defined labels, and no new errors (my container lacks the `lmodern` font package; commenting that line out was the only accommodation — not a manuscript issue). The notation A_c, U_c, V_ch^(L) in the inserted proposition matches the pre-existing charge section exactly.

## 3. Assessment of the scientific content

The central correction is the right one. The old Proposition (charge potential) said H_ch "is diagonal in composition space," which invites the invariant-subspace reading; the replacement makes it a compression statement and the new periodic proposition supplies exactly the missing object — the off-manifold block, in fully closed form, with singular values, Gram operator, norm, rank, and the complete zero set. The formulation via the composition-resolved centered action is elegant: it makes the SVD essentially one line once block-charge orthogonality is noted, and it produces the zero set by inspection of a manifestly nonnegative expression.

The boundary-scope split (T2.3a/T2.3b) is well judged. The structural observation driving it — that on the path the linear charge term stops being scalar because endpoint occupancy fluctuates, so the open operator picks up the −V₄B term and the leakage acquires V₄/A₄ dependence — is correct, clearly explained, and has the practically important consequence that open normalized leakage is meaningless without recording V₄/A₄. Flagging that the periodic L=4, n=3 accidental zeros do not survive the open truncation is exactly the kind of trap the Phase-3 numerics would otherwise have fallen into.

One edit deserves specific praise: the one-pair curvature paragraph. The old text said the order-four selector "splits the composition levels before the order-six mixing is applied," which is actually false in the one-pair sector (both compositions have n₁n₂ = 0, and I confirmed the leakage is exactly zero there too). The new text states the correct thing — zero selector shift and zero leakage at n=1, with the caveat deferred to higher fillings. That is a real correction, not just re-wording.

The two remaining occurrences of "selector-plus-mixing" in the revised manuscript are both deliberately qualified ("...only after an additional mechanism isolates the product-AGP manifold at order four or stronger"), which is consistent with the new claim discipline. The abstract, claim-hierarchy list, discussion, and conclusion edits are mutually consistent.

## 4. The Meissner paper: relevant context, no action required

The companion release (`microscopic_transfer_exact_meissner_revised_release`) is complementary rather than overlapping, and nothing in the leakage update depends on it or conflicts with it. It concerns the **matched/screened** family, where the order-four terms are canceled and the first manifold-breaking state source appears at weighted degree eight (the B_e⁴/commutator-square source, with the C₈ counterterm and the dark-branch construction as the two repair routes). The present update concerns the **unscreened** family, where the manifold already breaks at order λ⁴. Together they make a coherent pair of instances of the same phenomenon — exact compression without invariance, with an explicit state source at the leading available order in each family — and the leakage update's "valid routes" (present as compression / study the true charge-selected states / add an order-four isolation mechanism) are the unscreened analogue of the Meissner paper's routes A and B.

If you want to connect them, the natural touch point is language: the Meissner paper's two-sided row-ideal / state-source formalism (QRP) is exactly what R = (1−P)C₄P is an instance of, and a one-sentence cross-reference in the QGN discussion would unify the two papers' vocabulary. Purely optional; the leakage update stands alone.

## 5. Minor comments

1. **"order-one off-manifold singular values" (updated tex, after the periodic proposition).** Strictly the singular values are A_c·σ_L(c), i.e. order one only after normalizing by A_c (σ_L values are ~5.4–9.2 in the Phase-3 sectors, but A₄ = 0.098). Since the sentence says "before the microscopic factor λ⁴Δ is restored," consider "order-one per A_c" or "nonvanishing, of the same order as the selector shift" to avoid a literal misreading.

2. **The open corollary is never referenced from the main text.** `cor:charge-leakage-open` carries its label but no `\cref` points to it. A one-line pointer — most naturally in the numerical-certificates section or in the remark after the periodic proposition ("for the open truncation used by the DMRG plan see Appendix ...") — would keep open-chain readers from missing it.

3. **Certificate self-description.** The "explicit centered action/SVD norm" PASS lines verify moment identities (centered sum zero, centered square = N·variance); the operator-level SVD then follows analytically from diagonality + block-charge orthogonality. That inference is airtight, but a one-line note in the certificate header saying so would preempt a skeptical reader asking why no Fock-space SVD is computed.

4. **Enumeration range.** The default certificate covers L = 2,...,8 (one retained nondefault run reaches L = 9). Since the formulas are rational in n₁, n₂ per L but the derivations are for general L, extending the default range by one or two volumes is cheap insurance if you ever want it; not needed for correctness, as the analytic proofs carry the general-L claims.

5. **`RESULTS.md` decimal norms table** is labeled "At A₄ = 0.098, the open composition-resolved leakage norms" — these are τ values (A₄ included), which then multiply λ⁴Δ with A₄ *not* double-counted since the manuscript's eq. (unscreened-selector) keeps A₄ inside the bracket. Consistent, but worth a half-sentence in RESULTS.md stating that A₄ is included in the tabulated numbers, since the neighboring Gram table is normalized per A₄².

## 6. Bottom line

All 15 certificate checks reproduce; an independent re-implementation confirms every closed formula, zero set, and table exactly; the patch applies cleanly, compiles cleanly, and the prose edits are correct and internally consistent — including one place where the edit fixes a genuinely wrong statement (one-pair sector). The Meissner paper is complementary context only. I recommend applying the update as shipped, with the five minor comments above as optional polish.
