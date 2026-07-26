# Review: "Microscopic Origins of Multiplicative Null Interactions in Quantum-Geometric-Nesting Superconductors"

**Revised draft, July 26, 2026 (package `microscopic_qgn_revised_package_1`)**

---

## 1. Summary of the paper

The paper addresses a structural puzzle raised by the author's companion no-go theorem: an additive, fully Peierls-covariant one-body parent cannot generate off-diagonal composition hopping in a two-block QGN superconductor, so the minimal escape is the multiplicative null row K_x = Q_{a,x} B_x (odd-simultaneous-bond-parity projector times crossed transfer). The apparent artificiality of that row is the target. The paper shows, in increasing order of ambition:

1. **Local algebra.** K†K has an exact rank-three seniority-zero representation in conventional density/pair-hopping operators, and B² alone fails to annihilate the product-AGP zero manifold (so it must be compensated, not ignored).
2. **Sign obstruction.** Second-order passive elimination is negative semidefinite, so a bare +αK†K cannot arise passively; the sharp full-Fock bound K² ⪯ 9Q_a gives a clean positivity window for the sign-reversed passive route.
3. **Electron-only parity splitting.** A four-site plaquette with ordinary hopping, an onsite seniority gap, and intersite repulsion produces an even/odd simultaneous-parity doublet split at fourth order with an explicit closed-form coefficient J_parity.
4. **All-filling parity router (the central result).** The exact identities C_s² B P₀ = 4Q_s B P₀ and C_a² B P₀ = 4Q_a B P₀ show that parity projectors emerge from interference of ordinary orbital-current paths on the universal broken-pair source, at every local filling. A four-orbital control sector then yields the exact Schur complement −(t_B²/d_s)B² + αK†K under a full high-block gap.
5. **Lattice downfolding.** A weighted finite-order local Schrieffer–Wolff theorem for a periodic array of overlapping controls: the screened family gives the composition-degenerate parent at order λ⁶ with O(λ⁸) remainder and residual block coupling; the unscreened family retains an order-λ⁴ diagonal composition selector before the λ⁶ mixing.

Four Python verification programs with shipped certificates accompany the draft.

## 2. Verification performed for this review

**Reproduction.** All four `check_*.py` scripts were rerun (NumPy 2.4.4 / SciPy 1.17.1). All exit cleanly and report PASS. Outputs match the shipped certificates line by line except for last-digit differences (~10⁻¹⁵) in a handful of eigenvalue residuals — ordinary LAPACK platform noise, not a concern.

**Independent re-implementation.** I wrote a from-scratch 8-mode fermionic Fock-space verifier (Jordan–Wigner, explicit sign checks, no code shared with the package) and confirmed:

- spec(K²|_sz) = {0×13, 4×3} — the rank-three identity (Prop. 3.2 / Eq. rank-three);
- P₀BP₀ = 0 and BP₀ = P₂BP₀ — the universal source identity (Lemma 6.1);
- (J₁²−1)BP₀ = 0, (J₂²−1)BP₀ = 0, (J₁J₂−W)BP₀ = 0, and the routing identities C_s²BP₀ = 4Q_sBP₀, C_a²BP₀ = 4Q_aBP₀ — all exactly zero (Lemma 6.2);
- min spec(9Q_a − K²) = 0 to machine precision, with max|spec(B)| on the odd sector exactly 3 — Lemma 4.2 is sharp as claimed;
- P_Z₁B²P_Z₁ = 2[[1,1],[1,1]] (Eq. B2-compression), and K annihilates all nine local product-AGP components;
- independent ED of the plaquette: split/t⁴ → 0.0448122461 as t→0, matching the closed form 256(2Δ_p+3V_p)/[V_p(2Δ_p+V_p)(4Δ_p+3V_p)²] at Δ_p=10, V_p=3;
- the charge-potential formula (Prop. 7.2), including the accidental vanishing of the selector coefficient 16A_c(L−2)/[L(L−1)] at L=2.

**Internal numeric cross-checks.** All quoted numbers are mutually consistent: α = t_B²(1/9.5 − 1/9.875) = 1.958694×10⁻³; ᾱ = 4b²/Δ̄_m²(s²/Δ̄_s − a²/Δ̄_a) = 0.0018375; 4α = 7.83×10⁻³ against feedback 3.46×10⁻² at t_B = 0.7 and 1.55×10⁻⁵ vs 1.60×10⁻⁴ at t_B = 0.1 (consistent with the ~0.15·t_B⁴ law in the certificate table); the compensated one-pair curvature matrix in the certificate ([[2.1707, −0.1707], [−0.1707, 2.7707]]) matches the 2j_a + 8α structure of Eq. (one-pair-curvature); the perturbative expansion of the local Schur result reproduces Eqs. (virt-six) and (betas) term by term.

I found no mathematical errors. Every exactness claim I tested is exact, and every claim the paper labels as perturbative or conditional is labeled correctly.

## 3. Assessment

This is a careful, unusually honest piece of work. The claim hierarchy on page 2 and the "What is and is not exact" remark are exemplary — the revised draft consistently distinguishes finite-dimensional identities, gap-conditional exact statements, and finite-order perturbative statements, and the numerical section explicitly states what the L=2/L=3 one-pair audits do *not* certify. The central routing identity (Lemma 6.2) is a genuinely nice result: it is exact, filling-independent, and converts the seemingly ad hoc parity projector into current-path interference. The sharp bound K² ⪯ 9Q_a with an analytic proof, the closed-form plaquette coefficient verified two independent ways, and the L=2 selector blindness caveat all indicate the revision responded seriously to the earlier referee round.

The main scientific limitation — which the paper itself states — is that the construction demonstrates *existence within a conventional interaction class*, not naturalness: the control sector, coefficient matchings, and screening counterchannel are engineered. My major comments below are about making the cost of that engineering quantitatively explicit.

## 4. Major comments

**M1. Quantify the fine-tuning sensitivity of the screened family.** The screened theorem requires the exact matchings β₄ = b²/Δ̄_m and β₆ = 4b²s²/(Δ̄_m²Δ̄_s), plus a counterchannel that exactly cancels the monopole polynomials. Section 3.3 acknowledges qualitatively that imperfect cancellation is "a controlled QGN-breaking perturbation," but the lattice section never states the required accuracy. A mismatch δβ₄ leaves a term λ⁴Δ·δβ₄·ΣB_x², which is *two weighted orders larger* than the λ⁶ target; the composition-degenerate parent survives as the leading composition physics only if δβ₄ = O(λ²·ᾱ) (and δβ₆ = O(ᾱ), δA₄ likewise for the selector). One short subsection or remark stating these tolerances would materially sharpen the paper and preempt the obvious referee objection that the screened family is a measure-zero point. (The unscreened corollary already handles the generic case, which is a strength — say so more loudly.)

**M2. Dependence on unpublished companion manuscripts.** The motivating no-go theorem and the base two-block parent are cited as manuscripts (Refs. SledgianowskiParent/MultiBlock/NoGo) with no arXiv identifiers. The paper cannot be independently evaluated as a closed logical unit until those are available; in particular, the claim that K_x is "the minimal structural escape" is inherited, not proven here. Either add arXiv IDs/links, or include a short self-contained appendix stating the precise no-go hypotheses and conclusion being escaped. Section 2.2's one-paragraph summary is not quite enough for a referee to check that all no-go hypotheses except additivity really are retained.

**M3. Connect to the perturbative-gadget literature.** The four-orbital control sector is, structurally, a perturbative gadget (a mediator ancilla with engineered virtual denominators), and the weighted hierarchy with delayed cross-gadget corrections is closely analogous to gadget-overlap analysis in Hamiltonian complexity (Kempe–Kitaev–Regev; Oliveira–Terhal; Jordan–Farhi; Bravyi–DiVincenzo–Loss you already cite). Citing this literature would both situate the contribution and let you borrow standard language for the overlap-correction argument. It would also clarify what is genuinely new here: the *exactness* of the router identity on the full source space (gadget constructions are usually only perturbative), which deserves more emphasis.

**M4. The Kramers-doubled embedding is asserted, not constructed.** Section 6.3 states that a time-reversal-invariant realization exists by doubling the control sector and that "the two partners produce the same reduced coefficient." Since the current vertices are time-reversal odd and this is the step that makes the model physical, a short appendix writing the doubled Hamiltonian and verifying (or certificate-checking) the equality of the two partners' reduced coefficients would close a real gap. Relatedly, α > 0 requires s²/Δ̄_s > a²/Δ̄_a — same-chirality current scattering strictly stronger than opposite-chirality. Is there any physical mechanism (symmetry, geometry, screening) that generically produces this asymmetry, or is it another tuned dial? One sentence either way belongs in the Discussion.

**M5. The weighted regrading of Bravyi–DiVincenzo–Loss deserves slightly more detail.** Appendix E argues the linked-cluster machinery applies "coefficient by coefficient after this analytic regrading," and the fermionic issue is dispatched in one sentence (graded tensor product / sector-wise Jordan–Wigner in 1D). I believe both claims, but they carry the entire global theorem. In particular: (i) the BDL construction is formulated for a fixed perturbation V, whereas here different terms enter at different λ powers — a paragraph explaining why the recursion's convergence radius and locality estimates survive the regrading (e.g., by treating the family as analytic in λ and grading the Taylor coefficients) would help; (ii) the claim that shell/screening/base terms, being control-diagonal, "do not create connected multi-control virtual walks below that order" could state explicitly that they can still *extend the support* of a cluster and why that does not lower the degree-eight threshold.

## 5. Minor comments

1. **Notation drift Π⁻₁₂ vs Q_a.** The abstract and Eq. (intro-K) use Π⁻₁₂,ₓ; the body from Eq. (Qsa) onward uses Q_{a,x}; the certificates use `Pi12`. Pick one symbol (or state the identification once, immediately after Eq. (Qsa)).
2. **Abstract length.** At four dense paragraphs, the abstract reads like a second introduction. The first paragraph plus a compressed version of the router/downfolding results would serve better; the claim-hierarchy box already carries the fine print.
3. **Theorem 4.1** is the textbook second-order SW sign statement. Consider demoting to a lemma or remark, keeping "Theorem" for results with real content (4.2's sharp bound, 6.2, 6.3, 7.x).
4. **Lemma 4.2 proof:** "Both are fixed by the edge exchange and hence lie in the W = +1 sector" — for a four-fermion state, invariance under a mode permutation involves a sign; one clause noting that the sign works out (+1 for the filled ± shells) would make the proof airtight. (I verified numerically that it does.)
5. **"Weighted degree"** is used from Section 7.2 onward but only defined precisely in Appendix E. A one-line definition at first use would help.
6. **Section 6.3, Eq. (normal-four-fermion):** the minus sign on the right-hand side is convention-dependent (mode ordering); either state the ordering or drop the explicit sign.
7. **Table 2:** the row "Plaquette t=0.02: (E_a−E_s)/t⁴ = 0.044811672577" would be more useful with the analytic limit 0.044812246105 printed beside it, so the reader sees the convergence gap (~5.7×10⁻⁷) is O(t²)-consistent.
8. **Certificates:** the shipped certificate files differ from fresh reruns only in ~10⁻¹⁵ residuals, as expected across BLAS builds. Consider printing residuals with a PASS/FAIL threshold only (or rounding to 2 significant figures) so that certificate diffs are bitwise stable across platforms.
9. **References:** Peotta–Törmä (Ref. Peotta2015), Schrieffer–Wolff (Ref. Schrieffer1966), and Kato (Ref. Kato1976) appear in the bibliography but I could not find them cited in the text — either cite or prune.
10. **Figure 3** is helpful; consider annotating the shared active sites between adjacent controls, since "neighboring controls share active orbitals" is the entire reason the global theorem is nontrivial.
11. The phrase "minimal structural escape" is now used carefully (minimality in hypothesis space) — good; make sure the abstract's "the minimal structural escape" carries the same qualifier, since the abstract is where casual readers will over-read it.

## 6. Verdict

The mathematics checks out completely: I reproduced all four certificates and independently re-derived every central exact identity (rank-three reduction, source identity, current router, sharp bound, plaquette coefficient, charge projection), finding no errors. The revision notes' claimed fixes are all actually present in the draft. The remaining weaknesses are of positioning and completeness — fine-tuning tolerances left implicit (M1), reliance on unavailable companions (M2), missing gadget-literature context (M3), and two asserted-not-shown steps (M4, M5) — none of which threaten the core results.

**Recommendation: minor-to-moderate revision.** Address M1–M5 (M1 and M2 being the most important) and the minor list; no new results are required.

---

*Review performed against the July 26, 2026 package. Verification: all four `check_*.py` scripts rerun successfully (PASS, certificates reproduced modulo ~10⁻¹⁵ platform noise); independent 8-mode Fock-space re-implementation of the local identities, plaquette ED, and charge combinatorics all confirm the paper's claims.*
