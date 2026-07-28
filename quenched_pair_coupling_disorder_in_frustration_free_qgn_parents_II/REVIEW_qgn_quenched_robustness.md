# Review: "Quenched Pair-Coupling Disorder in Frustration-Free QGN Parents"

**Package reviewed:** `qgn_quenched_pair_coupling_robustness_fixed_cutoff_noise_package` (draft dated July 27, 2026; 3,043-line LaTeX source, ~65 pp. compiled, plus verification archive)

**Reviewer's overall assessment.** This is a strong, unusually honest research draft. The central exact result — that the many-body branch curvature of the weighted swap parent is *exactly* the classical weighted network cell energy, `C = γ_{N,n} min_φ Σ_e j_e(a_e + (B*φ)_e)²` — is clean, correct, and load-bearing for everything that follows. The reduction of the dynamical order-of-limits question to an H⁻¹ tail of the random-conductance drift is elegant, the macroscopic-modulation counterexample is genuinely clarifying, and the sharp-cutoff audit (Wegner, count-vs-current obstruction, randomization no-go, single-crossing Krein geometry, participation coarea, acoustic-CLT obstruction, Doob/noise-semigroup reformulation) is a coherent research program with its one remaining open input stated precisely. The claim hierarchy, CLAIM_STATUS.md, and source maps set a standard for proved/imported/open bookkeeping that most papers should envy.

The two most significant issues are self-containedness gaps, not errors: (1) the dynamical kernel 𝒦 is never defined in this paper, and (2) the proof of the periodic gradient-propagator lemma — the single most load-bearing analytic step in the quenched variance theorem — uses undefined symbols and unverifiable pointers into the imported literature. Details below.

---

## 1. What I verified

### Hand-checked algebra (all correct)

- **Tangent-space constants** (Lemma 3.1): ⟨f|g⟩ = α_{N,n} f*g with α = n(N−n)/[N(N−1)] follows correctly from Eν_x = n/N, Eν_xν_y = n(n−1)/[N(N−1)] and the zero-mean condition; ‖z_e‖² = 2α = γ via f·f = 2. ✓
- **Weighted Hodge proof** (Thm 4.1): the normal equations, residual orthogonality, and Moore–Penrose substitution are all right. ✓
- **Ring harmonic-mean law** (Prop 4.6) and the **modulated-ring counterexample** (Prop 6.1): I recomputed ∮dθ/(1+ε cos θ) = 2π/√(1−ε²), the first moment ⟨s, M s⟩ = ¼‖BJ1‖² = (Nε²/2)sin²(π/N), and the spectral-Markov step. All correct; the logic that all positive source mass collapses to vanishing target eigenvalues, so the fixed-ζ multiplier → 1, is sound. ✓
- **Site/target tail identity** (Prop 5.2), including the factor 4 between M_J = C*C/4 and L_J = CC*: correct via the singular-value correspondence ⟨v,g⟩ = √μ⟨u,s⟩. ✓
- **Massive-resolvent tightness** (Thm 7.5): the inequality 1/λ − 1/(λ+m) ≥ 1/(2λ) on λ ≤ m and the I(m) ↑ I monotone limit are fine. ✓
- **Heat-to-tail conversion** (Prop 8.1) and the **Abel mean bound**: ∫_{1/δ}^∞ (1+j₊t)^{−1−d/2} dt ≤ C j₊^{−1}(δ/j₊)^{d/2}; the Stieltjes integration by parts is right. ✓
- **Abel variance bookkeeping** (Thm 8.9): the τ^{−d/4} linear tail, the convolution tail (using integrability of h for d ≥ 1), the L^{−2d} → L^{−d} tensorization count over dL^d bonds, and the j₊ rescaling 𝔄^J(δ) = j₊𝔄^{J̃}(δ/j₊) all check out. ✓
- **Born sector** (Thm 8.13): ‖M_δ‖ ≤ |A|², rank ≤ Weyl count, the Hanson–Wright application, and the exact Fourier mean (I re-derived ⟨φ_q, B D_a² B* φ_q⟩ = Σ_μ A_μ² 4sin²(q_μ/2)). ✓
- **Relative Wegner** (Thm 9.1): the ramp-function/functional-derivative/rank-one-interlacing argument closes and the constant dL^dρ_∞ j₊ ℓ/E is right. ✓
- **Count-vs-drift obstruction** (Prop 9.2): total mass = arithmetic − harmonic → 1 − √(1−ε²); first moment O(N⁻²); Markov placement; O(1) clean low modes. ✓
- **Randomized-cutoff identities** (Thm 9.3) and the **heat-mixture no-go** (Prop 9.4): E_Θ T(Θ_δ) = 𝔄(δ) via survival function; the Cauchy–Schwarz CV bound Var Θ ≥ (EΘ)². ✓
- **Log-shell identity, localized bound, and pointwise maximalization** (Thm 9.5, Lem 9.7): interval-containment geometry correct. ✓
- **Krein crossing coordinates** (Prop 10.2): determinant lemma, Hellmann–Feynman speed m²/ṁ, and the cancellation h_e + t a_e m_e = h_e − a_e at t = −1/m_e. ✓
- **Dynamical-defect rates** (Cor 11.2): I redid the x = ru substitution; Ψ_d(r) = r^{d/2} (d<4), r²log(2/r) (d=4), r² (d>4) is exactly right, including the boundary term giving the +Cr². ✓
- **Doob/noise-semigroup representation** (eq. 10.31): standard Dirichlet-form identity for the product resampling semigroup. ✓
- **Gaussian half-space model** (eqs. 10.33–10.35): P_tF = qΦ(−a_t Z) with a_t = e^{−t}/√(1−e^{−2t}) is the standard OU computation. ✓

### Reproduction of the numerical audit

I ran two of the shipped scripts in a fresh Linux environment (Python 3 + numpy/scipy):

- `verify_abel_tail_identities.py`: **PASS** — Abel/heat identity to 4.3e−19, energy-excess to 5.5e−18, Duhamel derivative to 2.2e−19, both samplewise envelopes hold. Matches the shipped certificate.
- `verify_weighted_qgn_hodge.py`: **PASS** overall — all exact many-body/graph identities to ≤3e−15, gauge invariance, ellipticity, ring law, and the counterexample limit 0.0632503 all reproduce exactly. **However**, the random-torus Monte Carlo scan does *not* reproduce the paper's table: fresh run gives (L=4) mean 1.21287, std 0.1022 vs. the paper's 1.19523, std 0.1071. The deterministic checks are bit-stable; the MC draws are environment-dependent. See Minor comment M6.

---

## 2. Major comments

**MC1. 𝒦 is never defined.** The dynamical kernel 𝒦_{G,n}^J(iζ; a) — the object of Theorem 5.1, the counterexample, the Kohn–Abelian corollaries, and half the abstract — has no definition anywhere in the paper. The proof of Theorem 5.1 invokes "the positive-square dynamical least-squares formula" without statement, and Appendix A covers only the static curvature convention. A reader who does not have the companion zero-frequency paper cannot even state Theorem 5.1 precisely. This is the single most important fix: add the definition of 𝒦 and a derivation (or precise statement with proof) of the dynamical least-squares identity to Appendix A. As written the paper is fully self-contained *except* for this.

**MC2. Lemma "Periodic gradient-propagator estimate" is not checkable as written.** This lemma is the crucial analytic input for Theorem 8.9 (quenched Abel variance), which in turn drives Corollaries 8.10–8.12 and the entire fluctuation layer. Problems with the current proof:

- The weight ω_{L,s}(y) is **never defined** (presumably a Gaussian or polynomial parabolic weight, but the reader must guess, and the choice matters for the claimed cancellation).
- The parabolic Green function G_L(s, J, y, 0) is **never defined** on the torus (periodization convention, mean-zero projection, etc.).
- The imported inputs are cited as "Theorem 2 under condition (69b)," "Theorem 3(b)," and "Corollary 8 in the long version" of Gloria–Neukamm–Otto, with the long version identified only as "MiS Preprint 3/2013." Referees will want the precise statements, with their exact hypotheses, reproduced as numbered imported theorems (as you already do for Theorem 7.1). In particular the claim that both imported estimates are *uniform in the period L* is the entire point, and right now it rests on one sentence.
- The exponent bookkeeping does close (I checked: (1+s)^{−d/2−1/2+d/(2r)} × (1+s)^{d/(2r′)} = (1+s)^{−1/2} using 1/r + 1/r′ = 1), but the weighted-Hölder step and the role of α r′ > d deserve a full display.

Recommendation: give this lemma its own appendix with the weight defined, the two imported statements quoted verbatim with references to specific equation/theorem numbers in the published versions, and the Hölder/cancellation computation displayed. This is where any serious referee will spend most of their effort; make it easy for them.

**MC3. The periodization scheme should be fixed once, precisely.** Theorem 7.1 assumes "a standard stationary periodization scheme for which stochastic homogenization holds," and Lemma 7.2's proof says "away from a vanishing boundary fraction, periodization agrees with the infinite environment." Both are fine morally, but the scheme (restriction of the infinite sample with wraparound identification? independent resampling of the torus?) affects (i) whether Birkhoff applies as stated, (ii) whether the wrap bonds spoil stationarity of the torus ensemble, and (iii) what "i.i.d. periodic ensemble" means in Theorem 8.5 (i.i.d. *on the torus* vs. periodized infinite i.i.d. — these differ). Define the scheme once in Section 7 and use it consistently; state explicitly that for the i.i.d. case the torus law is the product law on the L^d fundamental-domain bonds.

**MC4. Consider splitting the paper.** Sections 1–7 plus the annealed Section 8 material form a complete, self-contained, publishable unit with a crisp story: exact kernel → weighted Hodge → ellipticity → counterexample → homogenization → almost-sure Kohn–Abelian equality → uniform annealed tail. The quenched-fluctuation and sharp-cutoff program (the Abel variance, Born sector, three routes, rank-one/coarea geometry, acoustic obstruction, noise-semigroup target) is a second paper's worth of material addressed to a different (probabilist) audience, and it ends at an open problem. At ~65 dense pages with a 1.5-page abstract, the current draft asks a single referee to be expert in frustration-free lattice models, stochastic homogenization, and random-operator spectral statistics simultaneously. If you keep it unified, at minimum: add a table of contents, and compress the abstract to ~250 words (the current abstract is a second introduction; the "Claim hierarchy" already does that job well).

**MC5. Two statements slightly oversell their content — one-sentence fixes.**
- Corollary 8.12 ("Every-cutoff centered bound in the bounded-mode window"): the bound Var T ≤ C(1+m_eff)V_L comes from the *uncentered* second moment (Var ≤ E T²). At bounded m_eff the scales coincide, so the conclusion is fine, but a sentence noting that no actual centering is used (and that the bound is vacuous as a "centered" statement when m_eff → ∞) would prevent a referee from reading it as a concentration result.
- The abstract sentence "the same second moment already proves the full centered sharp-cutoff scale at every prescribed cutoff whenever this effective mode count stays bounded" is accurate but should carry the same caveat.

---

## 3. Minor comments

**M1. Uncited bibliography entries.** `GloriaOtto2011`, `deBuyerMourrat2015`, and `DuerinckxGloria2025` appear in the bibliography but are never cited in the text. Either cite them (de Buyer–Mourrat seems intended for the environment-process heat decay discussion; Duerinckx–Gloria 2025 for the acoustic/dispersive discussion near Prop 10.5) or delete them.

**M2. Notation collisions.**
- `d` is the spatial dimension throughout, but Prop 8.1 says "for some d > 0" as the decay exponent — shadowing the ambient dimension in the very proposition that couples them. Rename the exponent.
- `z_e` is a many-body state |z_e⟩ in Lemma 3.1 and a scalar edge-gradient z_e(t) in eq. (8.24). Rename one (e.g. ζ_e(t) is taken; w_e(t) works).
- `q` is overloaded four ways: charge operator q_x, failure probability q, crossing jump q_e(r), participation weight q_k(J). At least rename the probability (use p or θ).

**M3. "Kohn–Abelian equality" is used repeatedly (including in two theorem names) but never defined or referenced.** One sentence at first use — and a citation for the terminology — would help readers outside the immediate subfield.

**M4. Theorem 9.1 (relative Wegner), proof.** The step from the functional-derivative identity to the conditional one-bond integral silently uses that s ↦ Tr F(L_{J^{(e,s)}}) is nondecreasing (true since F is nondecreasing and ∂_s L = b_e b_e* ⪰ 0, by monotonicity of eigenvalues); state it. Also "a smooth approximation removes the two irrelevant corner points" — say what is approximated (F, mollified at the corners, with monotone convergence of both sides).

**M5. Prop 10.2 (Krein coordinates).** Worth one sentence noting that a crossing at positive t requires m_e(r) < 0, i.e. r sits above at least one eigenvalue of the deleted-bond operator — this orients the reader on which spectral windows can produce crossings.

**M6. Numerical reproducibility.** The deterministic identity checks are bit-stable across environments, but the random-torus MC table in §12 does not reproduce from the shipped `verify_weighted_qgn_hodge.py` in a fresh environment (see §1 above). Since the paper explicitly frames these as "a representative seeded run," either (a) make the RNG fully deterministic (numpy `default_rng` with fixed per-case seeds, and record the numpy version in the manifest), or (b) add a sentence that the MC rows are environment-dependent while the identity rows are exact. Given how much care the archive otherwise takes with CHECKSUMS and certificates, (a) is worth doing.

**M7. Forward reference from front matter.** Claim-hierarchy item 5 cites \eqref{eq:integrated-noise-energy-goal}, an equation ~2,400 lines later. Consider stating the open problem's equation number in words ("the fixed-cutoff disorder-noise estimate of Problem 10.7") instead.

**M8. §8 assumes d ≥ 2 while all counterexamples live on rings (d = 1).** This is fine and even elegant, but say explicitly early in §8 that the quantitative theory is d ≥ 2 and that d = 1 is exactly solvable and used only for obstructions.

**M9. Prop 8.1 hypothesis.** H_L(t) ≤ C(1+t)^{−1−d/2} with F_L(0+) = 0 is used in the integration by parts (F_L(λ)/λ → 0); this is automatic from the bound but worth one clause.

**M10. Abstract/conclusion consistency.** Both are consistent with the claim list — I checked each abstract claim against a theorem. The only mismatch of tone: the abstract's "This removes the earlier multiple-crossing ambiguity" refers to a prior package version; a fresh reader has no "earlier" to compare against. Rephrase self-containedly ("a one-bond path crosses a prescribed cutoff at most once").

---

## 4. Assessment of the claim boundary

I checked CLAIM_STATUS.md against the paper. The proved/imported/open division is accurate; in particular:

- Items 1–11 (exact finite-graph layer): fully proved in the paper; I verified the key algebra independently. Agreed.
- Item 12 (uniform finite-torus tail): correctly labeled as resting on the imported GNO periodic divergence-semigroup estimate; the derivation from it is complete modulo MC2.
- Items 15, 29–38: the derivations are correct *given* Lemma "periodic gradient-propagator"; the lemma itself is the weak link (MC2).
- The "Not yet claimed" list is candid and matches the text, including the honest acknowledgment that the acoustic-CLT statement is a surrogate, not a conductance-model counterexample.

The one place the labels could mislead: CLAIM_STATUS item 32 calls the bounded-mode result a "centered theorem" — same caveat as MC5.

## 5. Recommendation

The exact layer (kernel, weighted Hodge, ellipticity, ring law, dynamical filter, counterexample, ergodic tightness) is correct and ready. Before circulating beyond a research-draft audience: fix MC1 (define 𝒦 — this is essential), rebuild the gradient-propagator proof per MC2, pin down the periodization scheme (MC3), and decide on the split (MC4). The minor list is mechanical. With MC1–MC3 addressed I would consider the annealed half of this paper solid enough for submission; the quenched half is publishable as a "program + partial results + precisely isolated open problem" paper, and Problem 10.7 (the fixed-cutoff disorder-noise energy bound) is a well-posed, attractive target.
