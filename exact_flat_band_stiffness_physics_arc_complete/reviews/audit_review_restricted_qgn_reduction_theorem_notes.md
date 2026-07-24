# Pre-expert audit: "A Restricted Reduction Theorem for Hermitian QGN Superconductors"

Reviewed: `restricted_qgn_reduction_theorem_package` (July 23, 2026) — manuscript source, algebra certificate, archived outputs, and the search memo — against the three prior packages and the GHK/QGN literature. This review is written for the stated purpose: hardening the draft before it goes to a PhD-level expert for (i) the least-squares/intertwiner proof and (ii) the H1–H4 audit of GHK's intended models.

## Verdict in one paragraph

I re-derived every lemma and theorem in the draft by hand and found **no mathematical errors**. Lemma 1 (least-squares curvature), Lemma 2 (the S⊕(S−1) decomposition including the norm ratios (23)–(24)), Theorem 1's intertwiner argument (in particular the two load-bearing steps: (30) — the longitudinal source lies in ker D† and escapes screening — and (31) — equivariance of the kernel projector with the ladder isometries), the translation-invariant Γ formula (36), the roots-of-unity filter, the locality bounds, and the normalization chain (50)–(58) all check out. Beyond reproducing their certificate (tensor error 1.25×10⁻¹⁸), I subjected Theorem 1 to two independent tests it passed: a **new stress-test model outside every prior test set** (non-UPC, asymmetric, time-reversal flat band; identity verified with Γ ≠ 0 to the finite-difference floor), and an **analytic cross-package check** (the draft's Γ formula reproduces the Phase-I PSD defect exactly). I also found a structural simplification the draft misses that materially changes the shape of the requested audit: **H1–H3 hold automatically for every projected PSD-Hubbard model with time-reversal-partner bands — UPC is not needed** — so the model-by-model audit of GHK's systems reduces to H4 (simple gapped AGP branch) plus the lift/locality assumption. The main conditional of the whole program is H4, and it is genuinely open at general size.

---

## 1. Line-by-line proof verification

### Lemma 1 (§3) — correct

Re-derived from second-order perturbation theory: with H = ½D†D and D|ψ⟩ = 0, the diamagnetic term is Re⟨b_i,b_j⟩, the ∂²D terms vanish on ψ, the first-order equation is exactly the normal equation (15), and the identity ‖Πb‖² = ‖b‖² + Re⟨b, Dχ⟩ (via ⟨Dχ, Dχ+b⟩ = ⟨χ, D†(Dχ+b)⟩ = 0) reconciles (13) with the Kohn dia-minus-para form. Solvability of (15) needs ⟨ψ|H′|ψ⟩ = 0, which holds since D|ψ⟩ = 0; the minimizer in (14) can be taken ⊥ψ because Dψ = 0. Two writing points for the expert version: state the solvability/uniqueness of χ_i on ψ⊥ explicitly, and note that multi-parameter analyticity of the simple branch follows from Kato (cited) — currently both are implicit.

### Lemma 2 (§4) — correct

All components verified: β_B = Tr b/L from the one-pair RDM JJ†/L = I/L (this is exactly where skew-unitarity/equal singular values enters — the draft's remark after (24) is accurate); the commutator-pull-through (21) is valid because quadratic pure-creation operators commute; every two-particle state ⊥ η⁺|0⟩ is annihilated by η⁻ (its image is ∝⟨η⁺0|ζ⟩|0⟩ = 0) and is therefore lowest weight of exactly S−1; the norms (23), (24) follow the standard k!(2S)!/(2S−k)! pattern with 2S = L and L−2; combining them gives the transverse norm ratio precisely ρ_n = n(L−n)/(L−1) — I recomputed this chain and it is exact. One rigor gap for the expert: the isometry U_n must be argued to be a *single* map on the whole S−1 lowest-weight multiplicity space (i.e., ζ-independent), which follows because (24) holds with the same constant for every lowest-weight vector — one sentence, currently omitted.

### Theorem 1 (§5) — correct; the two key steps hold

- **(30), longitudinal unscreened:** D_n†(|n⟩⊗v) = Σv_λS_λ†|n⟩ = Σv_λS_λ|n⟩ = 0 by Hermiticity + H2. So the longitudinal source sits inside ker D†, Π acts on it as the identity, and — the physical crux — the paramagnetic correction cannot touch it because H′|ψ⟩ = ½ΣS_λ(Ṡ_λ|ψ⟩) kills the component of Ṡ_λ|ψ⟩ along the ground manifold. Verified.
- **(31), equivariance:** D is an SU(2) intertwiner by H1, hence ker D† is invariant and Π commutes with the pseudospin action, hence with normalized ladder maps between weight spaces of the same isotypic component. Correct, but this is the step the expert will want expanded: write the target as ⊕ of isotypic components (S−1)⊗(multiplicity), show Π acts as id⊗P_mult there, and that U_n^⊕ acts as (ladder)⊗id. That's three lines of Schur and currently one sentence of prose.
- Sector-orthogonality of the two components after Π (needed for the cross terms to drop) follows from the same equivariance; worth one explicit sentence.
- The final assembly C⁽ⁿ⁾ = n²Γ + ρ_n(C⁽¹⁾ − Γ) and positivity (28) via n² − ρ_n = Ln(n−1)/(L−1) ≥ 0: verified.

### §§6–8 (translation invariance, winding, locality) — correct

(36) re-derived (the V factor and V̂(0) contraction are right; the "real Hermitian square root" caveat for keeping factors Hermitian is correctly flagged as part of the hypothesis). Proposition 2 is the standard multidimensional aliasing filter. Theorem 2's bounds check out, including the |ℓ_i|-weighted geometric sums; the polynomial-times-e^{−2μL_min} propagation into (58) works because the correction prefactor grows only polynomially in V while Γ decays exponentially, and (59) keeps the denominator alive. §8.1's gauge-invariant-trace argument is a clean and correct way to make the estimate topology-blind; citing Brouder et al. for why the *projected symbol* is only exponentially (not finitely) localized is apt, and the accompanying warning is important and correct.

### §9–10 (normalization and thermodynamic theorem) — correct with one presentational gap

The algebra of (54)–(58) checks out exactly (I verified (58)'s coefficient L²(n−1)/[(L−1)(L−n)] independently). The gap: **(53), C⁽¹⁾ = 4m_pair⁻¹, is presented as a convention but is actually a small lemma** — it needs E₁(A) = E_pair(Q₀+2A) along the analytic branch, which holds for translation-covariant lifts (both electrons shift by A) but should be stated as such with the assumption visible. An expert will poke exactly here.

---

## 2. Independent verification evidence (new in this review)

1. **Certificate reproduced.** Same seed: max tensor-identity error 1.25×10⁻¹⁸, fd cross-check 4.5×10⁻⁸, all assertions pass. The certificate design is good: random Haar-conjugated J (genuinely multiband), the *complete* 36-dimensional pseudospin commutant as square factors, generic one-body twist sources with no structure beyond H3 — so it tests the theorem's hypotheses and nothing narrower. (An attempted L=5 run exceeded my time budget — the dense 2^{2L} Fock construction scales badly; consider a fixed-number-sector implementation so referees can run L=5–6 cheaply.)
2. **New out-of-sample stress test.** I built (independent code, momentum-space route) a model in the theorem class that appears in none of the four packages: two-orbital TRS flat band with u↑(k) = (cos θ, sin θ), θ = 0.3 + 0.8cos k + 0.4sin 2k — deliberately **non-UPC** (orbital weights non-uniform) and k-asymmetric so Γ ≠ 0. Results at M=6: AGP is the exact zero-energy unique ground state in every sector (‖H|AGP⟩‖ ~ 10⁻¹⁶, gap 0.13); Γ from independent trace-symbol derivatives = 1.5758×10⁻³; predicted defects (n²−ρ_n)Γ = 3.781884×10⁻³ (n=2) and 1.134565×10⁻² (n=3) match the measured C_n − ρ_nC₁ to residuals ~2×10⁻⁹ (fd floor). The identity holds with nonzero Γ on a model the authors never saw.
3. **Analytic cross-package closure.** Applying the draft's (36)+(41) to the Phase-I winding model in the PSD convention: f₁ = −sin(2c sin qMk), f₂ = +sin(2c sin qMk), grid averages alias to ∓(1/M)sin(2c sin qMA), giving β = ∓2cq and Γ = V·U·Σβ² = 8Uc²q²M — which reproduces the Phase-I PSD defect 8Uc²q²M²n(n−1)/(M−1) = (n²−ρ_n)Γ **exactly**. The theorem's Γ formula quantitatively contains the original counterexample. Strongly recommend adding this one-paragraph check to the draft; it closes the loop across all four packages.

---

## 3. A structural finding that simplifies the requested audit

The draft treats H1 as a hypothesis to be checked model-by-model. In fact, for the projected PSD Hubbard interaction it is **automatic given time reversal**. Short derivation (worth adding as a lemma): for S = n̄_{Rα↑} − n̄_{Rα↓} with band-basis kernels b↑, b↓ and η⁺ = Σ_k d†_{k↑}d†_{−k↓},

[S, η⁺] = Σ_{pq} ( b↑_{p,−q} − b↓_{q,−p} ) d†_{p↑}d†_{q↓},

and for projected onsite densities b↓_{q,−p} = ū↓(q)u↓(−p)·(phase) = u↑(−q)ū↑(p)·(phase) = b↑_{p,−q} when u↓(k) = u↑*(−k). So H1 holds for **any** TRS pair of bands; H2 then follows from H1 + S_λ|0⟩ = 0 (so H2 need not be a separate hypothesis — a simplification worth making); H3 is immediate (∂_A of a projected density is a pure one-body operator); J from time reversal is automatically skew-unitary with equal singular values. My non-UPC numerical test above confirms all of this concretely (AGP exactly annihilated without UPC).

Consequences for the audit:

- **The H1–H3 half of the requested "line-by-line verification for GHK models" is a two-line corollary**, not a model-by-model computation. UPC appears nowhere in H1–H4; the theorem class is strictly larger than UPC-Hubbard. (UPC matters for GHK's *evaluation* of m_pair via the minimal quantum metric, not for the reduction identity itself.)
- **The entire model-dependent burden concentrates in H4** (simple gapped analytic AGP branch at every size) **and the lift/locality assumption (38)**. For GHK Models I and II, H4 is verified numerically at 2×2, 2×3, 3×3 (unique ground states, gaps 0.25–0.5, from the search package I audited yesterday), but there is no proof at general size — and Phase II is the standing demonstration that H4 can fail (reducibility). The natural companion result is a connectivity ⇒ uniqueness lemma in the Lieb / Mielke–Tasaki tradition; without it, Theorem 3 applied to any concrete GHK model is conditional on H4 at every size along the sequence.

I recommend the expert be pointed at three things in priority order: the isotypic-decomposition writing of (31); the H4 question for Models I/II (is there a Perron–Frobenius/connectivity route?); and the status of the non-Hermitian generalized-nesting model that violated the law by 0.53 in the search package — the draft excludes non-Hermitian constructions in §11 but does not cite that concrete example, which is the *evidence* that the Hermiticity hypothesis is essential rather than convenient. Its QGN-membership status (per the PRX definitions) should be settled in writing.

---
