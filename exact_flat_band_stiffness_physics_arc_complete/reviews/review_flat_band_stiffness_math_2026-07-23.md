# Review: "Exact Flat-Band Stiffness from Pair Mobility in a QGN Model"

**Reviewer:** Claude (independent math review, July 23, 2026)
**Version reviewed:** `exact_flat_band_stiffness_physics_arc` (md/docx/pdf, package of 2026-07-23)

## Summary verdict

The mathematics is sound. I hand-verified every central derivation, re-ran the package
certificates in a fresh environment (all pass), and — most importantly — wrote a fully
independent exact-diagonalization implementation of projected Model II from scratch,
sharing no code with the package. On the 3×3 torus, for ξ ∈ {0.7, 1.0} and n ∈ {1, 2},
it reproduces E(0) = 0 to 1e-16 and the headline curvature
C^(n) = [n(V−n)/(V−1)]·|U|ξ²/2 to ~4e-10 (the finite-difference floor), including the
V/(V−1) factor. I found no errors in any theorem, lemma, or coefficient. The issues below
are about completeness of two asserted claims, and about notation/presentation.

## What was checked and confirmed

**Hand-verified algebra.** Lemma 1 (least-squares curvature; standard second-order
perturbation theory for H = D†D/2 with D|ψ⟩ = 0, correctly reduced to the projected
residual). Lemma 2, including the norm-ratio computation: I independently reproduced
n²·N_n²·[(n−1)!(L−2)!/(L−n−1)!]·L = ρ_n from (23)–(24), the coefficient β_B = Tr b/L
(the 1-RDM of |1⟩ is I/L because J is skew-unitary), the fact that [B, η⁺] commutes with
η⁺ (pure creation bilinear), and the ΔS ≤ 1 selection rule (B is a rank ≤ 1 pseudospin
tensor). Theorem 1: the longitudinal source lies in ker D_n† exactly as claimed in (30)
(this is where Hermiticity of the square factors is load-bearing — correctly flagged),
the intertwining relation (31), and the coefficient algebra, including the n = 1
self-consistency check C^(1) = Γ + (C^(1) − Γ) and n² − ρ_n = nL(n−1)/(L−1) ≥ 0, which
gives (28). Equation (10c): I verified [dΓ↑(a) − dΓ↓(ā), η⁺] = a J − J āᵀ vanishes
exactly when a is Hermitian in the TR-conjugate pairing convention. Proposition 1: the
vectorization argument, the reduction of the commutant via irreducibility of Λⁿ under
U(V), and the zero-mode generating function (10j) all check out. Proposition 2 and
Theorem 2 (roots-of-unity filter and locality bounds, including the role of bounded
aspect ratio in Theorem 3's exponential case — that hypothesis is genuinely needed and
correctly stated). Lemma 3 (C^(1) = 4 m⁻¹ from Q = 2A). Lemma 4: δα_k(Q) in (M8), the
grid moment identities (V⁻¹Σ sin k_i = 0 for all N_i; V⁻¹Σ sin²k_i = 1/2 requiring
N_i ≥ 3, which is exactly the 2×2 anomaly), the rank-two Gram matrix (M7), reality of F
near Q = 0 (k → −k symmetry kills the imaginary part, so |F| is smooth at F(0) = 1), the
Hessian |U|ξ²/8, and the Bessel factorization (M9) with its consistency against (M4).
Corollary 4.1 (constant orbital weights ⇒ the trace symbol is A-independent ⇒ Γ = 0
exactly). The normalization chain (54)–(58): I verified both the identity
ρ_n/V = (N_flat/2)ν(1−ν)L/[(L−1)V]·… behind (56) and the coefficient
L²(n−1)/[(L−1)(L−n)] in (58) from scratch. The winding-obstruction algebra (O4)–(O5),
including n − ρ_n = n(n−1)/(M−1). The §2.2 counterexample logic (component-polarized zero
mode, E_pair = (U/2)sin²Q ⇒ m⁻¹ = U, R_M(1/2) = 0).

**Certificates re-run (fresh environment).** `model_II_pair_dispersion_certificate.py`:
PASS, with errors matching §14's table (direct Q-block vs formula 4.3e-16, twist vs
Q = 2A 5.9e-16, mass Hessian 3.6e-9, Bessel limit 3.2e-16, many-body curvature 5.1e-9).
`restricted_qgn_reduction_certificate.py` (L = 4, random multiband skew-unitary J): PASS,
tensor-identity residual 1.25e-18, matching the "< 1.3e-18" claim.

**Independent reproduction.** A from-scratch ED of the projected PSD-Hubbard Model II
(own conventions, own second quantization, no shared code) confirms the headline theorem
numerically; see summary verdict above. As a side observation, my first attempt had a
conjugation error in the spin-down projected orbitals and produced E(0) ≠ 0 — the
frustration-free zero mode is a sharp built-in test, which speaks well for the
robustness of the claimed structure.

**External checks.** Ref. [1] exists as cited, including the nonstandard-looking DOI
10.1103/gw85-5r92, which resolves at APS to PRL "Bootstrapping Flatband Superconductors:
Rigorous Lower Bounds on Superfluid Stiffness." The conjecture as stated in (1) matches
the source, with the same N_flat/2 and ν conventions. The fixed-layout docx (page images)
renders correctly and matches the markdown source.

## Findings

### F1 (main gap): the all-size connectivity claim for 0 < ξ < π/2 is asserted, not proven

Theorem 4's unconditional range rests on the §3.2 sentence: "For 0 < ξ < π/2, its
finite-grid Fourier coefficients at displacements 0, e_x, and e_y are nonzero on every
rectangular torus." No proof or citation is given in the paper, and this is the one
load-bearing step a referee cannot check from the text. The claim is true — I scanned the
1D alias coefficients c_r^(N)(ξ) = N⁻¹Σ_k e^{iξcos k}e^{−ikr} over N = 3…60 and
ξ ∈ (0, π/2) and found min |c_0| ≈ 0.471 (attained at ξ → π/2, N = 6) and c_1 ≈ ξ/2 for
small ξ, bounded away from zero — and it is easily provable: c_r^(N) is the alias sum
Σ_ℓ i^{r+ℓN} J_{r+ℓN}(ξ), and for ξ < π/2 the tail is dominated by
|J_n(ξ)| ≤ (ξ/2)ⁿ/n! with ξ/2 < π/4 < 1, so the leading Bessel term wins for every
N ≥ 3. Recommendation: display this two-line inequality (or cite the package note that
contains it), and add one sentence deriving graph connectivity from nonvanishing at
displacements {0, e_x, e_y} (the diagonal projector blocks of Model II are constant
= 1/2, so all edges run between the two orbital layers; links at 0 and e_x, e_y make the
bipartite graph connected). The parallel Model I assertion in the same paragraph
("nonzero nearest-neighbor diagonal coefficients on every isolated N×N torus") deserves
the same treatment. Relatedly, my scan shows why the restriction matters: near
ξ ≈ 2.405 (first zero of J₀), c_0 nearly vanishes at large N, so the certificate-only
status of ξ > π/2 cases is the right call.

### F2: boxed-but-unproven statements in §12

(61b) and (61c) are boxed, but (61c)'s Δ_{N,1} = O(N⁻²) is supported only by a one-line
trial-state sketch, and the representation-theoretic content of (61a)–(61b) is compressed
past the point of verifiability from the text (the objects w_x, q_x, r_x and the claim
that the conjugation representation contains the traceless one-pair representation are
not developed). Since §12 is explicitly non-load-bearing and deferred to a standalone
note, either unbox these or add a short appendix. Boxing elsewhere in the paper signals
"proven here," and these two are the exceptions.

### F3: H2 derivation wording in §3.1

"H2 follows because the factors annihilate the vacuum" — annihilating the vacuum alone
gives S_λ|0⟩ = 0; S_λ|n⟩ = 0 additionally needs [S_λ, η⁺] = 0, i.e. (10c)/H1. One-word
fix: "H2 follows from H1 because the factors annihilate the vacuum."

### F4: notation collisions

(a) L is the pair-orbital count in §3–6 and (58), but L_j, L_min, L_max are torus side
lengths in §7–9, and both appear in Theorem 3's statement. (b) V is a vector space in
§3.1–3.2 and the volume from §7 on (they coincide for Model II, which partially masks
the collision). (c) The filling is ν = n/V in the abstract/§1 and ν = n/L in (55) —
consistent for Model II but the reader meets both. (d) U conventions drift: §2.1 uses
attractive U with E_n(0) = −Un/2 < 0 (so the Hamiltonian there is not the PSD form used
everywhere else, and the sign convention of U in (O4)–(O5) is implicit), §2.2 uses
PSD U/2, and the Model II sections use |U|. A short notation paragraph, plus one
sentence in §2.1 fixing the sign convention, would remove all four ambiguities.

### F5: equation-numbering irregularities

(10) appears after (10a)–(10j); (61) in §11 is followed by (61a)–(61c) in §12; (1) and
(1a)–(1d) similarly interleave. Harmless but worth fixing in production, since cross
references like "(10)" currently point below "(10j)."

### F6: minor presentational points

(a) §2.1: the exact factorization (O3) is stated without derivation; one sentence (the
deformation angle becomes k-independent under the uniform shift, so the twisted
deformation acts as a global rotation of the pair basis) would make it self-contained,
with the certificate as backup. (b) Abstract and §1 say "In the connected Hermitian
positive-square class, the fixed-pair-number Hessian obeys (62)" — Theorem 1 needs only
H1–H4; connectivity is the sufficient condition for H4 in the Hubbard subclass. The
current phrasing slightly overstates what connectivity is needed for. (c) Lemma 4's
"lowest branch" claim implicitly uses (1−|F|)/4 ≤ 1/2 ≤ (all other sector energies)/|U|;
one clause would close it. (d) Ref. [3] (Kato) lacks year/edition; Ref. [11] lacks
volume/pages. (e) The docx is a fixed-layout page-image rendering — fine as a viewing
copy, but not usable for journal submission or text extraction; keep the md/pdf as the
canonical sources.

## Bottom line

The proof structure — obstructions → restricted reduction theorem → connectivity
criterion → exact Model II closure — is coherent, and every quantitative claim I tested
(by hand, by re-running certificates, and by independent reimplementation) is correct,
including the exact finite-size factor V/(V−1) and Γ = 0 for Model II. The single
substantive revision item is F1: put the alias-sum inequality behind the 0 < ξ < π/2
connectivity claim into the paper (or an explicit citation), since Theorem 4's
unconditional range depends on it. Everything else is polish.
