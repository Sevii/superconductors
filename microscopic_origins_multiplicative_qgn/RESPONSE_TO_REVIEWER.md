# Response to the additional review

We thank the reviewer for a careful and constructive report. The manuscript has been revised throughout. The exact local results are unchanged; the main additions clarify the tuning cost of the screened family, make the one-body obstruction self-contained, place the control construction in the perturbative-gadget literature, provide an explicit time-reversal-invariant Kramers implementation, and strengthen the weighted local Schrieffer-Wolff argument.

## Major comments

### 1. Fine-tuning sensitivity of the screened construction

A new subsection, **“Mismatch tolerances and the generic unscreened regime,”** introduces physical coefficients

\[
\beta_q^{\rm phys}=\beta_q+\delta\beta_q
\]

and residual monopole interactions \(\delta\mathcal C_{q,L}\). It states the residual low-block interaction explicitly and gives a sufficient local-norm tolerance:

\[
\|\mathcal E_4\|_\mu\le \varepsilon\lambda^2\bar\alpha,
\qquad
\|\mathcal E_6\|_\mu\le \varepsilon\bar\alpha.
\]

The manuscript then spells out the componentwise consequences

\[
\delta\beta_4,\delta A_4=O(\varepsilon\lambda^2\bar\alpha),
\qquad
\delta\beta_6,\delta A_6=O(\varepsilon\bar\alpha),
\]

up to fixed local-norm factors. It now says explicitly that exact asymptotic composition degeneracy is a tuned submanifold, while any fixed design value of \(\lambda\) has a nonzero tolerance tube. The unscreened family is presented more prominently as the generic, robust selector-plus-mixing model rather than as a failed screened parent.

### 2. Dependence on unpublished companion manuscripts

A new self-contained appendix, **“Self-contained additive one-body Peierls obstruction,”** gives:

- the composition-curvature definition;
- the one-pair rigidity lemma;
- the five precise hypotheses of the no-go theorem;
- the block-diagonal conclusion for the exact curvature map; and
- the proof using common Kato transport and frustration-free source projection.

The appendix also identifies exactly which hypothesis is changed by the present construction: \(K_x=Q_{a,x}B_x\) is multiplicative and many body rather than additive and one body. The companion manuscripts remain cited for context, but the logical boundary used by this paper can now be checked without access to them.

### 3. Relation to perturbative gadgets

The introduction and discussion now cite and compare the construction with the perturbative-gadget literature of Kempe-Kitaev-Regev, Oliveira-Terhal, and Jordan-Farhi. The revised text distinguishes the present result from a generic mediator gadget in three ways:

1. the parity router is an exact all-filling operator identity before perturbative elimination;
2. every microscopic monomial receives the complete Peierls lift; and
3. connected overlap is controlled in a volume-uniform local interaction norm.

The manuscript does not claim a Hamiltonian-complexity reduction or optimal gadget scaling.

### 4. Explicit time-reversal-invariant Kramers embedding

The control construction is now made explicit in the main text and in a new appendix. For Kramers partners \(\kappa=\pm\), time reversal exchanges the partners, while the active transfer \(B\) is even and the active currents \(C_s,C_a\) are odd. Opposite signs on the partner router vertices therefore make the microscopic Hamiltonian time-reversal invariant.

The appendix proves that:

- with one control electron, the two partner sectors have identical high operators and identical Schur coefficients;
- a unique local time-reversal-invariant control ground state is obtained by filling the \(g_+g_-\) Kramers shell;
- scaling each partner bridge by \(1/\sqrt2\) reproduces the single-flavor degree-four and degree-six coefficients exactly; and
- processes exciting both Kramers electrons require four bridge vertices and begin at weighted degree eight.

A new verifier, `check_kramers_closed_shell_embedding.py`, constructs the exact closed-shell one-bond Feshbach map and confirms an \(O(\lambda^8)\) difference from the single-flavor target. The manuscript also clarifies that

\[
\frac{s^2}{\bar\Delta_s}>\frac{a^2}{\bar\Delta_a}
\]

is an open microscopic inequality, not a symmetry identity. A plaquette-generated even/odd gap ordering and unequal orbital overlaps provide plausible mechanisms, but the inequality is retained as an explicit parameter condition.

### 5. Weighted regrading and locality

The revised global section defines weighted degree at first use. The finite-order Schrieffer-Wolff subsection now writes

\[
V(\lambda)=\lambda V_1+\lambda^2V_2+\lambda^4V_4+\lambda^6V_6
\]

and shows that, for \(|\lambda|<1\),

\[
\|V(\lambda)\|_\mu\le |\lambda|\sum_{q\in\{1,2,4,6\}}\|V_q\|_\mu.
\]

The weighted appendix explains the auxiliary multivariate expansion, coefficientwise regrouping by total \(\lambda\)-weight, and the geometric tail estimate giving an \(O(\lambda^8)\) remainder after truncation through weight seven. It also introduces a control-excitation parity: only bridge vertices change that parity, so a connected two-control process requires two bridges on each control. Control-diagonal shell, screening, and base terms can enlarge support but cannot reduce the four-bridge, weight-eight threshold.

## Minor comments

- The odd simultaneous-parity projector is denoted \(Q_{a,x}\) consistently.
- The abstract was reduced to approximately 255 words and qualifies minimality as minimality in the no-go theorem's hypothesis space.
- The elementary second-order sign statement is now a lemma.
- The proof of the sharp bound states explicitly why the filled bonding or antibonding shell has positive swap sign.
- Weighted degree is defined before its first use.
- The normal-ordering equation now specifies the active-before-control mode convention and explains the convention-dependent sign.
- The certificate table places the finite-\(t\) plaquette value beside the analytic limit.
- Raw floating-point logs are retained, while `PLATFORM_STABLE_CERTIFICATE_SUMMARY.txt` reports threshold-based PASS results and rounded coefficients for cross-platform comparison.
- Peotta-Torma, the original Schrieffer-Wolff paper, and Kato are now cited in the body.
- The periodic-control figure marks an active cell shared by neighboring controls.
- “Minimal” is consistently restricted to the structural hypothesis changed in the one-body theorem.

## Verification

All five verifier programs pass. The LaTeX source compiles to a 27-page PDF with no undefined references, undefined citations, or overfull boxes. The final PDF was rendered page by page and visually inspected.
