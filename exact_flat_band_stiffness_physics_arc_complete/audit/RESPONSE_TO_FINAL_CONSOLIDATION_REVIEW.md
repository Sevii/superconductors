# Response to Final Consolidation Review

**Manuscript:** *Exact Flat-Band Stiffness from Pair Mobility in a QGN Model*  
**Author:** Nicholas Sledgianowski  
**Date:** July 23, 2026

The review recommendations have been incorporated as follows.

## Four local corrections

1. **Ground branch in Lemma 3.** The proof now states that translation labels remain good under the flat-connection family and that H4 makes the zero-label ground state simple and isolated; its Kato continuation therefore remains in that block for sufficiently small twist.
2. **The `N_i >= 3` condition.** A dedicated remark explains that on a two-point cycle, `sin k_i` vanishes identically. On a `2 x 2` torus both mass components vanish, exactly accounting for the earlier zero-curvature anomaly.
3. **The role of `xi < pi/2`.** The pair-mass lemma is stated for every `xi != 0`. The interval `0 < xi < pi/2` is used only to prove connectivity on every rectangular torus. The finite-size theorem is stated for any additional value of `xi` whose projected-frame graph is independently certified as connected.
4. **Exact cancellation of the longitudinal term.** Equation (M11) is now Corollary 4.1 and is advertised explicitly as exact on every finite torus: there is no winding, locality, or finite-size correction to estimate.

## Consolidation and scope

- The winding obstruction, reducibility obstruction, general reduction theorem, H4 connectivity theorem, locality theorem, and Model-II pair-dispersion closure have been consolidated into one physics paper.
- The stronger filling-independent sector-gap question has been removed from the load-bearing proof. Only the proved comparison `Delta_{L,n} <= Delta_{L,1}` and the closing collective branch remain in the physics paper.
- The proposed rank-one Aldous problem is deferred to a standalone mathematical note; a separate outline is archived.
- The paper proves the conjecture in Gao-Han-Khalaf's literal thermodynamic flat-connection-curvature definition of `D_s`. A stricter dynamical zero-frequency identification is isolated in one future-work paragraph rather than added as an unnecessary hypothesis.

## Standing debts cleared before submission

- **Non-Hermitian boundary:** the paper distinguishes perfect QGN of the band geometry from the Hermitian positive-square interaction class used by the theorem. The published SSH-type non-Hermitian extension is identified as outside the present proof. The archived random non-Hermitian stress test is described only as evidence that Hermiticity is substantive, not as a counterexample to the published SSH construction.
- **Prior art:** a separate audit records the closest exact-ground-state, pair-mass, minimal-metric, QGN, stiffness-bound, and irreducibility results and states the novelty claim narrowly.
- **Alon-Puder:** the reference has been verified as Gil Alon and Doron Puder, *Aldous-type Spectral Gaps in Unitary Groups*, arXiv:2603.00353 (submitted February 27, 2026). Its proved families have not been identified with the rank-one Householder-reflection Laplacians here, so it is cited only in the deferred mathematics outline and scope audit.

## Headline statement after revision

For projected-Hubbard Model II, `0 < xi < pi/2`, rectangular tori with `N_x,N_y >= 3`, and every nonsingular filling,

\[
\kappa_{ij}^{(n)}(V)
=
\frac{V}{V-1}\nu(1-\nu)\frac{|U|\xi^2}{8}\delta_{ij},
\]

while

\[
(m_{\mathrm{pair}}^{-1})_{ij}
=
\frac{|U|\xi^2}{8}\delta_{ij}.
\]

Therefore, in the thermodynamic flat-connection definition used by Gao, Han, and Khalaf,

\[
D_{s,ij}
=
\nu(1-\nu)(m_{\mathrm{pair}}^{-1})_{ij}.
\]
