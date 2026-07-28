# Response to `REVIEW_qgn_quenched_robustness.md`

The reviewer found the exact finite-graph layer correct and identified self-containedness, periodization, organization, and reproducibility issues. The revision implements the recommended two-paper split and addresses every listed comment.

## Major comments

| Comment | Revision |
|---|---|
| MC1: the dynamical kernel was undefined | Paper I, Section 5 now defines the imaginary-frequency kernel from current and stress operators. Appendix A derives the positive-square identity `K(i zeta)=<S,zeta^2/(L^2+zeta^2)S>` and its static limit. |
| MC2: the periodic gradient-propagator lemma was not checkable | Paper II, Appendix A now defines `G_L`, `omega_{L,t}`, the horizontal generator, the vertical derivative, `SG_L(rho)`, and the source norm. It restates GNO Theorem 2 under (69b), Remark 11/equation (70), Theorem 3(b), and equation (74), with period-independent constants, then displays the weighted-Holder calculation and exact exponent cancellation. |
| MC3: periodization was not fixed | Paper I, Definition 7.1 distinguishes (i) space-periodization of one infinite realization, used for almost-sure statements, and (ii) the product law on periodic torus bonds, used for quantitative estimates. Paper II uses only the second ensemble and says so at the outset. |
| MC4: split the paper | Implemented. Paper I is the exact/static/annealed paper; Paper II is the quenched-fluctuation/sharp-cutoff paper. Both are 21 pages and have separate abstracts, conclusions, bibliographies, theorem boundaries, and verification manifests. |
| MC5: bounded-mode “centered” language oversold the mechanism | Paper II’s abstract, claim hierarchy, corollary, theorem-boundary section, and conclusion now state explicitly that the estimate uses only `Var T <= E T^2`, has the target scale only at bounded effective-mode count, and is not a mesoscopic concentration theorem. |

## Minor comments

| Comment | Revision |
|---|---|
| M1: uncited bibliography entries | Every retained entry is cited. Unneeded entries were removed from Paper II; the acoustic and environment-process references are now cited where used. |
| M2: notation collisions | The generic heat exponent is `sigma`; the parabolic edge gradient is `w_e(t)`; the probability parameter is `vartheta`. The many-body tangent state `z_e` remains confined to Paper I. |
| M3: define Kohn-Abelian equality | Paper I defines the Kohn and Abelian orders of limits immediately after the dynamical filter and cites Kohn and Scalapino-White-Zhang. |
| M4: Wegner proof details | Paper II states that `partial_s L=b_e b_e^* >= 0`, so all eigenvalues and `Tr F(L)` are monotone, and specifies monotone smooth approximation of the two ramp corners. |
| M5: orientation of the Krein crossing | Paper II states that an admissible positive crossing requires `m_e(r)<0`, hence the cutoff lies above at least one deleted-bond eigenvalue. |
| M6: Monte Carlo reproducibility | Paper I now uses independent deterministic seed sequences for exact checks and each lattice size, fixes 24 samples per size, and records Python/NumPy/SciPy versions. The regenerated certificate and CSV are included. |
| M7: long forward reference in front matter | Paper II’s claim hierarchy refers to the named fixed-cutoff disorder-noise problem, not a far-forward equation number. |
| M8: dimensional scope | Both papers state that quantitative theory is for `d>=2`; one-dimensional rings are used only as exactly solvable obstructions. |
| M9: lower boundary in Stieltjes integration | Paper I’s heat-to-tail proposition explicitly notes that the heat bound implies `F_L(lambda)/lambda -> 0`, eliminating the lower-end boundary term. |
| M10: package-history wording | Replaced with the self-contained statement that a one-bond path crosses a prescribed cutoff at most once. |

## Remaining theorem boundary

The split does not convert the still-open mesoscopic sharp-cutoff problem into a theorem. Paper II isolates the precise remaining estimate as the fixed-cutoff disorder-noise quadratic-variation bound. All statements depending on imported stochastic-homogenization estimates remain explicitly labeled as such.
