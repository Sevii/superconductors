# Revision notes: Exact Flat-Band Stiffness from Pair Mobility in a QGN Model

## Mathematical corrections applied

1. **Finite-volume pair momentum was redefined without an invalid continuous momentum-block identification.** The revised paper defines the continuous lifted family
   \[
   H_{\mathrm{pair},\mathbf N}^{\mathrm{lift}}(Q)
   :=H_1(Q/2)|_{K=0}
   \]
   and its Hessian. It explicitly states that an arbitrary continuous \(Q\) need not be an ordinary untwisted finite-torus momentum block because the momentum grid is discrete. The exact chain-rule identity is now
   \[
   C_{ij}^{(1)}=4(m_{\mathrm{pair},\mathbf N}^{-1})_{ij}.
   \]

2. **The common Fock-space hypothesis was made explicit.** All fixed-number Hamiltonians are now stated to be restrictions of one analytic, particle-number-preserving Fock-space family with the same square factors. This supplies the hypothesis used when forming \(D=\bigoplus_n D_n\).

3. **The Hermitian positive-kernel hypothesis was formalized.** The translation kernel is now assumed real symmetric and positive semidefinite on the combined real-space/channel index, so it admits a real square root and produces Hermitian square factors. The finite-torus Fourier convention is written explicitly, fixing the factor of \(V\) in the longitudinal tensor.

4. **The connectivity proof was expanded.** The spanning-tree induction now displays the commutators that generate the antisymmetric and symmetric cross couplings and the new diagonal rank-one operator. The paper also clarifies that \(p_x=|v_x\rangle\langle v_x|\) is a rank-one positive operator, not generally an orthogonal projector.

5. **Lemma 2's valid range was stated explicitly:** \(L\ge2\) and \(1\le n\le L-1\).

6. **The locality constants were corrected.** The constants \(C_1,C_2\) in the exponential-locality estimate are now stated to depend on \(C_0\) as well as \(d,\mu\), the channel number, and \(V_0\).

7. **The exact Model-II cancellation \(\Gamma=0\) was proved by a displayed trace calculation.** The revised proof introduces the finite-grid orbital weight \(w_a(A)\), uses time reversal to obtain the down-spin trace, and shows directly that pointwise constant orbital weights force every longitudinal derivative to vanish.

8. **The counterexample statements were made reproducible and precise.** The four-cell winding example now states \(M=4\), \(q=1\), \(c=1/8\), \(\phi(k)=(\pi/2)\sin k\), and \(n=2\), deriving the raw defect \(2g/3\) and normalized defect \(g/24\). The reducible example is described as the analytic one-pair branch through \(Q=0\), without making an unnecessary global lowest-branch claim.

9. **The reproducibility section was synchronized with the supplied archive.** Unsupported wording was removed and the reported numerical residuals were replaced by the supplied or rerun certificate values.

## Independent checks rerun

- Exact four-cell SymPy certificate: raw discrepancy \(2/3\) and normalized discrepancy \(1/24\) for \(g=1\).
- Resonant/nonresonant winding implementation: the \(M=8\), harmonic-four control has defect \(3.33\times10^{-9}\), while the resonant cases reproduce their analytic defects.
- Restricted reduction certificate: maximum tensor-identity error \(1.2487\times10^{-18}\); maximum direct finite-difference error \(5.1744\times10^{-8}\).
- H4 certificate: full Lie dimensions \(16\) and \(25\) for the connected \(L=4,5\) tests, exact connected/disconnected nullity counts, and all 55 audited Model-I/II frame graphs connected.
- Model-II certificate: direct lifted pair matrix error \(4.259\times10^{-16}\), electronic-twist identity error \(5.881\times10^{-16}\), finite-grid mass error \(3.558\times10^{-9}\), and many-body curvature error \(2.496\times10^{-9}\).
- Reducible fixed-local certificate: the half-filled even-torus state remains an exact zero mode for tested twists while the one-pair inverse mass is positive.
- Sector-gap evidence: 160 retained rank-one cases agree across fillings to worst relative deviation \(1.019\times10^{-13}\); the generic non-rank-one Hermitian control has ratio \(0.432639\).

The standalone longitudinal-defect script was not rerun because its imported module `qgn_search_v3` and input JSON file were not included in the supplied archive. The revised Model-II \(\Gamma=0\) result is proved analytically in the paper and does not depend on that script.

## Build verification

The revised Markdown compiles to a 20-page PDF. The PDF was rendered and visually inspected, passed the PDF preflight check, and produced no overfull or underfull box warnings in the final LaTeX build. The headline coefficient and finite-size stiffness formula are unchanged.
