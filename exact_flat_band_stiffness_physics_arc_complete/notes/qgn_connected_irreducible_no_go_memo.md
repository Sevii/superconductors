# Search for a fixed-local, connected, irreducible QGN counterexample

**Status:** No counterexample found. The requested hypotheses instead lead to a proof-ready no-go identity showing that the thermodynamic ratio must approach one in the standard Hermitian QGN class.

**Date:** July 23, 2026

## 1. Target

We searched for a fixed local QGN Hamiltonian with:

1. one model held fixed as the torus grows;
2. a connected/irreducible projected one-particle graph;
3. a unique, gapped analytic AGP ground branch in each fixed-pair sector;
4. the Gao–Han–Khalaf flat-connection prescription `k -> k+A`;
5. Hermitian local QGN generators and a positive interaction kernel;
6. a thermodynamic stiffness-to-pair-mass ratio

\[
R_M(\nu)=
\frac{D_s(M,\nu)}{(N_{\rm flat}/2)\nu(1-\nu)m_{\rm pair}^{-1}(M)}
\]

converging to a value different from one.

No such example was found.

## 2. Computational search

### 2.1 Random connected rank-one QGN models

A search over 1,200 random fixed-local Hermitian-QGN Hamiltonians used four connected projectors and one or two independently weighted positive squares built from onsite, nearest-neighbor, and range-two neutral bilinears. Every retained model had a unique one-pair and half-filled ground state and a positive finite-size gap.

The largest canonical deviation at `M=6` was

\[
\left|\frac{E_n''}{\rho_n E_1''}-1\right|
=1.9310809\times10^{-3},
\qquad
\rho_n=\frac{n(M-n)}{M-1}.
\]

For the same candidate at `M=8`, it fell to

\[
1.2470577\times10^{-4}.
\]

Across the top 35 candidates, the median attenuation factor from `M=6` to `M=8` was approximately `0.0681`.

For the strongest candidate, converged finite-difference calculations gave:

| M | n | canonical deviation |
|---:|---:|---:|
| 5 | 2 | 3.61817e-2 |
| 6 | 3 | 1.93109e-3 |
| 7 | 3 | 2.92151e-4 |
| 8 | 4 | 1.24712e-4 |
| 9 | 4 | 2.56e-6 |

The sequence does not indicate a nonunit limiting ratio.

### 2.2 Published GHK model families

Independent projected-Hubbard exact diagonalization was run for the two electronic structures used by Gao, Han, and Khalaf (their Models I and II), using the literal `k -> k+A` prescription, at `2x2`, `2x3`, and `3x3` tori. Whenever the one-pair curvature was nonzero, the canonical ratio equaled one to roughly `1e-9`–`1e-10`. Cases with nearly zero pair curvature had meaningless relative ratios but absolute defects of order `1e-10`.

### 2.3 Wider controls

Additional searches found:

- 300 random full-spin/spin-mixing Hermitian-QGN models at `M=4`: maximum deviation `4.3e-9`;
- positive finite-range translation kernels: no deviation beyond finite-difference error;
- irreducibility-restoring perturbations of the reducible model: the canonical law returned;
- large deviations only after leaving the strict class, for example by using non-Hermitian generalized-nesting generators or nodal pairing form factors.

## 3. Exact finite-size defect identity

Let `L` be the number of pair orbitals, and define the normalized AGP states

\[
|n\rangle =
\left[\frac{(L-n)!}{n!L!}\right]^{1/2}
(\eta^+)^n|0\rangle,
\qquad
\eta^+=\sum_{a=1}^{L}c_{a\uparrow}^\dagger c_{a\downarrow}^\dagger.
\]

After diagonalizing the positive interaction kernel, a standard Hermitian QGN Hamiltonian can be written

\[
H(A)=\frac12\sum_\lambda S_\lambda(A)^2,
\]

where at zero twist

\[
[S_\lambda,\eta^+]=[S_\lambda,\eta^-]=[S_\lambda,\eta^z]=0,
\qquad S_\lambda|n\rangle=0.
\]

For a unique analytic branch, second-order frustration-free perturbation theory gives a least-squares problem for the sources `dot S_lambda |n>`. Since `dot S_lambda` is a number-conserving one-body operator, its action on the maximal-pseudospin state decomposes only into total pseudospin `S=L/2` and `S-1` components.

The maximal-spin component scales as `n`; the orthogonal component has the exact AGP norm ratio

\[
\rho_n=\frac{n(L-n)}{L-1}.
\]

Because every zero-twist `S_lambda` is a pseudospin scalar, the least-squares problem in the `S-1` sector is unitarily equivalent at every filling. Therefore

\[
\boxed{
E_n''(0)=
\rho_n E_1''(0)
+\bigl(n^2-\rho_n\bigr)\Gamma_L,
}
\]

where `Gamma_L >= 0` is the one-pair longitudinal source weight,

\[
\Gamma_L
=\sum_\lambda
\left|\langle1|\dot S_\lambda|1\rangle\right|^2
\]

up to the common Hamiltonian normalization.

A numerical certificate independently computes `Gamma_L` from trace derivatives and verifies the identity at every filling for `L=5,6,7`; the largest residual is `1.4e-8`, set by finite-difference accuracy.

## 4. Why locality removes the only defect

For a translation-invariant local model, the longitudinal matrix element is a discrete Brillouin-zone average of a shifted periodic function. If

\[
f(k)=\sum_{r\in\mathbb Z}\widehat f_r e^{irk},
\]

then

\[
\frac1L\sum_{m=0}^{L-1}f(k_m+A)
=\sum_{\ell\in\mathbb Z}\widehat f_{\ell L}e^{i\ell LA}.
\]

Thus `Gamma_L` is generated only by winding Fourier coefficients at multiples of the circumference.

- For a fixed strictly finite-range lift, these coefficients vanish exactly once the torus is larger than the Fourier support, so `Gamma_L=0` exactly.
- For a fixed uniformly exponentially local lift, they decay exponentially, so `Gamma_L -> 0` up to polynomial volume factors.

The earlier winding counterexample evades this statement by moving a Fourier coefficient to displacement `r=L` as the system size changes. The reducible counterexample evades the unique/irreducible-ground-branch hypothesis.

## 5. Consequence for the GHK ratio

Let

\[
L=(N_{\rm flat}/2)V,
\qquad \nu=n/L.
\]

When `Gamma_L=0`, the canonical finite-size identity is exact:

\[
E_n''=\frac{n(L-n)}{L-1}E_1''.
\]

Using `Q=2A`, `D_s=E_n''/(4V)`, and `m_{pair}^{-1}=E_1''/4`, the ratio used in the thermodynamic conjecture is

\[
R_L(\nu)=\frac{L}{L-1},
\]

so the leading finite-size correction is universal and

\[
\boxed{R_L(\nu)\longrightarrow1.}
\]

For exponentially local models, the additional longitudinal term also vanishes. Hence a fixed-local connected irreducible Hermitian-QGN model with a unique analytic AGP branch and `R_L -> c != 1` appears to be ruled out.

## 6. What remains to make this a theorem

The argument is proof-ready but should still be written with full operator-domain and multiband detail. The remaining obligations are:

1. state the second-order frustration-free least-squares lemma for a unique analytic ground branch;
2. formulate the `S`/`S-1` decomposition for the full multiband QGN pairing matrix;
3. prove the locality bound on `Gamma_L` uniformly in dimension and for topological projector bundles;
4. fix the exact electronic-twist/pair-momentum/stiffness normalization;
5. clarify whether every model intended by the GHK conjecture lies in the Hermitian-square class used here.

## 7. Research conclusion

The requested counterexample was not found because the requested conditions appear to remove both known failure mechanisms. The search has instead produced a candidate restricted theorem:

> **Restricted QGN reduction theorem.** For a fixed uniformly local Hermitian QGN model with a full-rank pairing matrix, a unique analytic AGP ground branch, and a specified Peierls lift, the many-pair twist curvature equals the one-pair curvature times `n(L-n)/(L-1)` plus a nonnegative winding term. The winding term vanishes for sufficiently large strictly finite-range tori and vanishes asymptotically for uniformly exponentially local models. Consequently the Gao–Han–Khalaf stiffness/pair-mass ratio tends to one.

This is a stronger outcome than another finite-size counterexample if the remaining proof details close.
