# T2.3 boundary-scope decision

## Decision

T2.3 should be split into two statements rather than left ambiguous.

- **T2.3a — periodic theorem boundary.** The microscopic downfolding theorem places one control on every translated bond of a periodic cycle. The exact cycle leakage proposition is therefore the theorem-level result required by the paper.
- **T2.3b — open numerical boundary.** The DMRG/MPS plan also calls for open-chain spectra. If the implementation forms the open system by deleting the wrap bond and retaining the path edges `1--2,...,(L-1)--L`, a separate leakage formula is required before those spectra are compared with the periodic selector law.

Thus an open-boundary theorem is **not needed to validate the stated periodic microscopic theorem**, but it **is needed to interpret the planned open-chain numerical data**. The package derives T2.3b for the standard path truncation used in the numerical plan.

## Why the cycle formula cannot be reused

On the cycle every active site lies on two charge bonds, so at fixed total pair number

\[
\sum_{x\in\mathbb Z_L}N_x=4(n_1+n_2)
\]

is a scalar. The linear term proportional to \(U_4\) therefore creates no leakage.

On the path the two endpoints have degree one. With

\[
B_a=X_{a,1}+X_{a,L},\qquad B=B_1+B_2,
\]

one instead has

\[
\sum_{x=1}^{L-1}N_x=4(n_1+n_2)-2B.
\]

Endpoint occupancy fluctuates in a product AGP, so the linear charge term participates in the open leakage. This is the structural reason the open result is not obtained by replacing \(L\) with \(L-1\) in the cycle formula.

## Exact open operator reduction

Let

\[
E_a=\sum_{x=1}^{L-1}X_{a,x}X_{a,x+1},\qquad
C=\sum_{x=1}^{L-1}(X_{1,x}+X_{1,x+1})(X_{2,x}+X_{2,x+1}),
\]

and \(Z=E_1+E_2+C\). Writing the physical crossed-shell parameter as

\[
V_4=4A_4-U_4>0,
\]

the path charge operator obeys the configuration-wise identity

\[
\boxed{
\mathcal C_{4,L}^{\rm op}
=(8A_4-2U_4)n+8A_4Z-V_4B.
}
\]

Its compressed selector is

\[
\boxed{
\langle n_1,n_2|\mathcal C_{4,L}^{\rm op}|n_1,n_2\rangle
=C_L^{\rm op}(n)+16A_4\frac{L-2}{L^2}n_1n_2,
}
\]

so the open selector coefficient is exactly \((L-1)/L\) times the periodic coefficient.

The exact open squared singular value is

\[
\boxed{
(\tau_{L;n_1,n_2}^{\rm op})^2
=64A_4^2\Gamma_0-16A_4V_4\Gamma_1+V_4^2\Gamma_2,
}
\]

where \(\Gamma_0=\operatorname{Var}(Z)\), \(\Gamma_1=\operatorname{Cov}(Z,B)\), and \(\Gamma_2=\operatorname{Var}(B)\). Closed rational formulas for all three are given in `DERIVATION.md` and `paper_insert.tex`.

## Consequences for the planned sectors

For physical shell parameters \(A_4>0\), \(V_4>0\) and \(L\ge3\), the open leakage vanishes exactly at the four corners

\[
(n_1,n_2)\in\{0,L\}\times\{0,L\}.
\]

In particular, the periodic accidental zeros at \(L=4,n=3\), namely \((3,0)\) and \((0,3)\), disappear on the path. At the project default \(V_4/A_4=1\), every open Phase-3 composition at \(L=4,6,8\) and total filling \(n=2,3\) has nonzero order-four leakage.

## Normalization qualification

The periodic result satisfies

\[
(\tau^{\rm per})^2/A_4^2=\sigma_L^2
\]

independently of \(U_4\). On the path,

\[
\frac{(\tau^{\rm op})^2}{A_4^2}
=64\Gamma_0-16\frac{V_4}{A_4}\Gamma_1
+\left(\frac{V_4}{A_4}\right)^2\Gamma_2.
\]

Therefore open data can be normalized by \(A_4\) only after the dimensionless shell ratio \(V_4/A_4\) has been fixed and recorded. A common rescaling of all shell parameters leaves the normalized result unchanged; varying \(A_4\) at fixed matching parameter \(\beta_4\) generally changes the ratio and hence changes the open normalized leakage.
