# Exact unscreened leakage certificate, version 2: results

## 1. Explicit SVD and general norm

Let

\[
\mathcal Z_n=\operatorname{span}\{|n_1,n_2\rangle:n_1+n_2=n\},
\qquad
R=(1-P_{\mathcal Z_n})\mathcal C P_{\mathcal Z_n}.
\]

Because the charge operator preserves both block charges, distinct compositions have orthogonal images. For each composition \(c\),

\[
R|c\rangle=(\mathcal C-\langle\mathcal C\rangle_c)|c\rangle,
\qquad
\tau_c^2=\operatorname{Var}_c(\mathcal C).
\]

For \(\tau_c>0\), define

\[
|\ell_c\rangle=\frac{(\mathcal C-\langle\mathcal C\rangle_c)|c\rangle}{\tau_c}.
\]

Then the exact singular-value resolution is

\[
\boxed{R=\sum_{\tau_c>0}\tau_c|\ell_c\rangle\langle c|,}
\]

with

\[
\boxed{R^\dagger R=\sum_c\tau_c^2|c\rangle\langle c|,}
\]

\[
\boxed{\|R\|=\max_c\tau_c,\qquad
\operatorname{rank}R=\#\{c:\tau_c>0\}.}
\]

This completes the requested action, singular-value, norm, and rank statement. The composition basis is the canonical right-singular basis selected by the conserved block charges; the operator identities themselves are basis independent.

## 2. Periodic cycle

For \(L\ge3\),

\[
\tau_{L;n_1,n_2}^{2,{\rm per}}=A_4^2\sigma_L^2(n_1,n_2),
\]

where

\[
\boxed{
\begin{aligned}
\sigma_L^2=64\Bigg[&
\frac{n_1(n_1-1)(L-n_1)(L-n_1-1)+n_2(n_2-1)(L-n_2)(L-n_2-1)}{(L-2)(L-1)^2}\\
&+\frac{2n_1n_2(L-n_1)(L-n_2)(3L-8)}{L^2(L-1)^2}
\Bigg].
\end{aligned}}
\]

The periodic zero set for \(L\ge3\) consists exactly of sectors in which one block is empty or full and the other has filling \(0,1,L-1\), or \(L\). At \(L=2\), every sector has zero leakage.

The projected selector remains

\[
16A_4\frac{L-2}{L(L-1)}n_1n_2.
\]

## 3. Standard open path

For the path with \(L\) sites and \(L-1\) charge bonds, define endpoint occupancy \(B\) and the interior quadratic variable \(Z\) as in `DERIVATION.md`. The exact configuration-wise decomposition is

\[
\boxed{
\mathcal C_{4,L}^{\rm op}
=(8A_4-2U_4)n+8A_4Z-V_4B,
\qquad V_4=4A_4-U_4.
}
\]

The compressed open selector is

\[
\boxed{
\langle\mathcal C_{4,L}^{\rm op}\rangle_{n_1,n_2}
=C_L^{\rm op}(n)+16A_4\frac{L-2}{L^2}n_1n_2.
}
\]

Thus the open selector coefficient is exactly \((L-1)/L\) times the periodic coefficient.

The exact open squared singular value is

\[
\boxed{
(\tau_{L;n_1,n_2}^{\rm op})^2
=64A_4^2\Gamma_0-16A_4V_4\Gamma_1+V_4^2\Gamma_2,
}
\]

where \(\Gamma_0=\operatorname{Var}(Z)\), \(\Gamma_1=\operatorname{Cov}(Z,B)\), and \(\Gamma_2=\operatorname{Var}(B)\), with closed rational formulas in `DERIVATION.md`.

For physical \(A_4>0,V_4>0\) and \(L\ge3\), open leakage vanishes exactly at

\[
(n_1,n_2)\in\{0,L\}\times\{0,L\}.
\]

The periodic \(L=4,n=3\) endpoint zeros therefore disappear after deleting the wrap bond.

## 4. Exact Phase-3 Gram matrices

Entries below are the diagonal of \(R^\dagger R/A_4^2\), ordered as \((0,n),(1,n-1),\ldots,(n,0)\). Open values use the project default \(V_4/A_4=1\).

### Periodic

| \(L\) | \(n\) | exact diagonal | rank |
|---:|---:|---|---:|
| 4 | 2 | \([128/9,\ 32,\ 128/9]\) | 3/3 |
| 4 | 3 | \([0,\ 512/9,\ 512/9,\ 0]\) | 2/4 |
| 6 | 2 | \([384/25,\ 320/9,\ 384/25]\) | 3/3 |
| 6 | 3 | \([576/25,\ 16256/225,\ 16256/225,\ 576/25]\) | 4/4 |
| 8 | 2 | \([640/49,\ 32,\ 640/49]\) | 3/3 |
| 8 | 3 | \([1280/49,\ 3328/49,\ 3328/49,\ 1280/49]\) | 4/4 |

### Open path, \(V_4/A_4=1\)

| \(L\) | \(n\) | exact diagonal | rank |
|---:|---:|---|---:|
| 4 | 2 | \([19,\ 65/2,\ 19]\) | 3/3 |
| 4 | 3 | \([81/4,\ 1015/12,\ 1015/12,\ 81/4]\) | 4/4 |
| 6 | 2 | \([16,\ 2660/81,\ 16]\) | 3/3 |
| 6 | 3 | \([146/5,\ 30826/405,\ 30826/405,\ 146/5]\) | 4/4 |
| 8 | 2 | \([369/28,\ 237/8,\ 369/28]\) | 3/3 |
| 8 | 3 | \([3165/112,\ 7521/112,\ 7521/112,\ 3165/112]\) | 4/4 |

At \(A_4=0.098\), the open composition-resolved leakage norms are:  The tabulated values are the full singular values \(\tau^{\rm op}\), so the factor \(A_4\) is already included; only the microscopic prefactor \(\lambda^4\Delta\) remains to be restored.

| \(L\) | \(n\) | norms in composition order |
|---:|---:|---|
| 4 | 2 | \(0.427172,\ 0.558686,\ 0.427172\) |
| 4 | 3 | \(0.441000,\ 0.901298,\ 0.901298,\ 0.441000\) |
| 6 | 2 | \(0.392000,\ 0.561596,\ 0.392000\) |
| 6 | 3 | \(0.529563,\ 0.854982,\ 0.854982,\ 0.529563\) |
| 8 | 2 | \(0.355763,\ 0.533403,\ 0.355763\) |
| 8 | 3 | \(0.520959,\ 0.803073,\ 0.803073,\ 0.520959\) |

These dimensionless values multiply \(\lambda^4\Delta\) in the order-four effective Hamiltonian.

## 5. Boundary and normalization consequences

The parent theorem is periodic, so the cycle proposition is the theorem-level statement. Open analysis is nevertheless required for the planned DMRG calculations because the project explicitly uses open systems where possible. The path result is not a finite-size approximation to the cycle result; it is the exact result for a different boundary operator.

On the cycle,

\[
(\tau^{\rm per})^2/A_4^2=\sigma_L^2
\]

is independent of \(U_4\). On the path,

\[
\frac{(\tau^{\rm op})^2}{A_4^2}
=64\Gamma_0-16(V_4/A_4)\Gamma_1+(V_4/A_4)^2\Gamma_2.
\]

Thus the open normalized certificate requires the ratio \(V_4/A_4\) to be recorded. A common rescaling of the shell leaves it unchanged, but changing \(A_4\) at fixed \(\beta_4\) generally changes that ratio.

## 6. Phase-3 implication

Neither boundary supports the interpretation of the unmodified order-four shell as an autonomous product-AGP selector Hamiltonian in the planned interior sectors. On the open path the conclusion is even stronger: all listed Phase-3 compositions leak, including the polarized \(L=4,n=3\) endpoints that are accidental periodic zeros.

The valid routes remain:

1. present the selector as an exact compression;
2. study the true order-four charge-selected states in the full seniority-zero space;
3. add or derive an order-four isolation mechanism that preserves the product-AGP manifold.
