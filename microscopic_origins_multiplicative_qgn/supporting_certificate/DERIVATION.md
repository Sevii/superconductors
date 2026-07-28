# Exact unscreened leakage: explicit SVD and periodic/open-boundary formulas

## 1. Scope and boundary decision

The microscopic lattice theorem is stated for a periodic array of bond controls. Its theorem-level leakage certificate is therefore the periodic result. The numerical project, however, explicitly proposes both periodic and open chains. Under the standard open truncation—keep the same `L` active sites but remove the wrap bond, leaving the path edges `x=1,...,L-1`—the periodic formula cannot be reused: endpoint sites have degree one, so the linear charge sum is no longer constant inside a fixed-filling product AGP.

Accordingly, the task is split as follows.

- **T2.3a (periodic):** the exact theorem-level result.
- **T2.3b (open):** a separate path formula needed before open-chain Phase-3 spectra are interpreted.

The formulas below complete both parts for the standard path truncation.

---

## 2. General composition-resolved action and singular-value decomposition

Let

\[
\mathcal Z_n=\operatorname{span}\{\lvert n_1,n_2\rangle:n_1+n_2=n\},
\qquad
R_{4,L}^{(n)}=(1-P_{\mathcal Z_n})\mathcal C_{4,L}P_{\mathcal Z_n}.
\]

The charge operator preserves the two block charges. Each fixed composition sector contains one normalized product AGP, so for

\[
\mu_{L;n_1,n_2}
=\langle n_1,n_2\rvert\mathcal C_{4,L}\lvert n_1,n_2\rangle
\]

one has the exact action

\[
\boxed{
R_{4,L}^{(n)}\lvert n_1,n_2\rangle
=
\bigl(\mathcal C_{4,L}-\mu_{L;n_1,n_2}\bigr)
\lvert n_1,n_2\rangle .
}
\]

Define

\[
\tau_{L;n_1,n_2}^2
=
\left\langle
\bigl(\mathcal C_{4,L}-\mu_{L;n_1,n_2}\bigr)^2
\right\rangle_{n_1,n_2}
=
\operatorname{Var}_{n_1,n_2}(\mathcal C_{4,L}).
\]

When \(\tau_{L;n_1,n_2}>0\), set

\[
\lvert\ell_{L;n_1,n_2}\rangle
=
\frac{
(\mathcal C_{4,L}-\mu_{L;n_1,n_2})\lvert n_1,n_2\rangle
}{\tau_{L;n_1,n_2}}.
\]

This vector is orthogonal to the product AGP in its composition sector. Vectors associated with different compositions have orthogonal support because their block charges differ. Therefore the exact singular-value decomposition is

\[
\boxed{
R_{4,L}^{(n)}
=
\sum_{\substack{n_1+n_2=n\\ \tau_{L;n_1,n_2}>0}}
\tau_{L;n_1,n_2}
\lvert\ell_{L;n_1,n_2}\rangle
\langle n_1,n_2\rvert .
}
\]

Consequently,

\[
\boxed{
(R_{4,L}^{(n)})^\dagger R_{4,L}^{(n)}
=
\sum_{n_1+n_2=n}
\tau_{L;n_1,n_2}^2
\lvert n_1,n_2\rangle\langle n_1,n_2\rvert,
}
\]

\[
\boxed{
\|R_{4,L}^{(n)}\|
=
\max_{n_1+n_2=n}\tau_{L;n_1,n_2},
\qquad
\operatorname{rank}R_{4,L}^{(n)}
=
\#\{(n_1,n_2):\tau_{L;n_1,n_2}>0\}.
}
\]

This is basis-independent at the operator level; the product-AGP composition basis simply gives the canonical right singular vectors selected by the conserved block charges.

---

## 3. Periodic cycle

Let \(X_{a,x}\in\{0,1\}\) be the hard-core pair occupation on site \(x\) in block \(a=1,2\). On the cycle,

\[
N_x=2\bigl(X_{1,x}+X_{1,x+1}+X_{2,x}+X_{2,x+1}\bigr),
\qquad
D_L=\sum_{x\in\mathbb Z_L}N_x^2,
\]

and

\[
\mathcal C_{4,L}^{\rm per}
=A_4D_L-\frac{U_4}{2}\sum_xN_x.
\]

At fixed total pair number \(n=n_1+n_2\),

\[
\sum_xN_x=4n
\]

is a scalar. Hence

\[
R_{4,L}^{(n),\rm per}\lvert n_1,n_2\rangle
=A_4(D_L-\langle D_L\rangle_{n_1,n_2})\lvert n_1,n_2\rangle.
\]

For \(L\ge3\),

\[
\tau_{L;n_1,n_2}^{2,{\rm per}}
=A_4^2\sigma_L^2(n_1,n_2),
\]

where

\[
\boxed{
\begin{aligned}
\sigma_L^2(n_1,n_2)=64\Bigg[&
\frac{
 n_1(n_1-1)(L-n_1)(L-n_1-1)
+n_2(n_2-1)(L-n_2)(L-n_2-1)
}{(L-2)(L-1)^2}\\
&+\frac{
2n_1n_2(L-n_1)(L-n_2)(3L-8)
}{L^2(L-1)^2}
\Bigg].
\end{aligned}
}
\]

The projected mean is

\[
\langle D_L\rangle
=
8n+
\frac{8[n_1(n_1-1)+n_2(n_2-1)]}{L-1}
+
\frac{32n_1n_2}{L},
\]

or, at fixed \(n\),

\[
\langle D_L\rangle
=
8n+\frac{8n(n-1)}{L-1}
+
16\frac{L-2}{L(L-1)}n_1n_2.
\]

Thus the selector remains exact as a compression, while the variance gives its off-manifold companion.

For \(L\ge3\), the periodic leakage vanishes exactly when one block is empty or full and the other has occupation \(0,1,L-1\), or \(L\). For \(L=2\), every sector has zero leakage.

---

## 4. Open path: operator decomposition

Now retain only the path edges \(x=1,\ldots,L-1\). Define

\[
B_a=X_{a,1}+X_{a,L},
\qquad
E_a=\sum_{x=1}^{L-1}X_{a,x}X_{a,x+1},
\]

\[
t_{a,x}=X_{a,x}+X_{a,x+1},
\qquad
C=\sum_{x=1}^{L-1}t_{1,x}t_{2,x},
\]

and

\[
B=B_1+B_2,
\qquad
Z=E_1+E_2+C.
\]

Direct expansion gives the configuration-wise identities

\[
\boxed{
\sum_{x=1}^{L-1}N_x=4n-2B,
\qquad
\sum_{x=1}^{L-1}N_x^2=8n-4B+8Z.
}
\]

Write the shell parameters in their physical form

\[
A_4=\frac{U_4+V_4}{4},
\qquad
V_4=4A_4-U_4>0.
\]

Then the open charge operator is

\[
\boxed{
\mathcal C_{4,L}^{\rm op}
=(8A_4-2U_4)n+8A_4Z-V_4B.
}
\]

This displays the essential boundary effect: the periodic linear-charge term was scalar, but on a path it combines with the endpoint contribution from \(D_L\) into the fluctuating term \(-V_4B\).

---

## 5. Open projected selector

For a uniform \(n_a\)-subset of \(L\) sites,

\[
\langle B_a\rangle=\frac{2n_a}{L},
\qquad
\langle E_a\rangle=\frac{n_a(n_a-1)}{L},
\]

and block independence gives

\[
\langle C\rangle
=\frac{4(L-1)n_1n_2}{L^2}.
\]

Therefore

\[
\left\langle\sum_{x=1}^{L-1}N_x\right\rangle
=\frac{4(L-1)n}{L},
\]

\[
\boxed{
\left\langle\sum_{x=1}^{L-1}N_x^2\right\rangle
=
\frac{8n(L+n-2)}{L}
+
16\frac{L-2}{L^2}n_1n_2.
}
\]

The compressed open charge operator is consequently

\[
\boxed{
\mu_{L;n_1,n_2}^{\rm op}
=
C_L^{\rm op}(n)
+
16A_4\frac{L-2}{L^2}n_1n_2,
}
\]

with

\[
C_L^{\rm op}(n)
=
\frac{8A_4n(L+n-2)-2U_4n(L-1)}{L}.
\]

The open selector coefficient is exactly \((L-1)/L\) times the periodic coefficient. It still vanishes at \(L=2\).

---

## 6. Open variance ingredients

Introduce the fixed-size subset probabilities

\[
p_k(n)=\frac{(n)_k}{(L)_k},
\]

with the value zero when the event is impossible. For one block at filling \(n\), exact path counting gives

\[
\boxed{
b_L(n):=\operatorname{Var}(B_a)
=
\frac{2n(L-n)(L-2)}{L^2(L-1)},
}
\]

\[
\boxed{
h_L(n):=\operatorname{Cov}(E_a,B_a)
=
-\frac{2n(n-1)(L-n)}{L^2(L-1)},
}
\]

\[
\boxed{
e_L(n):=\operatorname{Var}(E_a)
=
\frac{n(n-1)(L-n)(L-n+1)}{L^2(L-1)}.
}
\]

For the edge variable \(t_{a,x}\), define

\[
u_0(n)=2p_1(n)+2p_2(n),
\qquad
u_1(n)=p_1(n)+3p_2(n),
\qquad
u_2(n)=4p_2(n).
\]

These are respectively \(\mathbb E[t_x^2]\), \(\mathbb E[t_xt_y]\) for adjacent distinct edges, and the same moment for disjoint edges. Since the path has \(L-1\) equal ordered edge pairs, \(2(L-2)\) adjacent ordered pairs, and \((L-2)(L-3)\) disjoint ordered pairs,

\[
\boxed{
\begin{aligned}
w_L(n_1,n_2):=\operatorname{Var}(C)
={}&(L-1)u_0(n_1)u_0(n_2)
+2(L-2)u_1(n_1)u_1(n_2)\\
&+(L-2)(L-3)u_2(n_1)u_2(n_2)
-
\left[\frac{4(L-1)n_1n_2}{L^2}\right]^2.
\end{aligned}
}
\]

The mixed covariances follow without further edge enumeration. Conditioned on block 1,

\[
\mathbb E[C\mid X_1]
=
\frac{2n_2}{L}\sum_xt_{1,x}
=
\frac{2n_2}{L}(2n_1-B_1).
\]

Hence

\[
\operatorname{Cov}(C,E_1)
=-\frac{2n_2}{L}h_L(n_1),
\qquad
\operatorname{Cov}(C,B_1)
=-\frac{2n_2}{L}b_L(n_1),
\]

with the analogous formulas after exchanging the blocks.

---

## 7. Exact open leakage formula

Define

\[
\boxed{
\begin{aligned}
\Gamma_0={}&e_L(n_1)+e_L(n_2)+w_L(n_1,n_2)
-\frac{4n_2}{L}h_L(n_1)-\frac{4n_1}{L}h_L(n_2),\\
\Gamma_1={}&h_L(n_1)+h_L(n_2)
-\frac{2n_2}{L}b_L(n_1)-\frac{2n_1}{L}b_L(n_2),\\
\Gamma_2={}&b_L(n_1)+b_L(n_2).
\end{aligned}
}
\]

By construction,

\[
\Gamma_0=\operatorname{Var}(Z),
\qquad
\Gamma_1=\operatorname{Cov}(Z,B),
\qquad
\Gamma_2=\operatorname{Var}(B).
\]

Therefore the exact open leakage singular value is

\[
\boxed{
(\tau_{L;n_1,n_2}^{\rm op})^2
=
64A_4^2\Gamma_0
-16A_4V_4\Gamma_1
+V_4^2\Gamma_2.
}
\]

Equivalently, with \(q=V_4/A_4\),

\[
\boxed{
\frac{(\tau_{L;n_1,n_2}^{\rm op})^2}{A_4^2}
=
64\Gamma_0-16q\Gamma_1+q^2\Gamma_2.
}
\]

This formula, inserted into the general SVD of Section 2, gives the exact open composition-resolved action, every singular value, the rank, and the operator norm.

---

## 8. Exact open zero set

The sign structure is decisive:

\[
\Gamma_0\ge0,
\qquad
\Gamma_2\ge0,
\qquad
\Gamma_1\le0.
\]

The first and third statements are variances. The middle covariance obeys \(\Gamma_1\le0\) because \(h_L(n)\le0\), \(b_L(n)\ge0\), and every coefficient subtracted in its definition is nonnegative.

For physical shell parameters \(A_4>0\) and \(V_4>0\), all three terms in

\[
64A_4^2\Gamma_0-16A_4V_4\Gamma_1+V_4^2\Gamma_2
\]

are therefore nonnegative. If the leakage vanishes, then \(\Gamma_2=0\). For \(L\ge3\), the formula for \(b_L(n)\) shows that this occurs only when each block is empty or full. Conversely, at those four corners the charge operator is constant. Thus

\[
\boxed{
L\ge3,\ V_4>0:
\quad
\tau_{L;n_1,n_2}^{\rm op}=0
\iff
(n_1,n_2)\in\{0,L\}\times\{0,L\}.
}
\]

For \(L=2\), the path has one bond containing all sites, so \(N_1=2n\) is constant and every composition has zero leakage.

In particular, the periodic accidental zeros at \(L=4,n=3\) do **not** survive the open truncation.

---

## 9. Normalization consequence

On the cycle,

\[
(\tau^{\rm per})^2=A_4^2\sigma_L^2,
\]

so the normalized leakage is independent of \(U_4\) and of the absolute choice of \(A_4\).

On the path,

\[
\frac{(\tau^{\rm op})^2}{A_4^2}
=
64\Gamma_0-16(V_4/A_4)\Gamma_1+(V_4/A_4)^2\Gamma_2.
\]

Thus the open certificate is invariant under a common rescaling of \((A_4,U_4,V_4)\), but it must hold the physical ratio \(V_4/A_4\) fixed. The periodic T2.5 statement remains exact; its unrestricted `A4` normalization should not be transferred to open chains without this ratio qualification.

---

## 10. Verification range

The accompanying verifier exhaustively checks both boundary conditions, all compositions, and all hard-core configurations for \(L=2,\ldots,9\):

- 760 boundary/composition sectors;
- 699,040 configuration evaluations;
- exact rational agreement with every projected mean and variance formula;
- exact centered-action and SVD normalization checks;
- exact periodic and open zero sets;
- exact Phase-3 Gram diagonals at \(L=4,6,8\), \(n=2,3\).
