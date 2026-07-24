---
title: "Exact Flat-Band Stiffness from Pair Mobility in a QGN Model"
subtitle: "Finite-size obstructions, a many-body reduction theorem, and a proof of the Gao-Han-Khalaf formula for Model II"
author: "Nicholas Sledgianowski"
date: "July 23, 2026"
geometry: margin=0.88in
fontsize: 10.5pt
header-includes:
  - \usepackage{amsmath,amssymb,mathtools,bm}
  - \usepackage{booktabs}
  - \usepackage{microtype}
  - \usepackage{array}
  - \allowdisplaybreaks
abstract: |
  Gao, Han, and Khalaf conjectured that the thermodynamic flat-connection stiffness of a quantum-geometric-nesting (QGN) superconductor is exactly the hard-core filling factor times the inverse two-particle pair mass. We prove this equality for their projected-Hubbard Model II and derive a broader restricted reduction theorem. The route to the theorem passes through two exact boundary cases: a period-circumference winding deformation shows that finite-torus twist curvature is underdetermined unless a displacement lift is specified, while a fixed-local reducible QGN model shows that locality and positive pair mass do not suffice without irreducibility of the paired ground branch. Under H1-H4 in the Hermitian positive-square class, the fixed-pair-number Hessian obeys
  $C^{(n)}=\frac{n(L-n)}{L-1}C^{(1)}+(n^2-\frac{n(L-n)}{L-1})\Gamma$,
  where the longitudinal tensor $\Gamma$ is positive semidefinite and is a pure winding observable under a translation-covariant flat connection. In the projected PSD-Hubbard subclass, connectedness of the projected-orbital frame is equivalent to uniqueness of the antisymmetrized-geminal-power ground state in every nonsingular filling sector. For Model II, the frame is connected on every rectangular torus for $0<\xi<\pi/2$, $\Gamma=0$ exactly with no finite-size correction, and
  $(m_{\mathrm{pair}}^{-1})_{ij}=|U|\xi^2\delta_{ij}/8$ for every $N_x,N_y\ge3$. Hence
  $\kappa_{ij}^{(n)}=[V/(V-1)]\nu(1-\nu)|U|\xi^2\delta_{ij}/8$
  at finite size and
  $D_{s,ij}=\nu(1-\nu)(m_{\mathrm{pair}}^{-1})_{ij}$
  in the thermodynamic flat-connection definition used by Gao, Han, and Khalaf.
---

# 1. Introduction and headline theorem

Gao, Han, and Khalaf used the reduced-density-matrix bootstrap to obtain rigorous lower bounds on the zero-temperature flat-connection response of frustration-free flat-band superconductors. In every QGN example they tested, the lower bound saturated a variational upper bound determined by the two-particle pair mass. They consequently conjectured

$$
D_{s,ij}
=
\frac{N_{\mathrm{flat}}}{2}\,
\nu(1-\nu)\,
(m_{\mathrm{pair}}^{-1})_{ij},
\tag{1}
$$

and explicitly left a rigorous proof to future work [1]. The content of (1) is unusually strong: it says that a finite-density many-body response is fixed by one bound pair, with no dynamical dressing beyond the hard-core filling factor.

This paper gives a complete proof for the projected positive-semidefinite Hubbard realization of their Model II and, at the same time, identifies the exact hypotheses that make the reduction work.

**Headline theorem.** Let

$$
H_{\mathrm{II}}(k)
=
-t[\sigma_x\sin\alpha_k+\sigma_y\cos\alpha_k],
\qquad
\alpha_k=\xi(\cos k_x+\cos k_y),
\tag{T1}
$$

and project the positive-semidefinite Hubbard interaction onto its lower time-reversal-related flat bands. For $0<\xi<\pi/2$, every rectangular torus with $N_x,N_y\ge3$, $V=N_xN_y$, and $1\le n\le V-1$ has a unique analytic AGP ground branch with exact curvature

$$
\boxed{
\kappa_{ij}^{(n)}(V)
=
\frac{V}{V-1}\,
\nu(1-\nu)\,
\frac{|U|\xi^2}{8}\delta_{ij},
\qquad
\nu=\frac nV.
}
\tag{T2}
$$

The two-particle bound branch has

$$
\boxed{
(m_{\mathrm{pair}}^{-1})_{ij}
=
\frac{|U|\xi^2}{8}\delta_{ij}>0,
}
\tag{T3}
$$

so along every regular sequence with $n_V/V\to\nu\in(0,1)$,

$$
\boxed{
D_{s,ij}
=
\nu(1-\nu)(m_{\mathrm{pair}}^{-1})_{ij}.
}
\tag{T4}
$$

Here $D_s$ is exactly the thermodynamic ground-energy curvature under the flat gauge connection used in Ref. [1]. The finite-size factor $V/(V-1)$ is explicit and universal.

The proof emerged through the following sequence.

| Stage | Exact conclusion | Role in the final theorem |
|---|---|---|
| Winding counterexample | Untwisted finite-torus samples do not determine twist curvature. | The displacement-resolved Peierls lift must be part of the model. |
| Reducible QGN counterexample | A fixed-local QGN model can have zero stiffness and positive pair mass. | Connectivity or branch uniqueness is essential. |
| Connected-model search | All surviving deviations reduce to one nonnegative longitudinal tensor. | Suggests an exact algebraic reduction rather than another counterexample. |
| General Hermitian-QGN theorem | $C^{(n)}=\rho_nC^{(1)}+(n^2-\rho_n)\Gamma$. | Isolates the only finite-size defect. |
| H4 connectivity theorem | Connected projected frame $\Longleftrightarrow$ unique AGP zero mode. | Removes the branch hypothesis in projected Hubbard models. |
| Model-II dispersion lemma | Pair mass is exact and positive; $\Gamma=0$ exactly. | Closes every remaining hypothesis and proves (T2)-(T4). |

The standard QGN ideal Hamiltonians use Hermitian local null operators and a positive interaction kernel [2]. The abstract theorem below is restricted to that positive-square class, a full-rank equal-singular-value pairing matrix, and a specified translation-covariant twist lift. The final Model-II result is unconditional in the range stated above.

**Notation.** Throughout, $L$ is the number of collective pair orbitals, $\mathbf N=(N_1,\ldots,N_d)$ are torus side lengths, and $V=\prod_jN_j$ is the number of unit cells. The one-spin flat-band space is denoted $\mathcal V$. The filling is always $\nu=n/L$; for Model II, $L=V$. In Section 2.1, $g>0$ denotes the magnitude of an attractive coupling. From Section 2.2 onward, $U>0$ denotes the coefficient of a positive-semidefinite square, while the Model-II formulas use $|U|$.

# 2. Two exact obstructions and the corrected theorem class

The two counterexamples are not discarded preliminaries. Together they identify the two independent ways in which an overbroad finite-size theorem can fail: the twist may be underdetermined, or the paired ground branch may be reducible.

## 2.1 Winding/lift obstruction

On an $M$-point one-dimensional torus, consider the smooth rank-one spinor

$$
u_{\uparrow,c}^{(M,q)}(k)=
\begin{pmatrix}
\cos[\pi/4+c\sin(qMk)]\\
e^{i\phi(k)}\sin[\pi/4+c\sin(qMk)]
\end{pmatrix},
\qquad
u_\downarrow(k)=u_\uparrow^*(-k).
\tag{O1}
$$

At every allowed momentum $k_m=2\pi m/M$, the deformation is invisible because $\sin(qMk_m)=0$. Under the continuous shift $k_m\mapsto k_m+A$, however,

$$
\sin[qM(k_m+A)]=\sin(qMA)
\tag{O2}
$$

is independent of $m$. For a projected onsite attraction of magnitude $g>0$, the shifted deformation angle $x_A=c\sin(qMA)$ is independent of $k_m$. The two orbital pair amplitudes obey $\cos(\pi/4+x_A)\cos(\pi/4-x_A)=\sin(\pi/4+x_A)\sin(\pi/4-x_A)=\tfrac12\cos(2x_A)$. Every projected interaction matrix element therefore acquires the same squared factor, giving the exact operator factorization

$$
H_c(A)=
\cos^2[2c\sin(qMA)]\,H_0(A).
\tag{O3}
$$

If $E_n(0)=-gn/2$, then

$$
E_{n,c}''(0)
=
E_{n,0}''(0)+4gc^2q^2M^2n.
\tag{O4}
$$

Consequently, with $\rho_n=n(M-n)/(M-1)$,

$$
\boxed{
E_{n,c}''-\rho_nE_{1,c}''
=
\bigl(E_{n,0}''-\rho_nE_{1,0}''\bigr)
+
4gc^2q^2M^2\frac{n(n-1)}{M-1}.
}
\tag{O5}
$$

The four-cell realization has a unique gapped branch and an exact normalized defect $g/24$ at $n=2$. A nonresonant $M=8$, harmonic-four control obeys the filling law, whereas the resonant harmonic-eight deformation violates it again. The failure is therefore not caused by a generic high harmonic; it is a period-circumference winding resonance.

The structural conclusion is more important than the example:

> The untwisted finite-torus Hamiltonian fixes Fourier alias sums, not the displacement moments needed for a continuous twist. A finite-size curvature theorem must specify the displacement lift or the complete family $P_A$.

This example is not a fixed-local counterexample along a thermodynamic sequence because its resonant displacement grows with $M$. Fixed finite range removes the resonance once the torus is sufficiently large, and uniform exponential locality suppresses it exponentially.

## 2.2 Connectivity/reducibility obstruction

A distinct obstruction survives fixed locality. Take

$$
u_\uparrow(k)=
\begin{pmatrix}\cos k\\ \sin k\end{pmatrix},
\qquad
u_\downarrow(k)=u_\uparrow^*(-k),
\tag{O6}
$$

whose projector has range two and obeys $P(k+\pi)=P(k)$. On every even torus, the band space decomposes into even and odd cell components. The twist-stable UPC holds exactly, and the interaction is the projected positive-semidefinite Hubbard Hamiltonian

$$
H(A)=\frac U2\sum_{R,\alpha}
[\bar n_{R\alpha\uparrow}(A)-\bar n_{R\alpha\downarrow}(A)]^2.
\tag{O7}
$$

At half filling, the state with one component completely paired and the other empty is annihilated by every square factor for every $A$. Therefore

$$
E_0(A)=0,
\qquad
D_s(M,\nu=1/2)=0.
\tag{O8}
$$

The one-pair branch is nevertheless

$$
E_{\mathrm{pair}}(Q)=\frac U2\sin^2Q,
\qquad
m_{\mathrm{pair}}^{-1}=U.
\tag{O9}
$$

Thus

$$
\boxed{
R_M(1/2)=0
}
\tag{O10}
$$

for every even $M\ge6$. This is a fixed, strictly local, thermodynamically persistent discrepancy, but it occurs because the ground space is reducible and degenerate. On odd tori, where the step-two graph is connected, the reduction law returns.

The corrected theorem class must therefore include both:

1. a specified fixed-local twist lift; and
2. an irreducible, uniquely continued AGP ground branch.

The remainder of the paper proves that, in the standard Hermitian positive-square class, these repairs remove both known failure mechanisms and leave one explicitly controlled longitudinal winding tensor.

# 3. Finite-dimensional QGN setup

Let $\mathfrak h$ be a one-particle Hilbert space of dimension $2L$. Let

$$
J^T=-J,\qquad J^\dagger J=I_{2L}
\tag{2}
$$

be a skew-unitary pairing matrix. Define

$$
\eta^+=\frac12\sum_{a,b=1}^{2L}J_{ab}c_a^\dagger c_b^\dagger,
\qquad
\eta^-=(\eta^+)^\dagger,
\qquad
\eta^z=\frac{\widehat N-L}{2}.
\tag{3}
$$

These operators generate an $\mathfrak{su}(2)$ algebra. By unitary congruence, $J$ can be put into the canonical form $\bigoplus_{a=1}^L i\sigma_y$. Thus (3) is the fully multiband, basis-independent version of $L$ collective pair orbitals; $J$ may mix bands, spin, orbital labels, momenta, and gauge patches.

For $0\le n\le L$, define the normalized AGP state

$$
|n\rangle=
\left[\frac{(L-n)!}{n!L!}\right]^{1/2}
(\eta^+)^n|0\rangle.
\tag{4}
$$

The states $|n\rangle$ form the spin-$S=L/2$ multiplet, with $m=n-L/2$.

Let $A=(A_1,\ldots,A_d)$ be an electronic flat connection. In the sector with $2n$ particles, assume

$$
H_n(A)=\frac12\sum_{\lambda=1}^{M_s}
S_\lambda(A)^\dagger S_\lambda(A),
\tag{5}
$$

where the $S_\lambda(A)$ are analytic, Hermitian, particle-number-preserving operators. The index $\lambda$ may include real-space position, QGN channel, and a factorization index of the positive interaction kernel. We impose:

**H1: QGN scalar condition at zero twist.** With $S_\lambda=S_\lambda(0)$,

$$
[S_\lambda,\eta^+]=[S_\lambda,\eta^-]=[S_\lambda,\eta^z]=0.
\tag{6}
$$

**H2: Frustration freeness.** For the filling sectors considered,

$$
S_\lambda|n\rangle=0\quad\text{for every }\lambda.
\tag{7}
$$

**H3: One-body twist source.** The first derivatives

$$
B_{\lambda i}=\left.\partial_{A_i}S_\lambda(A)\right|_{A=0}
\tag{8}
$$

are number-conserving one-body operators,

$$
B_{\lambda i}=d\Gamma(b_{\lambda i})
=\sum_{a,b}(b_{\lambda i})_{ab}c_a^\dagger c_b,
\tag{9}
$$

with no $A$-dependent scalar centering term. Because the square factors are Hermitian, the matrices $b_{\lambda i}$ are Hermitian. Requiring a Hermitian square-factor representation is part of the restricted theorem.

**H4: Simple analytic branch.** At $A=0$, zero is a simple eigenvalue of $H_n(0)$ in the fixed-$2n$ sector and is separated from the rest of that finite-dimensional sector by a positive gap. Hence the ground eigenvalue $E_n(A)$ and eigenvector are analytic near $A=0$. No gap uniform in system size is assumed here.

Define the exact electronic-twist Hessian

$$
C_{ij}^{(n)}=\left.\partial_{A_i}\partial_{A_j}E_n(A)\right|_{A=0}.
\tag{10}
$$

## 3.1 Automatic H1-H3 for projected PSD-Hubbard models

Let $\mathcal V$ be the spin-up flat-band space, let the spin-down space be its time-reversal conjugate, and define projected orbital vectors

$$
v_x=P|x\rangle,\qquad p_x=|v_x\rangle\langle v_x|,\qquad \sum_xp_x=I_{\mathcal V}.
\tag{10a}
$$

The projected PSD-Hubbard factors are

$$
S_x=d\Gamma_\uparrow(p_x)-d\Gamma_\downarrow(\bar p_x).
\tag{10b}
$$

For every Hermitian one-body matrix $a=a^\dagger$,

$$
[d\Gamma_\uparrow(a)-d\Gamma_\downarrow(\bar a),\eta^+]=0.
\tag{10c}
$$

Thus H1 is automatic; H2 follows from H1 because the factors annihilate the vacuum; and differentiation of a twisted projector keeps the source one body, proving H3. Uniform pairing of the orbital weights is not needed for these algebraic facts.

## 3.2 Connectivity theorem closing H4 in the Hubbard subclass

Define the projected-frame overlap graph with active microscopic labels as vertices and

$$
x\sim y\quad\Longleftrightarrow\quad
\langle v_x,v_y\rangle=\langle x|P|y\rangle\ne0.
\tag{10d}
$$

**Proposition 1 (H4 from frame connectivity).** For the projected PSD-Hubbard Hamiltonian

$$
H_0=\frac U2\sum_x
[d\Gamma_\uparrow(p_x)-d\Gamma_\downarrow(\bar p_x)]^2,
\qquad U>0,
\tag{10e}
$$

the AGP is the unique zero-energy state in every total-$2n$ sector, $1\le n\le L-1$, if and only if the overlap graph is connected. Consequently the zero eigenvalue has a positive finite-size sector gap and its continuation under an analytic twist is a simple analytic ground branch.

*Proof.* Any zero-energy state is annihilated by every square factor. Summing those null equations and using $\sum_xp_x=I_{\mathcal V}$ gives $N_\uparrow=N_\downarrow=n$. Identify

$$
\Lambda^n\mathcal V\otimes\overline{\Lambda^n\mathcal V}
\simeq\operatorname{End}(\Lambda^n\mathcal V)
\tag{10f}
$$

by vectorization. If $A_x=d\Gamma_n(p_x)$, then a factor acts as

$$
S_x\operatorname{vec}(T)=\operatorname{vec}([A_x,T]),
\tag{10g}
$$

and hence

$$
\langle\operatorname{vec}(T),H_0\operatorname{vec}(T)\rangle
=\frac U2\sum_x\|[A_x,T]\|_{\mathrm{HS}}^2.
\tag{10h}
$$

Connectedness implies

$$
\operatorname{Lie}_{\mathbb R}\{ip_x\}=\mathfrak u(\mathcal V).
\tag{10i}
$$

One proof orders the frame along a spanning tree. If $\mathfrak u(W)$ has already been generated and the next vector has nonzero components in both $W$ and a new orthogonal direction, two commutators with its rank-one projector generate the symmetric and antisymmetric couplings to that direction; their commutator generates its diagonal projector. Induction reaches all of $\mathcal V$.

The exterior-power representation of $U(\mathcal V)$ on $\Lambda^n\mathcal V$ is irreducible for $0<n<L$. Schur's lemma therefore reduces the commutant in (10h) to scalar multiples of the identity. Its vectorization is the AGP state. Conversely, if the graph has components of span dimensions $L_1,\ldots,L_c$, the zero-mode dimension is

$$
[z^n]\prod_{r=1}^c(1+z+\cdots+z^{L_r}),
\tag{10j}
$$

which exceeds one at every nonsingular filling when $c>1$. Finite-dimensional positivity then supplies a strictly positive next eigenvalue, and analytic perturbation theory supplies the local analytic branch. $\square$

For a positive interaction kernel $H_K=\frac12\sum_{x,y}K_{xy}S_xS_y$, the same result holds whenever $K$ is strictly positive on the span of the active density-difference constraints. If the kernel loses rank, the correct replacement is irreducibility of the one-body Lie algebra generated by the surviving square factors.

For Model II, the all-size connectivity claim can be proved directly. Define the one-dimensional finite-grid coefficient

$$
c_r^{(N)}(\xi)
=\frac1N\sum_{m=0}^{N-1}
e^{i\xi\cos k_m}e^{-irk_m},
\qquad k_m=\frac{2\pi m}{N}.
\tag{C1}
$$

For $0<\xi<\pi/2$,

$$
\operatorname{Re}c_0^{(N)}(\xi)
=\frac1N\sum_m\cos(\xi\cos k_m)>0,
\tag{C2}
$$

because every summand is strictly positive. Inversion symmetry of the grid cancels the term odd in $\sin k_m$ and gives

$$
c_1^{(N)}(\xi)
=\frac1N\sum_m \cos k_m\,e^{i\xi\cos k_m},
\qquad
\operatorname{Im}c_1^{(N)}(\xi)
=\frac1N\sum_m \cos k_m\sin(\xi\cos k_m)>0.
\tag{C3}
$$

The last inequality follows because $x\sin(\xi x)>0$ for every nonzero $x\in[-1,1]$, and the grid contains $k=0$. Hence $c_0^{(N)}$ and $c_1^{(N)}$ are nonzero for every $N\ge2$. The two-dimensional off-diagonal projector coefficient factorizes, up to a fixed nonzero phase, as

$$
\widehat P_{12}(r_x,r_y)
\propto c_{r_x}^{(N_x)}(\xi)c_{r_y}^{(N_y)}(\xi).
\tag{C4}
$$

It is therefore nonzero at displacements $(0,0)$, $(1,0)$, and $(0,1)$. The first edge joins the two orbital layers within a cell, while the latter two join neighboring cells through opposite layers; two-step paths connect equal-orbital vertices across both lattice directions. Thus the projected-frame graph is connected on every rectangular torus for $0<\xi<\pi/2$. The restriction is imposed by the $c_0$ positivity argument, not by pair-mass positivity.

The analogous all-size Model-I statement is not used below. A separate gauge-invariant audit found connected graphs in the tested Model-I tori and in additional Model-II cases through $\xi=2$, but those results are retained as finite certificates rather than promoted to all-size theorems.

# 4. Frustration-free least-squares curvature

**Lemma 1 (least-squares curvature).** Let $D(A):\mathcal H\to\mathcal K$ be analytic between finite-dimensional Hilbert spaces and let

$$
H(A)=\frac12D(A)^\dagger D(A).
\tag{11}
$$

Suppose $D(0)|\psi\rangle=0$ and zero is a simple gapped eigenvalue of $H(0)$. Put $D=D(0)$ and

$$
b_i=(\partial_{A_i}D)(0)|\psi\rangle.
\tag{12}
$$

Then

$$
C_{ij}=\operatorname{Re}\langle \Pi b_i,\Pi b_j\rangle,
\qquad
\Pi=P_{\ker D^\dagger},
\tag{13}
$$

or, directionally,

$$
v_iC_{ij}v_j
=
\min_{\chi\perp\psi}
\left\|D\chi+v_i b_i\right\|^2.
\tag{14}
$$

*Proof.* Choose the analytic gauge $\langle\psi|\partial_{A_i}\psi\rangle=0$ and write $|\chi_i\rangle=\partial_{A_i}|\psi(A)\rangle|_0$. Since $D|\psi\rangle=0$, the first-order eigenvalue equation is the normal equation

$$
D^\dagger(D|\chi_i\rangle+b_i)=0.
\tag{15}
$$

The gap makes $H|_{\psi^\perp}$ invertible, so the normal equation has a unique solution $|\chi_i\rangle\perp|\psi\rangle$. Equivalently, $D|\chi_i\rangle+b_i$ is the orthogonal projection of $b_i$ onto $\ker D^\dagger$. Direct differentiation of $H=D^\dagger D/2$ gives

$$
C_{ij}=\operatorname{Re}\langle D\chi_i+b_i,D\chi_j+b_j\rangle.
\tag{16}
$$

Equation (13) follows, and (14) is the associated normal-equation characterization. Terms containing $\partial_i\partial_jD$ vanish because $D|\psi\rangle=0$. $\square$

For (5), take

$$
D_n|\chi\rangle=(S_1|\chi\rangle,\ldots,S_{M_s}|\chi\rangle),
\qquad
b_i^{(n)}=(B_{1i}|n\rangle,\ldots,B_{M_si}|n\rangle).
\tag{17}
$$

Then Lemma 1 is an exact version of the diamagnetic-minus-paramagnetic Kohn formula, reorganized as a positive least-squares residual.

# 5. Multiband $S\oplus(S-1)$ source decomposition

The following lemma closes the multiband step. It uses only the skew-unitarity of $J$ and therefore does not require a band-diagonal, spin-diagonal, or topologically trivial representation.

**Lemma 2 (one-body action on an AGP multiplet).** Let $B=d\Gamma(b)$ be any number-conserving one-body operator and define

$$
\beta_B=\langle1|B|1\rangle=\frac1L\operatorname{Tr}_{\mathfrak h}b.
\tag{18}
$$

There is a lowest-weight pair state $|\zeta_B\rangle$, orthogonal to $|1\rangle$, in total pseudospin $S-1=L/2-1$, such that

$$
B|n\rangle
=n\beta_B|n\rangle
+\sqrt{\rho_n}\,U_n\bigl(B|1\rangle-\beta_B|1\rangle\bigr),
\tag{19}
$$

where

$$
\rho_n=\frac{n(L-n)}{L-1},
\tag{20}
$$

and $U_n$ is the normalized pseudospin-ladder isometry from the one-pair weight space to the $n$-pair weight space inside every copy of the $S-1$ representation.

*Proof.* Since $B|0\rangle=0$,

$$
B(\eta^+)^n|0\rangle
=n(\eta^+)^{n-1}[B,\eta^+]|0\rangle.
\tag{21}
$$

Decompose the two-particle state

$$
[B,\eta^+]|0\rangle
=\beta_B\eta^+|0\rangle+|\zeta_B\rangle.
\tag{22}
$$

The coefficient is (18), because the one-particle reduced density matrix of $|1\rangle$ is $I_{2L}/L$. The orthogonal complement of $\eta^+|0\rangle$ in the two-particle sector is a direct sum of lowest-weight copies of pseudospin $S-1$. No lower total pseudospin can occur after acting once with a fermion bilinear on the maximal-spin vacuum.

The AGP norm is

$$
\|(\eta^+)^n|0\rangle\|^2
=\frac{n!L!}{(L-n)!}.
\tag{23}
$$

For a lowest-weight $S-1$ state,

$$
\|(\eta^+)^{n-1}|\zeta_B\rangle\|^2
=\frac{(n-1)!(L-2)!}{(L-n-1)!}
\|\zeta_B\|^2.
\tag{24}
$$

Combining (21)-(24) with the normalization in (4) gives exactly the squared ratio $\rho_n$ in (20). Because the norm factor in (24) is the same for every lowest-weight vector, the normalized ladder operation defines one linear isometry $U_n$ on the complete $S-1$ lowest-weight multiplicity space, not a vector-dependent choice. $\square$

The identity (18) is the point at which the full multiband pairing matrix simplifies: equal singular values make the one-pair density matrix proportional to the identity on the entire paired flat-band space. If the pairing form has several inequivalent singular-value blocks, the density matrix is not proportional to a single identity and the argument becomes blockwise rather than scalar.

# 6. Exact finite-size reduction identity

For every factor and twist direction, define

$$
\beta_{\lambda i}
=\langle1|B_{\lambda i}|1\rangle
=\frac1L\operatorname{Tr}b_{\lambda i},
\tag{25}
$$

and the longitudinal tensor

$$
\Gamma_{ij}
=\operatorname{Re}\sum_{\lambda=1}^{M_s}
\beta_{\lambda i}^*\beta_{\lambda j}.
\tag{26}
$$

It is positive semidefinite.

**Theorem 1 (finite-size restricted QGN reduction).** Under H1-H4, for every $1\le n\le L-1$,

$$
\boxed{
C_{ij}^{(n)}
=
\rho_n C_{ij}^{(1)}
+\bigl(n^2-\rho_n\bigr)\Gamma_{ij},
\qquad
\rho_n=\frac{n(L-n)}{L-1}.
}
\tag{27}
$$

In particular,

$$
C^{(n)}-\rho_nC^{(1)}\succeq0.
\tag{28}
$$

*Proof.* Let $D=\bigoplus_{n=0}^L D_n$ act on the full even-particle Fock space, with the target carrying the diagonal pseudospin action. Because every $S_\lambda$ commutes with the pseudospin algebra, $D$ is an $SU(2)$ intertwiner. Decompose the target into isotypic components $\bigoplus_S(V_S\otimes M_S)$. On each component, $\ker D^\dagger$ is $V_S\otimes K_S$ for a multiplicity subspace $K_S\subseteq M_S$, and its orthogonal projector has the form $I_{V_S}\otimes P_{K_S}$. The normalized ladder maps act as the standard weight-space ladder on $V_S$ tensored with $I_{M_S}$. Therefore the kernel projector commutes with those ladder maps, and distinct total-pseudospin sectors remain orthogonal after projection.

Lemma 2 applied componentwise gives the direct-sum source decomposition

$$
b_i^{(n)}
=n|n\rangle\otimes\vec\beta_i
+\sqrt{\rho_n}\,U_n^\oplus b_{i,\perp}^{(1)},
\tag{29}
$$

where $\vec\beta_i=(\beta_{1i},\ldots,\beta_{M_si})$. The two terms lie in total pseudospin $S$ and $S-1$, respectively.

The longitudinal term cannot be removed by the least-squares correction. Indeed,

$$
D_n^\dagger(|n\rangle\otimes\vec v)
=\sum_\lambda v_\lambda S_\lambda^\dagger|n\rangle=0,
\tag{30}
$$

where Hermiticity and $S_\lambda|n\rangle=0$ were used. Hence $\Pi_n$ acts as the identity on the longitudinal source. On the $S-1$ component, intertwining gives

$$
\Pi_nU_n^\oplus=U_n^\oplus\Pi_1.
\tag{31}
$$

The two pseudospin sectors are orthogonal. Lemma 1 therefore yields

$$
C_{ij}^{(n)}
=n^2\Gamma_{ij}
+\rho_n\bigl(C_{ij}^{(1)}-\Gamma_{ij}\bigr),
\tag{32}
$$

which is (27). $\square$

**Interpretation.** The $S-1$ response channel is the genuine pair-mobility channel. Its reduced matrix element is fixed in the one-pair sector and acquires the hard-core factor $\rho_n$. The only possible violation of the exact finite-size filling law in this class is the maximal-pseudospin trace component. It is nonnegative and kinematic, not a new finite-density dressing channel.

# 7. Translation-invariant positive kernels

The factorized notation (5) is convenient for the proof, but locality is clearer before factorizing the interaction kernel. Let $\Lambda_{\mathbf N}$ be a $d$-dimensional rectangular torus with side lengths $\mathbf N=(N_1,\ldots,N_d)$ and volume

$$
V=\prod_{j=1}^dN_j.
\tag{33}
$$

Suppose there are finitely many translated Hermitian QGN null operators $X_{R,a}(A)$ and an $A$-independent positive translation kernel $\mathcal V_{ab}(R-R')$:

$$
H(A)=\frac12\sum_{R,R'\in\Lambda_{\mathbf N}}
\sum_{a,b}
X_{R,a}(A)^\dagger
\mathcal V_{ab}(R-R')
X_{R',b}(A).
\tag{34}
$$

Let $\widehat{\mathcal V}(q)$ be its Fourier transform. Positivity means $\widehat{\mathcal V}(q)\succeq0$ for every finite-torus momentum.

Translation invariance of $|1\rangle$ makes

$$
\beta_{a,i}^{(\mathbf N)}
=\left\langle1\left|
\partial_{A_i}X_{0,a}(A)|_0
\right|1\right\rangle
\tag{35}
$$

independent of the cell. The longitudinal tensor is then

$$
\boxed{
\Gamma_{ij}^{(\mathbf N)}
=V\,\operatorname{Re}
\left[(\vec\beta_i^{(\mathbf N)})^\dagger
\widehat{\mathcal V}(0)
\vec\beta_j^{(\mathbf N)}\right].
}
\tag{36}
$$

Equation (36) is invariant under any real Hermitian square-root factorization of the positive kernel. The existence of such a factorization is included in the Hermitian-kernel hypothesis.

# 8. Winding formula in arbitrary dimension

Let

$$
\mathcal K_{\mathbf N}
=\left\{
\left(\frac{2\pi m_1}{N_1},\ldots,\frac{2\pi m_d}{N_d}\right):
0\le m_j<N_j
\right\}.
\tag{37}
$$

Assume the specified flat-connection lift is translation covariant in the following sense: for each local channel, the normalized one-pair trace can be written as a discrete average of a globally defined periodic scalar symbol $f_a(k)$:

$$
\langle1|X_{0,a}(A)|1\rangle
=\frac1V\sum_{k\in\mathcal K_{\mathbf N}}f_a(k+A)+\text{constant}.
\tag{38}
$$

Here $f_a$ is a trace over the flat-band fiber, normalized by the number of pair orbitals per cell. It is gauge invariant.

Write the Fourier series

$$
f_a(k)=\sum_{R\in\mathbb Z^d}\widehat f_a(R)e^{ik\cdot R}.
\tag{39}
$$

**Proposition 2 (multidimensional roots-of-unity filter).** The finite-grid average obeys

$$
\frac1V\sum_{k\in\mathcal K_{\mathbf N}}f_a(k+A)
=
\sum_{\ell\in\mathbb Z^d}
\widehat f_a(\ell\odot\mathbf N)
\,e^{i(\ell\odot\mathbf N)\cdot A},
\tag{40}
$$

where $\ell\odot\mathbf N=(\ell_1N_1,\ldots,\ell_dN_d)$. Consequently,

$$
\boxed{
\beta_{a,i}^{(\mathbf N)}
=
\sum_{\ell\ne0}
 i\ell_iN_i\,
\widehat f_a(\ell\odot\mathbf N).
}
\tag{41}
$$

*Proof.* Substitute (39) into the left side of (40). The product of the $d$ discrete roots-of-unity sums is one exactly when every component of $R$ is divisible by the corresponding $N_j$, and zero otherwise. Differentiation gives (41). $\square$

Thus the longitudinal term is a pure winding observable. It is absent in the Brillouin-zone integral and can survive at finite size only through Fourier coefficients at displacement vectors that wrap around the torus.

# 9. Locality and topological projector bundles

**Theorem 2 (locality bound).** Assume a fixed number of local QGN channels and a uniform bound

$$
\|\widehat{\mathcal V}(0)\|\le V_0.
\tag{42}
$$

Let $N_{\min}=\min_jN_j$ and $N_{\max}=\max_jN_j$.

1. **Strict finite range.** If there is an $R_0$ independent of $\mathbf N$ such that

$$
   \widehat f_a(R)=0\quad\text{whenever }|R|_\infty>R_0,
   \tag{43}
$$

   then for every torus with $N_{\min}>R_0$,

$$
   \beta_{a,i}^{(\mathbf N)}=0,
   \qquad
   \Gamma_{ij}^{(\mathbf N)}=0.
   \tag{44}
$$

2. **Uniform exponential locality.** If

$$
   |\widehat f_a(R)|\le C_0e^{-\mu|R|_1}
   \tag{45}
$$

   with $C_0,\mu>0$ independent of $\mathbf N$, then there are constants $C_1,C_2$ depending only on $d,\mu$, the channel number, and $V_0$ such that, for all sufficiently large $N_{\min}$,

$$
   |\beta_{a,i}^{(\mathbf N)}|
   \le C_1N_{\max}e^{-\mu N_{\min}},
   \tag{46}
$$

   and

$$
   \|\Gamma^{(\mathbf N)}\|
   \le C_2V N_{\max}^2e^{-2\mu N_{\min}}.
   \tag{47}
$$

*Proof.* Part 1 follows immediately from (41): every nonzero $\ell\odot\mathbf N$ has $|\ell\odot\mathbf N|_\infty\ge N_{\min}$.

For Part 2, (41) gives

$$
|\beta_{a,i}^{(\mathbf N)}|
\le C_0N_i
\sum_{\ell\ne0}|\ell_i|
\exp\left[-\mu\sum_j|\ell_j|N_j\right].
\tag{48}
$$

The product geometric series is bounded by a constant times $e^{-\mu N_{\min}}$ for large $N_{\min}$, proving (46). Substituting (46) into (36) gives (47). $\square$

## 9.1 Why topology does not obstruct the estimate

No global Bloch frame is used in (38)-(47). On a gauge patch, the projected one-body symbol is a bundle endomorphism $x_a(k)$. On an overlap it transforms by conjugation,

$$
x_a^{(q)}(k)=g_{pq}(k)^{-1}x_a^{(p)}(k)g_{pq}(k).
\tag{49}
$$

The normalized trace $f_a(k)$ is therefore patch independent and defines a global periodic scalar even when the flat-band vector bundle is topologically nontrivial. The Fourier estimate applies directly to this scalar. A Chern obstruction to a global orthonormal Wannier frame is not an obstruction to a globally defined projector or to gauge-invariant trace symbols.

The locality hypothesis must nevertheless be imposed on the *projected symbol and its specified displacement lift*. A finite-range microscopic parent does not automatically make the spectrally flattened projector strictly finite range; in generic isolated bands one expects exponential rather than compact support.

# 10. Normalization: electronic twist, pair momentum, stiffness, and pair mass

We now freeze conventions. Set the lattice constants, electron charge magnitude, and $\hbar$ to one.

1. $A$ is a uniform electronic crystal-momentum shift, implemented as $k\mapsto k+A$ in the projected interaction.
2. The boundary phase is $\phi_i=N_iA_i$. Thus

$$
   \partial_{A_i}\partial_{A_j}
   =N_iN_j\partial_{\phi_i}\partial_{\phi_j}.
   \tag{50}
$$

3. A pair contains two electrons, so its center-of-mass momentum is

$$
   Q=2A.
   \tag{51}
$$

4. Define the one-pair dispersion by

$$
   E_{\mathrm{pair}}(Q)
   =E_{\mathrm{pair}}(0)
   +\frac12Q_i(m_{\mathrm{pair}}^{-1})_{ij}Q_j+o(|Q|^2).
   \tag{52}
$$

**Lemma 3 (flat connection is pair momentum).** For a translation-covariant projected interaction in which the electronic flat connection is inserted as $u_\sigma(k)\mapsto u_\sigma(k+A)$, the zero-label one-pair block is the untwisted two-particle block at total momentum $Q=2A$:

$$
H_1(A)|_{K=0}\cong H_{\mathrm{pair}}(Q=2A),
\qquad
E_1(A)=E_{\mathrm{pair}}(2A).
\tag{52a}
$$

Consequently,

$$
   \boxed{C_{ij}^{(1)}=4(m_{\mathrm{pair}}^{-1})_{ij}.}
   \tag{53}
$$

*Proof.* A zero-label pair under the electronic shift has physical constituent momenta $k+A$ and $-k+A$, whose sum is $2A$. Relabeling the relative momentum identifies the complete matrix block, not merely its lowest eigenvalue. Translation labels remain good for the flat-connection family, and H4 makes the zero-label ground state simple and isolated; its Kato continuation therefore stays in that block for sufficiently small $A$. Two derivatives then give (53). $\square$

## 10.1 Exact Model-II pair dispersion and positive mass

For Model II,

$$
H_{\mathrm{II}}(k)
=-t[\sigma_x\sin\alpha_k+\sigma_y\cos\alpha_k],
\qquad
\alpha_k=\xi(\cos k_x+\cos k_y).
\tag{M1}
$$

Choose the smooth lower-band gauge

$$
u_\uparrow(k)=\frac1{\sqrt2}
\begin{pmatrix}
e^{i\alpha_k/2}\\ i e^{-i\alpha_k/2}
\end{pmatrix},
\qquad
u_\downarrow(k)=u_\uparrow^*(-k).
\tag{M2}
$$

The two orbital weights are pointwise constant, $|u_{1\sigma}(k)|^2=|u_{2\sigma}(k)|^2=1/2$.

**Lemma 4 (two-line Model-II pair dispersion).** In the one-up/one-down sector on a rectangular torus, the lowest branch at total momentum $Q$ is

$$
\boxed{
E_{\mathrm{pair},\mathbf N}(Q)
=\frac{|U|}{4}[1-|F_{\mathbf N}(Q)|],
\qquad
F_{\mathbf N}(Q)=\frac1V\sum_{k\in\mathcal K_{\mathbf N}}
e^{i[\alpha_{k+Q/2}-\alpha_{k-Q/2}]}.
}
\tag{M3}
$$

For $N_x,N_y\ge3$ its exact Hessian is

$$
\boxed{
(m_{\mathrm{pair},\mathbf N}^{-1})_{ij}
=\left.\partial_{Q_i}\partial_{Q_j}E_{\mathrm{pair},\mathbf N}(Q)\right|_{Q=0}
=\frac{|U|\xi^2}{8}\delta_{ij}\succ0
\quad(\xi\ne0).
}
\tag{M4}
$$

*Proof.* In the relative-momentum basis $|k;Q\rangle$, the PSD-Hubbard one-pair block is

$$
H_Q=|U|\left[\frac12I-\frac1V\sum_{a=1}^2
|g_{a,Q}\rangle\langle g_{a,Q}|\right],
\tag{M5}
$$

where, with $\delta\alpha_k(Q)=\alpha_{k+Q/2}-\alpha_{k-Q/2}$,

$$
g_{1,Q}(k)=\frac12e^{-i\delta\alpha_k(Q)/2},
\qquad
g_{2,Q}(k)=\frac12e^{+i\delta\alpha_k(Q)/2}.
\tag{M6}
$$

The two nonzero eigenvalues of the rank-two kernel in (M5) are those of

$$
\frac14
\begin{pmatrix}
1&F_{\mathbf N}(Q)\\
F_{\mathbf N}(Q)^*&1
\end{pmatrix},
\tag{M7}
$$

namely $(1\pm|F_{\mathbf N}|)/4$. The corresponding energies are $|U|(1\mp|F_{\mathbf N}|)/4$, while every state orthogonal to the rank-two kernel has energy $|U|/2$. Since $0\le|F_{\mathbf N}|\le1$, the branch in (M3) is the lowest one. Moreover,

$$
\delta\alpha_k(Q)
=-2\xi\sum_{i=x,y}\sin k_i\sin(Q_i/2).
\tag{M8}
$$

On every $N_i\ge3$ momentum grid, $V^{-1}\sum_k\sin k_i=0$ and $V^{-1}\sum_k\sin k_i\sin k_j=\delta_{ij}/2$. Hence $F(0)=1$, $\partial_iF(0)=0$, and $\partial_i\partial_jF(0)=-\xi^2\delta_{ij}/2$. Differentiating $|F|$ at the nonzero point $F(0)=1$ gives (M4). $\square$

**Remark (the $2\times2$ anomaly).** The hypothesis $N_i\ge3$ is substantive. On a two-point cycle $k_i\in\{0,\pi\}$, one has $\sin k_i=0$ identically, so the corresponding component of the pair-mass Hessian vanishes. On a $2\times2$ torus both components vanish. This exactly explains the zero-curvature anomaly seen in the earliest finite-size Model-II calculations; it is a discrete-grid exception, not a failure of the formula above.

In the thermodynamic limit the form factor factorizes explicitly,

$$
F_\infty(Q)
=J_0[2\xi\sin(Q_x/2)]J_0[2\xi\sin(Q_y/2)],
\tag{M9}
$$

so

$$
E_{\mathrm{pair}}(Q)
=\frac{|U|\xi^2}{16}(Q_x^2+Q_y^2)+O(|Q|^4),
\tag{M10}
$$

which independently reproduces (M4). **Corollary 4.1 (exact absence of the longitudinal correction).** The pointwise constant orbital weights imply that every longitudinal trace derivative vanishes. Therefore

$$
\boxed{\Gamma_{ij}^{\mathrm{Model\ II}}=0}
\tag{M11}
$$

on every finite torus. Equation (M11) is exact: there is no winding, locality, or finite-size correction to estimate.

With these conventions, define the finite-volume canonical Drude curvature

$$
   \kappa_{ij}^{(n)}(V)=\frac1{4V}C_{ij}^{(n)}.
   \tag{54}
$$

Let $N_{\mathrm{flat}}$ count flat bands including spin. The number of pair orbitals is

$$
L=\frac{N_{\mathrm{flat}}}{2}V,
\qquad
\nu=\frac nL=\frac{2n}{N_{\mathrm{flat}}V}.
\tag{55}
$$

Substituting (27) and (53) into (54) gives the exact finite-size formula

$$
\boxed{
\kappa_{ij}^{(n)}
=
\frac{N_{\mathrm{flat}}}2\nu(1-\nu)
\frac{L}{L-1}
(m_{\mathrm{pair}}^{-1})_{ij}
+
\frac{n^2-\rho_n}{4V}\Gamma_{ij}.
}
\tag{56}
$$

When $\Gamma=0$, the only difference from the thermodynamic conjecture is the universal canonical factor $L/(L-1)$.

For an anisotropic tensor, compare along a fixed real direction $v$. Define

$$
R_{\mathbf N,v}
=
\frac{v_i\kappa_{ij}^{(n)}v_j}
{(N_{\mathrm{flat}}/2)\nu(1-\nu)
 v_i(m_{\mathrm{pair}}^{-1})_{ij}v_j}.
\tag{57}
$$

Using $C^{(1)}=4m_{\mathrm{pair}}^{-1}$,

$$
\boxed{
R_{\mathbf N,v}
=
\frac{L}{L-1}
+
\frac{L^2(n-1)}{(L-1)(L-n)}
\frac{v_i\Gamma_{ij}v_j}{v_iC_{ij}^{(1)}v_j}.
}
\tag{58}
$$

# 11. Thermodynamic restricted reduction theorem

**Theorem 3 (thermodynamic QGN reduction).** Consider a sequence of regular tori with $N_{\min}\to\infty$ and bounded aspect ratio. Assume H1-H4 at each size, a fixed finite-range or uniformly exponentially local lift as in Theorem 2, and fillings $n_{\mathbf N}/L\to\nu\in(0,1)$. Suppose that in a direction $v$ the one-pair inverse mass remains nonzero:

$$
\liminf_{\mathbf N\to\infty}
 v_i(m_{\mathrm{pair}}^{-1})_{ij}v_j>0.
\tag{59}
$$

Then

$$
\boxed{
R_{\mathbf N,v}\longrightarrow1.
}
\tag{60}
$$

For fixed strict finite range, $\Gamma^{(\mathbf N)}=0$ exactly once the torus exceeds the support, and

$$
R_{\mathbf N,v}=\frac{L}{L-1}
\tag{61}
$$

at every larger size. For uniform exponential locality, the second term in (58) is bounded by a polynomial in the side lengths times $e^{-2\mu N_{\min}}$, and therefore vanishes.

*Proof.* The finite-range statement follows from (44), (56), and (58). In the exponential case, (47) and $L=(N_{\mathrm{flat}}/2)V$ show that the longitudinal correction in (58) is bounded by a polynomial factor times $e^{-2\mu N_{\min}}$. Assumption (59) prevents division by a vanishing pair curvature. The first term in (58) tends to one. $\square$

Gao, Han, and Khalaf define their target response through the second derivative of the ground-state energy under a flat connection and denote the thermodynamic quantity $D_s$ [1]. With that literal definition, no additional charge-gap hypothesis is needed to state their conjecture. Their finite-size calculations are more precisely Drude curvatures, and a stricter dynamical identification with a zero-frequency transport coefficient may still require an order-of-limits argument.

## 11.1 Headline application: Model II Hubbard

**Theorem 4 (Gao-Han-Khalaf equality for connected Model-II tori).** Consider the lower spin-degenerate flat band of (M1), projected with the PSD-Hubbard interaction

$$
H(A)=\frac{|U|}{2}\sum_{R,a}
[\bar n_{Ra\uparrow}(A)-\bar n_{Ra\downarrow}(A)]^2,
\tag{M12}
$$

and the flat connection $k\mapsto k+A$. Let $\xi\ne0$ and let $N_x,N_y\ge3$. On any rectangular torus for which the projected-frame graph is connected, the fixed-number ground branch is unique and analytic and, for $1\le n\le V-1$, its exact curvature is

$$
\boxed{
\kappa_{ij}^{(n)}(V)
=\frac1{4V}C_{ij}^{(n)}
=\frac{V}{V-1}\nu(1-\nu)
\frac{|U|\xi^2}{8}\delta_{ij},
\qquad \nu=\frac nV.
}
\tag{M13}
$$

The mass formula (M4) holds for every $\xi\ne0$; the restriction $0<\xi<\pi/2$ enters only through the all-size connectivity proof. In particular, for $0<\xi<\pi/2$ the frame graph is connected on every rectangular torus, and for any regular thermodynamic sequence with $n_V/V\to\nu\in(0,1)$,

$$
\boxed{
D_{s,ij}
:=\lim_{V\to\infty}\kappa_{ij}^{(n_V)}(V)
=\nu(1-\nu)(m_{\mathrm{pair}}^{-1})_{ij}
=\frac{|U|\xi^2}{8}\nu(1-\nu)\delta_{ij}.
}
\tag{M14}
$$

Thus the finite-size theorem applies at any audited-connected value of $\xi\ne0$, while the unconditional all-torus thermodynamic statement is proven throughout $0<\xi<\pi/2$.

*Proof.* Time reversal gives H1-H3. Connectedness of the projected-frame graph gives H4 by Proposition 1. For $0<\xi<\pi/2$ that connectedness is proven uniformly for every rectangular torus; outside this interval the same finite-size conclusion holds whenever connectedness is independently certified. Equation (M11) removes the longitudinal term exactly. With $L=V$ and $N_{\mathrm{flat}}=2$, Theorem 1 and Lemmas 3-4 give

$$
C_{ij}^{(n)}
=\frac{n(V-n)}{V-1}C_{ij}^{(1)}
=\frac{n(V-n)}{V-1}\frac{|U|\xi^2}{2}\delta_{ij}.
\tag{M15}
$$

Dividing by $4V$ proves (M13); taking $V\to\infty$ proves (M14). $\square$

# 12. Finite-volume isolation and the deferred sector-gap problem

Proposition 1 gives a unique AGP zero mode and a strictly positive spectral gap on every fixed connected torus. That finite-volume isolation is the only gap input used by Lemma 1 and Theorem 1. Neither Theorem 3 nor Theorem 4 assumes a lower bound uniform in volume.

The representation-theoretic decomposition of the full sector spectrum, comparisons of gaps across fillings, and the rank-one Aldous question are mathematically interesting but are not load-bearing for the stiffness theorem. They are deferred to a standalone mathematical note. No finite-size sector-gap ordering or $O(N^{-2})$ scaling claim is asserted or used in this paper.

This separation is deliberate: the Gao-Han-Khalaf conjecture addressed by Theorem 4 is stated in terms of thermodynamic flat-connection ground-energy curvature. A stricter dynamical identification and its order of limits are scoped as a subsequent project in Section 13.3.

# 13. Relation to prior work and theorem boundary

## 13.1 What is new relative to earlier flat-band results

Several ingredients of the proof have important precedents.

Tovmasyan *et al.* established exact paired ground states and an emergent pseudospin $SU(2)$ structure for projected attractive Hubbard models under the uniform-pairing condition [5]. Huhtinen *et al.* derived an exact Cooper-pair mass in that setting and related it to the minimal quantum metric, resolving orbital-embedding dependence [10]. The QGN construction generalized the solvable-ground-state framework and computed few-body excitations and pseudospin-based stiffness information for ideal models [2]. More recently, Herzog-Arbeitman *et al.* analytically classified pair bound-state spectra in exact eta-pairing flat-band models and identified the minimal metric in the quadratic pair dispersion [11]. Gao, Han, and Khalaf then supplied rigorous many-body lower bounds, the pair-mass upper bound, and the equality conjecture addressed here [1].

A targeted prior-art sweep through July 23, 2026 found no earlier proof of the exact Model-II many-body stiffness/pair-mass identity. The closest results establish exact paired ground states, exact or geometric pair masses, solvable pair spectra, stiffness bounds, and numerical saturation, but not the finite-density equality proved here.

The present contribution is different in four specific respects:

1. it proves an exact fixed-particle-number twist-Hessian reduction, including the finite-size factor $L/(L-1)$ and the complete nonnegative defect tensor $\Gamma$;
2. it proves that the defect is a winding trace and therefore vanishes exactly for sufficiently large fixed finite-range lifts and asymptotically for uniformly exponentially local lifts;
3. it gives a necessary-and-sufficient connectivity criterion for uniqueness of the AGP zero mode in projected PSD-Hubbard models; and
4. it combines these results with an exact Model-II pair branch to prove the GHK equality, rather than infer it from mean field, a variational upper bound, or numerical saturation.

The connectivity theorem is a superconducting operator-space analogue of the irreducibility criteria in rigorous flat-band ferromagnetism. Mielke proved necessary-and-sufficient uniqueness conditions in terms of irreducibility of the one-particle density matrix [8,9]. Here vectorization converts the paired zero-mode problem into a commutant on $\Lambda^n\mathcal V$, and connected rank-one projected orbitals generate the full one-body unitary algebra.

## 13.2 Non-Hermitian QGN boundary

Perfect QGN is a property of the flat-band electronic structure; it should not be conflated with the narrower interaction class used in this theorem. The standard ideal Hamiltonians of Ref. [2] are built from Hermitian local operators and a positive kernel. In its “Possible Extensions” section, the same work explicitly introduces an SSH-type pairing construction with a non-Hermitian interaction coefficient/form, describes it as a generalization of the standard form, and notes that its excitations are not exactly solvable [2].

Accordingly:

- a band used in such an extension may still satisfy perfect QGN;
- the associated interaction need not admit the Hermitian positive-square least-squares and pseudospin-intertwiner structure used here;
- the large discrepancy found in our non-Hermitian generalized-nesting stress test is evidence that Hermiticity is a substantive hypothesis, but that random test was not a reproduction of the published SSH extension and is not advertised as a counterexample to it.

The restricted theorem therefore proves the Hubbard Model-II case and the stated Hermitian positive-square class. It does not adjudicate the equality for every non-Hermitian extension that might fall under the broadest phrase “QGN model.”

## 13.3 Remaining scope

The theorem does not presently cover degenerate fixed-number ground spaces, several inequivalent singular-value blocks of the pairing matrix, rank-deficient interaction kernels without a separate Lie-algebra audit, or general non-Hermitian generalized-QGN constructions. The reducible counterexample classifies the first boundary exactly. The natural replacement in the multi-block case is a block-resolved filling law with several pseudospins.

Ref. [1] calls the thermodynamic ground-energy curvature under a flat connection $D_s$, and that literal conjecture is what Theorem 4 proves for Model II. A stricter identification with a dynamical zero-frequency superfluid response can depend on the order of thermodynamic, zero-momentum, and zero-frequency limits. No uniform-in-volume gap is assumed or needed here, and no sector-gap scaling statement is used. We leave the dynamical order-of-limits problem as the next project rather than attaching it as an unnecessary hypothesis to the result proved here.

# 14. Reproducibility and cross-package verification

The consolidated archive contains independent certificates for every stage of the argument.

| Component | Principal check |
|---|---|
| Winding obstruction | Exact SymPy Kohn response; independent fixed-number projected-density implementation; resonant and nonresonant controls. |
| Reducible obstruction | Fixed-local range-two ED; exact zero-curvature component-polarized state; positive one-pair mass. |
| General reduction | Random full-rank multiband skew-unitary pairing matrix; tensor identity residual $<1.3\times10^{-18}$. |
| Longitudinal formula | Independent trace derivative reproduces all filling defects and the winding-model PSD defect. |
| H4 connectivity | Random path-like and block-linked Parseval frames; exact zero-mode counts; Model-I/II graph audit. |
| Model-II closure | Direct $Q$-block versus closed pair formula, $Q=2A$, mass Hessian, exact $\Gamma=0$, and many-body curvature. |
| Published models | Independent projected-Hubbard ED for Gao-Han-Khalaf Models I and II at several finite tori. |

The Model-II certificate agrees with the direct pair block to $4.3\times10^{-16}$, verifies $E_1(A)=E_{\mathrm{pair}}(2A)$ to $5.9\times10^{-16}$, and reproduces the many-body finite-size formula to approximately $2.5\times10^{-9}$, the finite-difference floor. The H4 theorem was independently stress-tested on minimal path graphs and on two large frame blocks joined by a $10^{-6}$ overlap; both generate the full $\mathfrak u(L)$ as predicted.

The original counterexamples are now boundary tests of the final theorem rather than disconnected side projects. The winding defect is exactly the longitudinal tensor in a nonlocal resonant lift, and the even/odd counterexample is exactly the disconnected ground-space degeneracy counted by the H4 theorem.

# 15. Conclusion

The research arc began with a broad attempt to prove an exact many-body/two-body reduction and initially produced two counterexamples. The first showed that continuous finite-torus twist curvature is not determined by untwisted samples unless the displacement lift is fixed. The second showed that even a fixed-local, twist-stable QGN model can evade the pair-mass value when its paired ground space is reducible. Those failures identified, rather than destroyed, the theorem.

Under H1-H4 in the Hermitian positive-square class, the exact finite-size identity is

$$
\boxed{
C_{ij}^{(n)}
=
\frac{n(L-n)}{L-1}C_{ij}^{(1)}
+
\left[
n^2-\frac{n(L-n)}{L-1}
\right]\Gamma_{ij},
\qquad
\Gamma\succeq0.
}
\tag{62}
$$

The longitudinal tensor is a winding observable. Fixed locality removes it. In the projected PSD-Hubbard subclass, projected-frame connectivity is exactly the condition that supplies H4 by removing competing AGP occupancy sectors. In Model II, the cancellation is stronger: $\Gamma=0$ exactly on every finite torus, with no correction, while the one-pair mass is strictly positive and known in closed form. The resulting finite-size formula (M13) and thermodynamic equality (M14) prove the restricted Gao-Han-Khalaf conjecture for the projected-Hubbard Model II throughout $0<\xi<\pi/2$, and on any additional finite torus whose frame connectivity is certified.

The next physics project is extension, not completion of this case: determine the largest Hermitian QGN interaction class for which the intertwiner proof survives, derive the multi-singular-value block law, and address the dynamical order of limits under transport definitions stricter than the flat-connection curvature used in Ref. [1].

# Acknowledgments

Nicholas Sledgianowski guided the research direction, model selection, independent-review process, and interpretation. ChatGPT Sol 5.6 provided computational, algebraic, drafting, and document-preparation assistance. Independent numerical and expert reviews materially improved the scope, proof presentation, H4 theorem, and Model-II closure.

# References

[1] Q. Gao, Z. Han, and E. Khalaf, “Bootstrapping Flatband Superconductors: Rigorous Lower Bounds on Superfluid Stiffness,” *Physical Review Letters* **136**, 076503 (2026), DOI 10.1103/gw85-5r92; arXiv:2506.18969v4.

[2] Z. Han, J. Herzog-Arbeitman, B. A. Bernevig, and S. A. Kivelson, “Quantum Geometric Nesting and Solvable Model Flat-Band Systems,” *Physical Review X* **14**, 041004 (2024); arXiv:2401.04163v3.

[3] T. Kato, *Perturbation Theory for Linear Operators*, 2nd ed., Springer-Verlag, Berlin-Heidelberg-New York (1976).

[4] D. J. Scalapino, S. R. White, and S. Zhang, “Insulator, Metal, or Superconductor: The Criteria,” *Physical Review B* **47**, 7995-8007 (1993).

[5] M. Tovmasyan, S. Peotta, P. Törmä, and S. D. Huber, “Effective Theory and Emergent SU(2) Symmetry in the Flat Bands of Attractive Hubbard Models,” *Physical Review B* **94**, 245149 (2016).

[6] C. Brouder, G. Panati, M. Calandra, C. Mourougane, and N. Marzari, “Exponential Localization of Wannier Functions in Insulators,” *Physical Review Letters* **98**, 046402 (2007).

[7] E. H. Lieb, “Two Theorems on the Hubbard Model,” *Physical Review Letters* **62**, 1201-1204 (1989).

[8] A. Mielke, “Ferromagnetism in the Hubbard Model and Hund’s Rule,” *Physics Letters A* **174**, 443-448 (1993), DOI 10.1016/0375-9601(93)90207-G.

[9] A. Mielke, “Stability of Ferromagnetism in Hubbard Models with Degenerate Single-Particle Ground States,” *Journal of Physics A: Mathematical and General* **32**, 8411-8418 (1999), DOI 10.1088/0305-4470/32/48/304.

[10] K.-E. Huhtinen, J. Herzog-Arbeitman, A. Chew, B. A. Bernevig, and P. Törmä, “Revisiting Flat Band Superconductivity: Dependence on Minimal Quantum Metric and Band Touchings,” *Physical Review B* **106**, 014518 (2022).

[11] J. Herzog-Arbeitman, A. Chew, K.-E. Huhtinen, P. Törmä, and B. A. Bernevig, “Many-Body Superconductivity in Topological Flat Bands,” *Communications Physics*, published online 4 July 2026, DOI 10.1038/s42005-026-02732-2.

[12] S. Peotta and P. Törmä, “Superfluidity in Topologically Nontrivial Flat Bands,” *Nature Communications* **6**, 8944 (2015).

[13] P. Törmä, S. Peotta, and B. A. Bernevig, “Superconductivity, Superfluidity and Quantum Geometry in Twisted Multilayer Systems,” *Nature Reviews Physics* **4**, 528-542 (2022).

[14] W. Kohn, “Theory of the Insulating State,” *Physical Review* **133**, A171-A181 (1964).

[15] J. S. Hofmann, E. Berg, and D. Chowdhury, “Superconductivity, Charge Density Wave, and Supersolidity in Flat Bands with a Tunable Quantum Metric,” *Physical Review Letters* **130**, 226001 (2023).

