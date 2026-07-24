---
title: "Zero-Twist Uniform Pairing Does Not Fix Finite-Torus Twist Curvature"
subtitle: "A resonant winding counterexample, a nonresonant control, and the role of the twist prescription"
author: "ChatGPT Sol 5.6 - guided by Nicholas Sledgianowski"
date: "July 2026"
lang: en-US
---

*Revised research note. The four-cell calculation has been independently reproduced in a separate fixed-number bit-basis implementation. This manuscript is not peer reviewed.*

# Abstract

A proposed exact finite-size extension of flat-band pair-mass relations asserts that the fixed-pair-number boundary-twist curvature is determined by the one-pair sector, with all filling dependence contained in the hard-core factor $\nu(1-\nu)$. We show that the zero-twist uniform-pairing condition (UPC) is not sufficient for that statement. On a four-cell ring, a smooth two-orbital rank-one flat-band projector contains a harmonic proportional to $\sin 4k$ that vanishes on the untwisted momentum grid but becomes momentum independent after a continuous twist. For the projected onsite attraction this gives the exact operator identity

$$
H_c(A)=\cos^2\!\left[2c\sin(4A)\right]H_0(A).
$$

At $c=1/8$, exact symbolic response theory yields

$$
E_n''(0)=U\left[\frac{\pi^2}{24}n(4-n)+n\right],
$$

so the two-pair raw-curvature defect is $2U/3$ and the normalized defect is $U/24$. A second fixed-number implementation confirms the ground energies, the uniform $U/4$ gap, the operator factorization, and the curvature values.

The same Bloch data at $M=8$ provide a decisive control. The hidden harmonic aliases to $(-1)^m\sin4A$, rather than a common scalar, and the filling law is restored numerically to approximately $10^{-9}$. Replacing $\sin4k$ by $\sin8k$ restores scalar factorization and the violation. Thus the load-bearing mechanism of this construction is an exact period-$M$ winding resonance, not the mere presence of a high Fourier harmonic.

We then prove two Phase-II results. First, for every $M\ge 3$, a period-$M$ deformation produces an exact UPC counterexample with an analytic, gapped but degenerate paired ground manifold and normalized defect

$$
K_n-\rho_nK_1
=
Uc^2q^2M\frac{n(n-1)}{M-1},
\qquad
\rho_n=\frac{n(M-n)}{M-1}.
$$

Second, finite-torus samples determine only alias sums of displacement-resolved Fourier coefficients; their twist derivatives depend on the unresolved winding moments. This proves that twist curvature is not determined by the untwisted finite-torus Hamiltonian alone unless a displacement-resolved lift, or an equivalent family $P_A$, is specified. A roots-of-unity identity also characterizes twist-stable UPC exactly: it fails precisely when diagonal projector weights contain Fourier modes at nonzero multiples of $M$.

Gao, Han, and Khalaf introduce the flat connection through displacement-resolved Peierls phases and derive $k\mapsto k+A$ for the Bloch Hamiltonian and projected band eigenvectors. The winding-aware construction therefore follows their literal prescription once the full infinite-lattice lift is supplied. A nonwinding lift gives the same $H(0)$ but removes the extra curvature. The result establishes a limitation of a UPC-only exact finite-size theorem; it does not refute the QGN conjecture for a fixed local model in the thermodynamic limit.

# 1. Claim under examination and scope

Let $E_n(A)$ be the lowest analytic energy branch in the sector with $n$ spin-singlet pairs on a one-dimensional torus with $M$ unit cells. In the normalization used in the finite-size proposal,

$$
K_n=\frac{1}{4M}E_n''(0),
\qquad
\nu_n=\frac{n}{M}.
$$

The proposed finite-size filling law is

$$
\frac{K_n}{\nu_n(1-\nu_n)}
=
\frac{K_1}{\nu_1(1-\nu_1)}.
$$

Equivalently, defining

$$
\rho_n
=
\frac{\nu_n(1-\nu_n)}{\nu_1(1-\nu_1)}
=
\frac{n(M-n)}{M-1},
$$

one asks whether

$$
E_n''(0)=\rho_n E_1''(0)
$$

holds at every nonsingular filling.

This note addresses that exact finite-size statement for a projected onsite attraction with the twist inserted before projection. It does **not** claim that the thermodynamic relation conjectured for quantum-geometric-nesting (QGN) models is false. The distinction is important: the published work defines stiffness from a boundary twist, studies finite systems, and uses finite-size scaling toward a thermodynamic value, whereas the claim tested here asks for a sharpened exact identity on every arbitrary finite torus.

# 2. Four-cell construction

## 2.1 Flat-band projector

Take a one-dimensional torus with

$$
M=4,
\qquad
k_m=\frac{\pi m}{2},
\qquad
m=0,1,2,3,
$$

and two microscopic orbitals per cell. Let

$$
b=\frac{\pi}{2},
\qquad
c=\frac18,
$$

and define

$$
 u_\uparrow(k)=
\begin{pmatrix}
\cos\!\left(\frac{\pi}{4}+c\sin4k\right)\\[2mm]
 e^{ib\sin k}\sin\!\left(\frac{\pi}{4}+c\sin4k\right)
\end{pmatrix},
$$

with spin down related by time reversal,

$$
 u_\downarrow(k)=u_\uparrow^*(-k).
$$

Set

$$
P_\sigma(k)=|u_\sigma(k)\rangle\langle u_\sigma(k)|,
\qquad
h_\sigma(k)=I_2-P_\sigma(k).
$$

The parent Hamiltonian has exactly flat eigenvalues $0$ and $1$ at every momentum. The active band is therefore isolated and rank one. The projector is smooth, periodic, and analytic; its real-space coefficients decay rapidly, although the model is not strictly finite range.

At every allowed momentum,

$$
\sin 4k_m=0,
$$

so the two orbital weights are exactly equal:

$$
|u_{\sigma,1}(k_m)|^2
=
|u_{\sigma,2}(k_m)|^2
=
\frac12.
$$

The finite-torus UPC therefore holds at zero twist. The continuous Brillouin-zone average also gives equal diagonal projector weights.

## 2.2 Projected interaction and analytic branch

Insert the electronic twist before projection,

$$
P_{\sigma,A}(k_m)=P_\sigma(k_m+A),
$$

and project the onsite attraction,

$$
H_c(A)
=-U\sum_{R,\alpha=1}^{2}
\bar n_{R\alpha\uparrow}(A)
\bar n_{R\alpha\downarrow}(A),
\qquad U>0.
$$

At $A=0$, exact diagonalization gives

$$
E_n(0)=-\frac{Un}{2}
$$

and a unique ground state with a uniform gap $U/4$ for $n=1,2,3$:

| Pair number $n$ | Sector dimension | Ground energy | Gap |
|---:|---:|---:|---:|
| 1 | 16 | $-U/2$ | $U/4$ |
| 2 | 36 | $-U$ | $U/4$ |
| 3 | 16 | $-3U/2$ | $U/4$ |

The finite-dimensional Hamiltonian is analytic in $A$, and the ground eigenvalue is simple at zero twist, so each branch is analytic in a neighborhood of $A=0$.

# 3. Exact period-four factorization

On the four-point grid,

$$
\sin\bigl[4(k_m+A)\bigr]
=
\sin(2\pi m+4A)
=
\sin4A.
$$

Set

$$
x_A=c\sin4A.
$$

The orbital mixing angles are $\pi/4+x_A$ for spin up and $\pi/4-x_A$ for spin down. The products entering every projected interaction matrix element satisfy

$$
\cos^2\!\left(\frac\pi4+x_A\right)
\cos^2\!\left(\frac\pi4-x_A\right)
=
\frac14\cos^2(2x_A),
$$

and

$$
\sin^2\!\left(\frac\pi4+x_A\right)
\sin^2\!\left(\frac\pi4-x_A\right)
=
\frac14\cos^2(2x_A).
$$

All remaining momentum-dependent phases are identical to those of the $c=0$ model. Hence the full many-body operator obeys

$$
\boxed{
H_c(A)=f_c(A)H_0(A),
\qquad
f_c(A)=\cos^2\!\left[2c\sin4A\right].
}
$$

For $c=1/8$,

$$
f_c(0)=1,
\qquad
f_c'(0)=0,
\qquad
f_c''(0)=-2.
$$

Every eigenvalue branch factorizes. Since $E_n(0)=-Un/2$,

$$
E_{n,c}''(0)
=
E_{n,0}''(0)+f_c''(0)E_n(0)
=
E_{n,0}''(0)+Un.
$$

The extra term is linear in $n$, rather than proportional to $n(4-n)$.

# 4. Exact curvature and independent verification

For the $c=0$ reference model, exact symbolic Kohn response gives

$$
E_{n,0}''(0)
=
\frac{U\pi^2}{24}n(4-n).
$$

The diamagnetic and paramagnetic contributions are:

| $n$ | Diamagnetic | Paramagnetic | Total $E_{n,0}''(0)$ |
|---:|---:|---:|---:|
| 1 | $U\pi^2/4$ | $-U\pi^2/8$ | $U\pi^2/8$ |
| 2 | $U\pi^2/3$ | $-U\pi^2/6$ | $U\pi^2/6$ |
| 3 | $U\pi^2/4$ | $-U\pi^2/8$ | $U\pi^2/8$ |

Thus

$$
\boxed{
E_n''(0)
=
U\left[\frac{\pi^2}{24}n(4-n)+n\right].
}
$$

At two pairs,

$$
E_2''(0)-\frac43 E_1''(0)=\frac{2U}{3}.
$$

Since $K_n=E_n''/(16)$ at $M=4$,

$$
\boxed{
K_2-\frac43K_1=\frac{U}{24}.
}
$$

The particle-hole-related comparison is also violated:

$$
K_3-K_1=\frac{U}{8}.
$$

A separate fixed-number bit-basis implementation, assembled from projected local densities rather than the symbolic response code, reproduced the following checks:

1. $E_n(0)=-Un/2$ and gap $U/4$ in all three sectors.
2. $H_c(0)=H_0(0)$ exactly up to machine precision.
3. At $A=0.137$, the maximum matrix-element error in the factorization identity is approximately $1.1\times10^{-16}$.
4. The curvature formula
   
   $$
   E_n''(0)=\frac{U\pi^2}{24}n(4-n)+64Uc^2n
   $$
   
   holds numerically for $c=0$, $0.05$, and $1/8$.
5. The exact defects $2U/3$, $U/24$, and $K_3-K_1=U/8$ are reproduced.

These checks certify the arithmetic and isolate the interpretation, rather than the computation, as the remaining issue.

# 5. The $M=8$ control: hidden is not enough

The same Bloch function with the same harmonic $\sin4k$ can be evaluated on an eight-cell torus,

$$
M=8,
\qquad
k_m=\frac{\pi m}{4}.
$$

Again,

$$
\sin4k_m=\sin(\pi m)=0,
$$

so the $c$ deformation is invisible at zero twist and UPC holds. Under a twist, however,

$$
\sin\bigl[4(k_m+A)\bigr]
=(-1)^m\sin4A.
$$

The deformation is now momentum dependent. The best scalar fit to $H_c(A)$ at $A=0.137$ leaves a relative Frobenius residual of approximately $7.4\%$ in the two-pair sector, so the four-cell factorization mechanism is absent.

The numerical curvatures for $c=1/8$ are:

| $n$ | $E_n''(0)/U$ | $\rho_n$ | $E_n''(0)/(U\rho_n)$ |
|---:|---:|---:|---:|
| 1 | $2.233700549$ | $1$ | $2.233700549$ |
| 2 | $3.829200945$ | $12/7$ | $2.233700551$ |
| 3 | $4.78650117$ | $15/7$ | $2.23370055$ |

Within the numerical precision of the independent implementations, the reduced curvature is constant. At $n=2$, the raw defect is below $4\times10^{-9}$.

This control sharpens the interpretation:

> A harmonic that is invisible on the untwisted grid does not automatically violate the filling law. The scalar-factorization counterexample constructed here requires a period-$M$ resonance that aliases the deformation to the same value at every momentum and thereby multiplies the entire interaction operator by a scalar.

This comparison does not prove that winding resonance is the only possible failure mechanism for a UPC-only theorem. It identifies the exact mechanism responsible for the present family and rules out the broader claim that any hidden high harmonic is sufficient.

To test this diagnosis, replace $\sin4k$ by $\sin8k$ on the same eight-cell torus. Then

$$
\sin\bigl[8(k_m+A)\bigr]=\sin8A
$$

is again momentum independent. The scalar factorization returns to machine precision, and the two-pair raw defect is

$$
E_2''(0)-\frac{12}{7}E_1''(0)=\frac{8U}{7}.
$$

The comparison is summarized in Figure 1. The plotted quantity is $E_n''/\rho_n$, which must be independent of $n$ if the filling law holds.

![Reduced curvatures for the four-cell resonance, the eight-cell nonresonant control using the same harmonic, and the restored eight-cell resonance. The $M=8$, $r=4$ control is flat within numerical accuracy; the period-$M$ cases are not.](resonance_comparison.png){width=88%}

# 6. Twist prescription and finite-torus lift

The counterexample is prescription dependent in a precise sense, and that dependence must be stated explicitly.

## 6.1 The prescription used in the motivating work

Gao, Han, and Khalaf define the stiffness from the curvature of the ground-state energy under a flat connection, written in the main text as $k\mapsto k+A$. Their supplemental material then begins with displacement-resolved real-space hopping amplitudes $t_{\alpha\beta,R-R'}$ and introduces the connection through the Peierls phase

$$
 e^{iA\cdot(R+r_\alpha-R'-r_\beta)}.
$$

It then derives

$$
t(k;A)=t(k+A;0)
$$

and uses band eigenvectors evaluated at $k+A$ in the projected interaction. Their explicit projected Hamiltonian likewise contains factors $U(k_i+A)$. Thus, when a model is specified by its full infinite-lattice displacement data, the winding number carried by each hopping is part of the model and the present $k\mapsto k+A$ construction follows the cited prescription.

## 6.2 What the finite-torus Hamiltonian does not determine

If one specifies only the untwisted finite-torus samples $P(k_m)$, the twisted family is not unique. In the present example,

$$
P_c(k_m)=P_0(k_m)
$$

for every allowed $k_m$. Therefore the winding lift $P_c(k)$ and the nonwinding lift $P_0(k)$ define exactly the same finite-torus Hamiltonian at $A=0$, with the same UPC, spectrum, gap, and ground state. Nevertheless,

$$
P_c(k_m+A)\neq P_0(k_m+A)
$$

and their curvatures differ by the linear-in-$n$ term.

This is the more fundamental finite-size statement:

> Boundary-twist curvature is not a function of the untwisted finite-torus Hamiltonian alone. It also depends on a displacement-resolved lift, or equivalently on a specified gauge-covariant family $P_A$.

The alternative is explicit in this example: choose the nonwinding lift $P_0(k)$ rather than $P_c(k)$. It gives the same finite-torus Hamiltonian at $A=0$ but removes the $c$-dependent curvature exactly. A minimal-image construction is another way of selecting such a nonwinding representative when the finite real-space data admit one. This does not invalidate the counterexample under the displacement-resolved winding prescription; it shows that the theorem must say which lift is being held fixed.

# 7. Phase-II obstruction theorem: period-$M$ resonance

The four-cell identity is one member of a general algebraic family. Let $M\ge3$, let $q$ be a positive integer, and define

$$
 u_{\uparrow,c}^{(M,q)}(k)=
\begin{pmatrix}
\cos\!\left(\frac\pi4+c\sin(qMk)\right)\\[2mm]
 e^{i\phi(k)}\sin\!\left(\frac\pi4+c\sin(qMk)\right)
\end{pmatrix},
$$

where $\phi(k)$ is a smooth odd periodic phase, and set $u_\downarrow(k)=u_\uparrow^*(-k)$. At the allowed momenta $k_m=2\pi m/M$,

$$
\sin(qMk_m)=0,
$$

whereas after the twist,

$$
\sin\bigl[qM(k_m+A)\bigr]=\sin(qMA)
$$

is independent of $m$.

## 7.1 Scalar-factorization proposition

**Proposition 1 (period-$M$ scalar factorization).** For the projected cross-attraction,

$$
H_c(A)=f_{M,q,c}(A)H_0(A),
\qquad
f_{M,q,c}(A)=\cos^2\!\left[2c\sin(qMA)\right].
$$

Consequently,

$$
f_{M,q,c}''(0)=-8c^2q^2M^2.
$$

If a zero-twist eigenbranch has energy $E_n(0)=-Un/2$, then

$$
\boxed{
E_{n,c}''(0)
=
E_{n,0}''(0)+4Uc^2q^2M^2n.
}
$$

*Proof.* Under the resonance, the mixing angles $\pi/4\pm c\sin(qMA)$ are momentum independent. For orbital 1 and orbital 2, respectively,

$$
\cos^2\!\left(\frac\pi4+x\right)
\cos^2\!\left(\frac\pi4-x\right)
=
\sin^2\!\left(\frac\pi4+x\right)
\sin^2\!\left(\frac\pi4-x\right)
=
\frac14\cos^2(2x).
$$

All phase-dependent matrix elements are those of the $c=0$ operator, proving the factorization. Differentiating the eigenvalue identity $E_{n,c}(A)=f_{M,q,c}(A)E_{n,0}(A)$ and using $f(0)=1$, $f'(0)=0$, and $E_n(0)=-Un/2$ gives the curvature formula. $\square$

Define

$$
\Delta_n(c)=E_{n,c}''(0)-\rho_nE_{1,c}''(0),
\qquad
\rho_n=\frac{n(M-n)}{M-1}.
$$

Proposition 1 immediately gives

$$
\boxed{
\Delta_n(c)
=
\Delta_n(0)
+
4Uc^2q^2M^2\frac{n(n-1)}{M-1}.
}
$$

Thus any reference family obeying the filling law is driven away from it by every nonzero resonant deformation.

## 7.2 An exact counterexample for every finite size

Set $\phi(k)=0$. At $c=0$, the two orbital components are constant and the projected interaction reduces to

$$
H_0=-\frac{U}{2}\sum_R n_{R\uparrow}n_{R\downarrow}.
$$

In the sector with $n$ particles of each spin, the ground energy is $-Un/2$, the ground-space dimension is $\binom{M}{n}$, and the excitation gap is $U/2$. The untwisted Hamiltonian is independent of $A$, so $E_{n,0}''(0)=0$. Proposition 1 therefore yields the following result.

**Theorem 2 (arbitrary-$M$ UPC obstruction).** For every $M\ge3$, every $q\ge1$, every $c\ne0$, and every $2\le n\le M-1$, the smooth rank-one flat-band family above satisfies zero-twist UPC and has an analytic, gapped paired ground manifold, but violates the filling law by

$$
\boxed{
E_{n,c}''(0)-\rho_nE_{1,c}''(0)
=
4Uc^2q^2M^2\frac{n(n-1)}{M-1}>0.
}
$$

In the normalization $K_n=E_n''/(4M)$,

$$
\boxed{
K_{n,c}-\rho_nK_{1,c}
=
Uc^2q^2M\frac{n(n-1)}{M-1}>0.
}
$$

The ground state is degenerate, but this does not create a branch ambiguity: the entire ground manifold is multiplied by the same positive analytic scalar $f_{M,q,c}(A)$ near $A=0$. The four-cell phase choice $\phi(k)=b\sin k$ shows separately that the obstruction can coexist with a unique gapped branch. Constructing an equally explicit unique-ground-state family for every $M$ remains open.

# 8. Finite sampling, twist-stable UPC, and locality repairs

The resonance theorem can be recast as two elementary aliasing statements. These statements separate the well-posedness issue from the many-body dynamics.

## 8.1 Roots-of-unity characterization of twist-stable UPC

Let a periodic diagonal projector weight have the Fourier series

$$
w_\alpha(k)=\sum_{r\in\mathbb Z}\widehat w_{\alpha,r}e^{irk},
$$

and define its finite-grid twisted average

$$
W_{\alpha,M}(A)
=
\frac1M\sum_{m=0}^{M-1}w_\alpha(k_m+A).
$$

Using

$$
\frac1M\sum_{m=0}^{M-1}e^{irk_m}
=
\begin{cases}
1,&r\in M\mathbb Z,\\
0,&r\notin M\mathbb Z,
\end{cases}
$$

one obtains the exact identity

$$
\boxed{
W_{\alpha,M}(A)
=
\sum_{\ell\in\mathbb Z}
\widehat w_{\alpha,\ell M}e^{i\ell MA}.
}
$$

Zero-twist UPC constrains only $W_{\alpha,M}(0)$. UPC is stable under a continuous twist in a neighborhood of zero precisely when all nonconstant alias modes $\widehat w_{\alpha,\ell M}$, $\ell\ne0$, vanish.

This identity explains all three numerical cases:

- $M=4,r=4$: diagonal weights contain nonzero modes at multiples of four, so finite-grid UPC changes under the twist.
- $M=8,r=4$: the relevant diagonal-weight harmonics are odd multiples of four and therefore are not multiples of eight; their finite-grid average remains exactly constant.
- $M=8,r=8$: nonzero modes at multiples of eight return, and UPC again changes under the twist.

Twist-stable UPC therefore excludes the present resonance family. It is a natural repair, but the current work does not prove that it is sufficient for the full many-body filling law.

## 8.2 Displacement-lift ambiguity

Let a matrix-valued periodic projector be written as

$$
P(k)=\sum_{R\in\mathbb Z}P_Re^{ikR}.
$$

Choose one representative $r$ for each residue class modulo $M$. The untwisted finite-grid samples depend only on the alias sums

$$
\widetilde P_r
=
\sum_{\ell\in\mathbb Z}P_{r+\ell M},
\qquad
P(k_m)=\sum_r e^{ik_mr}\widetilde P_r.
$$

By contrast, the twist derivatives are

$$
\boxed{
\left[\partial_A^jP(k_m+A)\right]_{A=0}
=
\sum_r e^{ik_mr}
\sum_{\ell\in\mathbb Z}
\left[i(r+\ell M)\right]^j P_{r+\ell M}.
}
$$

The $j=0$ data determine only alias sums. The $j=1,2$ responses require the first and second displacement moments within each alias class. Therefore the untwisted finite-torus Hamiltonian does not, by itself, determine its continuous boundary-twist curvature.

In the present construction, the lifts $P_c(k)$ and $P_0(k)$ satisfy

$$
P_c(k_m)=P_0(k_m)
$$

at every allowed momentum. The winding-aware lift retains the displacement-$M$ information and produces the counterexample. The nonwinding lift $P_0$ defines the same $H(0)$ and removes the extra curvature exactly.

## 8.3 What locality does and does not repair

Two exact consequences follow.

First, if all Fourier coefficients vanish outside a fixed representative interval $|R|<M/2$, then the discrete Fourier transform is injective and the finite samples determine the lift uniquely. This stronger minimal-image condition closes the prescription ambiguity.

Second, the specific scalar resonance is impossible whenever the lift has no nonzero Fourier support in the winding classes $R\in M\mathbb Z$. A simpler but weaker sufficient condition is a displacement range strictly smaller than the circumference.

For a sequence of lifts satisfying a uniform exponential bound

$$
\|P_R\|\le Ce^{-\mu|R|},
$$

the unresolved winding contribution to the $j$th projector derivative is bounded by a polynomial in $M$ times $e^{-\mu M/2}$. Thus this particular ambiguity is exponentially suppressed along a fixed exponentially local thermodynamic sequence. The counterexample family avoids that conclusion because its resonant Fourier weight is deliberately moved to displacement $R=M$ as the system size changes.

These locality statements remove the finite-torus winding pathology; they do not prove the desired $n(M-n)$ reduction. A genuinely dynamical condition is still needed to control the diamagnetic and reduced-resolvent terms.

# 9. Robustness to the positive-semidefinite QGN/Hubbard convention

The motivating QGN work writes the projected Hubbard interaction in the positive-semidefinite form

$$
H_{\mathrm{PSD}}(A)
=
\frac{U}{2}\sum_{R,\alpha}
\left[
\bar n_{R\alpha\uparrow}(A)
-
\bar n_{R\alpha\downarrow}(A)
\right]^2.
$$

At zero twist this differs from the cross-attraction by a fixed-number term under UPC. Away from zero twist, that term need not remain constant, so the convention must be tested explicitly.

For the four-cell resonant model at $c=1/8$, the independent bit-basis implementation gives, to numerical precision sufficient to recognize the exact values,

$$
E_{n,\mathrm{PSD}}''(0)
=
U\left[\frac{\pi^2}{24}n(4-n)+2n\right].
$$

Thus

$$
E_{2,\mathrm{PSD}}''-\frac43E_{1,\mathrm{PSD}}''
=
\frac{4U}{3},
$$

and

$$
K_{2,\mathrm{PSD}}-\frac43K_{1,\mathrm{PSD}}
=
\frac{U}{12}.
$$

The resonance is therefore not removed by the positive-semidefinite convention; its curvature contribution is doubled. The $M=8,r=4$ nonresonant control again obeys the filling law numerically, while $M=8,r=8$ violates it.

This robustness check does not convert the result into a refutation of the thermodynamic QGN conjecture. It shows that the finite-size winding obstruction is present under both natural projected-Hubbard conventions when the same displacement-resolved twist is used.

# 10. Phase-I and Phase-II research status

The broad UPC-only exact finite-size theorem should be replaced by a layered question: first make the twist response well posed, then identify the dynamical hypothesis that produces the many-body/two-body reduction.

| Research target | Status | Result |
|---|---|---|
| Exact four-cell arithmetic | Complete | Symbolic Kohn response gives the exact curvatures and defects. |
| Independent implementation | Complete | Fixed-number projected-density construction reproduces energies, gaps, factorization, and curvatures. |
| Nonresonant control | Complete | The same $r=4$ harmonic at $M=8$ obeys the law to about $10^{-9}$. |
| Restored resonance control | Complete | $M=8,r=8$ restores scalar factorization and gives defect $8U/7$ at $n=2$. |
| Twist-prescription audit | Complete | The motivating work uses displacement-resolved Peierls phases and $k\mapsto k+A$. |
| Interaction-convention check | Complete numerically | The positive-semidefinite convention retains the violation. |
| Arbitrary-$M$ obstruction theorem | Complete | Theorem 2 supplies an exact degenerate family for every $M\ge3$. |
| Alias/lift characterization | Complete | Finite samples fix alias sums, not winding moments; twist-stable UPC has an exact Fourier criterion. |
| Unique gapped family for all $M$ | Open | Proven at $M=4$; numerical evidence supports the phase-only reference family more broadly. |
| Sufficient repaired filling-law theorem | Open | Requires resolvent closure or another dynamical hypothesis beyond well-posedness. |

## 10.1 Phase I: certification package

The reproducibility package now contains:

1. an exact SymPy certificate for the four-cell Kohn response;
2. an independent NumPy/SciPy implementation constructed from projected local densities in fixed-number bit bases;
3. regression tests for the resonant and nonresonant cases;
4. CSV output for the cross-attraction and positive-semidefinite conventions;
5. the figure comparing reduced curvatures.

The remaining pre-submission task is external adversarial review of the revised scope and the arbitrary-$M$ proof, rather than further numerical validation of the four-cell arithmetic.

## 10.2 Phase II: next theorem targets

The obstruction and well-posedness parts are now established. The next analytic targets are:

1. prove the filling law for the nonwinding phase-only reference family
   
   $$
   u(k)=\frac{1}{\sqrt2}(1,e^{ib\sin k})^T
   $$
   
   at arbitrary $M$;
2. combine that proof with Proposition 1 to obtain a connected, unique-ground-state counterexample family for arbitrary $M$;
3. formulate the theorem on a fixed displacement-resolved infinite-lattice model, not on finite samples alone;
4. determine whether twist-stable UPC plus QGN representation/resolvent closure is sufficient for the filling law;
5. test fixed finite-range and uniformly exponentially local model sequences to separate finite-size winding artifacts from genuine many-body violations.

A repaired theorem should distinguish two levels:

- **Well-posedness:** the model includes a displacement-resolved Peierls lift or an explicitly specified family $P_A$, with a nonwinding/locality condition if the claim is intended to depend only on finite-torus data.
- **Dynamics:** the twist derivatives and reduced resolvent remain in representation channels whose matrix elements scale with the common factor $n(M-n)$.

The present work proves that the first level cannot be omitted. It does not yet prove that any proposed condition at the second level is sufficient.

# 11. Conclusion

The four-cell calculation is correct, but its scientific meaning is narrower and sharper than an unrestricted claim that flat-band filling laws fail. The exact finite-size result is:

> Zero-twist UPC, even together with exact flatness, time reversal, a projected onsite attraction, and a unique gapped paired branch, does not determine finite-torus twist curvature under a displacement-resolved $k\mapsto k+A$ prescription.

The failure occurs because a period-$M$ harmonic is invisible at zero twist yet becomes momentum independent under the twist, multiplying the entire projected interaction by a scalar. The $M=8,r=4$ control shows that a hidden high harmonic does not generically cause the violation. The arbitrary-$M$ theorem proves that the UPC-only hypothesis set is incomplete, while the finite-sampling identities identify the more fundamental ambiguity: finite-torus data determine alias sums but not winding moments.

For a fixed uniformly local model along a thermodynamic sequence, this exact resonance is absent or exponentially suppressed once the circumference exceeds the relevant displacement scale. The QGN conjecture for that setting therefore remains open. The revised research objective is to prove a restricted filling law after the twist lift, locality, and twist-stability assumptions are explicit, or to locate a second obstruction that survives those repairs and is genuinely dynamical.

# Acknowledgment

We thank an independent expert reviewer for reproducing the four-cell calculation in a separate fixed-number implementation and for proposing the eight-cell control that clarified the resonance mechanism.

# Appendix A. Reproducibility package

The computational package contains two code paths:

1. an exact SymPy implementation of the $M=4$ Kohn response, returning expressions in $\mathbb{Q}[\pi^2]$;
2. an independent numerical implementation that constructs projected local densities in a fixed-number bit basis and forms either the cross-attraction or positive-semidefinite Hamiltonian.

The numerical verifier reports the factorization residual, zero-twist energies and gaps, five-point energy curvatures, reduced curvature $E_n''/\rho_n$, and the analytic period-$M$ defect. Its default test suite covers three cases: $M=4,r=4$; $M=8,r=4$; and $M=8,r=8$.

# Appendix B. Short proof of the exponential-tail estimate

Choose residue representatives $r$ with $|r|\le M/2$. The difference between two lifts that agree in their alias sums enters the $j$th twist derivative through terms with $\ell\ne0$ and weights $(r+\ell M)^j$. Under $\|P_R\|\le Ce^{-\mu|R|}$,

$$
\sum_{\ell\ne0}|r+\ell M|^j\|P_{r+\ell M}\|
\le
C\sum_{\ell\ne0}|r+\ell M|^j e^{-\mu|r+\ell M|}.
$$

Since $|r+\ell M|\ge(|\ell|-1/2)M$, the right-hand side is bounded by

$$
C_j M^j e^{-\mu M/2}
$$

for a constant $C_j$ independent of $M$. This is a statement about the ambiguity in projector derivatives. Passing to a many-body curvature bound additionally requires norm estimates for the projected interaction and its resolvent.

# References

1. Q. Gao, Z. Han, and E. Khalaf, "Bootstrapping Flatband Superconductors: Rigorous Lower Bounds on Superfluid Stiffness," *Physical Review Letters* **136**, 076503 (2026), DOI: 10.1103/gw85-5r92; arXiv:2506.18969v4. In particular, see the supplemental derivation of the Peierls substitution and $k\mapsto k+A$, Eqs. (S83)-(S90).
2. *Finite-Size Flat-Band Stiffness Research Proposal*, revised 22 July 2026, unpublished manuscript.
3. M. Tovmasyan, S. Peotta, P. Torma, and S. D. Huber, "Effective Theory and Emergent SU(2) Symmetry in the Flat Bands of Attractive Hubbard Models," *Physical Review B* **94**, 245149 (2016).
4. W. Kohn, "Theory of the Insulating State," *Physical Review* **133**, A171-A181 (1964).
