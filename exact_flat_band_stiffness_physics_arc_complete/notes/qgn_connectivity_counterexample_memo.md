# Fixed-local QGN counterexample candidate from a reducible paired ground space

## Status

This memo records a second counterexample candidate to the broad statement

\[
D_s=\frac{N_{\rm flat}}2\,\nu(1-\nu)m_{\rm pair}^{-1}
\]

for all QGN models.  Unlike the earlier winding construction, the model below is fixed as the system size changes, has strictly finite real-space range, uses the `k -> k+A` flat-connection prescription, satisfies twist-stable UPC, and has a discrepancy that survives along an infinite thermodynamic sequence.

The caveat is equally sharp: the flat-band projector is **reducible**.  On every even torus it splits into two disconnected components, and the fixed-number ground space is degenerate across different pair distributions between them.  Thus the example refutes the broad hypothesis set only if connectivity/irreducibility and a unique AGP ground branch are not assumed.  It does not refute a corrected conjecture for connected QGN models with a unique paired ground state.

## 1. Fixed local band

Take a one-dimensional periodic lattice with two microscopic orbitals per cell and one rank-one flat band per spin.  Choose

\[
u_\uparrow(k)=
\begin{pmatrix}
\cos k\\
\sin k
\end{pmatrix},
\qquad
u_\downarrow(k)=u_\uparrow^*(-k)=
\begin{pmatrix}
\cos k\\
-\sin k
\end{pmatrix}.
\]

The spin-up projector is

\[
P_\uparrow(k)=\frac12
\begin{pmatrix}
1+\cos 2k & \sin 2k\\
\sin 2k & 1-\cos 2k
\end{pmatrix},
\]

and spin down is its time-reversal partner.  A parent Hamiltonian is

\[
h_\sigma(k)=I_2-P_\sigma(k).
\]

Its bands have exactly flat eigenvalues 0 and 1.  Because the projector contains only Fourier harmonics 0 and ±2, the parent is strictly finite range two, with no parameter depending on the circumference.

For every `M >= 4` relevant here,

\[
\frac1M\sum_{k_m}[P_\sigma(k_m+A)]_{11}
=
\frac1M\sum_{k_m}[P_\sigma(k_m+A)]_{22}
=\frac12
\]

for every continuous twist `A`.  The model therefore satisfies not merely zero-twist UPC but twist-stable UPC.  Spin conservation and time reversal put it in the perfect uniform-pairing particle-particle QGN class.

## 2. Interaction and twist

Use the projected positive-semidefinite Hubbard interaction

\[
H_M(A)=\frac U2\sum_{R,\alpha}
\left[\bar n_{R\alpha\uparrow}(A)-
      \bar n_{R\alpha\downarrow}(A)\right]^2,
\qquad U>0,
\]

with the flat connection inserted through `P_sigma(k+A)`.

In the flat-band real-space basis `d_{R sigma}`, the projected annihilation operators are

\[
\bar c_{R1\sigma}(A)
=\frac12\left(e^{iA}d_{R+1,\sigma}+e^{-iA}d_{R-1,\sigma}\right),
\]

and, up to a spin-dependent overall sign that drops out of the density,

\[
\bar c_{R2\sigma}(A)
=\frac1{2i}\left(e^{iA}d_{R+1,\sigma}-e^{-iA}d_{R-1,\sigma}\right).
\]

Every projected local density therefore preserves cell parity.

## 3. Exact thermodynamic obstruction

Let `M=2L`.  The band Hilbert space decomposes into two independent components,

\[
\mathcal H_{\rm band}=\mathcal H_{\rm even}\oplus\mathcal H_{\rm odd},
\]

each containing `L` one-particle states per spin.

At `n=L` pairs, hence `nu=n/M=1/2`, consider

\[
|\Phi_{\rm even}\rangle=
\prod_{R\ {\rm even}}
 d^\dagger_{R\uparrow}d^\dagger_{R\downarrow}|0\rangle.
\]

The even component is completely filled for both spins and the odd component is empty.  If `R` is odd, each projected density in the corresponding PSD factor acts on the full even component by the same scalar for up and down.  If `R` is even, it acts on the empty odd component and gives zero.  Hence

\[
\left[\bar n_{R\alpha\uparrow}(A)-
      \bar n_{R\alpha\downarrow}(A)\right]
|\Phi_{\rm even}\rangle=0
\]

for every `R`, `alpha`, and `A`.  Since the Hamiltonian is positive semidefinite,

\[
E_{0,M}(A)=0
\]

exactly for every even `M` and every sufficiently small twist.  Therefore

\[
D_s(M,\nu=1/2)=0.
\]

This is not a finite-size aliasing effect and not a choice-of-lift ambiguity: the same range-two Hamiltonian is used for all `M`.

## 4. Exact one-pair mass

In the one-pair sector with center-of-mass momentum `Q`, use relative-momentum states

\[
|k;Q\rangle=d^\dagger_{k+Q/2,\uparrow}
              d^\dagger_{-k+Q/2,\downarrow}|0\rangle.
\]

The two orbital form factors are

\[
g_1(k,Q)=\frac12(\cos Q+\cos 2k),
\qquad
g_2(k,Q)=\frac12(\cos Q-\cos 2k).
\]

For `M >= 6`, the finite momentum grid obeys

\[
\langle\cos 2k\rangle=0,
\qquad
\langle\cos^2 2k\rangle=\frac12.
\]

The nonzero eigenvalues of the attractive rank-two kernel are therefore the eigenvalues of

\[
G(Q)=\frac14
\begin{pmatrix}
\cos^2Q+\tfrac12 & \cos^2Q-\tfrac12\\
\cos^2Q-\tfrac12 & \cos^2Q+\tfrac12
\end{pmatrix},
\]

namely `cos^2(Q)/2` and `1/4`.  Near `Q=0` the lowest bound-state branch is

\[
E_{\rm pair}(Q)
=\frac U2-U\frac{\cos^2Q}{2}
=\frac U2\sin^2Q.
\]

Thus

\[
m_{\rm pair}^{-1}=E_{\rm pair}''(0)=U>0.
\]

For one flat band per spin, `N_flat=2`, and at `nu=1/2` the conjectured value is

\[
D_s^{\rm pred}
=\frac{2}{2}\frac14 U
=\frac U4.
\]

The exact ratio is consequently

\[
\boxed{R_M(1/2)=0}
\]

for every even `M >= 6`.  In particular, `R_M-1=-1` does not vanish in the thermodynamic limit.

## 5. Independent finite-size checks

The accompanying fixed-number exact-diagonalization code gives the following raw electronic-twist curvatures.  Ratios use the equivalent finite-size comparison `E_n''/[rho_n E_1'']`, with `rho_n=n(M-n)/(M-1)`.

| M | n | connectivity of step 2 mod M | E1-based prediction | exact curvature | ratio |
|---:|---:|:---:|---:|---:|---:|
| 5 | 2 | connected | 6 | 6 | 1 |
| 6 | 2 | two components | 6.4 | 4 | 5/8 |
| 6 | 3 | two components | 7.2 | 0 | 0 |
| 8 | 2 | two components | 48/7 | 16/3 | 7/9 |

For `M=6`, exact diagonalization finds a gap `U/4` above the sampled ground manifold in the `n=1,2,3` sectors.  The component-polarized half-filled state is annihilated channel by channel at nonzero twists to machine precision for both `M=6` and `M=8`.

The odd `M=5` result is an instructive control: because `gcd(2,5)=1`, displacement by two generates the whole finite ring and the filling law is recovered.  The contrast isolates reducibility rather than locality, UPC, QGN, or the twist convention.

## 6. Interpretation

The example reveals a second missing hypothesis, distinct from the winding issue:

> QGN, frustration freeness, twist-stable UPC, locality, and a positive pair mass do not force the exact ground-energy stiffness to equal the pair-mass upper bound when the paired ground space is reducible across disconnected components.

The ground manifold contains the uniform AGP superconducting state, but it also contains component-polarized full/empty states with zero stiffness.  Under a twist, the exact canonical ground energy follows the lowest-curvature branch.  Therefore a theorem must either select a particular AGP branch or prevent such competing ground states.

Natural repairs are:

1. require the projector/Wannier connectivity graph to be connected or primitive;
2. require irreducibility of the pair algebra;
3. require uniqueness of the fixed-number AGP ground state, up to unavoidable global symmetries;
4. formulate stiffness for a specified analytic ground-state branch rather than the lower envelope of a degenerate manifold.

## 7. What remains open

This is a stronger counterexample than the winding example in three respects: it is fixed local, twist-stable, and thermodynamically persistent.  It is weaker in one decisive respect: it uses a disconnected projector and a degenerate paired ground manifold.

A counterexample to the likely corrected conjecture would still need a connected/irreducible QGN model with a controlled paired branch and a nonzero limiting discrepancy.  Preliminary connected controls tested so far recover the filling law to numerical precision, so connectivity is currently the leading candidate for the missing assumption rather than a cosmetic repair.

## Reproducibility files

- `qgn_connectivity_counterexample.py`: independent fixed-number construction, spectral checks, curvature tests, structural diagnostics, and exact null-state checks.
- `qgn_connectivity_counterexample_results.csv`: archived numerical table.
