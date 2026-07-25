# Multi-Block Reduction and Irreducibility Laws for QGN Superconductors

* A review-integrated addendum to **Exact Flat-Band Stiffness from Pair Mobility in a QGN Model**
* Nicholas Sledgianowski
* Review-integrated draft, July 24, 2026

## Abstract

The companion paper proved an exact finite-size filling law for a Hermitian positive-square quantum-geometric-nesting (QGN) superconductor when the full-rank pairing form has one nonzero singular value. This addendum treats several inequivalent nonzero singular values and closes the corresponding irreducibility problem in the blockwise-tight rank-one projected-Hubbard subclass.

A Hermitian one-body null factor that commutes with the weighted pairing operator automatically preserves every singular-value eigenspace and therefore commutes with a product pseudospin group $\prod_a SU(2)_a$. At fixed total pair number the zero space consequently contains one product-AGP state for every allowed block filling. For a block-preserving twist and a fixed composition $\mathbf n=(n_1,\ldots,n_q)$, the exact positive reduction is

$$
C_{ij}^{(\mathbf n)}
=
\operatorname{Re}\left\langle\sum_a n_a\boldsymbol\beta_i^a,\sum_b n_b\boldsymbol\beta_j^b\right\rangle
+
\sum_a\frac{n_a(L_a-n_a)}{L_a-1}\Theta_{ij}^a,
$$

with the transverse term set to zero for $L_a=1$. If every block trace source vanishes, the response is a sum of independent hard-core laws, and in pair-mass normalization

$$
D_{s,ij}^{(\boldsymbol\nu)}
=
\sum_a\ell_a\nu_a(1-\nu_a)(m_{\mathrm{pair},a}^{-1})_{ij}.
$$

The microscopic result separates the semisimple and central parts of multi-block connectivity. Blockwise frame connectivity and absence of unitary synchronization are necessary and sufficient for uniqueness of the product AGP in every fixed balanced composition sector. They do not, by themselves, exclude block-imbalanced determinant zero modes. Those are governed by the normalized trace-profile matrix $T$: the fixed-even-total-number kernel is exactly the product-AGP span if and only if there is no nonzero signed capacity vector $\boldsymbol\delta$, with $\delta_a\in\{0,\pm L_a\}$ and even occupied capacity, satisfying $T\boldsymbol\delta=0$. The reason is that $d\Gamma$ is a Lie-algebra representation rather than an associative-algebra homomorphism, so associative fullness need not generate the relative block charges.

At unresolved total filling, a block-preserving twist gives a lower envelope over allowed zero-mode sectors. A generic block-mixing twist instead gives a positive semidefinite operator on the complete zero manifold, with universal intrablock and transfer Clebsch factors. Thus the multi-singular-value problem has an exact block-resolved branch law, an exact microscopic irreducibility criterion, and, when block charge is not conserved, a composition-space matrix law rather than one scalar filling polynomial.

# Introduction and headline statements

The parent paper established a finite-dimensional reduction theorem for Hermitian positive-square QGN Hamiltonians with one full-rank skew-unitary pairing block. In that setting a one-body source acting on the maximal pseudospin multiplet decomposes into a longitudinal $S$ channel and a transverse $S-1$ channel, and the fixed-pair-number Hessian obeys $$C_{ij}^{(n)}
 =
 \rho_n C_{ij}^{(1)}
 +(n^2-\rho_n)\Gamma_{ij},
 \qquad
 \rho_n=\frac{n(L-n)}{L-1}.
 \label{eq:parent-law}$$ When the trace tensor $\Gamma$ vanishes, the many-body response is exactly the one-pair response multiplied by the hard-core factor. The equal-singular-value assumption is the point at which the one-pair density matrix becomes proportional to one identity on the full paired space. The parent paper identified several inequivalent singular-value blocks as the natural next extension but left both the filling law and its irreducibility hypothesis open .

This addendum closes both parts in the Hermitian one-body class.

#### Product pseudospins and branch reduction.

Distinct nonzero singular values force a symmetry $$SU(2)_1\times\cdots\times SU(2)_q.
 \label{eq:product-group-intro}$$ The unresolved fixed-total-number zero space is therefore degenerate: it contains one product AGP for every allowed occupation vector $\mathbf n=(n_1,\ldots,n_q)$. For a twist that preserves the singular blocks, each fixed composition is an invariant branch. If all block trace sources vanish, its exact curvature is $$\boxed{
 C_{ij}^{(\mathbf n)}
 =
 \sum_{a=1}^q
 \frac{n_a(L_a-n_a)}{L_a-1}
 C_{ij}^{(\mathbf e_a)}.}
 \label{eq:headline-block-law}$$ The convention is that the $a$th term is zero when $L_a=1$. Here $\mathbf e_a$ denotes one pair in block $a$ and no pairs elsewhere.

#### Microscopic irreducibility.

For blockwise-tight rank-one projected-Hubbard factors, the one-block connectivity theorem has a precise multi-block replacement. Geometry and traces control different pieces of the constraint algebra: $$\boxed{
 \text{frame geometry controls the semisimple Lie algebra,
 while trace profiles control its center}.}
 \label{eq:geometry-center-slogan}$$ Blockwise connectivity and absence of synchronized equal-dimensional labeled frames are necessary and sufficient for uniqueness of the product AGP in every fixed balanced composition sector. A second, independent finite test is needed at fixed total particle number: no signed collection of completely filled blocks may have vanishing trace profile. This distinction corrects the tempting but false assertion that fullness of the generated associative $*$-algebra automatically supplies all relative block-charge constraints. The one-line reason is that $x\mapsto d\Gamma(x)$ is a Lie-algebra representation, not an associative-algebra homomorphism: block idempotents obtained only by multiplying the $x_\lambda$ need not appear as one-body charge constraints.

#### Unresolved filling and generic twists.

If the zero manifold is exactly the product-AGP span and the full twist preserves block charges, the physical fixed-total-number response is the lower envelope over compositions. If the source mixes singular blocks, degenerate least-squares perturbation theory produces a positive semidefinite matrix $$Q_v^{(n)}=T_v^\dagger T_v
 \label{eq:headline-Q}$$ on the finite zero manifold. Its diagonal entries contain the intrablock factors $$\rho_a(n_a)=\frac{n_a(L_a-n_a)}{L_a-1}
 \label{eq:rho-intro}$$ and the fermion-transfer factors $$\tau_{a\leftarrow b}(\mathbf n)
 =
 \frac{n_b(L_a-n_a)}{L_a},
 \label{eq:tau-intro}$$ while its only possible off-diagonal entries move one pair between two blocks with exact product-pseudospin Clebsch weights.

The weighted AGP generated by the original pairing form is only one superposition inside this zero manifold. Unequal singular values change that superposition but do not split its components. A block selector, a microscopic irreducibility certificate, or the complete degenerate response operator must therefore be specified before an unresolved stiffness can be assigned.

# Finite-dimensional setting

Let the paired one-particle Hilbert space be $$\mathfrak h=\bigoplus_{a=1}^q\mathfrak h_a,
 \qquad
 \dim\mathfrak h_a=2L_a,
 \qquad
 L=\sum_{a=1}^qL_a.
 \label{eq:h-decomp}$$ We exclude unpaired zero-singular-value modes. Let $F^T=-F$ be a full-rank pairing form with skew-Takagi decomposition $$F=U\left(\bigoplus_{a=1}^q \sigma_aJ_a\right)U^T,
 \qquad
 \sigma_a>0,
 \qquad
 \sigma_a\ne\sigma_b\quad(a\ne b),
 \label{eq:takagi}$$ where $$J_a^T=-J_a,
 \qquad
 J_a^\dagger J_a=I_{2L_a}.
 \label{eq:Ja}$$ After applying $U$ to the fermion modes, define the block pseudospins $$\eta_a^+=\frac12 c^\dagger J_ac^\dagger,
 \qquad
 \eta_a^-=(\eta_a^+)^\dagger,
 \qquad
 \eta_a^z=\frac{N_a-L_a}{2}.
 \label{eq:block-eta}$$ Different blocks commute, and each triple generates $\mathfrak{su}(2)$. This product-pseudospin organization is the natural singular-value refinement of the QGN and flat-band $\eta$-pairing constructions . The weighted pair operator is $$\eta_F^+=\frac12c^\dagger Fc^\dagger
 =\sum_{a=1}^q\sigma_a\eta_a^+.
 \label{eq:weighted-eta}$$

For a block occupation vector $$\mathbf n=(n_1,\ldots,n_q),
 \qquad
 0\le n_a\le L_a,
 \label{eq:composition}$$ define the normalized product AGP $$|\mathbf n\rangle
 =\bigotimes_{a=1}^q
 \left[\frac{(L_a-n_a)!}{n_a!L_a!}\right]^{1/2}
 (\eta_a^+)^{n_a}|0\rangle.
 \label{eq:product-agp}$$ It has $2n_a$ fermions in block $a$ and total pair number $|\mathbf n|=\sum_an_a$.

Consider an analytic Hermitian positive-square family $$H(A)=\frac12\sum_{\lambda=1}^{M_s}S_\lambda(A)^2
 =\frac12D(A)^\dagger D(A),
 \label{eq:H-family}$$ where $A=(A_1,\ldots,A_d)$ and $$D(A)|\psi\rangle
 =\bigl(S_1(A)|\psi\rangle,\ldots,S_{M_s}(A)|\psi\rangle\bigr).
 \label{eq:D-def}$$ At zero twist assume:

1.  **One-body weighted QGN symmetry.** Each $S_\lambda=S_\lambda(0)$ is an uncentered Hermitian one-body operator, $$S_\lambda=d\Gamma(s_\lambda),
     \qquad
     s_\lambda=s_\lambda^\dagger,
     \qquad
     [S_\lambda,\eta_F^+]=0,
     \qquad
     S_\lambda|0\rangle=0.
     \label{eq:MB1}$$

2.  **One-body twist source.** The first derivatives are Hermitian number-conserving one-body operators, $$B_{\lambda i}
     =\left.\partial_{A_i}S_\lambda(A)\right|_{A=0}
     =d\Gamma(b_{\lambda i}),
     \qquad
     b_{\lambda i}=b_{\lambda i}^\dagger.
     \label{eq:MB2}$$

3.  **Block-resolved irreducibility when invoked.** In each fixed block-occupation sector under consideration, $|\mathbf n\rangle$ is the unique zero mode of $H(0)$. Since the sector is finite dimensional, the next eigenvalue is then strictly positive. <span id="item:MB3" label="item:MB3"></span>

For the branch theorem we additionally assume that the full family preserves every block charge near $A=0$, $$=0
 \qquad
 \text{for every }\lambda,a.
 \label{eq:block-preserving-family}$$ For the generic-twist theorem we drop <a href="#eq:block-preserving-family" data-reference-type="eqref" data-reference="eq:block-preserving-family">[eq:block-preserving-family]</a> and instead require the complete zero space in the fixed total-pair-number sector. In the abstract QGN class this is an assumption. Section <a href="#sec:projected-irred" data-reference-type="ref" data-reference="sec:projected-irred">4</a> proves necessary-and-sufficient finite certificates for both MB3 and product-AGP completeness in the blockwise-tight rank-one projected-Hubbard subclass.

<div class="remark">

**Remark 1** (Scope of the splitting result). *The automatic splitting proved below uses the one-body form of the unperturbed null factors. For more general many-body null operators, product-pseudospin symmetry may be imposed directly as a hypothesis, but it does not follow solely from commutation with the weighted operator $\eta_F^+$. This is the same Hermitian one-body QGN class in which the singular-value decomposition has its sharpest algebraic meaning.*

</div>

# Distinct singular values force product pseudospins

<div id="prop:block-splitting" class="proposition">

**Proposition 2** (Automatic singular-block splitting). *Let $S=d\Gamma(s)$ with $s=s^\dagger$. If $$=0,
 \label{eq:weighted-commute}$$ then $s$ is block diagonal in the decomposition <a href="#eq:h-decomp" data-reference-type="eqref" data-reference="eq:h-decomp">[eq:h-decomp]</a>, and $$=[S,\eta_a^-]=[S,\eta_a^z]=0
 \qquad
 (a=1,\ldots,q).
 \label{eq:all-block-commute}$$*

</div>

<div class="proof">

*Proof.* The standard bilinear commutator identity converts <a href="#eq:weighted-commute" data-reference-type="eqref" data-reference="eq:weighted-commute">[eq:weighted-commute]</a> into $$sF+Fs^T=0.
 \label{eq:sF}$$ Taking the adjoint and using $s=s^\dagger$ gives $$F^\dagger s+s^TF^\dagger=0.
 \label{eq:adj-sF}$$ Multiplying <a href="#eq:sF" data-reference-type="eqref" data-reference="eq:sF">[eq:sF]</a> on the right by $F^\dagger$, multiplying <a href="#eq:adj-sF" data-reference-type="eqref" data-reference="eq:adj-sF">[eq:adj-sF]</a> on the left by $F$, and subtracting yields $$sFF^\dagger=FF^\dagger s.
 \label{eq:commute-FF}$$ In the skew-Takagi basis, $$FF^\dagger=\bigoplus_{a=1}^q\sigma_a^2I_{2L_a}.
 \label{eq:FF-spectrum}$$ The $\sigma_a^2$ are distinct, so <a href="#eq:commute-FF" data-reference-type="eqref" data-reference="eq:commute-FF">[eq:commute-FF]</a> forces $s=\bigoplus_as_a$. Restricting <a href="#eq:sF" data-reference-type="eqref" data-reference="eq:sF">[eq:sF]</a> to block $a$ gives $$s_aJ_a+J_as_a^T=0,
 \label{eq:symplectic-block}$$ which is equivalent to $[d\Gamma(s_a),\eta_a^+]=0$. Hermiticity gives the lowering commutator by adjunction, and block diagonality gives $[S,N_a]=0$, hence the $\eta_a^z$ commutator. ◻

</div>

<div id="cor:degeneracy" class="corollary">

**Corollary 3** (Automatic total-sector degeneracy). *Under MB1, every product AGP <a href="#eq:product-agp" data-reference-type="eqref" data-reference="eq:product-agp">[eq:product-agp]</a> is annihilated by every unperturbed square factor. At fixed total pair number $n$, $$\dim\ker H(0)\big|_{2n}
 \ge [z^n]\prod_{a=1}^q(1+z+\cdots+z^{L_a}).
 \label{eq:zero-count}$$ If $q>1$, the right-hand side exceeds one for every $$1\le n\le L-1.
 \label{eq:nonsaturated}$$ Thus a simple zero branch in the unresolved fixed-total-number sector is impossible throughout the nonsaturated range.*

</div>

<div class="proof">

*Proof.* Proposition <a href="#prop:block-splitting" data-reference-type="ref" data-reference="prop:block-splitting">2</a>, the vacuum condition, and repeated commutation give $$S_\lambda|\mathbf n\rangle=0
 \label{eq:all-product-zero}$$ for every allowed $\mathbf n$. States with different compositions are orthogonal because they have different eigenvalues of the commuting block charges $N_a$. Counting the allowed compositions gives <a href="#eq:zero-count" data-reference-type="eqref" data-reference="eq:zero-count">[eq:zero-count]</a>. With at least two positive capacities, every nonsaturated total occupation admits at least two allocations, so the coefficient is greater than one. ◻

</div>

The weighted AGP generated by $F$ is one vector in this zero manifold. Expanding <a href="#eq:weighted-eta" data-reference-type="eqref" data-reference="eq:weighted-eta">[eq:weighted-eta]</a> gives $$(\eta_F^+)^n|0\rangle
 \propto
 \sum_{|\mathbf n|=n}
 \left[
 \prod_{a=1}^q
 \sigma_a^{n_a}\sqrt{\binom{L_a}{n_a}}
 \right]
 |\mathbf n\rangle.
 \label{eq:weighted-expansion}$$ The singular values determine the amplitudes in <a href="#eq:weighted-expansion" data-reference-type="eqref" data-reference="eq:weighted-expansion">[eq:weighted-expansion]</a>; they do not lift the exact zero-energy degeneracy.

<div class="remark">

**Remark 4** (Equal singular values form one block). *If two subspaces carry the same singular value, the argument above does not force a null factor to preserve them separately. They must be combined into one maximal singular-value block. The decomposition in <a href="#eq:takagi" data-reference-type="eqref" data-reference="eq:takagi">[eq:takagi]</a> is therefore indexed by distinct nonzero singular values, not by an arbitrary refinement of the one-particle space.*

</div>

# Projected-Hubbard irreducibility: geometry, synchronization, and center

We now specialize the abstract product-pseudospin setting to the rank-one projected-Hubbard class that supplied the parent paper’s connectivity theorem. This section discharges MB3 and the fixed-total-number zero-manifold hypothesis by two logically distinct tests.

## Blockwise-tight rank-one factors and vectorization

Write the one-spin flat-band space as $$\mathcal V=\bigoplus_{a=1}^q\mathcal V_a,
 \qquad
 L_a=\dim\mathcal V_a.
 \label{eq:V-blocks}$$ The paired block in Section 2 is then $\mathfrak h_a=\mathcal V_a\oplus\overline{\mathcal V_a}$ with the canonical time-reversal pairing form. For each common factor label $\lambda$, let $$p_\lambda^a=|v_\lambda^a\rangle\langle v_\lambda^a|\succeq0
 \qquad\text{on }\mathcal V_a,
 \label{eq:p-lambda-a}$$ where zero vectors are allowed on labels inactive in a block. Assume the blockwise tight-frame identities $$\sum_\lambda p_\lambda^a=c_aI_{\mathcal V_a},
 \qquad
 c_a>0.
 \label{eq:blockwise-tight}$$ Put $$x_\lambda=\bigoplus_{a=1}^q p_\lambda^a
 \label{eq:x-lambda}$$ and define $$S_\lambda
 =
 d\Gamma_\uparrow(x_\lambda)
 -
 d\Gamma_\downarrow(\overline{x_\lambda}),
 \qquad
 H_0=\frac12\sum_\lambda S_\lambda^2.
 \label{eq:multiblock-Hubbard}$$ A strictly positive interaction kernel on the active constraint span gives the same null space after factorization.

The relevant one-body object is the compact real Lie algebra $$\mathfrak g
 =
 \operatorname{Lie}_{\mathbb R}\{ix_\lambda\}
 \subseteq
 \bigoplus_{a=1}^q\mathfrak u(\mathcal V_a).
 \label{eq:joint-lie}$$ Define the normalized trace-profile matrix $$T_{\lambda a}
 =
 \frac{\operatorname{Tr}p_\lambda^a}{L_a},
 \qquad
 Z=\operatorname{rowspan}_{\mathbb R}T\subseteq\mathbb R^q.
 \label{eq:trace-profile}$$

Let $\mathbf r=(r_1,\ldots,r_q)$ and $\mathbf s=(s_1,\ldots,s_q)$ be the spin-up and spin-down occupations of the blocks and set $$K_{\mathbf r}
 =
 \bigotimes_{a=1}^q\Lambda^{r_a}\mathcal V_a.
 \label{eq:K-r}$$ The charge sector is naturally $$K_{\mathbf r}\otimes\overline{K_{\mathbf s}}
 \simeq
 \operatorname{Hom}(K_{\mathbf s},K_{\mathbf r}).
 \label{eq:charge-vectorization}$$ If $$A_\lambda^{(\mathbf r)}
 =
 \sum_a d\Gamma_{r_a}(p_\lambda^a),
 \label{eq:A-lambda-r}$$ then $$S_\lambda\operatorname{vec}(X)
 =
 \operatorname{vec}\!\left(
 A_\lambda^{(\mathbf r)}X
 -
 XA_\lambda^{(\mathbf s)}
 \right).
 \label{eq:block-vectorization}$$ Positivity of <a href="#eq:multiblock-Hubbard" data-reference-type="eqref" data-reference="eq:multiblock-Hubbard">[eq:multiblock-Hubbard]</a> therefore gives the master identity $$\boxed{
 \ker H_0\big|_{(\mathbf r,\mathbf s)}
 \simeq
 \operatorname{Hom}_{\mathfrak g}(K_{\mathbf s},K_{\mathbf r}).}
 \label{eq:master-kernel}$$ The zero-mode problem sees the Lie algebra generated by the one-body constraints. It does not see arbitrary multiplicative combinations in the associative algebra generated by the $x_\lambda$.

## Associative fullness is not enough

The distinction in the previous paragraph is substantive.

<div id="ex:center-counterexample" class="example">

**Example 5** (Full associative algebra but a missing relative center). *Take $\mathcal V_A=\mathcal V_B=\mathbb C^2$ with basis $e_1,e_2$, and put $$f_\theta=\cos\theta\,e_1+\sin\theta\,e_2,
 \qquad
 g_\theta=-\sin\theta\,e_1+\cos\theta\,e_2.
 \label{eq:rotated-basis}$$ For $a=A,B$, define $$\begin{aligned}
 p_1^a&=\frac12|e_1\rangle\langle e_1|,
 &
 p_2^a&=\frac12|e_2\rangle\langle e_2|,
 \nonumber\\
 p_3^a&=\frac12|f_{\theta_a}\rangle\langle f_{\theta_a}|,
 &
 p_4^a&=\frac12|g_{\theta_a}\rangle\langle g_{\theta_a}|,
 \label{eq:center-counterexample-frames}
\end{aligned}$$ with $$\theta_A=\frac\pi4,
 \qquad
 \theta_B=\frac\pi6.
 \label{eq:center-counterexample-angles}$$ Each block is Parseval and connected. The labeled frames are not unitarily equivalent because $$\operatorname{Tr}(p_1^Ap_3^A)=\frac18,
 \qquad
 \operatorname{Tr}(p_1^Bp_3^B)=\frac{3}{16}.
 \label{eq:center-counterexample-mismatch}$$ Hence the generated complex $*$-algebra is $$\mathcal A=M_2(\mathbb C)\oplus M_2(\mathbb C).
 \label{eq:center-counterexample-assoc}$$ Nevertheless every $p_\lambda^a$ has trace $1/2$, so $$\operatorname{rank}T=1
 \label{eq:center-counterexample-rank}$$ and the real Lie algebra is only $$\mathfrak g
 =
 \mathfrak{su}(2)_A
 \oplus
 \mathfrak{su}(2)_B
 \oplus
 i\mathbb R(I_A\oplus I_B).
 \label{eq:center-counterexample-lie}$$ It lacks $i(I_A\oplus-I_B)$. The states with block $A$ filled in the up sector and block $B$ filled in the down sector, and the spin-reversed state, are annihilated by every square factor. At four particles the kernel has dimension five, whereas the balanced product-AGP compositions $(2,0)$, $(1,1)$, and $(0,2)$ account for only three states.*

</div>

The center of <a href="#eq:center-counterexample-assoc" data-reference-type="eqref" data-reference="eq:center-counterexample-assoc">[eq:center-counterexample-assoc]</a> contains both block idempotents, but $x\mapsto d\Gamma(x)$ is a Lie representation rather than an associative-algebra homomorphism: $$d\Gamma(xy)\ne d\Gamma(x)d\Gamma(y)
 \label{eq:dGamma-not-assoc}$$ in general. Relative block identities obtained as noncommutative polynomials in the $x_\lambda$ therefore need not be available as one-body charge constraints.

## Joint Lie algebra and complete charge-sector kernel

For each block, form the overlap graph on active labels, $$\lambda\sim_a\mu
 \quad\Longleftrightarrow\quad
 \langle v_\lambda^a,v_\mu^a\rangle\ne0.
 \label{eq:block-overlap-graph}$$ Assume in this subsection that every block frame spans and that every graph is connected. A disconnected block may first be refined into the spans of its connected components.

For equal-dimensional blocks define $a\sim b$ when there is a unitary $$W_{ba}:\mathcal V_a\longrightarrow\mathcal V_b
 \label{eq:sync-unitary}$$ such that $$W_{ba}p_\lambda^aW_{ba}^\dagger=p_\lambda^b
 \qquad
 \text{for every }\lambda.
 \label{eq:synchronization}$$ The equivalence classes are called synchronization classes and are denoted by $C$; their common dimension is $L_C$.

<div id="thm:joint-lie" class="theorem">

**Theorem 6** (Joint one-body Lie algebra). *Under blockwise connectivity, $$\boxed{
 \mathfrak g
 =
 \left(
 \bigoplus_C\Delta_C\mathfrak{su}(L_C)
 \right)
 \oplus iZ.}
 \label{eq:joint-lie-structure}$$ Here $\Delta_C$ is the diagonal embedding, up to the unitaries <a href="#eq:sync-unitary" data-reference-type="eqref" data-reference="eq:sync-unitary">[eq:sync-unitary]</a>, across the blocks in synchronization class $C$, and $$iZ
 =
 \left\{
 i\bigoplus_a z_aI_{\mathcal V_a}:
 (z_1,\ldots,z_q)\in Z
 \right\}.
 \label{eq:iZ}$$*

</div>

<div class="proof">

*Proof.* The parent paper’s spanning-tree induction applies to each coordinate family and gives $$\operatorname{Lie}_{\mathbb R}\{ip_\lambda^a\}
 =
 \mathfrak u(\mathcal V_a).
 \label{eq:coordinate-full}$$ Thus every coordinate projection of $\mathfrak g$ is surjective. The derived algebra $[\mathfrak g,\mathfrak g]$ is a compact semisimple subdirect product of the simple factors $\mathfrak{su}(\mathcal V_a)$, so its simple ideals are diagonally embedded across classes of isomorphic coordinate representations .

It remains to identify those classes. Write $$\widetilde p_\lambda^a
 =
 p_\lambda^a-T_{\lambda a}I_{\mathcal V_a}.
 \label{eq:centered-p}$$ A diagonal identification between two coordinate factors induces an automorphism of $\mathfrak{su}(L)$ carrying $i\widetilde p_\lambda^a$ to $i\widetilde p_\lambda^b$ for every label. For $L\ge3$, an outer automorphism negates the transpose of a traceless Hermitian matrix. The centered part of a nonzero positive rank-one matrix has one positive eigenvalue and $L-1$ equal negative eigenvalues; negation reverses these multiplicities and cannot be another centered positive rank-one matrix. For $L=2$, every automorphism is inner. The same spectral comparison forces $\operatorname{Tr}p_\lambda^a=\operatorname{Tr}p_\lambda^b$ whenever the centered projectors are identified, so the implementing unitary intertwines the full projectors rather than only their traceless parts. Hence every identification is inner and is exactly the unitary synchronization <a href="#eq:synchronization" data-reference-type="eqref" data-reference="eq:synchronization">[eq:synchronization]</a>. The active label sets also coincide: for $L\ge2$ the centered part of a nonzero rank-one positive matrix is nonzero, and an automorphism cannot map it to zero.

Finally, commutators have zero trace in every coordinate. The central part of $ix_\lambda$ is $$i\bigoplus_aT_{\lambda a}I_{\mathcal V_a}.
 \label{eq:central-generator}$$ The corresponding traceless tuple lies in the semisimple part, so the central tuple itself lies in $\mathfrak g$. Their span is $iZ$, and no further central directions can arise from Lie brackets. This proves <a href="#eq:joint-lie-structure" data-reference-type="eqref" data-reference="eq:joint-lie-structure">[eq:joint-lie-structure]</a>. ◻

</div>

For a synchronization class $C$, identify its blocks with $\mathbb C^{L_C}$ and define $$R_C(\mathbf r)
 =
 \bigotimes_{a\in C}\Lambda^{r_a}\mathbb C^{L_C}.
 \label{eq:R-C-r}$$ Decompose $$R_C(\mathbf r)
 \simeq
 \bigoplus_\mu
 V_{C,\mu}\otimes\mathbb C^{m_{C,\mu}(\mathbf r)}.
 \label{eq:RC-decomposition}$$

<div id="thm:complete-kernel" class="theorem">

**Theorem 7** (Complete charge-sector kernel formula). *Under blockwise connectivity, $$\boxed{
 \dim\ker H_0\big|_{(\mathbf r,\mathbf s)}
 =
 \mathbf 1_{\mathbf r-\mathbf s\in Z^\perp}
 \prod_C
 \left[
 \sum_\mu
 m_{C,\mu}(\mathbf r)m_{C,\mu}(\mathbf s)
 \right].}
 \label{eq:complete-kernel}$$ The multiplicities in <a href="#eq:RC-decomposition" data-reference-type="eqref" data-reference="eq:RC-decomposition">[eq:RC-decomposition]</a> are the usual Littlewood–Richardson multiplicities.*

</div>

<div class="proof">

*Proof.* By <a href="#eq:master-kernel" data-reference-type="eqref" data-reference="eq:master-kernel">[eq:master-kernel]</a>, the kernel is the intertwiner space between $K_{\mathbf s}$ and $K_{\mathbf r}$. The semisimple part in Theorem <a href="#thm:joint-lie" data-reference-type="ref" data-reference="thm:joint-lie">6</a> acts independently by the synchronized class representations <a href="#eq:R-C-r" data-reference-type="eqref" data-reference="eq:R-C-r">[eq:R-C-r]</a>. Schur decomposition therefore gives the product of multiplicity inner products in <a href="#eq:complete-kernel" data-reference-type="eqref" data-reference="eq:complete-kernel">[eq:complete-kernel]</a>. A center element associated with $\mathbf z\in Z$ acts on $K_{\mathbf r}$ with character $\sum_az_ar_a$. An intertwiner survives all central constraints exactly when $$\sum_az_a(r_a-s_a)=0
 \qquad
 \text{for every }\mathbf z\in Z,
 \label{eq:center-character-match}$$ which is $\mathbf r-\mathbf s\in Z^\perp$. ◻

</div>

## Corrected multi-block connectivity theorem

<div id="thm:MB3-criterion" class="theorem">

**Theorem 8** (Exact criterion for fixed-composition uniqueness (MB3)). *The product AGP is the unique zero mode in every fixed balanced composition sector if and only if:*

1.  *every block frame spans and its active-label overlap graph is connected; and*

2.  *no synchronization class with $L_C\ge2$ contains more than one block.*

*No trace-rank condition is needed for this fixed-composition statement. We call conditions (i)–(ii) the MB3 criterion.*

</div>

<div class="proof">

*Proof.* If every nontrivial synchronization class is a singleton, then $$K_{\mathbf n}=\bigotimes_a\Lambda^{n_a}\mathcal V_a
 \label{eq:K-n-irred}$$ is an external tensor product of irreducible representations of independent simple factors. Its commutant is scalar, and vectorization of that scalar is the product AGP. If a block graph is disconnected, the one-block geometric-series degeneracy appears by occupying that block alone. If two blocks of common dimension $L\ge2$ are synchronized, choose one pair in each. The diagonal algebra acts on $$\mathbb C^L\otimes\mathbb C^L
 =
 \operatorname{Sym}^2\mathbb C^L
 \oplus
 \Lambda^2\mathbb C^L,
 \label{eq:sync-two-fund}$$ so the commutant has dimension at least two. ◻

</div>

For two synchronized blocks of common dimension $L$, the Pieri rule gives a multiplicity-free decomposition of $\Lambda^r\mathbb C^L\otimes\Lambda^s\mathbb C^L$ . Consequently, $$\boxed{
 \dim\ker H_0\big|_{(r,s)\ {\rm balanced}}
 =
 \min(r,s,L-r,L-s)+1.}
 \label{eq:synchronized-count}$$ This is the exact interblock dark-state count associated with synchronization.

One-dimensional blocks carry no semisimple factor. Two $L=1$ blocks with identical labeled weight profiles are therefore trivially synchronized, but this is harmless for fixed-composition uniqueness; any relative-charge obstruction between them is detected entirely by the center test below.

The MB3 criterion of Theorem <a href="#thm:MB3-criterion" data-reference-type="ref" data-reference="thm:MB3-criterion">8</a> does not exclude block-imbalanced determinant modes. Assume that criterion holds and define $$\mathcal D_L
 =
 \left\{
 \boldsymbol\delta\in\mathbb Z^q:
 \delta_a\in\{0,+L_a,-L_a\}
 \right\},
 \label{eq:det-difference-set}$$ and $$F(\boldsymbol\delta)
 =
 \sum_{\delta_a\ne0}L_a.
 \label{eq:F-delta}$$

<div id="thm:det-resonance" class="theorem">

**Theorem 9** (Determinant-resonance criterion). *Assume the MB3 criterion of Theorem <a href="#thm:MB3-criterion" data-reference-type="ref" data-reference="thm:MB3-criterion">8</a>. A block-imbalanced zero mode exists in charge sector $(\mathbf r,\mathbf s)$ if and only if $$\boldsymbol\delta=\mathbf r-\mathbf s
 \in
 \mathcal D_L\cap Z^\perp.
 \label{eq:det-resonance-sector}$$ Every such allowed charge sector contributes a one-dimensional zero space. Hence the kernel in every even total-particle sector is exactly the span of the product AGPs if and only if $$\boxed{
 \left\{
 \boldsymbol\delta\in\mathcal D_L\cap Z^\perp:
 F(\boldsymbol\delta)\ {\rm even}
 \right\}
 =
 \{0\}.}
 \label{eq:det-resonance-even}$$*

</div>

<div class="proof">

*Proof.* With singleton nontrivial synchronization classes, the semisimple intertwiner problem factorizes by block. For $SU(L_a)$, the exterior powers $\Lambda^{r_a}\mathcal V_a$ and $\Lambda^{s_a}\mathcal V_a$ are isomorphic exactly when $r_a=s_a$ or $\{r_a,s_a\}=\{0,L_a\}$. In each allowed case the intertwiner is one dimensional. Their tensor product survives the Lie center precisely when $T(\mathbf r-\mathbf s)=0$, giving <a href="#eq:det-resonance-sector" data-reference-type="eqref" data-reference="eq:det-resonance-sector">[eq:det-resonance-sector]</a>.

A nonzero $\boldsymbol\delta$ contributes to an even total sector exactly when its forced filled-block particle number $F(\boldsymbol\delta)$ is even. Balanced occupations add particles in pairs and do not change this parity, proving <a href="#eq:det-resonance-even" data-reference-type="eqref" data-reference="eq:det-resonance-even">[eq:det-resonance-even]</a>. ◻

</div>

#### Odd-capacity resonances.

The parity qualifier in <a href="#eq:det-resonance-even" data-reference-type="eqref" data-reference="eq:det-resonance-even">[eq:det-resonance-even]</a> has independent physical content. If $F(\boldsymbol\delta)$ is odd, the resonance cannot enter any even pairing sector, but it produces exact zero-energy states in odd total-particle sectors. These are fermionic determinant dark states rather than competing paired AGPs. The smallest example has capacities $(L_1,L_2)=(1,2)$, center $Z=\operatorname{span}\{(2,1)\}$, and $\boldsymbol\delta=(1,-2)$: every even sector is product-AGP complete, while the three-particle sector contains the two spin-reversed determinant zero modes.

Still under the MB3 criterion of Theorem <a href="#thm:MB3-criterion" data-reference-type="ref" data-reference="thm:MB3-criterion">8</a>, the exact total-$2n$ nullity is $$\boxed{
 \begin{aligned}
 \dim\ker H_0\big|_{2n}
 ={}&
 \sum_{\substack{
 \boldsymbol\delta\in\mathcal D_L\cap Z^\perp\\
 F(\boldsymbol\delta)\le2n\\
 F(\boldsymbol\delta)\equiv2n\ ({\rm mod}\ 2)
 }}
 [z^{\,n-F(\boldsymbol\delta)/2}]
 \prod_{\delta_a=0}
 (1+z+\cdots+z^{L_a}).
 \end{aligned}}
 \label{eq:total-nullity-count}$$ The $\boldsymbol\delta=0$ term is the product-AGP composition count. Every nonzero term is a determinant-transfer dark sector.

<div id="cor:full-trace-rank" class="corollary">

**Corollary 10** (Simple sufficient center certificate). *If $$\operatorname{rank}T=q,
 \label{eq:full-trace-rank}$$ then $Z^\perp=\{0\}$ and no determinant resonance is possible. Together with connectivity and nonsynchronization, $$\mathfrak g=\bigoplus_a\mathfrak u(\mathcal V_a)
 \label{eq:full-joint-u}$$ and both fixed-composition uniqueness and product-AGP completeness hold.*

</div>

Full trace rank is sufficient but not necessary. For example, two nonsynchronized connected blocks with capacities $L_1=2$, $L_2=3$ and $Z=\operatorname{span}\{(1,1)\}$ have a one-dimensional center deficit, but no nonzero vector in $$\{0,\pm2\}\times\{0,\pm3\}
 \label{eq:unequal-capacity-set}$$ is orthogonal to $(1,1)$. Their product-AGP zero manifold is therefore exact in every sector.

The three failure mechanisms are now distinct:

<div class="center">

| Mechanism             | Missing algebraic structure | Zero-mode signature                |
|:----------------------|:----------------------------|:-----------------------------------|
| Frame disconnection   | reducible coordinate action | geometric-series component count   |
| Frame synchronization | independent simple factor   | Littlewood–Richardson multiplicity |
| Determinant resonance | relative block center       | signed full-determinant sector     |

</div>

## Finite certificate and explicit multi-singular-value models

The microscopic audit is finite:

1.  check spanning and connectivity of every block frame;

2.  for each equal-dimensional block pair, solve $$Wp_\lambda^a=p_\lambda^bW
     \qquad\text{for all }\lambda;
     \label{eq:sylvester-sync}$$ a nonzero solution is, after rescaling, unitary;

3.  form $T$ and enumerate at most $3^q$ signed capacity vectors in $\mathcal D_L$.

A mismatch in either $\operatorname{Tr}p_\lambda^a$ or $$\operatorname{Tr}(p_\lambda^ap_\mu^a)
 \label{eq:gram-sync-certificate}$$ already excludes synchronization.

<div id="ex:repaired-model" class="example">

**Example 11** (Minimal repaired two-block model). *Keep the angles <a href="#eq:center-counterexample-angles" data-reference-type="eqref" data-reference="eq:center-counterexample-angles">[eq:center-counterexample-angles]</a>, but use $$\begin{aligned}
 p_1^A&=\frac13|e_1\rangle\langle e_1|,
 &
 p_2^A&=\frac13|e_2\rangle\langle e_2|,
 &
 p_3^A&=\frac23|f_{\pi/4}\rangle\langle f_{\pi/4}|,
 &
 p_4^A&=\frac23|g_{\pi/4}\rangle\langle g_{\pi/4}|,
 \nonumber\\
 p_1^B&=\frac23|e_1\rangle\langle e_1|,
 &
 p_2^B&=\frac23|e_2\rangle\langle e_2|,
 &
 p_3^B&=\frac13|f_{\pi/6}\rangle\langle f_{\pi/6}|,
 &
 p_4^B&=\frac13|g_{\pi/6}\rangle\langle g_{\pi/6}|.
 \label{eq:repaired-frames}
\end{aligned}$$ Both blocks remain Parseval, connected, and nonsynchronized, while $$T=
 \begin{pmatrix}
 1/6&1/3\\
 1/6&1/3\\
 1/3&1/6\\
 1/3&1/6
 \end{pmatrix},
 \qquad
 \operatorname{rank}T=2.
 \label{eq:repaired-T}$$ Choose distinct pairing singular values, for example $$\sigma_A=1,
 \qquad
 \sigma_B=2.
 \label{eq:repaired-singular-values}$$ The resulting eight-mode positive-square QGN model has $$\mathfrak g=\mathfrak u(2)_A\oplus\mathfrak u(2)_B
 \label{eq:repaired-lie}$$ and exact even-sector nullities $$(1,2,3,2,1),
 \label{eq:repaired-nullities}$$ the coefficients of $(1+z+z^2)^2$.*

</div>

<div id="prop:lattice-ready" class="proposition">

**Proposition 12** (Lattice-ready multi-channel construction). *Let $\{p_x^a\}_x$ be connected one-block projected-Hubbard frames, padded by zero matrices to a common label set. Choose a real invertible $q\times q$ matrix $W$ with positive entries and introduce factor channels $$x_{x\mu}
 =
 \bigoplus_{a=1}^qW_{\mu a}p_x^a,
 \qquad
 \mu=1,\ldots,q.
 \label{eq:mixed-channel-factors}$$ Then $$\operatorname{Lie}_{\mathbb R}\{ix_{x\mu}\}
 =
 \bigoplus_a\mathfrak u(\mathcal V_a).
 \label{eq:mixed-channel-full}$$ For arbitrary distinct positive $\sigma_a$, the pairing form $$F=\bigoplus_a\sigma_aJ_a
 \label{eq:mixed-channel-F}$$ therefore gives a local Hermitian positive-square multi-singular-value QGN model satisfying both uniqueness hypotheses.*

</div>

<div class="proof">

*Proof.* For fixed $x$, invertibility of $W$ gives $$\sum_\mu(W^{-1})_{a\mu}x_{x\mu}
 =
 0\oplus\cdots\oplus p_x^a\oplus\cdots\oplus0.
 \label{eq:invert-W}$$ Thus every block generator lies independently in the joint Lie algebra. Connectivity supplies $\mathfrak u(\mathcal V_a)$ in each block. Every factor is block diagonal and hence commutes with every block pseudospin; expanding the squares may still produce interblock density couplings, so the Hamiltonian need not be a direct sum. ◻

</div>

The representation-theoretic open problem is therefore closed within this projected rank-one class. What remains model dependent is not the criterion but the physical realization: which singular blocks occur in a given band structure, how the experimental flat connection acts on them, and what block pair masses and trace sources it produces.

# The multi-block one-body source decomposition

Let $P_a$ be the one-particle projector onto $\mathfrak h_a$. For a number-conserving one-body operator $$B=d\Gamma(b),
 \label{eq:B-generic}$$ write $$b=\sum_{a,b}b^{ab},
 \qquad
 b^{ab}=P_abP_b,
 \qquad
 B^{ab}=d\Gamma(b^{ab}).
 \label{eq:block-B}$$ Let $\mathbf e_a$ denote the composition with one pair in block $a$ and zero elsewhere. For the diagonal blocks define $$\beta_a(B)
 =\langle\mathbf e_a|B|\mathbf e_a\rangle
 =\frac1{L_a}\operatorname{Tr}_{\mathfrak h_a}b^{aa},
 \label{eq:beta-def}$$ and the transverse one-pair seed $$|\zeta_a(B)\rangle
 =B^{aa}|\mathbf e_a\rangle-\beta_a(B)|\mathbf e_a\rangle.
 \label{eq:zeta-diag}$$ For $a\ne b$, define the transfer seed $$|\zeta_{a\leftarrow b}(B)\rangle
 =B^{ab}|\mathbf e_b\rangle.
 \label{eq:zeta-transfer}$$ We use the filling factors $$\rho_a(n_a)=
 \begin{cases}
 \dfrac{n_a(L_a-n_a)}{L_a-1},&L_a\ge2,\\[0.8em]
 0,&L_a=1,
 \end{cases}
 \label{eq:rho-block}$$ and, for ordered distinct blocks, $$\tau_{a\leftarrow b}(\mathbf n)
 =\frac{n_b(L_a-n_a)}{L_a}.
 \label{eq:tau-block}$$

<div id="lem:multiblock-action" class="lemma">

**Lemma 13** (Multi-block action of a one-body operator). *For every allowed composition $\mathbf n$, there are normalized product-pseudospin ladder isometries $U_{a,\mathbf n}$ and $U_{a\leftarrow b,\mathbf n}$, independent of the particular matrix $b$, such that $$\begin{aligned}
 B|\mathbf n\rangle
 ={}&\left(\sum_{a=1}^q n_a\beta_a(B)\right)|\mathbf n\rangle
 \nonumber\\
 &+\sum_{a=1}^q\sqrt{\rho_a(n_a)}\,
 U_{a,\mathbf n}|\zeta_a(B)\rangle
 \nonumber\\
 &+\sum_{a\ne b}\sqrt{\tau_{a\leftarrow b}(\mathbf n)}\,
 U_{a\leftarrow b,\mathbf n}|\zeta_{a\leftarrow b}(B)\rangle.
 \label{eq:multiblock-action}
\end{aligned}$$ The displayed channels are mutually orthogonal. The first lies in the product of maximal-spin representations. The second lowers block $a$ from spin $S_a=L_a/2$ to $S_a-1$. The third lowers blocks $a$ and $b$ to $S_a-\tfrac12$ and $S_b-\tfrac12$, respectively.*

</div>

<div class="proof">

*Proof.* For $B^{aa}$, the one-block source lemma applies inside $\mathfrak h_a$ and is tensored with the maximal-spin AGPs in the other blocks. Its longitudinal coefficient is $n_a\beta_a(B)$ and its transverse norm ratio is $\rho_a(n_a)$. Summing over $a$ gives the first two lines of <a href="#eq:multiblock-action" data-reference-type="eqref" data-reference="eq:multiblock-action">[eq:multiblock-action]</a>.

It remains to identify the transfer coefficient. Let $x:\mathfrak h_b\to\mathfrak h_a$ and $$B_x=c_a^\dagger x c_b.
 \label{eq:Bx}$$ The one-body density matrix of a block AGP is $$\langle\mathbf n|c_{a,i}^\dagger c_{a,j}|\mathbf n\rangle
 =\frac{n_a}{L_a}\delta_{ij}.
 \label{eq:block-1rdm}$$ Because the state is a product over blocks, for $x,y:\mathfrak h_b\to\mathfrak h_a$ one obtains $$\begin{aligned}
 \langle B_y\mathbf n|B_x\mathbf n\rangle
 &={\frac{n_b}{L_b}}
 \left(1-\frac{n_a}{L_a}\right)
 \operatorname{Tr}(y^\dagger x),
 \label{eq:transfer-inner-n}\\
 \langle B_y\mathbf e_b|B_x\mathbf e_b\rangle
 &={\frac1{L_b}}\operatorname{Tr}(y^\dagger x).
 \label{eq:transfer-inner-one}
\end{aligned}$$ The ratio of <a href="#eq:transfer-inner-n" data-reference-type="eqref" data-reference="eq:transfer-inner-n">[eq:transfer-inner-n]</a> and <a href="#eq:transfer-inner-one" data-reference-type="eqref" data-reference="eq:transfer-inner-one">[eq:transfer-inner-one]</a> is exactly <a href="#eq:tau-block" data-reference-type="eqref" data-reference="eq:tau-block">[eq:tau-block]</a>. Polarization shows that normalized pseudospin raising defines one isometry on the complete transfer-seed multiplicity space, independent of $x$. This is $U_{a\leftarrow b,\mathbf n}$.

The three channel types have different product-pseudospin labels or different block-charge weights. Distinct diagonal blocks lower different factors of the product representation, while distinct ordered transfer pairs produce different odd block-charge patterns. They are therefore mutually orthogonal. ◻

</div>

<div class="remark">

**Remark 14** (Meaning of the transfer factor). *The factor $n_b/L_b$ in <a href="#eq:transfer-inner-n" data-reference-type="eqref" data-reference="eq:transfer-inner-n">[eq:transfer-inner-n]</a> is the occupation probability in the source block, while $1-n_a/L_a$ is the hole probability in the target block. Relative to the one-pair seed normalization, their product becomes $\tau_{a\leftarrow b}=n_b(1-n_a/L_a)$. The new channel is therefore the exact hard-core phase-space factor for transferring one fermion between two paired pseudospins.*

</div>

# Block-resolved finite-size curvature

We first recall the least-squares form of frustration-free curvature. Let a finite-dimensional invariant sector have a unique normalized zero mode $|\psi\rangle$ at $A=0$. Put $$b_i=\left.\partial_{A_i}D(A)\right|_0|\psi\rangle,
 \qquad
 \Pi=P_{\ker D(0)^\dagger}.
 \label{eq:LS-data}$$ Then the exact branch Hessian is $$C_{ij}=\operatorname{Re}\langle\Pi b_i,\Pi b_j\rangle.
 \label{eq:LS}$$ Equivalently, in a real direction $v$, $$v_iC_{ij}v_j
 =\min_{\chi\perp\psi}
 \left\|D(0)\chi+v_i b_i\right\|^2.
 \label{eq:LS-min}$$ This is the least-squares lemma used in the parent paper; it follows by differentiating $H=D^\dagger D/2$ and solving the first-order normal equation.

Assume now the block-preserving condition <a href="#eq:block-preserving-family" data-reference-type="eqref" data-reference="eq:block-preserving-family">[eq:block-preserving-family]</a> and MB3 in a fixed composition sector $\mathbf n$. For each factor and twist direction define $$\beta_{\lambda i}^a
 =\frac1{L_a}\operatorname{Tr}_{\mathfrak h_a}b_{\lambda i}^{aa},
 \qquad
 \boldsymbol\beta_i^a=(\beta_{1i}^a,\ldots,\beta_{M_s i}^a)\in\mathbb C^{M_s}.
 \label{eq:beta-vector}$$ Let $$\boldsymbol\zeta_i^a
 =\bigl(\zeta_a(B_{1i}),\ldots,\zeta_a(B_{M_si})\bigr)
 \label{eq:zeta-vector}$$ be the target-space one-pair seed, and let $\Pi_a$ be the target-kernel projector in the sector $\mathbf e_a$. Define $$\Theta_{ij}^a
 =\operatorname{Re}\langle\Pi_a\boldsymbol\zeta_i^a,\Pi_a\boldsymbol\zeta_j^a\rangle,
 \label{eq:Theta}$$ so each $\Theta^a$ is positive semidefinite. Finally set $$\Gamma_{ij}^{ab}
 =\operatorname{Re}\langle\boldsymbol\beta_i^a,\boldsymbol\beta_j^b\rangle.
 \label{eq:Gamma-ab}$$

<div id="thm:block-resolved" class="theorem">

**Theorem 15** (Exact block-resolved reduction). *Under MB1–MB3 and the block-preserving twist condition <a href="#eq:block-preserving-family" data-reference-type="eqref" data-reference="eq:block-preserving-family">[eq:block-preserving-family]</a>, the exact Hessian of the simple branch in sector $\mathbf n$ is $$\boxed{
 C_{ij}^{(\mathbf n)}
 =\operatorname{Re}\left\langle
 \sum_an_a\boldsymbol\beta_i^a,
 \sum_bn_b\boldsymbol\beta_j^b
 \right\rangle
 +\sum_a\rho_a(n_a)\Theta_{ij}^a.}
 \label{eq:MB-positive}$$ Equivalently, because $$C_{ij}^{(\mathbf e_a)}=\Gamma_{ij}^{aa}+\Theta_{ij}^a,
 \label{eq:onepair-block-C}$$ one has $$\boxed{
 C_{ij}^{(\mathbf n)}
 =\sum_a\rho_a(n_a)C_{ij}^{(\mathbf e_a)}
 +\sum_{a,b}n_an_b\Gamma_{ij}^{ab}
 -\sum_a\rho_a(n_a)\Gamma_{ij}^{aa}.}
 \label{eq:MB-rearranged}$$*

</div>

<div class="proof">

*Proof.* Let $D_{\mathbf n}$ be the unperturbed null map restricted to the composition sector $\mathbf n$. By Proposition <a href="#prop:block-splitting" data-reference-type="ref" data-reference="prop:block-splitting">2</a>, the full null map intertwines the product group $\prod_aSU(2)_a$. Its target-kernel projector therefore commutes with the normalized ladder maps in Lemma <a href="#lem:multiblock-action" data-reference-type="ref" data-reference="lem:multiblock-action">13</a> and preserves the mutually orthogonal product-spin channels.

Apply Lemma <a href="#lem:multiblock-action" data-reference-type="ref" data-reference="lem:multiblock-action">13</a> to every $B_{\lambda i}$. Because the source is block diagonal, the transfer line is absent. The target source is the orthogonal sum of $$|\mathbf n\rangle\otimes\sum_an_a\boldsymbol\beta_i^a
 \label{eq:longitudinal-target}$$ and the block-transverse terms $\sqrt{\rho_a(n_a)}U_{a,\mathbf n}^{\oplus}\boldsymbol\zeta_i^a$. The longitudinal vector cannot be removed by a least-squares correction: $$D_{\mathbf n}^\dagger\bigl(|\mathbf n\rangle\otimes\mathbf v\bigr)
 =\sum_\lambda v_\lambda S_\lambda|\mathbf n\rangle=0.
 \label{eq:longitudinal-kernel}$$ On the block-$a$ transverse channel, intertwining gives $$\Pi_{\mathbf n}U_{a,\mathbf n}^{\oplus}
 =U_{a,\mathbf n}^{\oplus}\Pi_a.
 \label{eq:block-intertwine-Pi}$$ The least-squares formula <a href="#eq:LS" data-reference-type="eqref" data-reference="eq:LS">[eq:LS]</a> therefore reduces to the sum of the longitudinal Gram form and the $q$ orthogonal transverse Gram forms, proving <a href="#eq:MB-positive" data-reference-type="eqref" data-reference="eq:MB-positive">[eq:MB-positive]</a>. Substituting <a href="#eq:onepair-block-C" data-reference-type="eqref" data-reference="eq:onepair-block-C">[eq:onepair-block-C]</a> gives <a href="#eq:MB-rearranged" data-reference-type="eqref" data-reference="eq:MB-rearranged">[eq:MB-rearranged]</a>. ◻

</div>

<div class="remark">

**Remark 16** (The manifestly positive form is canonical). *Equation <a href="#eq:MB-positive" data-reference-type="eqref" data-reference="eq:MB-positive">[eq:MB-positive]</a> is positive in every real twist direction. In contrast with the one-block law, the difference $$C^{(\mathbf n)}-\sum_a\rho_a(n_a)C^{(\mathbf e_a)}
 \label{eq:not-positive-difference}$$ need not be positive semidefinite. The block trace vectors can interfere destructively in the first term of <a href="#eq:MB-positive" data-reference-type="eqref" data-reference="eq:MB-positive">[eq:MB-positive]</a>. The rearranged form <a href="#eq:MB-rearranged" data-reference-type="eqref" data-reference="eq:MB-rearranged">[eq:MB-rearranged]</a> is useful for comparison with one-pair data, but it should not be read as a positive correction theorem.*

</div>

<div id="cor:trace-free" class="corollary">

**Corollary 17** (Blockwise trace-free filling law). *If $$\beta_{\lambda i}^a=0
 \qquad
 \text{for every }a,\lambda,i,
 \label{eq:block-trace-free}$$ then $$\boxed{
 C_{ij}^{(\mathbf n)}
 =\sum_a
 \frac{n_a(L_a-n_a)}{L_a-1}
 C_{ij}^{(\mathbf e_a)},}
 \label{eq:trace-free-law}$$ with the $L_a=1$ terms understood as zero.*

</div>

<div class="remark">

**Remark 18** (Blockwise, not merely total, cancellation). *The condition $\sum_a\boldsymbol\beta_i^a=0$ only cancels the longitudinal source for the special composition with all $n_a$ equal to one. A uniform block law requires the stronger statement <a href="#eq:block-trace-free" data-reference-type="eqref" data-reference="eq:block-trace-free">[eq:block-trace-free]</a>, or at least the composition-dependent cancellation $\sum_an_a\boldsymbol\beta_i^a=0$ on the particular branch of interest.*

</div>

## Pair-mass normalization

Let the system live on a torus of volume $V$. As in the parent paper, define the canonical flat-connection curvature $$\kappa_{ij}^{(\mathbf n)}(V)=\frac1{4V}C_{ij}^{(\mathbf n)}.
 \label{eq:kappa-def}$$ Suppose the invariant one-pair branch in block $a$ has lifted center-of-mass dispersion $$E_{\mathrm{pair},a}(Q)
 =E_{\mathrm{pair},a}(0)
 +\frac12Q_i(m_{\mathrm{pair},a}^{-1})_{ij}Q_j
 +o(|Q|^2),
 \label{eq:block-pair-disp}$$ and that the electronic connection obeys $Q=2A$. Then $$C_{ij}^{(\mathbf e_a)}=4(m_{\mathrm{pair},a}^{-1})_{ij}.
 \label{eq:block-C-mass}$$ The mass in <a href="#eq:block-pair-disp" data-reference-type="eqref" data-reference="eq:block-pair-disp">[eq:block-pair-disp]</a> is block resolved: it is the curvature of the branch in the invariant $\mathbf e_a$ sector, not necessarily the smallest eigenvalue of the unresolved one-pair problem.

<div id="cor:mass-law" class="corollary">

**Corollary 19** (Finite-size block pair-mass law). *Under the blockwise trace-free condition, $$\boxed{
 \kappa_{ij}^{(\mathbf n)}(V)
 =\frac1V\sum_a
 \frac{n_a(L_a-n_a)}{L_a-1}
 (m_{\mathrm{pair},a}^{-1})_{ij}.}
 \label{eq:finite-mass-law}$$ Writing $\nu_a=n_a/L_a$ gives $$\kappa_{ij}^{(\mathbf n)}(V)
 =\sum_a\frac{L_a}{V}\frac{L_a}{L_a-1}
 \nu_a(1-\nu_a)(m_{\mathrm{pair},a}^{-1})_{ij}.
 \label{eq:finite-mass-law-nu}$$*

</div>

# Blockwise locality and the thermodynamic law

The longitudinal term in the parent paper is a winding trace. The same roots-of-unity argument applies separately to every singular block, provided the projected trace symbol is defined blockwise. For a translation-covariant local channel $\alpha$, let the normalized block trace be represented by a global periodic scalar $f_\alpha^a(k)$ and write $$f_\alpha^a(k)=\sum_{R\in\mathbb Z^d}\widehat f_\alpha^a(R)e^{ik\cdot R}.
 \label{eq:block-symbol-fourier}$$ On a rectangular torus with side lengths $\mathbf N=(N_1,\ldots,N_d)$, the finite-grid trace derivative is $$\beta_{\alpha i}^{a,(\mathbf N)}
 =\sum_{\ell\ne0}
 i\ell_iN_i\widehat f_\alpha^a(\ell\odot\mathbf N).
 \label{eq:block-winding}$$ Thus only block-resolved Fourier coefficients that wind around the torus contribute.

If every $\widehat f_\alpha^a$ has fixed finite range, then all $\boldsymbol\beta_i^a$ vanish exactly once the torus exceeds that range. If $$|\widehat f_\alpha^a(R)|\le C e^{-\mu|R|_1}
 \label{eq:block-exp-local}$$ uniformly in size and block, then every entry of the complete block Gram matrix $\Gamma^{ab}$ is exponentially suppressed up to polynomial side-length factors. Locality of only the unresolved sum $\sum_af_\alpha^a$ is insufficient: different block traces can cancel before the individual winding amplitudes are controlled.

<div id="thm:thermo-block" class="theorem">

**Theorem 20** (Thermodynamic block-resolved reduction). *Consider a regular sequence of tori with volume $V\to\infty$, a fixed finite number of singular blocks, and $$\frac{L_a}{V}\longrightarrow\ell_a,
 \qquad
 \frac{n_a}{L_a}\longrightarrow\nu_a.
 \label{eq:block-density-limits}$$ Assume block-resolved irreducibility at every size, a block-preserving twist, convergent block pair masses, and either fixed finite-range or uniformly exponentially local block trace symbols as above. Then $$\boxed{
 D_{s,ij}^{(\boldsymbol\nu)}
 =\sum_{a=1}^q
 \ell_a\nu_a(1-\nu_a)
 (m_{\mathrm{pair},a}^{-1})_{ij}.}
 \label{eq:thermo-block-law}$$ Blocks with $\ell_a=0$ do not contribute.*

</div>

<div class="proof">

*Proof.* The locality hypothesis removes the block longitudinal Gram term exactly at sufficiently large size in the finite-range case and asymptotically in the exponential case. Equation <a href="#eq:finite-mass-law-nu" data-reference-type="eqref" data-reference="eq:finite-mass-law-nu">[eq:finite-mass-law-nu]</a> then applies. Since $L_a/(L_a-1)\to1$ for every extensive block, taking the limit gives <a href="#eq:thermo-block-law" data-reference-type="eqref" data-reference="eq:thermo-block-law">[eq:thermo-block-law]</a>. ◻

</div>

Equation <a href="#eq:thermo-block-law" data-reference-type="eqref" data-reference="eq:thermo-block-law">[eq:thermo-block-law]</a> is the direct multi-pseudospin analogue of the Gao–Han–Khalaf hard-core factor : each independently conserved pseudospin contributes its own filling factor and its own block-resolved pair mobility.

# Unresolved total filling: a lower envelope, not one branch

The product AGPs are exactly degenerate at $A=0$. When the full zero manifold is exactly their span, even a twist that preserves every block charge requires comparison of the invariant composition sectors. In the projected rank-one subclass, this completeness is equivalent to the determinant-resonance test <a href="#eq:det-resonance-even" data-reference-type="eqref" data-reference="eq:det-resonance-even">[eq:det-resonance-even]</a>. If that test fails, the additional determinant sectors must also be included.

Let $$\mathcal C_n=\{\mathbf n:0\le n_a\le L_a,\ |\mathbf n|=n\}.
 \label{eq:composition-set}$$ For a real twist direction $v$, put $$C_v^{(\mathbf n)}=v_iC_{ij}^{(\mathbf n)}v_j.
 \label{eq:Cv}$$

<div id="prop:lower-envelope" class="proposition">

**Proposition 21** (Physical directional curvature for a block-preserving twist). *Assume the full family preserves all spin-resolved block charges, every sector $\mathbf n\in\mathcal C_n$ has the simple branch described above, and the complete zero space at total pair number $n$ is $\operatorname{span}\{|\mathbf n\rangle:\mathbf n\in\mathcal C_n\}$. Then the fixed-total-number ground energy satisfies $$E_0(tv)=\frac{t^2}{2}
 \min_{\mathbf n\in\mathcal C_n}C_v^{(\mathbf n)}+o(t^2),
 \label{eq:envelope-expansion}$$ and the physical directional curvature is $$\boxed{
 C_v^{\mathrm{phys}}(n)
 =\min_{\mathbf n\in\mathcal C_n}C_v^{(\mathbf n)}.}
 \label{eq:physical-min}$$ A single Hessian tensor need not exist if the minimizing composition depends on $v$.*

</div>

<div class="proof">

*Proof.* Block conservation makes the total sector a direct sum of invariant composition sectors for every sufficiently small $A$. Each lowest branch has zero constant and linear terms and expansion $E_{\mathbf n}(tv)=t^2C_v^{(\mathbf n)}/2+o(t^2)$. The total ground energy is their finite minimum. Taking the minimum of the leading quadratic coefficients gives <a href="#eq:envelope-expansion" data-reference-type="eqref" data-reference="eq:envelope-expansion">[eq:envelope-expansion]</a>. ◻

</div>

<div class="remark">

**Remark 22** (When a determinant resonance is present). *If <a href="#eq:det-resonance-even" data-reference-type="eqref" data-reference="eq:det-resonance-even">[eq:det-resonance-even]</a> fails, formula <a href="#eq:physical-min" data-reference-type="eqref" data-reference="eq:physical-min">[eq:physical-min]</a> is only the envelope inside the balanced product-AGP subspace. The exact physical curvature is obtained by applying the degenerate least-squares construction of Theorem <a href="#thm:degenerate-LS" data-reference-type="ref" data-reference="thm:degenerate-LS">24</a> to the complete kernel described by Theorem <a href="#thm:complete-kernel" data-reference-type="ref" data-reference="thm:complete-kernel">7</a>. A resonant determinant state is exactly flat only when its trace cancellation persists under the chosen twisted family.*

</div>

The trace-free case has a sharp saturation principle. Define the nonnegative directional block masses $$\mu_a(v)=v_i(m_{\mathrm{pair},a}^{-1})_{ij}v_j\ge0.
 \label{eq:mu-a}$$ Then, up to the common normalization, the branch objective is $$\Phi_v(\mathbf n)
 =\sum_a\frac{n_a(L_a-n_a)}{L_a-1}\mu_a(v).
 \label{eq:Phi}$$

<div id="prop:block-polarization" class="proposition">

**Proposition 23** (Block polarization of the minimizing composition). *The function <a href="#eq:Phi" data-reference-type="eqref" data-reference="eq:Phi">[eq:Phi]</a> is separable and concave in the block occupations. At fixed total pair number, a minimizer can be chosen with all but at most one block either empty or full. In particular, if $$n=\sum_{a\in A}L_a
 \label{eq:commensurate-n}$$ for some subset $A$ of blocks, then $$\boxed{C_v^{\mathrm{phys}}(n)=0}
 \label{eq:commensurate-zero}$$ for every direction $v$.*

</div>

<div class="proof">

*Proof.* If two occupations $n_a,n_b$ are both strictly between their bounds, restrict <a href="#eq:Phi" data-reference-type="eqref" data-reference="eq:Phi">[eq:Phi]</a> to the line $(n_a+t,n_b-t)$ while keeping the total fixed. This is a concave function of the allowed integer variable $t$, so one of the endpoints does not increase the objective. Moving to that endpoint saturates at least one of the two blocks. Repeating the exchange leaves at most one partially filled block. For <a href="#eq:commensurate-n" data-reference-type="eqref" data-reference="eq:commensurate-n">[eq:commensurate-n]</a>, fill the blocks in $A$ completely and leave the others empty. Every factor $n_a(L_a-n_a)$ vanishes, so the branch curvature is zero; positivity then makes it the physical minimum. ◻

</div>

This zero-curvature phenomenon does not contradict positive pair mobility. A full block and an empty block have no available hard-core phase space, even if a partially filled version of the same block has positive inverse pair mass. The total-number degeneracy allows the system to choose such a block-polarized composition.

In the thermodynamic trace-free limit, let the total pair density per cell be $$p=\frac nV=\sum_a\ell_a\nu_a.
 \label{eq:p-density}$$ Propositions <a href="#prop:lower-envelope" data-reference-type="ref" data-reference="prop:lower-envelope">21</a> and <a href="#prop:block-polarization" data-reference-type="ref" data-reference="prop:block-polarization">23</a> give the directional envelope $$\boxed{
 D_v^{\mathrm{phys}}(p)
 =\min_{\substack{0\le\nu_a\le1\\
 \sum_a\ell_a\nu_a=p}}
 \sum_a\ell_a\nu_a(1-\nu_a)\mu_a(v).}
 \label{eq:thermo-envelope}$$ A minimizer lies at a vertex of the filling polytope, so at most one block is partially occupied. At densities equal to a sum of complete block capacities, $p=\sum_{a\in A}\ell_a$, the directional stiffness vanishes exactly.

# Generic block-mixing twists

We now drop the block-preserving condition. The unperturbed Hamiltonian still has product pseudospin symmetry, but the first twist source may contain $b_{\lambda i}^{ab}$ with $a\ne b$. Fixed compositions then mix already at second order, and the correct object is a matrix on the complete zero manifold. The formulas below are written first for the product-AGP-complete case certified by Theorem <a href="#thm:det-resonance" data-reference-type="ref" data-reference="thm:det-resonance">9</a>; the same least-squares theorem applies without change to the larger kernel when a synchronization or determinant resonance is present.

Fix the total pair number $n$ and assume, as certified for example by Theorems <a href="#thm:MB3-criterion" data-reference-type="ref" data-reference="thm:MB3-criterion">8</a> and <a href="#thm:det-resonance" data-reference-type="ref" data-reference="thm:det-resonance">9</a>, that the complete zero space is $$\mathcal Z_n=\operatorname{span}\{|\mathbf n\rangle:\mathbf n\in\mathcal C_n\}.
 \label{eq:Zn}$$ Assume this zero manifold is separated from the positive spectrum at the given finite size. Standard finite-dimensional analytic cluster perturbation theory then applies . Let $$\Pi=P_{\ker D(0)^\dagger}
 \label{eq:global-Pi}$$ act on the target of the full null map. For a real direction $v$, put $$B_{\lambda,v}=v_iB_{\lambda i},
 \qquad
 \mathcal B_v|\psi\rangle
 =\bigl(B_{1,v}|\psi\rangle,\ldots,B_{M_s,v}|\psi\rangle\bigr),
 \label{eq:B-map}$$ and define $$T_v=\Pi\mathcal B_v\big|_{\mathcal Z_n}.
 \label{eq:T-v}$$

<div id="thm:degenerate-LS" class="theorem">

**Theorem 24** (Degenerate least-squares curvature). *The exact second-order curvature operator on the zero manifold is $$\boxed{Q_v^{(n)}=T_v^\dagger T_v\succeq0.}
 \label{eq:Q-def}$$ The lowest fixed-total-number branch obeys $$E_0(tv)
 =\frac{t^2}{2}\lambda_{\min}(Q_v^{(n)})+O(t^3),
 \label{eq:E-Q}$$ so its physical directional curvature is $$\boxed{C_v^{\mathrm{phys}}(n)=\lambda_{\min}(Q_v^{(n)}).}
 \label{eq:C-lambda-min}$$*

</div>

<div class="proof">

*Proof.* For a normalized $|\psi\rangle\in\mathcal Z_n$, allow a first-order correction $|\chi\rangle$ orthogonal to $\mathcal Z_n$. To first order in $t$, the null-map residual is $$D(0)|\chi\rangle+\mathcal B_v|\psi\rangle.
 \label{eq:degenerate-residual}$$ Minimization over $\chi$ projects the source onto $\ker D(0)^\dagger$, giving the quadratic form $$\min_\chi\|D(0)\chi+\mathcal B_v\psi\|^2
 =\|T_v\psi\|^2
 =\langle\psi|Q_v^{(n)}|\psi\rangle.
 \label{eq:degenerate-min}$$ Diagonalizing this finite positive form on $\mathcal Z_n$ gives the second-order branch curvatures, and the smallest eigenvalue gives the ground branch. ◻

</div>

## Universal composition-space matrix elements

For the direction $v$, define $$\boldsymbol\beta_v^a=v_i\boldsymbol\beta_i^a.
 \label{eq:beta-v}$$ Let $$\boldsymbol\zeta_v^a
 =\bigl(\zeta_a(B_{1,v}),\ldots,\zeta_a(B_{M_s,v})\bigr),
 \label{eq:zeta-v-a}$$ and, for $a\ne b$, $$\boldsymbol\zeta_v^{a\leftarrow b}
 =\bigl(\zeta_{a\leftarrow b}(B_{1,v}),\ldots,
 \zeta_{a\leftarrow b}(B_{M_s,v})\bigr).
 \label{eq:zeta-v-transfer}$$ Let $\Pi_a$ be the target-kernel projector in the two-particle sector with one pair in block $a$, and let $\Pi_{ab}$ be the target-kernel projector in the two-particle sector with one fermion in each of blocks $a,b$. Set $$\begin{aligned}
 k_v^a&=\|\Pi_a\boldsymbol\zeta_v^a\|^2,
 \label{eq:k-a}\\
 k_v^{a\leftarrow b}&=\|\Pi_{ab}\boldsymbol\zeta_v^{a\leftarrow b}\|^2,
 \label{eq:k-transfer}\\
 \omega_v^{ab}&=
 \left\langle
 \Pi_{ab}\boldsymbol\zeta_v^{b\leftarrow a},
 \Pi_{ab}\boldsymbol\zeta_v^{a\leftarrow b}
 \right\rangle.
 \label{eq:omega}
\end{aligned}$$ The phases of the composition basis and ladder maps are chosen by the natural convention generated by the ordered products of the $\eta_a^+$.

<div id="thm:composition-matrix" class="theorem">

**Theorem 25** (Composition-space matrix law). *In the occupation basis $\{|\mathbf n\rangle:\mathbf n\in\mathcal C_n\}$, the diagonal elements of <a href="#eq:Q-def" data-reference-type="eqref" data-reference="eq:Q-def">[eq:Q-def]</a> are $$\begin{aligned}
 (Q_v^{(n)})_{\mathbf n\mathbf n}
 ={}&\left\|\sum_an_a\boldsymbol\beta_v^a\right\|^2
 +\sum_a\rho_a(n_a)k_v^a
 \nonumber\\
 &+\sum_{a\ne b}
 \tau_{a\leftarrow b}(\mathbf n)k_v^{a\leftarrow b}.
 \label{eq:Q-diagonal}
\end{aligned}$$ All off-diagonal entries vanish except between neighboring compositions. If $$\mathbf m=\mathbf n+\mathbf e_a-\mathbf e_b,
 \qquad
 n_b\ge1,
 \qquad
 n_a\le L_a-1,
 \label{eq:neighbor-composition}$$ then $$\boxed{
 (Q_v^{(n)})_{\mathbf m\mathbf n}
 =\sqrt{
 \frac{
 n_b(L_a-n_a)(n_a+1)(L_b-n_b+1)
 }{L_aL_b}}
 \,\omega_v^{ab}.}
 \label{eq:Q-offdiag}$$ The reverse matrix element is its complex conjugate.*

</div>

<div class="proof">

*Proof.* Apply Lemma <a href="#lem:multiblock-action" data-reference-type="ref" data-reference="lem:multiblock-action">13</a> to the target source $\mathcal B_v|\mathbf n\rangle$. The longitudinal, intrablock, and ordered transfer channels are mutually orthogonal. Product-pseudospin intertwining allows the target-kernel projection to be evaluated on the corresponding one-pair seed. Taking norms gives <a href="#eq:Q-diagonal" data-reference-type="eqref" data-reference="eq:Q-diagonal">[eq:Q-diagonal]</a>.

For two different compositions, diagonal and intrablock channels have different block charges and cannot overlap. A transfer source from $\mathbf n$ through $a\leftarrow b$ has odd block occupations $$N_a=2n_a+1,
 \qquad
 N_b=2n_b-1.
 \label{eq:odd-intermediate}$$ The only source from another composition with the same charge pattern begins at $\mathbf m=\mathbf n+\mathbf e_a-\mathbf e_b$ and transfers one fermion back through $b\leftarrow a$. The overlap coefficient is $$\begin{aligned}
 &\sqrt{
 \tau_{a\leftarrow b}(\mathbf n)
 \tau_{b\leftarrow a}(\mathbf m)}
 \nonumber\\
 &\quad=
 \sqrt{
 \frac{
 n_b(L_a-n_a)(n_a+1)(L_b-n_b+1)
 }{L_aL_b}},
 \label{eq:clebsch-product}
\end{aligned}$$ while the reduced overlap is $\omega_v^{ab}$. This gives <a href="#eq:Q-offdiag" data-reference-type="eqref" data-reference="eq:Q-offdiag">[eq:Q-offdiag]</a>; no other pair of source channels can land in the same product-spin and block-charge sector. ◻

</div>

The square-root factor in <a href="#eq:Q-offdiag" data-reference-type="eqref" data-reference="eq:Q-offdiag">[eq:Q-offdiag]</a> is the normalized product-pseudospin matrix element of pair transfer, $$\frac{
 \langle\mathbf n+\mathbf e_a-\mathbf e_b|
 \eta_a^+\eta_b^-|\mathbf n\rangle
 }{\sqrt{L_aL_b}}.
 \label{eq:eta-transfer-matrix}$$ Thus $Q_v^{(n)}$ is a positive semidefinite tight-binding operator on the bounded composition lattice, with exact Clebsch weights and model-dependent one-pair Gram data as onsite and hopping amplitudes.

<div class="remark">

**Remark 26** (One-pair masses are not complete many-body data). *The eigenvalues of $Q_v^{(1)}/4$ describe the unresolved one-pair curvature eigenmodes. At higher filling they do not determine $Q_v^{(n)}$ by themselves. The separate transfer norms $k_v^{a\leftarrow b}$ and relative overlaps $\omega_v^{ab}$ are also required. This is the precise sense in which the multi-block extension is matrix valued rather than a scalar function of total filling and one-pair eigenmasses.*

</div>

# Minimal two-block interference certificate

The smallest example already rules out a universal scalar total-filling law. Take two blocks with $$L_1=L_2=1,
 \qquad
 F=\sigma_1J_1\oplus\sigma_2J_2,
 \qquad
 \sigma_1\ne\sigma_2,
 \label{eq:minimal-F}$$ and fermion modes $a_\pm,b_\pm$. Define the unperturbed Hermitian null factors $$S_a=n_{a+}-n_{a-},
 \qquad
 S_b=n_{b+}-n_{b-}.
 \label{eq:minimal-S}$$ At one pair, the zero space is $$|a\rangle=a_+^\dagger a_-^\dagger|0\rangle,
 \qquad
 |b\rangle=b_+^\dagger b_-^\dagger|0\rangle.
 \label{eq:minimal-basis}$$ Let $$B=b_+^\dagger a_++b_-^\dagger a_-
 +a_+^\dagger b_++a_-^\dagger b_-
 \label{eq:minimal-B}$$ and consider the analytic positive-square family $$H(A)=\frac12S_a^2+\frac12S_b^2+\frac12(AB)^2.
 \label{eq:minimal-H}$$ Direct fermionic algebra gives $$B|a\rangle=B|b\rangle
 =a_+^\dagger b_-^\dagger|0\rangle
 -a_-^\dagger b_+^\dagger|0\rangle,
 \qquad
 \|B|a\rangle\|^2=2.
 \label{eq:minimal-action}$$ Therefore the curvature operator on $\operatorname{span}\{|a\rangle,|b\rangle\}$ is $$Q=2
 \begin{pmatrix}
 1&1\\
 1&1
 \end{pmatrix},
 \qquad
 \operatorname{spec}Q=\{0,4\}.
 \label{eq:minimal-Q}$$ Each composition has positive diagonal curvature, but the antisymmetric superposition $$\frac{|a\rangle-|b\rangle}{\sqrt2}
 \label{eq:minimal-dark}$$ is annihilated by $B$ and remains an exact zero mode of <a href="#eq:minimal-H" data-reference-type="eqref" data-reference="eq:minimal-H">[eq:minimal-H]</a> for every $A$. The unresolved physical curvature is therefore zero. No scalar law built only from the total filling and the diagonal block responses can reproduce both the individual branches and their interference.

# Relation to the parent theorem and physical interpretation

For $q=1$, there are no transfer channels and no composition degeneracy. Equation <a href="#eq:MB-positive" data-reference-type="eqref" data-reference="eq:MB-positive">[eq:MB-positive]</a> becomes $$C_{ij}^{(n)}
 =n^2\Gamma_{ij}+\rho_n\Theta_{ij}
 =\rho_n C_{ij}^{(1)}+(n^2-\rho_n)\Gamma_{ij},
 \label{eq:recover-parent}$$ which is exactly the parent law <a href="#eq:parent-law" data-reference-type="eqref" data-reference="eq:parent-law">[eq:parent-law]</a>. The present result is therefore a genuine extension, but it also identifies why the extension cannot retain the same scalar form.

The hierarchy is:

1.  A single nonzero singular value gives one maximal pseudospin, one AGP state at fixed pair number, and the scalar reduction <a href="#eq:parent-law" data-reference-type="eqref" data-reference="eq:parent-law">[eq:parent-law]</a>.

2.  Several distinct singular values give a product pseudospin, a composition-degenerate zero manifold, and the block-resolved branch law <a href="#eq:MB-positive" data-reference-type="eqref" data-reference="eq:MB-positive">[eq:MB-positive]</a> when the twist preserves the blocks.

3.  A block-mixing source turns the same zero manifold into the composition-space eigenvalue problem <a href="#eq:Q-def" data-reference-type="eqref" data-reference="eq:Q-def">[eq:Q-def]</a>–<a href="#eq:Q-offdiag" data-reference-type="eqref" data-reference="eq:Q-offdiag">[eq:Q-offdiag]</a>.

This has two complementary physical readings.

First, if the singular blocks label genuinely conserved pairing channels, then $\nu_a$ are thermodynamic control parameters and <a href="#eq:thermo-block-law" data-reference-type="eqref" data-reference="eq:thermo-block-law">[eq:thermo-block-law]</a> is the natural stiffness law. Each channel contributes its own hard-core factor and pair mass, with no additional finite-density dressing in the trace-free local limit.

Second, if only the total particle number is fixed, the exact QGN degeneracy allows the system to redistribute pairs between blocks. The physical response is then an envelope or a matrix eigenvalue, not the curvature of the weighted AGP chosen by $F$. The block-polarized zero-curvature points in Proposition <a href="#prop:block-polarization" data-reference-type="ref" data-reference="prop:block-polarization">23</a> are the multi-pseudospin analogue of the reducibility obstruction in the parent paper: positive local pair mobility does not by itself guarantee a positive unresolved many-body stiffness when the paired ground space is reducible.

The result does not alter the one-block Model-II theorem of the parent paper. It closes one of that paper’s explicitly deferred theorem boundaries and, in the projected rank-one subclass, replaces the former irreducibility assumptions by the finite tests of Section <a href="#sec:projected-irred" data-reference-type="ref" data-reference="sec:projected-irred">4</a>. The remaining model-dependent questions are whether the physical twist preserves the singular blocks, whether additional terms select a composition, whether a concrete material realization avoids synchronization and determinant resonances, and whether its block one-pair branches admit simple mass formulas.

# Reproducibility certificates

Two independent finite-Fock-space certificates accompany this draft.

The source-reduction certificate uses canonical skew blocks and checks:

1.  the longitudinal, intrablock, and transfer decomposition in Lemma <a href="#lem:multiblock-action" data-reference-type="ref" data-reference="lem:multiblock-action">13</a>;

2.  the transfer inner-product identity <a href="#eq:transfer-inner-n" data-reference-type="eqref" data-reference="eq:transfer-inner-n">[eq:transfer-inner-n]</a>;

3.  the diagonal and nearest-neighbor matrix elements <a href="#eq:Q-diagonal" data-reference-type="eqref" data-reference="eq:Q-diagonal">[eq:Q-diagonal]</a> and <a href="#eq:Q-offdiag" data-reference-type="eqref" data-reference="eq:Q-offdiag">[eq:Q-offdiag]</a>; and

4.  the exact spectrum $\{0,4\}$ in the minimal interference example.

Across the retained random trials, its largest discrepancy is $$3.553\times10^{-14}.
 \label{eq:numerical-error}$$

The connectivity certificate independently reconstructs the associative and Lie closures, trace-profile ranks, and many-body nullities. In the center-resonant counterexample it also tests the two center candidates directly: the real-span residual is $1.59\times10^{-15}$ for $i(I_A\oplus I_B)$ and $2$ for $i(I_A\oplus-I_B)$. The two explicit determinant states are annihilated by every factor to residual at most $1.11\times10^{-16}$. The basic cases give

<div class="center">

| Case                           | $\dim_{\mathbb C}\mathcal A$ | $\dim_{\mathbb R}\mathfrak g$ | Even-sector nullities |
|:-------------------------------|:----------------------------:|:-----------------------------:|:---------------------:|
| Center-resonant counterexample |             $8$              |              $7$              |     $(1,2,5,2,1)$     |
| Full-center repaired model     |             $8$              |              $8$              |     $(1,2,3,2,1)$     |
| Unequal capacities $2,3$       |             $13$             |             $12$              |    $(1,2,3,3,2,1)$    |
| Synchronized control           |             $4$              |              $4$              |    $(1,4,10,4,1)$     |

</div>

Two adversarial cases exercise clauses not covered by those even-sector totals:

<div class="center">

| Case                                                    | $\dim_{\mathbb C}\mathcal A$ | $\dim_{\mathbb R}\mathfrak g$ | Nullities for every total particle number |
|:--------------------------------------------------------|:----------------------------:|:-----------------------------:|:------------------------------------------|
| Odd resonance $(L_1,L_2)=(1,2)$                         |             $5$              |              $4$              | $(1,0,2,2,2,0,1)$                         |
| One-dimensional block plus synchronized pair, $(1,2,2)$ |             $5$              |              $4$              | $(1,0,5,6,14,8,14,6,5,0,1)$               |

</div>

In the first case $Z=\operatorname{span}\{(2,1)\}$ and $\boldsymbol\delta=(1,-2)$ has odd capacity, so the even sectors remain exactly the product-AGP span while the three-particle sector contains two determinant dark states. The second case has synchronization class $\{B,C\}$, center $Z=\operatorname{span}\{(2,1,1)\}$, the even resonance $(0,2,-2)$, and the odd resonances $(\pm1,0,\mp2)$ and $(\pm1,\mp2,0)$. In both models the complete charge-sector formula <a href="#eq:complete-kernel" data-reference-type="eqref" data-reference="eq:complete-kernel">[eq:complete-kernel]</a> agrees with exact diagonalization in every charge sector, with zero mismatches.

The computations are not used as proofs; they check the exact algebraic dimensions, Lie-center membership, parity conditions, synchronization multiplicities, and counting formulas.

# Conclusion

Several inequivalent pairing singular values change both the response law and the zero-mode problem. In the Hermitian one-body QGN class, distinct singular values force independently conserved pseudospins. At fixed composition the one-body source decomposes into longitudinal, intrablock, and interblock-transfer channels with universal filling coefficients. In the blockwise trace-free local limit, $$D_{s,ij}^{(\boldsymbol\nu)}
 =
 \sum_a\ell_a\nu_a(1-\nu_a)
 (m_{\mathrm{pair},a}^{-1})_{ij}.
 \label{eq:conclusion-law}$$

The projected-Hubbard irreducibility problem has a separate two-layer answer. Connectivity and nonsynchronization determine the semisimple constraint algebra and are exactly what is required for a unique product AGP in each balanced composition. The trace-profile matrix determines the Lie center and decides whether block-imbalanced full determinants survive at fixed total particle number. Full associative-algebra generation does not replace this center test. Odd-capacity resonances leave every even pairing sector intact but create exact fermionic determinant zero modes in odd sectors.

Without block resolution, the same algebra predicts a lower envelope of branches for block-preserving twists and a positive composition-space matrix for generic twists. If a determinant resonance is present, the complete-kernel version of that matrix must include the additional full-determinant sectors. Explicit finite and lattice-ready constructions show that distinct-singular-value models satisfying all uniqueness conditions are nonempty and structurally broad.

The theorem architecture is therefore complete in finite dimension for the blockwise-tight rank-one Hermitian projected-Hubbard class. The next physics task is to identify natural band structures and interactions realizing these blocks, certify their trace profiles under the physical twist, and compute their block pair masses and locality data.

<div class="thebibliography">

99

N. Sledgianowski, “Exact Flat-Band Stiffness from Pair Mobility in a QGN Model: Finite-size obstructions, a many-body reduction theorem, and a proof of the Gao–Han–Khalaf formula for Model II,” manuscript (2026), <https://github.com/Sevii/superconductors/tree/master/exact_flat_band_stiffness_physics_arc_complete>.

Z. Han, J. Herzog-Arbeitman, B. A. Bernevig, and S. A. Kivelson, “Quantum Geometric Nesting and Solvable Model Flat-Band Systems,” *Phys. Rev. X* **14**, 041004 (2024), arXiv:2401.04163.

Q. Gao, Z. Han, and E. Khalaf, “Bootstrapping Flatband Superconductors: Rigorous Lower Bounds on Superfluid Stiffness,” *Phys. Rev. Lett.* **136**, 076503 (2026), arXiv:2506.18969.

T. Kato, *Perturbation Theory for Linear Operators*, 2nd ed. (Springer, Berlin, 1976).

M. Tovmasyan, S. Peotta, P. Törmä, and S. D. Huber, “Effective Theory and Emergent $SU(2)$ Symmetry in the Flat Bands of Attractive Hubbard Models,” *Phys. Rev. B* **94**, 245149 (2016).

B. C. Hall, *Lie Groups, Lie Algebras, and Representations: An Elementary Introduction*, 2nd ed. (Springer, Cham, 2015).

W. Fulton and J. Harris, *Representation Theory: A First Course* (Springer, New York, 1991).

</div>
