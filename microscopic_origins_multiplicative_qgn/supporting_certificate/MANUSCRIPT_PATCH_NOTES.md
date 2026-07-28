# Manuscript patch notes: explicit SVD and boundary scope

## 1. Correct the existing charge proposition’s wording

The current sentence

> “For \(L\ge2\) and the product AGP \(|n_1,n_2\rangle\), \(H_{\rm ch}\) is diagonal in composition space and …”

should be replaced by an explicit compression statement:

> “For fixed total pair number \(n\), the compression of \(H_{\rm ch}\) to
> \(\mathcal Z_n=\operatorname{span}\{|n_1,n_2\rangle:n_1+n_2=n\}\) is diagonal:
> \[
> P_{\mathcal Z_n}H_{\rm ch}P_{\mathcal Z_n}
> =\sum_{n_1+n_2=n}V_{\rm ch}^{(L)}(n_1,n_2)
> |n_1,n_2\rangle\langle n_1,n_2|.
> \]”

The existing formula for \(V_{\rm ch}^{(L)}\) remains exact. The correction prevents “diagonal compression” from being read as the invariant-subspace claim \(H_{\rm ch}\mathcal Z_n\subseteq\mathcal Z_n\).

## 2. Insert the periodic SVD proposition

Insert the periodic part of `paper_insert.tex` immediately after the projected charge proposition. It supplies:

- the exact action \(R|c\rangle=(H_{\rm ch}-\mu_c)|c\rangle\);
- the closed singular values \(A_c\sigma_L(c)\);
- normalized left singular vectors;
- the exact SVD and Gram operator;
- the general norm and rank formulas;
- the complete periodic zero set.

This closes T2.1 without changing the projected selector formula.

## 3. Place the open-path result in the numerical appendix

The microscopic theorem itself is periodic. The open corollary is needed because the DMRG/MPS plan removes the wrap bond for open calculations. It is best placed in a numerical/certificate appendix rather than presented as part of the periodic downfolding theorem.

The open path obeys

\[
H_{{\rm ch},L}^{\rm op}=(8A_c-2U_c)n+8A_cZ-V_cB,
\]

and its selector coefficient is

\[
16A_c\frac{L-2}{L^2}n_1n_2.
\]

The open singular values depend on \(V_c/A_c\); the manuscript or numerical report must state that ratio whenever open normalized leakage is quoted.

## 4. Narrow the unscreened corollary’s interpretation

The operator expansion for the unscreened family is unchanged. The composition formula should be described as the compression of the order-four term, not automatically as the spectrum of an isolated low multiplet.

A safe replacement is:

> “The compression of the order-four charge operator to the product-AGP composition space contains the displayed positive selector. Proposition~`prop:charge-leakage-periodic` gives its simultaneous off-manifold singular values. Therefore, without an additional order-four isolation mechanism, the selector is a projected law rather than an autonomous composition-space low-energy Hamiltonian.”

## 5. Revise discussion and conclusion language

Statements that the generic unscreened family directly realizes a robust low-energy “selector-plus-mixing problem” should distinguish three possibilities:

1. the exact projected selector law;
2. the full charge-selection problem in the seniority-zero Hilbert space;
3. a modified model with order-four or stronger isolation of the product-AGP composition manifold.

For open chains, emphasize that endpoint terms change both the selector coefficient and the leakage singular values. The periodic \(L=4,n=3\) endpoint zeros are not open-boundary zeros.

## 6. Ready-to-apply files

The package includes both a complete revised TeX source and a unified diff:

- `microscopic_origins_multiplicative_qgn_leakage_updated.tex`;
- `microscopic_origins_multiplicative_qgn_leakage.patch`.

The revised source integrates the periodic SVD proposition, adds the open-path result to the charge-combinatorics appendix, and narrows the abstract, theorem interpretation, discussion, and conclusion wherever “selector-plus-mixing” previously implied an isolated product-AGP low-energy manifold.
