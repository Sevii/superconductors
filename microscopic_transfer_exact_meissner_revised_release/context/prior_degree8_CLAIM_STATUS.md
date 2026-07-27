# Claim status: degree-eight microscopic transfer and dark-branch completion

## Analytically established in the paper

1. **Complete degree-eight resolvent class list for the single-flavor screened hierarchy.** After the two outer weight-two bridges, the only nonzero high-resolvent classes of remaining weight four are \(V_1^4\), \(V_{2,Q}^2\), and one matched \(W_4\) insertion. Every \(V_1V_1V_{2,Q}\) ordering vanishes by control-toggle parity. Degree-six vertices cannot enter at total weight eight.

2. **Closed source formula.** On the complete product-AGP manifold,
   \[
   R^{(8)}P_{\mathcal Z}
   =\frac{b^2}{\bar\Delta_m^3}
   \left[b^2\mathcal M-\frac{16s^4}{\bar\Delta_s^2}\mathcal B_2\right]P_{\mathcal Z},
   \]
   where
   \[
   \mathcal B_2=\sum_eB_e^2,
   \qquad
   \mathcal M=\sum_eB_e^4+\frac12\sum_{e<f}[B_e,B_f]^\dagger[B_e,B_f].
   \]
   The commutator form is linked because disjoint bridge operators commute.

3. **Closed-shell Kramers correction.** The time-reversal-invariant duplicate channel replaces \(\mathcal M\) by
   \[
   \mathcal M_\Theta=\frac12\sum_eB_e^4+\frac12\sum_{e<f}[B_e,B_f]^\dagger[B_e,B_f].
   \]
   It halves only the local fourth-power coefficient; the inter-bond coefficient is unchanged.

4. **All-filling obstruction.** On the three-cell two-pair sector there are explicit states with
   \[
   \langle\Phi|\mathcal B_2|0,2\rangle=0,
   \quad
   \langle\Phi|\mathcal M|0,2\rangle=8\sqrt3,
   \quad
   \langle\Phi|\mathcal M_\Theta|0,2\rangle=4\sqrt3.
   \]
   Therefore no scalar retuning compatible with the required degree-four cancellation can preserve the complete product-AGP manifold for nonzero bridge amplitude. The one-pair cancellation is accidental and does not extend to higher filling.

5. **Explicit degree-eight completion.** A finite-range active counterinteraction cancels the complete unsafe degree-eight source and leaves \((c_s-c_a)\sum_eK_e^\dagger K_e\), positive under the same-chirality dominance condition.

6. **Exact coherent dark branch.** For
   \[
   \eta_-^+=\eta_1^+-\eta_2^+,
   \qquad
   |\Omega_n^-\rangle\propto(\eta_-^+)^n|0\rangle,
   \]
   every bridge obeys \([B_e,\eta_-^+]=0\) and \(B_e|\Omega_n^-\rangle=0\). Inside the product-composition manifold the common bridge kernel is one dimensional. The exact compressed \(\sum_eB_e^2\) matrix is an irreducible Jacobi matrix with the entries displayed in the paper.

7. **All-orders zero-field eigenbranch.** With every control in its ground orbital, the dark AGP is an exact zero-energy eigenstate of the screened microscopic Hamiltonian at \(A=0\) for every finite coupling for which the Hamiltonian is defined.

8. **Positive-semidefinite finite-coupling completion.** If the direct bridge-square coefficient satisfies \(\beta>t_B^2/d_*\), each local control-edge block is positive semidefinite and every local zero vector has its control in \(g\) and obeys \(B_e\psi=0\). Adding the exact swap base parent makes the dark AGP the unique fixed-number zero state on a connected lattice.

9. **Microscopic Meissner floor.** In a complete \(d\ge2\) Peierls-covariant extension, the swap-row source cannot be canceled by the gadget rows. For all allowed transverse momenta and all imaginary frequencies,
   \[
   \mathcal K_{L,T}(q,i\zeta)\ge\mathcal C_{L,T}(q)\ge
   \lambda^6\Delta\|\mathcal S_{\rm sw,T}(q)\|^2.
   \]
   On the dark branch,
   \[
   D_{{\rm sw},L}^{(-)}=(j_1+j_2)\frac{n(2L-n)}{L(2L-1)},
   \]
   and every thermodynamic Kohn, Abelian, and transverse liminf has the positive floor
   \[
   2\lambda^6\Delta(j_1+j_2)\rho(1-\rho).
   \]

## Computationally certified

- Every degree-eight source class on all product-AGP fillings of the three-cell ring.
- Exact vanishing of the mixed \(V_1V_1V_2\) class.
- The linked commutator-square identity and explicit source completion.
- Failure of scalar tuning at two and three pairs for both single-flavor and Kramers variants.
- Independent small-coupling extraction of the single-flavor and Kramers Feshbach coefficients using random noncommuting active matrices.
- Bridge darkness, one-dimensional composition-space kernel, exact Jacobi entries, and the finite-volume swap-floor formula for \(L=2,3\).

Every retained certificate ends in `PASS`.

## Not established

- The original exactly matched finite-coupling microscopic Hamiltonian is not proved to have the dark eigenbranch as its ground state.
- The full product-AGP manifold is not preserved by the unmodified model beyond weighted degree six; the paper proves its failure at degree eight.
- The explicit order-eight counterinteraction has not been synthesized solely from the original restricted two-electron exchange/correlated-hopping vertices. Its normal ordering contains higher-body active terms.
- No convergent all-orders recursive completion preserving the full composition manifold is proved.
- The positive-semidefinite dark-branch family is a modified overcompensated microscopic model, not the original matched Schrieffer-Wolff realization.
- The theorem gives a strictly positive floor for the total Kohn, Abelian, and Meissner kernels. It does not yet prove equality of the additional gadget contribution in all three limits.
- No material parameter extraction, finite-temperature transition theorem, real-axis conductivity decomposition, or dynamical Maxwell-field screening calculation is claimed.

## Reviewer-safe headline

> The complete weighted-degree-eight audit rules out an all-composition transfer of the unmodified matched microscopic family, supplies an explicit finite-order counterterm completion, and reveals an exact coherent dark AGP. A strictly positive local completion turns that branch into a finite-coupling microscopic ground state with a rigorous zero-frequency transverse Meissner floor.
