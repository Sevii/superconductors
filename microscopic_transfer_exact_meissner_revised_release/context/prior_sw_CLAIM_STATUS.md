# Claim status: infrared-stable microscopic Meissner transfer

## Established in this package

1. **Target-metric response theorem.** For
   \[
   H_G(A)=\frac{\varepsilon}{2}D(A)^\dagger G(A)D(A),
   \qquad g_-I\le G(0)\le g_+I,
   \]
   with exact transverse source \(D^\dagger S=0\), the static curvature and every imaginary-frequency kernel satisfy
   \[
   \varepsilon g_-S^\dagger S\le C_G\le K_G(i\zeta)\le \varepsilon g_+S^\dagger S.
   \]
   The estimate is independent of the smallest positive singular value of \(D\).

2. **Current-source \(H_{-1}\) bound.** For \(G=I+2\eta X\),
   \[
   0\le K_G(i\zeta)-C_G
   \le \frac{4\varepsilon\eta^2}{g_-}S^\dagger X^\dagger X S.
   \]
   With \(\varepsilon=\lambda^6\Delta\) and \(\eta=\lambda^2\), the entire static response changes at order \(\lambda^8\Delta\), while the possible temporal defect starts at order \(\lambda^{10}\Delta\).

3. **Necessary refinement of the previous roadmap.** A generic SW remainder has both a zero-field state source \(QRP\) and a current source \(QJ_RP\). Current-source control alone is insufficient. A two-sided row-ideal remainder removes the state source exactly.

4. **Sharp anisotropic counterexample.** A bounded local row-metric perturbation can have a nonzero constant current \(H_{-1}\) weight and a finite Drude-Meissner mismatch. Thus boundedness is not tightness, and row-relative form alone does not prove the exact order-of-limits identification.

5. **Isotropic transverse theorem.** Under hydrodynamic closure, quasi-locality, and a point-group scalar zero-momentum vector symbol, transverse-longitudinal mixing is \(O(q^2)\), the soft \(H_{-1}\) weight is \(O(q^2)\), and
   \[
   D_K=D_A=D_M
   =\lambda^6\Delta D_{\rm parent}+O(\lambda^8\Delta).
   \]

6. **Finite-cluster row-ideal criterion.** A Hermitian cluster term is state-source safe exactly when
   \[
   r_X=q_Xr_Xq_X=D_X^\dagger X_XD_X,
   \qquad X_X=(D_X^+)^\dagger r_XD_X^+.
   \]
   For the swap parent, the finite-cluster pseudoinverse grows only polynomially, so exponential locality survives with a small loss of exponent.

7. **First degree-eight microscopic audit.** Two neighboring equal-gap bridge controls give the connected coefficient
   \[
   -\frac{\lambda^8\Delta b^4}{2\bar\Delta_m^3}\{B_e,B_f\}^2.
   \]
   On a three-cell one-pair product-AGP space,
   \[
   P_{\mathcal Z_1}\{B_0,B_1\}^2P_{\mathcal Z_1}
   =\frac{16}{3}\begin{pmatrix}1&1\\1&1\end{pmatrix},
   \]
   so the bridge-only overlap is not in the local two-sided parent-row ideal.

8. **Order-eight response-safe completion.** A Peierls-covariant quasi-local active counterinteraction can cancel the unsafe degree-eight projection and replace the safe part by a point-group-averaged geometric row metric. This proves the desired transfer for the completed effective expansion through weighted degree eight.

## Not yet established

- The complete sum of every weighted-degree-eight cluster in the **unmodified** microscopic family has not been enumerated. Other clusters may cancel the unsafe bridge-only compression.
- The required counterinteraction has not been synthesized using only the original restricted two-electron exchange/correlated-hopping channels.
- Hydrodynamic closure of the full operator-valued degree-eight row metric has not yet been verified.
- The existing microscopic construction is one-dimensional; the point-group theorem applies to a symmetric \(d\ge2\) extension whose bond controls still need to be written explicitly.
- No convergent all-orders local SW/Feshbach map has been proved. Therefore the package does not claim an exact, unmodified, finite-\(\lambda\) microscopic Meissner theorem.
- No material parameter extraction, finite-temperature transition theorem, or Maxwell-field screening calculation is claimed.

## Reviewer-safe headline

> We prove an infrared-stable response transfer theorem for row-relative, point-group-isotropic Schrieffer-Wolff remainders; identify a sharp anisotropic failure mode; and reduce the unmodified microscopic problem to a finite degree-eight cluster projection. The first bridge-overlap cluster is explicitly unsafe, so cancellation or a counterchannel must be demonstrated rather than assumed.
