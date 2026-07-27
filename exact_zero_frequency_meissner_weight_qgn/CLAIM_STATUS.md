# Claim status and theorem boundary

## Unconditional theorem proved in this paper

For the explicitly defined local gauge-covariant swap-QGN/AGP parent on periodic lattices of dimension `d >= 2`, every finite-volume lattice-transverse field satisfies

\[
\frac{1}{N}\mathcal K_T(q,i\zeta)
=
\frac{1}{N}\mathcal C_T(q)
=
j\frac{2n(N-n)}{N(N-1)}
\]

for every allowed transverse wave vector and every `zeta > 0`. Therefore, along any fixed-density sequence `n/N -> rho` with `0 < rho < 1`,

\[
\boxed{\kappa_{\rm K}=D_{\rm A}=D_{\rm M}=2j\rho(1-\rho)>0.}
\]

In the pair-phase convention `Q=2A`,

\[
\boxed{D_s=\frac{j}{2}\rho(1-\rho).}
\]

This is an unconditional zero-temperature transverse matter-Meissner/superfluid-weight theorem for the exact parent. The result is not inferred from ODLRO or a uniform gap: the finite-volume response is evaluated directly, while the fixed-number gap closes as `L^{-2}`.

## Definitions now fixed explicitly

- `K_{L,T}(q,0)` means the finite-volume static curvature evaluated on the transverse profile: the frequency limit is taken before the thermodynamic limit.
- The finite-`q` Abelian kernel `K_{L,T}(q,i zeta)` is separately defined.
- Because allowed momenta depend on finite volume, `q -> 0` is taken through admissible momentum sequences and cofinal size subsequences.
- A thermodynamic-first frequency limit at fixed nonzero `q` is a distinct order unless uniform temporal control is also proved there.
- For the exact swap parent the distinction disappears because the finite-size transverse kernel is independent of both `q` and `zeta`.

## General exact theorems

1. **Finite-volume Kubo–Kohn identity.** The imaginary-frequency current kernel differs from the Kohn curvature by an explicit positive operator.
2. **Sharp temporal criterion.** Under a uniform bound on the total current `H^{-1}` norm, temporal Kohn–Drude equality is equivalent to uniform `H^{-1}` tightness/equi-integrability at zero energy.
3. **Positive-square dynamical least squares.** For `H(A)=D(A)^†D(A)/2`, the complete imaginary-frequency response is a Tikhonov-regularized target-space residual.
4. **Conditional transverse completion.** Temporal equality plus continuity of the finite-volume-static transverse kernel at zero momentum identifies Kohn, Abelian Drude, and Meissner weights in a selected thermodynamic phase.
5. **Space-time summability criterion.** An absolutely summable Euclidean current-correlation envelope gives both the required temporal tightness and static transverse continuity.
6. **Exact graph-Hodge completion.** For the swap parent, the temporal and transverse problems are solved algebraically at finite size.
7. **Protected multiblock floor.** A swap-row source lying in the complete target kernel supplies a positive transverse operator floor even when the target Laplacian has nonzero swap/residual off-diagonal blocks.
8. **Independent failure mechanisms.** A soft positive target singular mode violates the temporal order; a clean normal gas has Drude response without a transverse Meissner response.

## Exact model properties

- The fixed-`2n` zero space is one dimensional and spanned by the paired Dicke/AGP state.
- The paired compression is the symmetric-exclusion generator, equivalently the isotropic spin-1/2 ferromagnet at fixed magnetization.
- The static response is the graph cycle-space norm.
- The finite-frequency response is an exact edge-Laplacian filter.
- Gauge-gradient fields have zero static curvature.
- Every divergence-free field has a purely diamagnetic, frequency-independent positive response at finite size.
- The pair correlator is exactly `n(N-n)/[N(N-1)]` off diagonal.
- The electronic-link Meissner kernel is `2j` times that off-diagonal correlator.
- The onsite `U > 0` selects seniority zero but drops out of the paired response formulas.
- The proof permits a closing complete fixed-number gap.

## Multiblock consequence

For a selected product composition with block densities `rho_b`, the swap backbone gives

\[
D_{\rm M}\ge 2\sum_b j_b\rho_b(1-\rho_b).
\]

The operator inequality does not assume that the complete target Laplacian is block diagonal. A scalar response for an unresolved degenerate composition manifold still requires a branch selector, coherent state, or density matrix.

## What is not claimed

- No finite-temperature phase theorem.
- No theorem for every projected QGN parent.
- No unrestricted real-axis conductivity decomposition beyond the Abelian zero-frequency kernel.
- No scalar response for a degenerate zero manifold without a physical selection rule.
- No one-dimensional magnetic-screening claim; the transverse Meissner definition requires `d >= 2`.
- No claim that ODLRO alone implies a Meissner effect.
- No inference from an addition/removal charge gap to current-sector regularity.
- No coupling to a dynamical Maxwell field or material penetration-depth prediction.
- No unconditional transfer to the original finite-coupling microscopic Schrieffer–Wolff realization.

## Separate microscopic transfer obligation

The microscopic construction reproduces the exact parent at leading weighted order and carries a quasi-local higher-order remainder. Because the exact parent has soft longitudinal modes, local operator-norm smallness alone does not control the response. A finite-coupling transfer theorem must additionally control the selected branch, the zero-field state source, the physical current source in the `H^{-1}` topology, the finite-wave-vector transverse window, and residual low/high block coupling. Those questions belong to the separate microscopic-transfer paper and do not weaken the exact-parent theorem proved here.
