# Claim status: weighted row-space Hodge theorem

## Closed exactly at finite volume

1. **Algebraic existence.** For `H_s = D^† G D / 2`, with complete zero projector `Z = P_ker(D)`, a factorization

   `J_R P = D^† eta P`

   exists if and only if `Z J_R P = 0`.

2. **Full-zero-space simplification.** When `P = Z`, the assumed condition `P J_R P = 0` already implies existence. When `P < Z`, it does not: the residual current can mix the selected branch into `Z-P`.

3. **Exact weighted representative.**

   `eta_* P = G D (D^† G D)^+ J_R P = (1/2) G D H_s^+ J_R P`.

4. **Loewner-order variational identity.** For every feasible `eta = eta_* + kappa`,

   `P eta^† G^{-1} eta P = (1/2) P J_R H_s^+ J_R P + P kappa^† G^{-1} kappa P`.

   The Loewner minimizer is unique and the minimum has no duality gap.

5. **Operator-norm equality condition.** Let `C_*` and `C_kappa` denote the two positive terms above and `lambda_* = ||C_*||`. A feasible witness attains the same operator norm as `eta_*` exactly when

   `C_kappa <= lambda_* P - C_*`.

   Hodge projection is sufficient, but is not necessary in general when `rank(P) >= 2`. The projected witness is the unique norm minimizer only when the spectral slack vanishes, in particular for `rank(P)=1`.

6. **Approximate-factorization fallback.** If `cP = Q_D(J_RP-D^†etaP)`, then for every `theta>0`,

   `b <= (1+theta)/V ||P eta^†G^{-1}eta P|| + (1+theta^{-1})/(2V) ||P c^† H_s^+ c P||`.

## Quasi-local existence proved under explicit sufficient hypotheses

1. **Buffered cluster-current criterion.** For every cluster, the theorem now lists as separate hypotheses:

   - convergence of the current-cluster decomposition;
   - frustration freeness `p_X P = P`;
   - local current safety `p_X j_X p_X = 0`.

   Then `eta_X = (D_X^+)^† j_X` gives the local factorization. Its support is contained in `X` enlarged by the fixed parent-row buffer `R_D`.

2. **Explicit exponential decay with support bookkeeping.** If `||D_X^+|| <= C_D(1+diam X)^nu` and the current interaction decays as `exp(-mu diam X)`, then for every `mu'<mu`,

   `||eta||_(mu';R_D) <= N_(R_D) exp(2 mu' R_D) C_D M_nu(mu-mu') ||J_R||_mu`.

   This also supplies convergence of the infinite-volume sum.

3. **Covariant row-metric residual.** If

   `R(A) = D(A)^† Y(A) D(A) / 2`,

   then `eta_0 = Y(0) dot(D)(0)P / 2` is a closed-form local witness.

4. **Uniformly reduced-gapped analytic route.** Exponential decay follows if the active row range has an analytic constant-rank projector and a uniformly invertible reduced matrix throughout a complex strip. These hypotheses exclude a genuine real-axis soft rank drop unless its singularity has already been analytically removed.

## Finite-range Fourier factorization is decidable

1. **One dimension.** Smith-normal-form divisibility is necessary and sufficient.
2. **All dimensions.** Over an effective coefficient field, Laurent-module membership is reduced to polynomial-module saturation `N:h^infinity` and decided by a terminating Gröbner-basis/syzygy computation. Positive output constructs a finite-range witness; negative output gives an algebraic nonmembership obstruction.
3. **Complexity boundary.** No polynomial worst-case Gröbner complexity is claimed. Floating-point material data require exact reconstruction or interval enclosure for a rigorous nonmembership claim.

## Gapless soft modes: characterized but not generally solved

With a controlled spectral decomposition into soft projectors `E_a`, bounded Hodge cost is equivalent to channelwise source vanishing at least as fast as the associated soft singular value. This is a theorem about the inverse-energy cost. It is not a general positive theorem for exponentially local, non-polynomial factorization across rank drops; analytic divisibility remains additional.

## Local computability established in three regimes

1. Sparse local-row bounds, polynomial in the number of retained local monomials.
2. Chebyshev-localized inverse-metric bounds with explicit range and error.
3. Fixed-fiber momentum-symbol bounds with rigorous grid-plus-Lipschitz enclosure.

## Not claimed

- No universal polynomial-time exact algorithm for arbitrary dense many-body cluster matrices.
- No general exponentially local soft-rank-drop division theorem.
- No proof yet that the complete unmodified microscopic degree-eight current satisfies the cluster criterion.
- No automatic Meissner theorem from a bounded current `H^{-1}` norm alone; low-energy tightness and transverse continuity remain separate.
