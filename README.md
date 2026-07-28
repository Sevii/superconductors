
*Welcome to the Sledgeworx.dev Lab Superconductors repository.*

This is where we post our research papers on superconductors.

## What is Sledgeworx.dev Labs?

It is a LLM powered mathematical and compuational research group.

## What is the hypothesis / alpha?

My theory is that science as a whole has been under math'ed its entire history. Mathematics is hard. Not every Scientist is Terrance Tao.
Not enough energy has been put into mathematics of the sciences.

Now that LLMs can perform top level mathematics we have the oppertunity to try to do more math.

## What is our plan?

The plan is to use AI to help/do mathematical and computational science on superconductors. To start we are examing proofs around superconductivity.

## How are our papers verified?

We use LLM adversarial review then post them in public.

**If you have a paper that you think needs more math effort invested please contact us.**

nick@sledgeworx.dev

## Our Papers

### Exact Flat-Band Stiffness from Pair Mobility in a QGN Model

We prove the Gao–Han–Khalaf conjecture that the superfluid stiffness of a quantum-geometric-nesting superconductor is fixed exactly by the two-particle pair mass and the hard-core filling factor, with no finite-size correction for their Model II. Along the way we establish a general many-body reduction theorem and pin down the precise conditions — irreducibility of the paired ground branch and a properly lifted flat connection — under which this striking single-pair reduction holds.

📄 [Read the paper (PDF)](exact_flat_band_stiffness_physics_arc_complete/paper/exact_flat_band_stiffness_physics_arc.pdf)

### Multi-Block Reduction and Irreducibility Laws for QGN Superconductors

*An addendum to the paper above.* The original result assumed the pairing form has a single nonzero singular value. Here we handle several inequivalent ones, where the symmetry splits into a product of independent pseudospin blocks. We derive an exact block-resolved filling law — the stiffness becomes a sum of independent hard-core contributions, one per block — and give a sharp microscopic criterion for when the paired ground state is unique: blockwise frame connectivity controls the semisimple part of the constraint algebra, while trace profiles control its center. When block charge is not conserved, the single filling polynomial is replaced by a composition-space matrix law.

📄 [Read the paper (PDF)](exact_flat_band_stiffness_multiblock_qgn_addendum/multiblock_qgn_addendum_reviewed.pdf)

### Exact Two-Block QGN Stiffness and Twist-Covariance Obstructions

A two-block QGN framework that separates exact blockwise stiffness, physical electromagnetic twisting, and composition-space mixing — with an exact one-pair reduction in each block and an exact middle-filling zero despite positive pair mobility in both. We then audit two proposed mechanisms for composition hopping and show both fail: once the active–remote bridge is itself Peierls twisted, its projected null row stays identically zero, so the apparent source was a gauge artifact. What survives is sharper — a fermionic monomial is twist invariant exactly when its net displacement vanishes, which permits position-balanced pair-hop vertices. We also reclassify our earlier twisted WSe₂ calculations: the geometric diagnostics remain useful screening tools, but they do not establish a controlled material realization.

📄 [Read the paper (PDF)](exact_two_block_qgn_stiffness/two_block_qgn_stiffness.pdf)

### Peierls-Transport Obstruction and the Minimal Many-Body Escape

We prove a no-go theorem: in exact multi-block projected-Hubbard and QGN parents, if the zero manifold contains the separate one-pair AGP states, then *every* number-conserving one-body null row preserves the blocks already at zero twist, so the composition-space curvature has no off-diagonal entries at all. Finite range or uniform exponential locality closes the last finite-torus loophole. We then construct the minimal escape — a strictly nearest-neighbor, translation-invariant many-body model that keeps every hypothesis of the theorem except one, and does produce genuine non-diagonal composition transfer.

📄 [Read the paper (PDF)](peierls_transport_obstruction_minimal_many_body_escape/peierls_transport_no_go_many_body_escape_reviewed.pdf)

### Microscopic Origins of Multiplicative Null Interactions in QGN Superconductors

The escape route from the no-go theorem is a multiplicative many-body row — a bond operator that transfers electrons between blocks *and* selects odd endpoint parity. On its face that looks engineered, so here we derive it from an ordinary multiorbital two-electron interaction. The key is an exact identity we call the all-filling current router: ordinary orbital-current path interference squares into the emergent parity projector, on the whole broken-pair source space and at every filling. From there, a four-orbital control sector with a gapped high block yields the desired term as an exact local Schur complement, and a closed-shell Kramers doubling reproduces the same coefficients. Assembling the bond gadgets into a periodic array gives two lattice theories. The matched screened family produces the composition-degenerate target at order λ⁶, with overlap, residual block coupling, and the first two Peierls derivatives of the remainder all pushed to λ⁸. The generic unscreened family is more delicate than we first reported: its order-λ⁴ charge operator is diagonal only as a *compression* to the product-AGP manifold, and the same order-four term leaks off that manifold. We give the exact composition-resolved action, singular values, norm, rank, and complete zero set of that leakage, for both the periodic lattice and the open path used by finite-system DMRG, plus an isolation gate showing the leakage is not suppressed by taking λ → 0 — so an autonomous composition Hamiltonian needs projection, exact cancellation, or order-four protection. Six finite-Fock-space verifiers audit the exact identities, the lattice scaling, and the leakage formulas; the exact scan passes 15/15 theorem-level checks over 760 boundary/composition sectors.

📄 [Read the paper (PDF)](microscopic_origins_multiplicative_qgn/paper/microscopic_origins_multiplicative_qgn.pdf)

### Exact Zero-Frequency Meissner Weight in a Gapless Frustration-Free QGN Parent

A finite-volume Peierls curvature is a Kohn response, not automatically a bulk Meissner weight — the thermodynamic, frequency, and transverse-momentum limits need not commute, and the earlier papers deliberately claimed no dynamical Meissner theorem. Here we close that gap. For any positive-square parent the entire finite-frequency response is exactly a regularized least-squares residual, and the normalized Abelian kernel matches the Kohn curvature precisely when the current-generated spectral measure is uniformly tight at zero energy. We then settle the transverse limit for an explicit local paired parent: on a connected graph, gauge-covariant site-swap rows turn the static subtraction into a graph Hodge decomposition, and because every lattice-transverse field is annihilated by the incidence matrix, its static and dynamical kernels coincide at every finite size, wave vector, and regulator. In d ≥ 2 all three weights collapse to the same positive value 2jρ(1−ρ). The state has exact pair off-diagonal long-range order while its fixed-number gap closes as L⁻², so no uniform many-body gap is used anywhere — a genuine zero-temperature Meissner theorem for a gapless parent.

📄 [Read the paper (PDF)](exact_zero_frequency_meissner_weight_qgn/exact_zero_frequency_meissner_weight_qgn.pdf)

### Microscopic Transfer of Exact Meissner Weight

The exact parent has a positive Meissner weight; the microscopic construction reproduces that parent only to order λ⁶, leaving a quasi-local λ⁸ remainder. Does the physics survive the transfer? Three results. First, target-metric deformations admit response bounds independent of the closing parent gap — but a bounded anisotropic metric still produces a finite Drude–Meissner mismatch, so row-relative form alone is not enough; adding hydrodynamic closure and point-group structure restores tightness and gives D_K = D_A = D_M. Second, a complete weighted-degree-eight audit yields a no-go: a two-pair witness shows the unmodified matched and Kramers families cannot preserve the full composition manifold at that order. An explicit finite-range higher-body counterterm cancels the unsafe source and restores the row ideal through degree eight. Third, the primitive bridge turns out to have an exact coherent dark channel, and local Schur overcompensation makes it the unique fixed-number ground state of a positive-semidefinite finite-coupling family — giving a rigorous transverse floor D₊ > 0. That is a finite-coupling Meissner lower bound for one selected branch, not an all-composition transfer.

📄 [Read the paper (PDF)](microscopic_transfer_exact_meissner_revised_release/paper/microscopic_transfer_exact_meissner.pdf)

### Weighted Row-Space Hodge Theorem for QGN Response Certificates

The response bounds in the papers above rest on a variational functional whose sharpness was never settled. Here we solve it. Among all target fields with the required divergence, an explicit weighted Hodge representative minimizes the target energy in Loewner order, and the resulting minimum is exactly the Schur complement — so the proposed operator-norm upper bound is not merely an upper bound, it is attained, with no duality gap and no uniform many-body gap needed. When the selected branch is a proper subspace of the zero manifold, the exact obstruction is the component that mixes it with the rest. Locality is a genuinely separate problem, and we treat it as such: a cluster-current sufficient theorem with exponential decay, a finite-fiber theorem, and — for finite-range translation-invariant factorization — a characterization as Laurent-module membership, which makes it decidable by Smith divisibility in 1D and by Gröbner saturation in any dimension. We do *not* claim a general positive theorem for exponentially local factorization through soft singular modes; analytic divisibility across the rank-drop set stays open. Three volume-independent certificates make the bound usable in practice.

📄 [Read the paper (PDF)](weighted_row_space_hodge_qgn/weighted_row_space_hodge_qgn.pdf)

### Quenched Pair-Coupling Disorder in Frustration-Free QGN Parents I

Everything so far assumed spatially uniform pair exchange. It turns out not to be required: on any finite connected graph with bounded conductances, the fixed-pair-number Dicke/AGP state is still the unique zero state, and its Peierls curvature is exactly the weighted network cell energy — so the many-body stiffness inherits sharp deterministic ellipticity bounds. Disorder does change the finite-size source condition, and a uniformly elliptic modulated ring shows ellipticity alone is *not* enough to make the limits commute. For a stationary ergodic conductance field, the curvature density converges almost surely to a homogenized tensor, and cell-energy convergence plus resolvent locality give almost-sure source tightness — hence equality between the thermodynamic Kohn curvature and the Abel-regularized zero-frequency response. On the i.i.d. torus we get explicit dynamical-defect rates, an exact subgap cutoff, and a diffusive diagonal limit.

📄 [Read the paper (PDF)](quenched_pair_coupling_disorder_in_frustration_free_qgn_parents_I/quenched_pair_coupling_disorder_in_frustration_free_qgn_parents_I.pdf)

### Quenched Pair-Coupling Disorder in Frustration-Free QGN Parents II

*Companion to Paper I.* Paper I's disorder results are annealed; this one asks what happens sample by sample. The obstacle is that the sharp spectral projector is discontinuous in each bond variable, so we work with a positive Abel envelope instead, derive an exact bond-sensitivity identity, and get a variance bound that yields a one-sided quenched bound on the sharp tail. We are explicit that this is not yet a mesoscopic concentration theorem. In the frozen-spectrum Born sector the sharp tail *is* an exact quadratic form and does concentrate. The rest of the paper is an honest audit of three routes to a nonlinear fixed-cutoff theorem and why each falls short: a modulated cycle shows eigenvalue counts don't control current weight, positive heat-mixture thresholds can't be relatively narrow, and log-shell averaging removes exceptional cutoffs only at the uncentered floor. Along the way, rank-one interlacing proves a single bond crosses a given cutoff at most once, and a coarea identity pins the exact current-weighted crossing surface. What remains is isolated as one clean open problem.

📄 [Read the paper (PDF)](quenched_pair_coupling_disorder_in_frustration_free_qgn_parents_II/quenched_pair_coupling_disorder_in_frustration_free_qgn_parents_II.pdf)