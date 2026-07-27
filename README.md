
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

The escape route from the no-go theorem is a multiplicative many-body row — a bond operator that transfers electrons between blocks *and* selects odd endpoint parity. On its face that looks engineered, so here we derive it from an ordinary multiorbital two-electron interaction. The key is an exact identity we call the all-filling current router: ordinary orbital-current path interference squares into the emergent parity projector, on the whole broken-pair source space and at every filling. From there, a four-orbital control sector with a gapped high block yields the desired term as an exact local Schur complement, and a closed-shell Kramers doubling reproduces the same coefficients. Assembling the bond gadgets into a periodic array gives two lattice theories: a screened family that produces the target at order λ⁶ with overlap corrections pushed to λ⁸, and the generic unscreened family, which keeps a diagonal charge selector at order λ⁴ before the composition mixing appears. Five finite-Fock-space verifiers audit the exact identities and the predicted scaling.

📄 [Read the paper (PDF)](microscopic_origins_multiplicative_qgn/microscopic_origins_multiplicative_qgn_draft.pdf)

### Exact Zero-Frequency Meissner Weight in a Gapless Frustration-Free QGN Parent

A finite-volume Peierls curvature is a Kohn response, not automatically a bulk Meissner weight — the thermodynamic, frequency, and transverse-momentum limits need not commute, and the earlier papers deliberately claimed no dynamical Meissner theorem. Here we close that gap. For any positive-square parent the entire finite-frequency response is exactly a regularized least-squares residual, and the normalized Abelian kernel matches the Kohn curvature precisely when the current-generated spectral measure is uniformly tight at zero energy. We then settle the transverse limit for an explicit local paired parent: on a connected graph, gauge-covariant site-swap rows turn the static subtraction into a graph Hodge decomposition, and because every lattice-transverse field is annihilated by the incidence matrix, its static and dynamical kernels coincide at every finite size, wave vector, and regulator. In d ≥ 2 all three weights collapse to the same positive value 2jρ(1−ρ). The state has exact pair off-diagonal long-range order while its fixed-number gap closes as L⁻², so no uniform many-body gap is used anywhere — a genuine zero-temperature Meissner theorem for a gapless parent.

📄 [Read the paper (PDF)](exact_zero_frequency_meissner_weight_qgn/exact_zero_frequency_meissner_weight_qgn.pdf)

### Microscopic Transfer of Exact Meissner Weight

The exact parent has a positive Meissner weight; the microscopic construction reproduces that parent only to order λ⁶, leaving a quasi-local λ⁸ remainder. Does the physics survive the transfer? Three results. First, target-metric deformations admit response bounds independent of the closing parent gap — but a bounded anisotropic metric still produces a finite Drude–Meissner mismatch, so row-relative form alone is not enough; adding hydrodynamic closure and point-group structure restores tightness and gives D_K = D_A = D_M. Second, a complete weighted-degree-eight audit yields a no-go: a two-pair witness shows the unmodified matched and Kramers families cannot preserve the full composition manifold at that order. An explicit finite-range higher-body counterterm cancels the unsafe source and restores the row ideal through degree eight. Third, the primitive bridge turns out to have an exact coherent dark channel, and local Schur overcompensation makes it the unique fixed-number ground state of a positive-semidefinite finite-coupling family — giving a rigorous transverse floor D₊ > 0. That is a finite-coupling Meissner lower bound for one selected branch, not an all-composition transfer.

📄 [Read the paper (PDF)](microscopic_transfer_exact_meissner_revised_release/paper/microscopic_transfer_exact_meissner.pdf)