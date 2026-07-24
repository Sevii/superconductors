---
title: "Prior-Art and Scope Audit for the QGN Stiffness Reduction Paper"
author: "Nicholas Sledgianowski"
date: "July 23, 2026"
geometry: margin=0.9in
fontsize: 10.5pt
---

# 1. Purpose and search scope

This audit accompanies *Exact Flat-Band Stiffness from Pair Mobility in a QGN Model*. It records the literature and theorem-boundary checks completed before submission. The sweep was targeted rather than bibliographically exhaustive. It searched for:

- an exact finite-density equality between many-body flat-connection stiffness and two-particle pair mass;
- exact paired ground states and pseudospin reductions in flat-band Hubbard models;
- exact pair spectra and minimal-quantum-metric formulas;
- uniqueness and irreducibility criteria for flat-band many-body ground states;
- non-Hermitian extensions of the quantum-geometric-nesting construction;
- spectral-gap results analogous to the deferred rank-one Aldous problem.

Searches used arXiv, APS journal pages, publisher records, and title/DOI queries through July 23, 2026.

# 2. Closest prior results and the remaining novelty

| Work | Established result relevant here | What it does not establish |
|---|---|---|
| Tovmasyan *et al.*, PRB **94**, 245149 (2016) | Exact BCS/AGP ground state for the projected attractive-Hubbard model under uniform pairing; emergent pseudospin $SU(2)$. | No exact finite-density stiffness-equals-pair-mass theorem of the form proved for Model II here. |
| Huhtinen *et al.*, PRB **106**, 014518 (2022) | Exact Cooper-pair mass in the UPC Hubbard setting and its relation to the minimal quantum metric; clarifies embedding dependence. | Does not reduce the exact many-body twist Hessian at arbitrary filling to that pair mass. |
| Han *et al.*, PRX **14**, 041004 (2024) | Defines QGN, constructs Hermitian positive-kernel ideal Hamiltonians, exact ordered ground states, and solvable parts of the excitation spectrum. | Does not state the exact finite-size canonical reduction theorem or prove the GHK equality for Model II. |
| Herzog-Arbeitman *et al.*, *Commun. Phys.* (2026), DOI 10.1038/s42005-026-02732-2 | Exact eta-pairing ground states and analytic Cooper-pair bound-state spectra in a broad family; quadratic pair dispersion tied to the minimal metric. | Does not prove that the full finite-density ground-energy curvature is exactly the hard-core filling factor times the pair mass. |
| Gao, Han, Khalaf, PRL **136**, 076503 (2026) | Rigorous RDM-bootstrap lower bounds; variational pair-mass upper bound; numerical saturation; explicitly conjectures exact equality and leaves proof open. | This is the conjecture closed here for projected-Hubbard Model II. |
| Mielke, PLA **174**, 443 (1993); J. Phys. A **32**, 8411 (1999) | Necessary-and-sufficient uniqueness/irreducibility criteria in flat-band ferromagnetic Hubbard problems. | Different particle channel and Hamiltonian; provides structural precedent, not the superconducting commutant theorem used here. |

The targeted sweep found no earlier proof of the following combined statement:

$$
\kappa_{ij}^{(n)}(V)
=
\frac{V}{V-1}\nu(1-\nu)
(m_{\mathrm{pair}}^{-1})_{ij}
$$

for every nonsingular filling of Gao-Han-Khalaf Model II at finite volume, together with its thermodynamic equality in their flat-connection definition of $D_s$. This is a statement about the result of the targeted sweep, not a claim that no unpublished or differently phrased result exists.

# 3. Exact relation to Gao-Han-Khalaf

The GHK paper makes four points that define the target precisely:

1. it obtains rigorous lower bounds on stiffness in frustration-free superconductors;
2. for the QGN examples tested, those lower bounds numerically meet a rigorous variational upper bound;
3. the upper bound is expressed through the inverse pair mass evaluated at the same finite size;
4. the authors conjecture equality for QGN models and explicitly leave a rigorous proof to future work.

The consolidated paper proves that equality for the projected PSD-Hubbard realization of Model II in the connected range $0<\xi<\pi/2$, with the exact finite-size correction $V/(V-1)$. It does not claim to prove the conjecture for every QGN interaction.

# 4. Non-Hermitian QGN-boundary adjudication

This boundary must be phrased carefully because “QGN” describes the flat-band electronic geometry, while the reduction theorem also assumes a particular interaction structure.

## 4.1 Standard ideal QGN Hamiltonians

Section II of Han *et al.* constructs ideal interactions using:

- a positive-semidefinite kernel;
- an infinite family of Hermitian local operators;
- projected commutation with the order parameter.

That is the Hermitian positive-square class used by the least-squares/intertwiner proof.

## 4.2 Published non-Hermitian extension

Section VI (“Possible Extensions”) of the same work explicitly presents an SSH-chain construction in which an interaction coefficient/form is non-Hermitian. The authors describe it as a generalization of the standard form, retain exact zero-energy ground states, and state that the excitations are not exactly solvable.

Therefore:

- the underlying band can still have perfect QGN;
- the interaction need not belong to the Hermitian positive-square class;
- perfect QGN alone is not equivalent to the hypotheses of the restricted reduction theorem.

## 4.3 Status of the computational stress test

The archived search package contains a random non-Hermitian generalized-nesting stress test with a large violation of the canonical law. It establishes only that the Hermitian-square hypothesis is substantive. It is **not** a reproduction of the published SSH extension and is **not** presented as a counterexample to that model.

The equality for published or future non-Hermitian QGN extensions remains open. The consolidated paper states this boundary explicitly and makes no broader claim.

# 5. Alon-Puder verification and the deferred gap note

The reference is real and current:

> Gil Alon and Doron Puder, “Aldous-type Spectral Gaps in Unitary Groups,” arXiv:2603.00353, submitted February 27, 2026.

The paper formulates unitary-group analogues of Aldous-type spectral-gap phenomena and proves several nontrivial cases. It is related in spirit to the observed filling-independent gaps in rank-one reflection Laplacians, but it is not presently a theorem that can be imported into the QGN proof:

- Alon-Puder study hypergraph-induced probability measures and associated $U(n)$ representation Laplacians.
- The QGN sector-gap operator is a weighted sum of conjugations by specific rank-one Householder reflections.
- No identification has been established that places the latter inside one of the cases proved by Alon-Puder.
- The stronger equality $\Delta_{L,n}=\Delta_{L,1}$ is not load-bearing for the stiffness theorem; the proved comparison $\Delta_{L,n}\leq\Delta_{L,1}$ is sufficient to show that a uniform global gap cannot be assumed.

For these reasons, the “rank-one Aldous” problem is spun out into a later standalone mathematical note and omitted from the physics arc paper except for a brief statement that the stronger gap-ordering question is deferred.

# 6. Embedding, lift, and finite-size prior art

Huhtinen *et al.* show that conventional geometric quantities and superfluid formulas can depend on orbital embedding and motivate the minimal quantum metric. The winding counterexample in the present project is parallel in spirit but distinct:

- embedding fixes how microscopic orbitals are represented inside a unit cell;
- the finite-torus displacement lift fixes which real-space winding representative is used when inserting a continuous flat connection;
- the new obstruction is an exact roots-of-unity aliasing mechanism at displacement equal to the circumference.

The consolidated paper cites the minimal-metric literature and treats the displacement lift as part of the finite-size model definition.

# 7. Submission claims authorized by this audit

The following claims are supported by the present proof and literature sweep:

1. **Exact restricted theorem:** the Hermitian positive-square QGN class satisfying the stated full-rank pairing, analytic-branch, and locality hypotheses obeys the finite-size reduction with nonnegative longitudinal tensor.
2. **Hubbard uniqueness theorem:** projected-frame connectedness is necessary and sufficient for a unique AGP zero mode at every nonsingular filling.
3. **Model-II finite-size theorem:** for $\xi\neq0$, $N_x,N_y\ge3$, and any finite torus whose projected-frame graph is connected, the exact formula holds.
4. **Model-II all-size theorem:** $0<\xi<\pi/2$ guarantees the required connectivity on all rectangular tori in the stated sequence.
5. **Model-II thermodynamic GHK equality:** in GHK’s literal flat-connection-curvature definition of $D_s$, the equality follows directly from the exact finite-size identity.
6. **Exact absence of the Model-II defect:** $\Gamma=0$ on every finite torus; this is not merely asymptotic.
7. **Scope boundary:** no claim is made for general non-Hermitian QGN extensions or for a stricter dynamical order-of-limits definition of superfluid stiffness.

# 8. Claims deliberately withheld

The paper does not claim:

- proof of the GHK equality for all QGN interactions;
- proof for disconnected or degenerate paired ground spaces;
- proof for several inequivalent pairing singular-value blocks;
- a uniform-in-volume many-body excitation gap;
- equality of all filling-sector gaps with the one-pair gap;
- adjudication of the published non-Hermitian SSH extension;
- equivalence between flat-connection curvature and every dynamical transport definition.

# 9. Bibliographic records checked

1. Q. Gao, Z. Han, E. Khalaf, PRL **136**, 076503 (2026), DOI 10.1103/gw85-5r92; arXiv:2506.18969v4.
2. Z. Han, J. Herzog-Arbeitman, B. A. Bernevig, S. A. Kivelson, PRX **14**, 041004 (2024); arXiv:2401.04163v3.
3. M. Tovmasyan, S. Peotta, P. Törmä, S. D. Huber, PRB **94**, 245149 (2016), DOI 10.1103/PhysRevB.94.245149.
4. K.-E. Huhtinen *et al.*, PRB **106**, 014518 (2022), DOI 10.1103/PhysRevB.106.014518.
5. J. Herzog-Arbeitman *et al.*, *Communications Physics* (2026), DOI 10.1038/s42005-026-02732-2.
6. A. Mielke, *Physics Letters A* **174**, 443-448 (1993), DOI 10.1016/0375-9601(93)90207-G.
7. A. Mielke, *Journal of Physics A* **32**, 8411-8418 (1999), DOI 10.1088/0305-4470/32/48/304.
8. G. Alon, D. Puder, “Aldous-type Spectral Gaps in Unitary Groups,” arXiv:2603.00353 (2026).

