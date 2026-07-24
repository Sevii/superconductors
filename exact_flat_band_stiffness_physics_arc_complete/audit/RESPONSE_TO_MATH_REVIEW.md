# Response to mathematical review of “Exact Flat-Band Stiffness from Pair Mobility in a QGN Model”

**Revision date:** July 23, 2026

The reviewer found the central mathematics and all numerical certificates correct. The revision addresses every finding as follows.

## F1 — all-size Model-II connectivity

Added a self-contained proof. For the one-dimensional finite-grid coefficients

$$
c_r^{(N)}(\xi)=N^{-1}\sum_m e^{i\xi\cos k_m}e^{-irk_m},
$$

$\operatorname{Re}c_0>0$ follows termwise for $0<\xi<\pi/2$, while inversion symmetry gives

$$
\operatorname{Im}c_1=N^{-1}\sum_m \cos k_m\sin(\xi\cos k_m)>0.
$$

The two-dimensional off-diagonal projector coefficient factorizes, so the coefficients at $(0,0)$, $(1,0)$, and $(0,1)$ are nonzero. These edges make the two-layer frame graph connected on every rectangular torus. The unused all-size Model-I assertion was withdrawn; only its finite audit remains.

## F2 — non-load-bearing sector-gap statements

Removed the boxed sector-gap equations and the compressed representation-theoretic assertions from the physics paper. Section 12 now states only the finite-volume gap input actually used and defers the rank-one Aldous and gap-scaling questions to a standalone mathematical note.

## F3 — H2 wording

Changed the sentence to: “H2 follows from H1 because the factors annihilate the vacuum.”

## F4 — notation

Added a notation paragraph, renamed the one-particle vector space $\mathcal V$, and renamed torus side lengths $\mathbf N=(N_1,\ldots,N_d)$, leaving $L$ exclusively for the number of pair orbitals. The filling is defined globally as $\nu=n/L$, with $L=V$ only in Model II. The winding example now uses $g>0$ for the attraction magnitude, separating it from the PSD-square and Model-II coupling conventions.

## F5 — equation numbering

Moved the Hessian definition (10) before (10a)-(10j), replaced the headline sublabels (1a)-(1d) by (T1)-(T4), and removed the obsolete (61a)-(61c) sequence with the deferred sector-gap material.

## F6 — presentation

Added the one-line derivation of the winding factorization; corrected the abstract and conclusion so connectivity is identified as an H4 criterion only in the projected-Hubbard subclass; added the energy ordering that makes Lemma 4’s branch the lowest one; and completed the Kato citation. Reference 11 is cited as an online publication because the publisher’s current early-access record does not yet display a volume/article number.

No theorem coefficient or load-bearing proof step changed in this revision.
