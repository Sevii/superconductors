# Exact unscreened product-AGP leakage certificate, version 2.1

This package completes the two requested extensions of the unscreened leakage analysis:

1. an explicit composition-resolved action, singular-value decomposition, Gram operator, general norm, and rank formula for the periodic charge leakage;
2. a separate exact derivation for the standard open path used by the DMRG/MPS plan.

The periodic result is the theorem-level certificate for the microscopic paper. The open result is a numerical-boundary certificate: it is required when the wrap bond is removed before interpreting open-chain spectra.

## Reproduce

```bash
bash run_certificate.sh
```

The verifier uses only the Python standard library. Every PASS/FAIL decision is made with integers and `fractions.Fraction`; decimal values are display-only.
The centered-action PASS lines check exact first and second moments; the displayed operator SVD is then an analytic consequence of occupation-basis diagonality and block-charge orthogonality, not a floating-point matrix decomposition.

## Central formulas

For either boundary, if `c=(n1,n2)` and

\[
R=(1-P_{\mathcal Z_n})\mathcal C P_{\mathcal Z_n},
\qquad
\mu_c=\langle c|\mathcal C|c\rangle,
\]

then

\[
R|c\rangle=(\mathcal C-\mu_c)|c\rangle.
\]

Writing \(\tau_c^2=\operatorname{Var}_c(\mathcal C)\) and normalizing the nonzero centered vectors gives

\[
R=\sum_{\tau_c>0}\tau_c|\ell_c\rangle\langle c|,
\qquad
R^\dagger R=\sum_c\tau_c^2|c\rangle\langle c|,
\]

\[
\|R\|=\max_c\tau_c,
\qquad
\operatorname{rank}R=\#\{c:\tau_c>0\}.
\]

On the periodic cycle, \(\tau_c=A_4\sigma_L(c)\) with the closed formula in `DERIVATION.md`.

On the open path,

\[
\mathcal C_{4,L}^{\rm op}=(8A_4-2U_4)n+8A_4Z-V_4B,
\qquad V_4=4A_4-U_4,
\]

and

\[
(\tau_c^{\rm op})^2
=64A_4^2\Gamma_0-16A_4V_4\Gamma_1+V_4^2\Gamma_2.
\]

The endpoint term \(-V_4B\) is why the open formula is not obtained by replacing \(L\) with \(L-1\) in the periodic result.

## Scope decision

- The parent microscopic theorem is periodic, so the periodic proposition closes the paper-level leakage question.
- The numerical plan also uses open chains. Under the standard truncation to \(L-1\) path bonds, open spectra require the separate formula in this package.
- If a future implementation adds edge counterterms or retains a special wrap-shell term, that defines another boundary operator and must be checked separately.

See `BOUNDARY_SCOPE_DECISION.md` for the full reasoning.

## Files

- `check_unscreened_agp_leakage.py` — exact cycle/path exhaustive verifier.
- `unscreened_agp_leakage_certificate.txt` — full exact certificate and focused matrices.
- `unscreened_agp_leakage_periodic_all_sectors.csv` — every periodic composition for \(L=2,\ldots,9\).
- `unscreened_agp_leakage_periodic_phase3_sectors.csv` — periodic \(L=4,6,8\), \(n=2,3\).
- `unscreened_agp_leakage_open_all_sectors.csv` — every open composition for \(L=2,\ldots,9\).
- `unscreened_agp_leakage_open_phase3_sectors.csv` — open \(L=4,6,8\), \(n=2,3\).
- `unscreened_agp_leakage_all_sectors.csv` and `unscreened_agp_leakage_phase3_sectors.csv` — backward-compatible periodic aliases.
- `DERIVATION.md` — full analytic derivation.
- `BOUNDARY_SCOPE_DECISION.md` — precise T2.3 scope decision.
- `T2_1_T2_3_COMPLETION.md` — task-completion record.
- `RESULTS.md` — principal exact results and implications.
- `paper_insert.tex` — manuscript-ready periodic proposition and open corollary.
- `microscopic_origins_multiplicative_qgn_leakage_updated.tex` — full uploaded manuscript with the leakage corrections integrated.
- `microscopic_origins_multiplicative_qgn_leakage.patch` — unified diff against the uploaded draft.
- `MANUSCRIPT_PATCH_NOTES.md` — recommended text changes.
- `PLATFORM_STABLE_CERTIFICATE_SUMMARY.txt` — portable exact summary.
- `BUILD_AND_TEST_RESULTS.txt` — retained validation runs.
- `ENVIRONMENT.txt`, `MANIFEST.txt`, `SHA256SUMS` — reproducibility metadata.

## Principal conclusion

The projected selector formulas are exact for both boundaries, but neither generally implies invariance of the product-AGP manifold. On the cycle there are a few saturation zeros; on the physical open path, for \(L\ge3\), leakage vanishes only when each block is completely empty or completely full. At the project default \(V_4/A_4=1\), every planned open Phase-3 composition at \(L=4,6,8\), \(n=2,3\) leaks at order \(\lambda^4\).
