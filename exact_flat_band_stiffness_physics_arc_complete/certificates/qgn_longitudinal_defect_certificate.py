#!/usr/bin/env python3
"""Certificate for the Hermitian-QGN longitudinal/AGP decomposition.

For a representative fixed-local connected Hermitian-QGN model, this script
checks across every filling (M=5,6,7) that

    E_n'' - rho_n E_1'' = (n^2-rho_n) Gamma_M,
    rho_n = n(M-n)/(M-1),

where Gamma_M is computed independently from the ground-multiplet (trace)
component of the first twist derivative of the local QGN square factors.
It also prints the symbolic AGP norm ratio rho_n.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import sympy as sp

# The helper module and input records ship alongside this script
# (extracted from legacy_packages/qgn_connected_irreducible_search_package.zip).
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT))
import qgn_search_v3 as q

RECORDS = ROOT / "search_connected_candidates_m6.json"
OUT = ROOT / "qgn_no_counterexample_evidence.csv"


def build_case():
    rec = json.loads(RECORDS.read_text())[0]
    base = q.standard_channels()
    names = (
        "P0", "P1", "I", "sx", "sy", "sz", "bI1c", "bI1s",
        "bsx1c", "bsx1s", "bsz1c", "bI2c",
    )
    pool = [base[k] for k in names]
    channels = [q.combine_channels(pool, np.asarray(rec["coeff1"]), "representative")]
    return rec, channels


def trace_gamma(M: int, model, params: dict, channels, h: float = 1e-5) -> float:
    """One-pair longitudinal curvature Gamma_M from trace derivatives.

    The Hamiltonian convention in qgn_search_v3 is H=(1/2) sum w S_R^2, so
    the second-derivative contribution is sum w |<1|dot S_R|1>|^2.
    """
    k = q.grid_1d(M)
    gamma = 0.0
    for channel in channels:
        samples = []
        for A in (-2*h, -h, h, 2*h):
            Avec = np.array([A])
            up, down = q.spinors(model, k, Avec, params)
            ups = q._projected_components(k, up, 1, channel.components, Avec)
            dns = q._projected_components(
                k, down, 1, q.qgn_down_components(channel.components), Avec
            )
            # The normalized one-pair AGP has one-particle density I/M per spin.
            samples.append([
                (np.trace(Su) + np.trace(Sd)) / M for Su, Sd in zip(ups, dns)
            ])
        for R in range(M):
            fm2, fm1, fp1, fp2 = (samples[j][R] for j in range(4))
            deriv = (fm2 - 8*fm1 + 8*fp1 - fp2) / (12*h)
            gamma += float(channel.weight) * abs(deriv)**2
    return float(gamma.real)


def symbolic_ratio() -> sp.Expr:
    n, L = sp.symbols("n L", integer=True, positive=True)
    # ||(eta+)^n |0>||^2 = n! L!/(L-n)!
    Nn = sp.factorial(n) * sp.factorial(L) / sp.factorial(L-n)
    # Lowest-weight S=L/2-1 source raised n-1 times.
    Cn = sp.factorial(n-1) * sp.factorial(L-2) / sp.factorial(L-n-1)
    ratio = sp.simplify(n**2 * Cn * L / Nn)
    return sp.factor(ratio)


def main() -> None:
    ratio = symbolic_ratio()
    expected = sp.symbols("n") * (sp.symbols("L") - sp.symbols("n")) / (sp.symbols("L") - 1)
    print("AGP source norm ratio:", ratio)

    rec, channels = build_case()
    rows_out = []
    worst = 0.0
    for M in (5, 6, 7):
        gamma = trace_gamma(M, q.model_phase_1d, rec["params"], channels)
        rows = q.evaluate_case(
            "representative", M, q.model_phase_1d, rec["params"], channels,
            range(1, M), h=0.002,
        )
        c1 = rows[0]["curvature"]
        print(f"M={M}: Gamma_M(trace)={gamma:.15g}")
        for row in rows[1:]:
            n = row["n"]
            rho = row["rho"]
            defect = row["defect"]
            predicted = (n*n-rho)*gamma
            residual = defect-predicted
            worst = max(worst, abs(residual))
            print(
                f"  n={n}: defect={defect:.15g}, predicted={predicted:.15g}, "
                f"residual={residual:.3e}"
            )
            rows_out.append({
                "M": M,
                "n": n,
                "rho_n": rho,
                "E1pp": c1,
                "Enpp": row["curvature"],
                "ratio_canonical": row["ratio"],
                "defect": defect,
                "Gamma_trace": gamma,
                "predicted_defect": predicted,
                "residual": residual,
                "gap": row["gap"],
                "degeneracy": row["deg"],
            })

    import csv
    with OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows_out[0]))
        writer.writeheader(); writer.writerows(rows_out)

    # Finite-difference errors dominate beyond this threshold.
    assert worst < 2e-8, worst
    print("worst residual:", worst)
    print("PASS")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
