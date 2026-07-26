#!/usr/bin/env python3
"""Step-size convergence for gauge-invariant singular-projector geometry."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from interaction_singular_subbundle_scan_revised import (
    ARCHIVED_CANDIDATES,
    coarse_metrics,
    fine_geometry,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("projector_geometry_convergence.json"))
    parser.add_argument("--shell", type=int, default=3)
    parser.add_argument("--nk", type=int, default=7)
    parser.add_argument("--steps", type=float, nargs="+", default=[1.0e-3, 2.0e-3, 4.0e-3])
    args = parser.parse_args()

    payload: dict[str, object] = {
        "method": "central finite differences of rank-one projectors",
        "delta_frac_values": args.steps,
        "shell": args.shell,
        "nk": args.nk,
        "candidates": {},
    }

    for name, archived in ARCHIVED_CANDIDATES.items():
        params = {key: archived[key] for key in ("theta", "V", "psi", "w", "mstar", "vz")}
        eta = float(archived["eta"])
        density_asym = float(archived["density_asym"])
        metrics = coarse_metrics(params, nk=9, eta=eta, density_asym=density_asym, shell=args.shell)
        active_pair = metrics["active_pair"]
        rows = []
        for step in args.steps:
            geom = fine_geometry(
                params,
                eta,
                density_asym,
                active_pair,
                shell=args.shell,
                nk=args.nk,
                delta_frac=float(step),
            )
            rows.append(
                {
                    "delta_frac": float(step),
                    "Gamma_proxy_sum_aM2": np.asarray(geom["Gamma_proxy_sum_aM2"]).tolist(),
                    "ideal_gamma_sum_aM2": np.asarray(geom["ideal_gamma_sum_aM2"]).tolist(),
                    "metric_active_aM2": np.asarray(geom["metric_aM2"][:2]).tolist(),
                    "chern_approx": np.asarray(geom["chern_approx"]).tolist(),
                }
            )
        reference = np.asarray(rows[1]["Gamma_proxy_sum_aM2"], dtype=float)
        spread = max(
            float(np.max(np.abs(np.asarray(row["Gamma_proxy_sum_aM2"], dtype=float) - reference)))
            for row in rows
        )
        payload["candidates"][name] = {
            "active_pair": list(active_pair),
            "singular_gap_min": float(metrics["singular_gap_min"]),
            "rows": rows,
            "maximum_Gamma_spread_from_middle_step": spread,
        }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(args.output)
    for name, item in payload["candidates"].items():
        print(name, item["maximum_Gamma_spread_from_middle_step"])


if __name__ == "__main__":
    main()
