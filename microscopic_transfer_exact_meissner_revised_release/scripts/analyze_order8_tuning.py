#!/usr/bin/env python3
"""Test whether a scalar bridge-amplitude tuning can remove the degree-8 source."""
from __future__ import annotations
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
import scipy
import scipy.linalg as la
import scipy.optimize as opt

MODULE = Path(__file__).with_name("audit_order8_all_fillings.py")
spec = importlib.util.spec_from_file_location("order8_all", MODULE)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
assert spec.loader is not None
spec.loader.exec_module(mod)


def raw_sources(L, n, ctrl, dm, s, a):
    active = mod.active_sector(L, 2*n)
    # b=1: A4 is U, B2+C is V in R8 Z = b^2 U + b^4 V.
    C2, A, Bq, Cq, r0 = mod.build_operators(active, ctrl, 1.0, s, a, dm)
    Z = active.Z
    U = -mod.word_action(C2, (A, A, A, A), r0, Z)
    Vb = -mod.word_action(C2, (Bq, Bq), r0, Z)
    Vc = +mod.word_action(C2, (Cq,), r0, Z)
    V = Vb + Vc
    mixed = sum((mod.word_action(C2, w, r0, Z) for w in ((A,A,Bq),(A,Bq,A),(Bq,A,A))), start=np.zeros_like(Z))
    return active, U, V, mixed


def norms(active, S):
    Z = active.Z
    C = mod.herm(Z.conj().T @ S)
    X = S - Z @ C
    return {
        "source_op": mod.opnorm(S),
        "source_fro": float(la.norm(S, "fro")),
        "compression_op": mod.opnorm(C),
        "cross_op": mod.opnorm(X),
        "compression": C,
        "compression_eigs": la.eigvalsh(C),
    }


def serial(x):
    if isinstance(x, np.ndarray):
        if x.ndim == 1:
            return [float(v) for v in x]
        return [[[float(z.real), float(z.imag)] for z in row] for row in x]
    if isinstance(x, (np.floating, np.integer)):
        return x.item()
    raise TypeError(type(x).__name__)


def main():
    L, dm, ds, da = 3, 10.0, 8.0, 8.0
    s, a = 0.9, 0.45
    ctrl = mod.control_space(L, dm, ds, da)
    out = []
    print("SCALAR-TUNING OBSTRUCTION FOR THE COMPLETE DEGREE-EIGHT SOURCE")
    print("="*78)
    print(f"L={L}, dm={dm}, ds={ds}, da={da}, s={s}, a={a}")
    print("R8 P_Z = b^2 U + b^4 V = b^2(U+xV), x=b^2")
    print()

    for n in (1,2,3):
        active, U, V, mixed = raw_sources(L,n,ctrl,dm,s,a)
        uf = float(la.norm(U,"fro")); vf=float(la.norm(V,"fro"))
        ip = float(np.real(np.vdot(U,V)))
        cosang = ip/(uf*vf) if uf*vf else 0.0
        x_fro = max(0.0, -ip/(vf*vf)) if vf else 0.0
        f = lambda x: mod.opnorm(U+x*V)
        res = opt.minimize_scalar(f, bounds=(0.0, 1.0), method="bounded", options={"xatol":1e-14})
        x_op=float(res.x)
        data = {
            "pairs":n,
            "compositions":active.compositions,
            "frobenius_cosine":cosang,
            "best_x_fro":x_fro,
            "best_b_fro":math.sqrt(x_fro),
            "best_x_op":x_op,
            "best_b_op":math.sqrt(x_op),
            "at_fro":norms(active,U+x_fro*V),
            "at_op":norms(active,U+x_op*V),
            "U":norms(active,U),
            "V":norms(active,V),
            "mixed_norm":mod.opnorm(mixed),
        }
        out.append(data)
        print(f"n={n}, compositions={active.compositions}")
        print(f"  cos_F(U,V)                    = {cosang:+.12f}")
        print(f"  best b (Frobenius)            = {math.sqrt(x_fro):.12e}")
        print(f"  min ||U+b^2V||_F              = {data['at_fro']['source_fro']:.12e}")
        print(f"  residual / ||U||_F            = {data['at_fro']['source_fro']/uf:.12e}")
        print(f"  best b (operator source norm) = {math.sqrt(x_op):.12e}")
        print(f"  min source op norm            = {data['at_op']['source_op']:.12e}")
        print(f"  compression op norm there     = {data['at_op']['compression_op']:.12e}")
        print(f"  cross-source op norm there    = {data['at_op']['cross_op']:.12e}")
        print(f"  mixed AAB norm                = {data['mixed_norm']:.12e}")
        print()

    x1 = 4*s**4/ds**2
    print(f"Analytic one-pair cancellation: b*=sqrt(4s^4/ds^2)={math.sqrt(x1):.12e}")
    exact_one = abs(out[0]["best_x_fro"]-x1) < 1e-12 and out[0]["at_fro"]["source_fro"] < 1e-12
    higher_obstruction = all(d["at_op"]["source_op"] > 1e-6 for d in out[1:])
    print("one-pair formula verified:", "YES" if exact_one else "NO")
    print("higher-filling scalar-tuning obstruction:", "YES" if higher_obstruction else "NO")
    print("OVERALL:", "PASS" if exact_one and higher_obstruction else "FAIL")

    Path(__file__).with_name("order8_tuning_analysis.json").write_text(json.dumps(out,default=serial,indent=2)+"\n")
    return 0 if exact_one and higher_obstruction else 1

if __name__ == "__main__":
    raise SystemExit(main())
