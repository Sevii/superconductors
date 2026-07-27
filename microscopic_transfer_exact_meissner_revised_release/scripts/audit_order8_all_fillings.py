#!/usr/bin/env python3
"""All-filling degree-eight Feshbach audit for the screened microscopic ring.

The calculation works directly on the product-AGP zero manifold, sector by
sector, and therefore avoids constructing the full degree-eight active
operator.  It retains every single-flavor weighted-degree-eight resolvent
class of the unmodified screened model:

  * V1^4              (two outer bridges and four routers),
  * V2,Q^2            (four bridge toggles),
  * one W4 insertion  (matched B^2 feedback),
  * the A,A,B mixed words (verified to vanish by control parity).

The output separates ground-manifold compression from the orthogonal state
source.  A nonzero compression is a genuine order-eight branch selector even
when the cross source vanishes.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
import scipy
import scipy.linalg as la
import scipy.sparse as sp


def annihilate(state: int, i: int):
    if not ((state >> i) & 1):
        return None
    sign = -1 if (state & ((1 << i) - 1)).bit_count() & 1 else 1
    return state ^ (1 << i), sign


def create(state: int, i: int):
    if (state >> i) & 1:
        return None
    sign = -1 if (state & ((1 << i) - 1)).bit_count() & 1 else 1
    return state | (1 << i), sign


def fixed_basis(nmodes: int, npart: int) -> list[int]:
    return [sum(1 << i for i in occ) for occ in itertools.combinations(range(nmodes), npart)]


def cdagc(basis: Sequence[int], dst: int, src: int) -> sp.csr_matrix:
    index = {s: k for k, s in enumerate(basis)}
    rows: list[int] = []
    cols: list[int] = []
    vals: list[complex] = []
    for col, state in enumerate(basis):
        x = annihilate(state, src)
        if x is None:
            continue
        state1, s1 = x
        x = create(state1, dst)
        if x is None:
            continue
        state2, s2 = x
        rows.append(index[state2])
        cols.append(col)
        vals.append(complex(s1 * s2))
    n = len(basis)
    return sp.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=complex)


def mode(L: int, block: int, cell: int, spin: int) -> int:
    return 2 * (block * L + cell) + spin


@dataclass
class Active:
    basis: list[int]
    eye: sp.csr_matrix
    B: list[sp.csr_matrix]
    Cs: list[sp.csr_matrix]
    Ca: list[sp.csr_matrix]
    Z: np.ndarray
    compositions: list[tuple[int, int]]


def active_sector(L: int, npart: int) -> Active:
    basis = fixed_basis(4 * L, npart)
    lookup = {state: k for k, state in enumerate(basis)}
    dim = len(basis)
    eye = sp.eye(dim, format="csr", dtype=complex)
    B: list[sp.csr_matrix] = []
    Cs: list[sp.csr_matrix] = []
    Ca: list[sp.csr_matrix] = []

    for x in range(L):
        y = (x + 1) % L
        bx = sp.csr_matrix((dim, dim), dtype=complex)
        for spin in (0, 1):
            for db, dc, sb, sc in ((0, y, 1, x), (0, x, 1, y)):
                hop = cdagc(basis, mode(L, db, dc, spin), mode(L, sb, sc, spin))
                bx += hop + hop.getH()
        js = []
        for block in (0, 1):
            j = sp.csr_matrix((dim, dim), dtype=complex)
            for spin in (0, 1):
                hop = cdagc(basis, mode(L, block, x, spin), mode(L, block, y, spin))
                j += 1j * hop - 1j * hop.getH()
            js.append(j.tocsr())
        B.append(bx.tocsr())
        Cs.append((js[0] + js[1]).tocsr())
        Ca.append((js[0] - js[1]).tocsr())

    if npart % 2:
        return Active(basis, eye, B, Cs, Ca, np.zeros((dim, 0), complex), [])
    n = npart // 2
    comps = [(n1, n - n1) for n1 in range(max(0, n - L), min(L, n) + 1)]
    cols: list[np.ndarray] = []
    for n1, n2 in comps:
        v = np.zeros(dim, dtype=complex)
        subsets1 = list(itertools.combinations(range(L), n1))
        subsets2 = list(itertools.combinations(range(L), n2))
        norm = math.sqrt(len(subsets1) * len(subsets2))
        for S1 in subsets1:
            for S2 in subsets2:
                state = 0
                for x in S1:
                    state |= 1 << mode(L, 0, x, 0)
                    state |= 1 << mode(L, 0, x, 1)
                for x in S2:
                    state |= 1 << mode(L, 1, x, 0)
                    state |= 1 << mode(L, 1, x, 1)
                v[lookup[state]] = 1.0 / norm
        cols.append(v)
    Z = np.column_stack(cols) if cols else np.zeros((dim, 0), complex)
    return Active(basis, eye, B, Cs, Ca, Z, comps)


@dataclass
class Control:
    states: list[tuple[int, ...]]
    index: dict[tuple[int, ...], int]
    energies: np.ndarray
    mvec: list[sp.csr_matrix]
    Rs: list[sp.csr_matrix]
    Ra: list[sp.csr_matrix]
    V2: list[sp.csr_matrix]


def transition(states, index, site, dst, src):
    rows, cols = [], []
    for col, state in enumerate(states):
        if state[site] != src:
            continue
        out = list(state)
        out[site] = dst
        row = index.get(tuple(out))
        if row is not None:
            rows.append(row)
            cols.append(col)
    return sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(states), len(states)), dtype=complex)


def control_space(L: int, dm: float, ds: float, da: float) -> Control:
    ground = (0,) * L
    states = [x for x in itertools.product(range(4), repeat=L) if x != ground]
    index = {x: k for k, x in enumerate(states)}
    gaps = [0.0, dm, ds, da]
    energies = np.array([sum(gaps[q] for q in x) for x in states], float)
    mvec, Rs, Ra, V2 = [], [], [], []
    for site in range(L):
        st = [0] * L
        st[site] = 1
        row = index[tuple(st)]
        mvec.append(sp.csr_matrix(([1.0], ([row], [0])), shape=(len(states), 1), dtype=complex))
        Rs.append((transition(states, index, site, 2, 1) + transition(states, index, site, 1, 2)).tocsr())
        Ra.append((transition(states, index, site, 3, 1) + transition(states, index, site, 1, 3)).tocsr())
        V2.append((transition(states, index, site, 1, 0) + transition(states, index, site, 0, 1)).tocsr())
    return Control(states, index, energies, mvec, Rs, Ra, V2)


def herm(x: np.ndarray) -> np.ndarray:
    return (x + x.conj().T) / 2


def opnorm(x: np.ndarray) -> float:
    if x.size == 0:
        return 0.0
    return float(la.norm(x, 2))


def build_operators(active: Active, ctrl: Control, b: float, s: float, a: float, dm: float):
    adim = active.eye.shape[0]
    qdim = len(ctrl.states)
    C2 = sp.csr_matrix((qdim * adim, adim), dtype=complex)
    A = sp.csr_matrix((qdim * adim, qdim * adim), dtype=complex)
    Bq = sp.csr_matrix((qdim * adim, qdim * adim), dtype=complex)
    for x in range(len(active.B)):
        C2 += b * sp.kron(ctrl.mvec[x], active.B[x], format="csr")
        A += s * sp.kron(ctrl.Rs[x], active.Cs[x], format="csr")
        A += a * sp.kron(ctrl.Ra[x], active.Ca[x], format="csr")
        Bq += b * sp.kron(ctrl.V2[x], active.B[x], format="csr")
    W4 = (b * b / dm) * sum((x @ x for x in active.B), start=sp.csr_matrix((adim, adim), dtype=complex))
    Cq = sp.kron(sp.eye(qdim, format="csr"), W4, format="csr")
    r0 = np.repeat(1.0 / ctrl.energies, adim)
    return C2, A, Bq, Cq, r0


def r0_apply(y: np.ndarray, r0: np.ndarray) -> np.ndarray:
    return r0[:, None] * y


def word_action(C2, word, r0, Z):
    y = r0_apply(C2 @ Z, r0)
    for op in reversed(word):
        y = r0_apply(op @ y, r0)
    return C2.getH() @ y


def sector_audit(L, npart, ctrl, dm, b, s, a):
    active = active_sector(L, npart)
    Z = active.Z
    C2, A, Bq, Cq, r0 = build_operators(active, ctrl, b, s, a, dm)
    src_a4 = -word_action(C2, (A, A, A, A), r0, Z)
    src_b2 = -word_action(C2, (Bq, Bq), r0, Z)
    src_c = +word_action(C2, (Cq,), r0, Z)
    src_aab = sum((word_action(C2, w, r0, Z) for w in ((A, A, Bq), (A, Bq, A), (Bq, A, A))), start=np.zeros_like(Z))
    classes = {"A4": src_a4, "B2": src_b2, "C": src_c, "AAB": src_aab}
    total = sum(classes.values(), start=np.zeros_like(Z))
    P = Z @ Z.conj().T

    def summarize(source):
        comp = herm(Z.conj().T @ source)
        cross = source - Z @ comp
        return {
            "source_norm": opnorm(source),
            "compression_norm": opnorm(comp),
            "cross_norm": opnorm(cross),
            "compression": comp,
            "eigenvalues": la.eigvalsh(comp) if comp.size else np.array([]),
        }

    return active, {k: summarize(v) for k, v in classes.items()}, summarize(total)


def matrix_text(M: np.ndarray) -> list[str]:
    return ["[" + ", ".join(f"{z.real:+.8e}" if abs(z.imag) < 1e-11 else f"{z.real:+.8e}{z.imag:+.2e}j" for z in row) + "]" for row in M]


def serial(obj):
    if isinstance(obj, np.ndarray):
        if np.iscomplexobj(obj):
            return [[[float(z.real), float(z.imag)] for z in row] for row in obj]
        return [float(x) for x in obj]
    if isinstance(obj, (np.floating, np.integer)):
        return obj.item()
    raise TypeError(type(obj).__name__)


def run_case(label, L, dm, ds, da, b, s, a):
    ctrl = control_space(L, dm, ds, da)
    print(f"CASE {label}: b={b:.12g}, s={s:.12g}, a={a:.12g}")
    print("-" * 78)
    records = []
    for n in range(0, 2 * L + 1):
        npart = 2 * n
        active, classes, total = sector_audit(L, npart, ctrl, dm, b, s, a)
        rec = {
            "pairs": n,
            "particles": npart,
            "active_dim": len(active.basis),
            "compositions": active.compositions,
            "classes": classes,
            "total": total,
        }
        records.append(rec)
        print(f"n={n:2d}  dim={len(active.basis):4d}  comps={active.compositions}")
        print(f"     total: compression={total['compression_norm']:.10e}  cross={total['cross_norm']:.10e}  source={total['source_norm']:.10e}")
        print(f"     class: A4={classes['A4']['compression_norm']:.3e}  B2={classes['B2']['compression_norm']:.3e}  C={classes['C']['compression_norm']:.3e}  AAB={classes['AAB']['source_norm']:.3e}")
        for line in matrix_text(total["compression"]):
            print("       " + line)
        print("     eig = [" + ", ".join(f"{x:+.8e}" for x in total["eigenvalues"]) + "]")
    print()
    return records


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, default=3)
    ap.add_argument("--dm", type=float, default=10.0)
    ap.add_argument("--ds", type=float, default=8.0)
    ap.add_argument("--da", type=float, default=8.0)
    ap.add_argument("--b", type=float, default=0.35)
    ap.add_argument("--s", type=float, default=0.90)
    ap.add_argument("--a", type=float, default=0.45)
    ap.add_argument("--json", type=Path)
    args = ap.parse_args()

    print("ALL-FILLING WEIGHTED-DEGREE-EIGHT MICROSCOPIC AUDIT")
    print("=" * 78)
    print("numpy", np.__version__, "scipy", scipy.__version__)
    print(f"L={args.L}, dm={args.dm}, ds={args.ds}, da={args.da}")
    print()

    generic = run_case("generic matched model", args.L, args.dm, args.ds, args.da, args.b, args.s, args.a)
    bstar = 2.0 * args.s**2 / args.ds
    tuned = run_case("one-pair cancellation b*=2 s^2/ds", args.L, args.dm, args.ds, args.da, bstar, args.s, args.a)

    max_aab = max(x["classes"]["AAB"]["source_norm"] for x in generic + tuned)
    generic_cross = max(x["total"]["cross_norm"] for x in generic)
    tuned_comp = max(x["total"]["compression_norm"] for x in tuned)
    tuned_cross = max(x["total"]["cross_norm"] for x in tuned)
    print("SUMMARY")
    print("-" * 78)
    print(f"max mixed-AAB source norm            = {max_aab:.12e}")
    print(f"max generic cross-source norm        = {generic_cross:.12e}")
    print(f"max tuned ground-compression norm    = {tuned_comp:.12e}")
    print(f"max tuned cross-source norm          = {tuned_cross:.12e}")
    print("generic model preserves zero manifold:", "YES" if max(x['total']['compression_norm'] for x in generic) < 1e-9 and generic_cross < 1e-9 else "NO")
    print("one-pair tuning is all-filling safe:", "YES" if tuned_comp < 1e-9 and tuned_cross < 1e-9 else "NO")

    checks = {
        "mixed_AAB_vanishes": max_aab < 1e-9,
        "particle_hole_endpoints_zero": generic[0]["total"]["source_norm"] < 1e-12 and generic[-1]["total"]["source_norm"] < 1e-12,
        "one_pair_cross_zero": generic[1]["total"]["cross_norm"] < 1e-9,
    }
    print("AUDIT IDENTITIES:", "PASS" if all(checks.values()) else "FAIL")

    if args.json:
        payload = {
            "parameters": vars(args) | {"json": str(args.json), "bstar": bstar},
            "generic": generic,
            "tuned": tuned,
            "summary": {
                "max_mixed_aab": max_aab,
                "max_generic_cross": generic_cross,
                "max_tuned_compression": tuned_comp,
                "max_tuned_cross": tuned_cross,
            },
            "checks": checks,
        }
        args.json.write_text(json.dumps(payload, default=serial, indent=2) + "\n")
    return 0 if all(checks.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
