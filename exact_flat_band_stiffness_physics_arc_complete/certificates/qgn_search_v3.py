#!/usr/bin/env python3
"""Fixed-local QGN counterexample search under the Gao-Han-Khalaf k->k+A prescription.

The code works in a single rank-one flat band per spin, with time reversal and
fixed pair number N_up=N_down=n.  It constructs Hermitian local neutral
bilinears S_R satisfying the singlet QGN condition at A=0 and positive
translation-invariant interactions H=sum_{RR'} V_{R-R'} S_R S_R'.

Important convention: both the band spinors and the microscopic neutral
bilinears are Peierls-substituted, so their momentum-space form factors are
evaluated at p+A,q+A, as in GHK Eqs. (S88)-(S90).
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np
from scipy.linalg import eigh

Array = np.ndarray


@lru_cache(None)
def masks(M: int, n: int) -> tuple[int, ...]:
    return tuple(sum(1 << i for i in c) for c in itertools.combinations(range(M), n))


@lru_cache(None)
def bilinear_dense(M: int, n: int) -> Array:
    """Matrices of a_p^dag a_q in the fixed-n momentum basis."""
    basis = masks(M, n)
    index = {m: i for i, m in enumerate(basis)}
    d = len(basis)
    out = np.zeros((M, M, d, d), dtype=np.complex128)
    for p in range(M):
        for q in range(M):
            mat = out[p, q]
            for col, mask in enumerate(basis):
                if not ((mask >> q) & 1):
                    continue
                sign = -1 if (mask & ((1 << q) - 1)).bit_count() % 2 else 1
                after = mask ^ (1 << q)
                if (after >> p) & 1:
                    continue
                if (after & ((1 << p) - 1)).bit_count() % 2:
                    sign *= -1
                final = after | (1 << p)
                mat[index[final], col] = sign
    return out


@lru_cache(None)
def mask_momentum_1d(M: int, n: int) -> Array:
    return np.array(
        [sum(i for i in range(M) if (mask >> i) & 1) % M for mask in masks(M, n)],
        dtype=np.int16,
    )


@lru_cache(None)
def sector_pairs_1d(M: int, n: int, Q: int = 0) -> tuple[tuple[int, int], ...]:
    mom = mask_momentum_1d(M, n)
    return tuple(
        (iu, idn)
        for iu, mu in enumerate(mom)
        for idn, md in enumerate(mom)
        if (int(mu) + int(md) - Q) % M == 0
    )


def grid_1d(M: int) -> Array:
    return (2.0 * np.pi * np.arange(M) / M)[:, None]


def lower_spinor(d: Array) -> Array:
    sx = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    sy = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    sz = np.array([[1, 0], [0, -1]], dtype=np.complex128)
    out = np.empty((len(d), 2), dtype=np.complex128)
    for i, (x, y, z) in enumerate(d):
        _, v = eigh(x * sx + y * sy + z * sz)
        out[i] = v[:, 0]
    return out


def model_phase_1d(k: Array, A: Array, p: dict) -> Array:
    """1D restriction of GHK Model II type; exactly finite-range only after spectral flattening."""
    x = k[:, 0] + A[0]
    alpha = p.get("xi", 1.1) * np.cos(x) + p.get("zeta", 0.0) * np.cos(2 * x)
    return np.column_stack((-np.sin(alpha), -np.cos(alpha), np.zeros_like(alpha)))


def model_flattened_1d(k: Array, A: Array, p: dict) -> Array:
    """Fixed trigonometric parent vector, spectrally flattened to an exact band projector."""
    x = k[:, 0] + A[0]
    v = np.column_stack(
        (
            p.get("a", 0.3) + p.get("bx", 0.8) * np.cos(x) + p.get("cx", 0.0) * np.cos(2 * x),
            p.get("by", 0.7) * np.sin(x) + p.get("cy", 0.0) * np.sin(2 * x),
            p.get("m", 1.4) + p.get("bz", 0.5) * np.cos(x) + p.get("cz", 0.0) * np.cos(2 * x),
        )
    )
    norm = np.linalg.norm(v, axis=1)
    if float(norm.min()) < 1e-9:
        raise ValueError("parent gap closes")
    return v / norm[:, None]


def spinors(model: Callable[[Array, Array, dict], Array], k: Array, A: Array, p: dict) -> tuple[Array, Array]:
    d_up = model(k, A, p)
    up = lower_spinor(d_up)
    d_minus = model(-k, -A, p)
    d_down = d_minus.copy()
    d_down[:, 1] *= -1  # complex conjugation sends sigma_y -> -sigma_y
    down = lower_spinor(d_down)
    return up, down


@dataclass(frozen=True)
class DirectedComponent:
    """One microscopic term coeff c^dag_{R+left} Aorb c_{R+right}."""

    Aorb: tuple[tuple[complex, ...], ...]
    left: int = 0
    right: int = 0
    coeff: complex = 1.0

    def matrix(self) -> Array:
        return np.asarray(self.Aorb, dtype=np.complex128)


@dataclass(frozen=True)
class HermitianChannel:
    name: str
    components: tuple[DirectedComponent, ...]
    weight: float = 1.0


def onsite_channel(Aorb: Array, name: str = "onsite", weight: float = 1.0) -> HermitianChannel:
    Aorb = np.asarray(Aorb, dtype=np.complex128)
    if np.linalg.norm(Aorb - Aorb.conj().T) > 1e-12:
        raise ValueError("onsite Aorb must be Hermitian")
    return HermitianChannel(
        name,
        (DirectedComponent(tuple(map(tuple, Aorb)), 0, 0, 1.0),),
        weight,
    )


def bond_channel(Aorb: Array, delta: int, phase: float = 0.0, name: str | None = None, weight: float = 1.0) -> HermitianChannel:
    """Hermitian e^{i phase} c_R^dag A c_{R+delta} + h.c."""
    Aorb = np.asarray(Aorb, dtype=np.complex128)
    c = np.exp(1j * phase)
    comps = (
        DirectedComponent(tuple(map(tuple, Aorb)), 0, delta, c),
        DirectedComponent(tuple(map(tuple, Aorb.conj().T)), delta, 0, np.conjugate(c)),
    )
    return HermitianChannel(name or f"bond{delta}", comps, weight)


def _projected_components(
    k: Array,
    band_spinor: Array,
    n: int,
    components: Sequence[DirectedComponent],
    twist: Array,
) -> list[Array]:
    """Projected translations of a microscopic channel for one spin.

    Following GHK Eqs. (S88)-(S90), the microscopic neutral bilinear is also
    Peierls-substituted, so displacement phases use p+A and q+A in addition
    to the band spinors u(k+A).
    """
    M = len(k)
    B = bilinear_dense(M, n)
    base = np.zeros((M, M), dtype=np.complex128)
    x = k[:, 0]
    xt = x + float(np.asarray(twist)[0])
    for comp in components:
        Aorb = comp.matrix()
        orb = np.einsum("pa,ab,qb->pq", np.conjugate(band_spinor), Aorb, band_spinor, optimize=True)
        left_phase = np.exp(-1j * xt * comp.left)
        right_phase = np.exp(1j * xt * comp.right)
        base += comp.coeff * np.einsum("p,pq,q->pq", left_phase, orb, right_phase)
    out: list[Array] = []
    for R in range(M):
        phase = np.exp(1j * x * R)
        one = np.einsum("p,pq,q->pq", np.conjugate(phase), base, phase) / M
        mat = np.einsum("pq,pqij->ij", one, B, optimize=True)
        out.append((mat + mat.conj().T) / 2)
    return out


def qgn_down_components(up_components: Sequence[DirectedComponent]) -> tuple[DirectedComponent, ...]:
    """Microscopic time-reversed QGN partner yielding B_pq=-A_{-q,-p}."""
    out = []
    for comp in up_components:
        A = comp.matrix()
        out.append(
            DirectedComponent(
                tuple(map(tuple, -A.T)),
                left=comp.right,
                right=comp.left,
                coeff=comp.coeff,
            )
        )
    return tuple(out)


class HermitianQGNSectorHamiltonian:
    """Hamiltonian in total momentum Q=0 sector for sum_R,ch w_ch S_R,ch^2.

    Each microscopic channel S_R is Hermitian and its down-spin block is the
    QGN partner of the up-spin block.  Positive combinations of channels can be
    added.  Cross-channel positive kernels can be represented by first taking
    linear combinations at the microscopic-channel level.
    """

    def __init__(self, k: Array, model, params: dict, n: int, A: Array, channels: Sequence[HermitianChannel], U: float = 1.0, Q: int = 0):
        self.k = np.asarray(k)
        self.M = len(k)
        self.n = n
        self.U = U
        self.pairs = sector_pairs_1d(self.M, n, Q)
        self.iu = np.array([p[0] for p in self.pairs], dtype=np.int32)
        self.idn = np.array([p[1] for p in self.pairs], dtype=np.int32)
        self.eq_u = self.iu[:, None] == self.iu[None, :]
        self.eq_d = self.idn[:, None] == self.idn[None, :]
        up, down = spinors(model, self.k, np.asarray(A, dtype=float), params)
        self.ops: list[tuple[float, Array, Array]] = []
        for channel in channels:
            ups = _projected_components(self.k, up, n, channel.components, A)
            downs = _projected_components(self.k, down, n, qgn_down_components(channel.components), A)
            for Su, Sd in zip(ups, downs):
                self.ops.append((float(channel.weight), Su, Sd))

    def dense(self) -> Array:
        iu, idn = self.iu, self.idn
        out = np.zeros((len(iu), len(iu)), dtype=np.complex128)
        for w, Su, Sd in self.ops:
            Au = Su @ Su
            Ad = Sd @ Sd
            # (Su⊗I + I⊗Sd)^2 = Su^2⊗I + I⊗Sd^2 + 2 Su⊗Sd.
            out += (self.U / 2) * w * (
                Au[iu[:, None], iu[None, :]] * self.eq_d
                + self.eq_u * Ad[idn[:, None], idn[None, :]]
                + 2.0 * Su[iu[:, None], iu[None, :]] * Sd[idn[:, None], idn[None, :]]
            )
        return (out + out.conj().T) / 2


def low_spectrum(k, model, params, n, A, channels, count: int = 5) -> Array:
    H = HermitianQGNSectorHamiltonian(k, model, params, n, A, channels).dense()
    hi = min(count - 1, len(H) - 1)
    return eigh(H, eigvals_only=True, subset_by_index=[0, hi], driver="evr")


def curvature(k, model, params, n, channels, h: float = 1.5e-3) -> dict:
    A0 = np.zeros(1)
    spec = low_spectrum(k, model, params, n, A0, channels, 8)
    e0 = float(spec[0])
    distinct = [float(x - spec[0]) for x in spec[1:] if x - spec[0] > 1e-8]
    gap = distinct[0] if distinct else float("nan")
    deg = int(np.sum(np.abs(spec - spec[0]) < 1e-8))
    es = []
    for s in (-2, -1, 1, 2):
        A = np.array([s * h])
        es.append(float(low_spectrum(k, model, params, n, A, channels, 1)[0]))
    em2, em1, ep1, ep2 = es
    c = (-ep2 + 16 * ep1 - 30 * e0 + 16 * em1 - em2) / (12 * h * h)
    return {"n": n, "E0": e0, "gap": gap, "deg": deg, "curvature": float(c), "sector_dim": len(sector_pairs_1d(len(k), n, 0))}


def evaluate_case(name: str, M: int, model, params: dict, channels: Sequence[HermitianChannel], ns: Iterable[int], h: float = 1.5e-3) -> list[dict]:
    k = grid_1d(M)
    rows = [curvature(k, model, params, n, channels, h) for n in ns]
    c1 = rows[0]["curvature"]
    for row in rows:
        n = row["n"]
        rho = n * (M - n) / (M - 1)
        row.update(
            {
                "case": name,
                "M": M,
                "rho": rho,
                "ratio": row["curvature"] / (rho * c1),
                "defect": row["curvature"] - rho * c1,
            }
        )
    return rows


def combine_channels(channels: Sequence[HermitianChannel], coeffs: Sequence[float], name: str) -> HermitianChannel:
    comps: list[DirectedComponent] = []
    for ch, coeff in zip(channels, coeffs):
        for c in ch.components:
            comps.append(DirectedComponent(c.Aorb, c.left, c.right, c.coeff * coeff))
    return HermitianChannel(name, tuple(comps), 1.0)


def standard_channels() -> dict[str, HermitianChannel]:
    I = np.eye(2, dtype=np.complex128)
    sx = np.array([[0, 1], [1, 0]], dtype=np.complex128)
    sy = np.array([[0, -1j], [1j, 0]], dtype=np.complex128)
    sz = np.array([[1, 0], [0, -1]], dtype=np.complex128)
    P0 = np.array([[1, 0], [0, 0]], dtype=np.complex128)
    P1 = np.array([[0, 0], [0, 1]], dtype=np.complex128)
    return {
        "P0": onsite_channel(P0, "P0"),
        "P1": onsite_channel(P1, "P1"),
        "I": onsite_channel(I, "I"),
        "sx": onsite_channel(sx, "sx"),
        "sy": onsite_channel(sy, "sy"),
        "sz": onsite_channel(sz, "sz"),
        "bI1c": bond_channel(I, 1, 0.0, "bond-I-cos"),
        "bI1s": bond_channel(I, 1, np.pi / 2, "bond-I-sin"),
        "bsx1c": bond_channel(sx, 1, 0.0, "bond-sx-cos"),
        "bsx1s": bond_channel(sx, 1, np.pi / 2, "bond-sx-sin"),
        "bsz1c": bond_channel(sz, 1, 0.0, "bond-sz-cos"),
        "bI2c": bond_channel(I, 2, 0.0, "bond-I2-cos"),
    }


def deterministic_suite() -> list[dict]:
    ch = standard_channels()
    model = model_flattened_1d
    params = {"a": 0.3, "m": 1.4, "bz": 0.5, "bx": 0.8, "by": 0.7}
    cases = {
        "hubbard-orbital": [ch["P0"], ch["P1"]],
        "onsite-sx": [ch["sx"]],
        "bond-I-cos": [ch["bI1c"]],
        "bond-I-sin": [ch["bI1s"]],
        "bond-sx-cos": [ch["bsx1c"]],
        "mixed-onsite-bond": [combine_channels([ch["P0"], ch["bsx1c"]], [1.0, 0.7], "mix")],
        "two-independent": [ch["P0"], ch["bI1c"]],
        "range-two-mix": [combine_channels([ch["sx"], ch["bI2c"]], [1.0, 0.45], "r2mix")],
    }
    all_rows: list[dict] = []
    for case_name, channels in cases.items():
        for M in (4, 6, 8):
            ns = [1, M // 2] if M >= 6 else [1, 2, 3]
            rows = evaluate_case(case_name, M, model, params, channels, ns)
            all_rows.extend(rows)
            print(case_name, "M", M, [(r["n"], r["ratio"], r["gap"]) for r in rows], flush=True)
    return all_rows


def random_search(seed: int = 20260722, trials: int = 30, Ms: Sequence[int] = (4, 6, 8)) -> list[dict]:
    rng = np.random.default_rng(seed)
    base = standard_channels()
    pool = [base[k] for k in ("P0", "P1", "sx", "sy", "sz", "bI1c", "bI1s", "bsx1c", "bsx1s", "bsz1c", "bI2c")]
    model = model_flattened_1d
    params = {"a": 0.3, "m": 1.4, "bz": 0.5, "bx": 0.8, "by": 0.7, "cx": 0.12, "cy": -0.08, "cz": 0.11}
    summary = []
    for t in range(trials):
        coeff = rng.normal(size=len(pool))
        coeff /= np.linalg.norm(coeff)
        channel = combine_channels(pool, coeff, f"random-{t:03d}")
        trial_rows = []
        for M in Ms:
            rows = evaluate_case(channel.name, M, model, params, [channel], [1, M // 2])
            trial_rows.extend(rows)
        max_dev = max(abs(r["ratio"] - 1) for r in trial_rows if r["n"] != 1)
        summary.append({"trial": t, "coeffs": coeff.tolist(), "max_dev": max_dev, "rows": trial_rows})
        print("trial", t, "max_dev", max_dev, flush=True)
    summary.sort(key=lambda z: z["max_dev"], reverse=True)
    return summary


def write_csv(rows: Sequence[dict], path: Path) -> None:
    if not rows:
        return
    keys = list(rows[0].keys())
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--random", type=int, default=0, help="number of random Hermitian QGN channels")
    ap.add_argument("--out", type=Path, default=Path("/mnt/data/qgn_search_v2_results"))
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = deterministic_suite()
    write_csv(rows, args.out / "deterministic.csv")
    if args.random:
        res = random_search(trials=args.random)
        (args.out / "random.json").write_text(json.dumps(res, indent=2))
        best_rows = [row for z in res[:5] for row in z["rows"]]
        write_csv(best_rows, args.out / "random_top5.csv")


if __name__ == "__main__":
    main()
