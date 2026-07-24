#!/usr/bin/env python3
"""Exact-diagonalization search in the actual GHK Model-I/Model-II Hubbard class.

This script uses a fixed 2D projector as the torus size changes, the projected
positive-semidefinite Hubbard interaction
    H(A)=U/2 sum_{R,alpha} (nbar_{R alpha up}(A)-nbar_{R alpha down}(A))^2,
and the literal k -> k+A prescription of Gao-Han-Khalaf.

It works in the N_up=N_down=n and total crystal momentum (0,0) sector.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
from functools import lru_cache
from pathlib import Path
from typing import Callable, Iterable

import numpy as np
from scipy.linalg import eigh
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigsh

Array = np.ndarray


@lru_cache(None)
def masks(M: int, n: int) -> tuple[int, ...]:
    return tuple(sum(1 << i for i in c) for c in itertools.combinations(range(M), n))


@lru_cache(None)
def bilinear_dense(M: int, n: int) -> Array:
    basis = masks(M, n)
    index = {m: i for i, m in enumerate(basis)}
    d = len(basis)
    out = np.zeros((M, M, d, d), np.complex128)
    for p in range(M):
        for q in range(M):
            mat = out[p, q]
            for col, mask in enumerate(basis):
                if not ((mask >> q) & 1):
                    continue
                sign = -1 if (mask & ((1 << q) - 1)).bit_count() % 2 else 1
                ma = mask ^ (1 << q)
                if (ma >> p) & 1:
                    continue
                if (ma & ((1 << p) - 1)).bit_count() % 2:
                    sign *= -1
                mat[index[ma | (1 << p)], col] = sign
    return out


def grid_2d(Nx: int, Ny: int) -> tuple[Array, Array]:
    labels = np.array([(ix, iy) for iy in range(Ny) for ix in range(Nx)], dtype=np.int16)
    k = np.column_stack((2 * np.pi * labels[:, 0] / Nx, 2 * np.pi * labels[:, 1] / Ny))
    return k, labels


def lower_spinor(d: Array) -> Array:
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]], complex)
    sz = np.array([[1, 0], [0, -1]], complex)
    out = np.empty((len(d), 2), complex)
    for i, (x, y, z) in enumerate(d):
        _, v = eigh(x * sx + y * sy + z * sz)
        out[i] = v[:, 0]
    return out


def model_I_d(k: Array, A: Array, p: dict) -> Array:
    x = k[:, 0] + A[0]
    y = k[:, 1] + A[1]
    t1 = p.get("t1", 1.0)
    t2 = p.get("t2", 1 / np.sqrt(2))
    z = -t1 * (
        np.exp(1j * np.pi / 4)
        + np.exp(1j * (x + y + np.pi / 4))
        + np.exp(1j * (x - np.pi / 4))
        + np.exp(1j * (y - np.pi / 4))
    )
    bx = z.real
    by = z.imag
    bz = -2 * t2 * (np.cos(y) - np.cos(x))
    return np.column_stack((bx, by, bz))


def model_II_d(k: Array, A: Array, p: dict) -> Array:
    x = k[:, 0] + A[0]
    y = k[:, 1] + A[1]
    xi = p.get("xi", 1.0)
    alpha = xi * (np.cos(x) + np.cos(y))
    return np.column_stack((-np.sin(alpha), -np.cos(alpha), np.zeros_like(alpha)))


def spinors(model: Callable, k: Array, A: Array, p: dict) -> tuple[Array, Array]:
    up = lower_spinor(model(k, A, p))
    dm = model(-k, -A, p)
    dm = dm.copy(); dm[:, 1] *= -1
    down = lower_spinor(dm)
    return up, down


@lru_cache(None)
def species_momenta(Nx: int, Ny: int, n: int) -> Array:
    M = Nx * Ny
    labels = np.array([(ix, iy) for iy in range(Ny) for ix in range(Nx)], dtype=np.int16)
    out = []
    for mask in masks(M, n):
        occ = [i for i in range(M) if (mask >> i) & 1]
        s = labels[occ].sum(axis=0)
        out.append((int(s[0] % Nx), int(s[1] % Ny)))
    return np.asarray(out, dtype=np.int16)


@lru_cache(None)
def sector_pairs(Nx: int, Ny: int, n: int, Qx: int = 0, Qy: int = 0) -> tuple[tuple[int, int], ...]:
    mom = species_momenta(Nx, Ny, n)
    buckets: dict[tuple[int, int], list[int]] = {}
    for i, q in enumerate(mom):
        buckets.setdefault((int(q[0]), int(q[1])), []).append(i)
    out = []
    for iu, q in enumerate(mom):
        target = ((Qx - int(q[0])) % Nx, (Qy - int(q[1])) % Ny)
        out.extend((iu, idn) for idn in buckets.get(target, []))
    return tuple(out)


def local_orbital_densities(k: Array, labels: Array, Nx: int, Ny: int, spinor: Array, n: int) -> list[Array]:
    M = len(k)
    B = bilinear_dense(M, n)
    out = []
    for Rx, Ry in labels:
        phase = np.exp(1j * (k[:, 0] * Rx + k[:, 1] * Ry))
        for alpha in range(spinor.shape[1]):
            v = phase * spinor[:, alpha]
            one = np.outer(np.conjugate(v), v) / M
            mat = np.einsum("pq,pqij->ij", one, B, optimize=True)
            out.append((mat + mat.conj().T) / 2)
    return out


class HubbardSector:
    def __init__(self, Nx: int, Ny: int, n: int, A: Array, model: Callable, params: dict, U: float = 1.0):
        self.Nx, self.Ny, self.n, self.U = Nx, Ny, n, U
        self.k, self.labels = grid_2d(Nx, Ny)
        self.M = len(self.k)
        up, down = spinors(model, self.k, np.asarray(A, float), params)
        self.nu = local_orbital_densities(self.k, self.labels, Nx, Ny, up, n)
        self.nd = local_orbital_densities(self.k, self.labels, Nx, Ny, down, n)
        self.pairs = sector_pairs(Nx, Ny, n)
        self.iu = np.array([x[0] for x in self.pairs], dtype=np.int32)
        self.idn = np.array([x[1] for x in self.pairs], dtype=np.int32)
        self.eq_u = self.iu[:, None] == self.iu[None, :]
        self.eq_d = self.idn[:, None] == self.idn[None, :]

    def dense(self) -> Array:
        iu, idn = self.iu, self.idn
        H = np.zeros((len(iu), len(iu)), complex)
        for Nu, Nd in zip(self.nu, self.nd):
            Au = Nu @ Nu
            Ad = Nd @ Nd
            H += (self.U / 2) * (
                Au[iu[:, None], iu[None, :]] * self.eq_d
                + self.eq_u * Ad[idn[:, None], idn[None, :]]
                - 2 * Nu[iu[:, None], iu[None, :]] * Nd[idn[:, None], idn[None, :]]
            )
        return (H + H.conj().T) / 2


def low_eigs(Nx: int, Ny: int, n: int, A: Array, model: Callable, params: dict, count: int = 6) -> tuple[Array, int]:
    obj = HubbardSector(Nx, Ny, n, A, model, params)
    H = obj.dense()
    if len(H) <= 5000:
        vals = eigh(H, eigvals_only=True, subset_by_index=[0, min(count - 1, len(H) - 1)], driver="evr")
    else:
        vals = np.sort(eigsh(csr_matrix(H), k=min(count, len(H)-1), which="SA", tol=2e-11, maxiter=10000, return_eigenvectors=False).real)
    return vals, len(H)


def curvature(Nx: int, Ny: int, n: int, model: Callable, params: dict, h: float = 1.5e-3) -> dict:
    A0 = np.zeros(2)
    vals, dim = low_eigs(Nx, Ny, n, A0, model, params, 8)
    e0 = float(vals[0])
    gaps = [float(v - vals[0]) for v in vals[1:] if v - vals[0] > 1e-8]
    gap = gaps[0] if gaps else float("nan")
    deg = int(np.sum(np.abs(vals - vals[0]) < 1e-8))
    es = []
    for s in (-2, -1, 1, 2):
        A = np.array([s * h, 0.0])
        ev, _ = low_eigs(Nx, Ny, n, A, model, params, 1)
        es.append(float(ev[0]))
    em2, em1, ep1, ep2 = es
    curv = (-ep2 + 16 * ep1 - 30 * e0 + 16 * em1 - em2) / (12 * h * h)
    return {"Nx": Nx, "Ny": Ny, "M": Nx * Ny, "n": n, "E0": e0, "gap": gap, "deg": deg, "sector_dim": dim, "curvature": float(curv)}


def evaluate(model_name: str, model: Callable, params: dict, sizes: Iterable[tuple[int,int]], ns_for_size: Callable[[int], list[int]], h: float) -> list[dict]:
    rows = []
    for Nx, Ny in sizes:
        M = Nx * Ny
        local = []
        for n in ns_for_size(M):
            r = curvature(Nx, Ny, n, model, params, h)
            local.append(r)
            print(model_name, Nx, Ny, "n", n, "E", r["E0"], "gap", r["gap"], "c", r["curvature"], "dim", r["sector_dim"], flush=True)
        c1 = local[0]["curvature"]
        for r in local:
            rho = r["n"] * (M - r["n"]) / (M - 1)
            r.update({"model": model_name, "rho": rho, "ratio": r["curvature"] / (rho * c1), "defect": r["curvature"] - rho * c1})
            print("  ratio", r["n"], r["ratio"], "defect", r["defect"], flush=True)
        rows.extend(local)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path(__file__).resolve().parent / "ghk_models_results")
    ap.add_argument("--h", type=float, default=1.5e-3)
    ap.add_argument("--extended", action="store_true")
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    sizes = [(2,2),(2,3),(3,3)]
    if args.extended:
        sizes.append((3,4))
    ns = lambda M: [1, 2] + ([3] if M >= 6 else []) + ([M//2] if M % 2 == 0 and M//2 not in (1,2,3) and M <= 8 else [])
    rows = []
    rows += evaluate("GHK-Model-II-xi1", model_II_d, {"xi":1.0}, sizes, ns, args.h)
    rows += evaluate("GHK-Model-II-xi2", model_II_d, {"xi":2.0}, sizes, ns, args.h)
    rows += evaluate("GHK-Model-I", model_I_d, {}, sizes, ns, args.h)
    keys = list(rows[0])
    with (args.out / "results.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)
    (args.out / "results.json").write_text(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
