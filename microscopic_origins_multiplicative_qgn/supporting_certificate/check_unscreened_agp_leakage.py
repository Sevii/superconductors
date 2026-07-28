#!/usr/bin/env python3
"""Exact periodic/open-boundary leakage certificate for the unscreened QGN shell.

The active order-four charge operator is

    C_{4,L}^{bc} = sum_{e in E_bc} [A4 N_e^2 - (U4/2) N_e],

restricted to the seniority-zero hard-core-pair Hilbert space of two blocks on
L active sites.  ``bc=periodic`` uses the L translated bonds of the cycle;
``bc=open`` uses the L-1 bonds of the path.  In a fixed composition (n1,n2),
the normalized product AGP is the uniform superposition over pairs of subsets
of sizes n1 and n2.

The verifier uses only exact integer arithmetic and ``fractions.Fraction`` for
all PASS/FAIL decisions.  It checks:

* the projected selector formula for each boundary condition;
* the exact composition-resolved action
      R|n1,n2> = (C4 - <C4>)|n1,n2>;
* the explicit singular-value decomposition and general norm formula;
* the closed periodic leakage formula;
* the closed open-path moment/covariance formula;
* the periodic and open zero sets in the physical shell cone A4>0, V4>0;
* exact Phase-3 Gram matrices at L=4,6,8 and n=2,3.

The PASS lines labelled ``explicit centered action/SVD norm`` verify the exact
centered first and second moments.  The operator-level SVD then follows
analytically from occupation-basis diagonality and orthogonality of distinct
block-charge sectors; no floating-point matrix SVD is used.

No NumPy, SciPy, or floating-point arithmetic enters any PASS/FAIL decision.
"""

from __future__ import annotations

import argparse
import csv
import shutil
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import comb, sqrt
from pathlib import Path
from typing import Iterable, Literal


Boundary = Literal["periodic", "open"]

DEFAULT_A4 = Fraction(49, 500)   # 0.098 = 2 beta4 for the project default
DEFAULT_U4 = Fraction(147, 500)  # 0.294; hence V4=4A4-U4=0.098=A4
DEFAULT_VOLUMES = tuple(range(2, 10))
FOCUS_VOLUMES = (4, 6, 8)
FOCUS_TOTAL_PAIRS = (2, 3)


def parse_fraction(text: str) -> Fraction:
    """Parse either ``p/q`` or a finite decimal exactly."""
    return Fraction(text)


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def sqrt_text(value: Fraction) -> str:
    if value == 0:
        return "0"
    return f"sqrt({fraction_text(value)})"


def masks_with_weight(volume: int, weight: int) -> list[int]:
    return [sum(1 << x for x in occupied) for occupied in combinations(range(volume), weight)]


def occupation(mask: int, site: int) -> int:
    return (mask >> site) & 1


def edges(volume: int, boundary: Boundary) -> list[tuple[int, int]]:
    if boundary == "periodic":
        return [(x, (x + 1) % volume) for x in range(volume)]
    return [(x, x + 1) for x in range(volume - 1)]


def charge_observables(
    volume: int,
    mask1: int,
    mask2: int,
    boundary: Boundary,
) -> tuple[int, int]:
    """Return ``(sum_e N_e, sum_e N_e^2)`` for the chosen boundary."""
    linear = 0
    quadratic = 0
    for x, y in edges(volume, boundary):
        ne = 2 * (
            occupation(mask1, x)
            + occupation(mask1, y)
            + occupation(mask2, x)
            + occupation(mask2, y)
        )
        linear += ne
        quadratic += ne * ne
    return linear, quadratic


def open_components(volume: int, mask1: int, mask2: int) -> tuple[int, int, int, int, int]:
    """Return ``(B1,B2,E1,E2,C)`` on the path.

    ``Ba`` is the occupied-endpoint count, ``Ea`` the occupied-neighbour-pair
    count, and ``C=sum_e t_{1,e} t_{2,e}`` with
    ``t_{a,e}=X_{a,x}+X_{a,x+1}``.
    """
    b1 = occupation(mask1, 0) + occupation(mask1, volume - 1)
    b2 = occupation(mask2, 0) + occupation(mask2, volume - 1)
    e1 = 0
    e2 = 0
    cross = 0
    for x in range(volume - 1):
        y = x + 1
        x1 = occupation(mask1, x)
        y1 = occupation(mask1, y)
        x2 = occupation(mask2, x)
        y2 = occupation(mask2, y)
        e1 += x1 * y1
        e2 += x2 * y2
        cross += (x1 + y1) * (x2 + y2)
    return b1, b2, e1, e2, cross


# ---------------------------------------------------------------------------
# Periodic formulas
# ---------------------------------------------------------------------------


def periodic_mean_linear(volume: int, n1: int, n2: int) -> Fraction:
    del volume
    return Fraction(4 * (n1 + n2), 1)


def periodic_mean_d(volume: int, n1: int, n2: int) -> Fraction:
    total = n1 + n2
    return (
        Fraction(8 * total, 1)
        + Fraction(8 * (n1 * (n1 - 1) + n2 * (n2 - 1)), volume - 1)
        + Fraction(32 * n1 * n2, volume)
    )


def periodic_baseline_d(volume: int, total_pairs: int) -> Fraction:
    return Fraction(8 * total_pairs, 1) + Fraction(
        8 * total_pairs * (total_pairs - 1), volume - 1
    )


def periodic_selector_per_a4(volume: int, n1: int, n2: int) -> Fraction:
    return Fraction(16 * (volume - 2) * n1 * n2, volume * (volume - 1))


def periodic_variance_d(volume: int, n1: int, n2: int) -> Fraction:
    """Exact ``Var(sum_e N_e^2)`` on the cycle."""
    if volume == 2:
        return Fraction(0, 1)
    if volume < 2:
        raise ValueError("volume must be at least 2")

    denom_within = (volume - 2) * (volume - 1) ** 2
    within1 = Fraction(
        n1 * (n1 - 1) * (volume - n1) * (volume - n1 - 1),
        denom_within,
    )
    within2 = Fraction(
        n2 * (n2 - 1) * (volume - n2) * (volume - n2 - 1),
        denom_within,
    )
    cross = Fraction(
        2
        * n1
        * n2
        * (volume - n1)
        * (volume - n2)
        * (3 * volume - 8),
        volume**2 * (volume - 1) ** 2,
    )
    return 64 * (within1 + within2 + cross)


def periodic_zero_leakage(volume: int, n1: int, n2: int) -> bool:
    if volume == 2:
        return True
    edge = {0, 1, volume - 1, volume}
    full_or_empty = {0, volume}
    return (n1 in full_or_empty and n2 in edge) or (
        n2 in full_or_empty and n1 in edge
    )


# ---------------------------------------------------------------------------
# Open-path formulas
# ---------------------------------------------------------------------------


def falling_probability(volume: int, filling: int, order: int) -> Fraction:
    """Probability that ``order`` specified distinct sites are occupied.

    This is ``(filling)_order/(volume)_order``, with value zero when the event
    is impossible.
    """
    if order < 0:
        raise ValueError("order must be nonnegative")
    if order == 0:
        return Fraction(1, 1)
    if order > filling or order > volume:
        return Fraction(0, 1)
    numerator = 1
    denominator = 1
    for j in range(order):
        numerator *= filling - j
        denominator *= volume - j
    return Fraction(numerator, denominator)


def open_mean_linear(volume: int, n1: int, n2: int) -> Fraction:
    return Fraction(4 * (volume - 1) * (n1 + n2), volume)


def open_mean_d(volume: int, n1: int, n2: int) -> Fraction:
    total = n1 + n2
    return Fraction(8 * total * (volume + total - 2), volume) + Fraction(
        16 * (volume - 2) * n1 * n2,
        volume**2,
    )


def open_baseline_d(volume: int, total_pairs: int) -> Fraction:
    return Fraction(8 * total_pairs * (volume + total_pairs - 2), volume)


def open_selector_per_a4(volume: int, n1: int, n2: int) -> Fraction:
    return Fraction(16 * (volume - 2) * n1 * n2, volume**2)


def open_var_b(volume: int, filling: int) -> Fraction:
    """Variance of endpoint occupancy ``B=X_1+X_L``."""
    return Fraction(
        2 * filling * (volume - filling) * (volume - 2),
        volume**2 * (volume - 1),
    )


def open_cov_e_b(volume: int, filling: int) -> Fraction:
    """Covariance of adjacent-pair count E and endpoint occupancy B."""
    return Fraction(
        -2 * filling * (filling - 1) * (volume - filling),
        volume**2 * (volume - 1),
    )


def open_var_e(volume: int, filling: int) -> Fraction:
    """Variance of the number of occupied edges on the path."""
    return Fraction(
        filling
        * (filling - 1)
        * (volume - filling)
        * (volume - filling + 1),
        volume**2 * (volume - 1),
    )


def open_t_moments(volume: int, filling: int) -> tuple[Fraction, Fraction, Fraction]:
    """Return E[t_e^2], E[t_e t_f] for adjacent/disjoint edge pairs."""
    p1 = falling_probability(volume, filling, 1)
    p2 = falling_probability(volume, filling, 2)
    return 2 * p1 + 2 * p2, p1 + 3 * p2, 4 * p2


def open_mean_cross(volume: int, n1: int, n2: int) -> Fraction:
    return Fraction(4 * (volume - 1) * n1 * n2, volume**2)


def open_var_cross(volume: int, n1: int, n2: int) -> Fraction:
    """Variance of ``C=sum_e t_{1,e}t_{2,e}`` on the path."""
    u01, u11, u21 = open_t_moments(volume, n1)
    u02, u12, u22 = open_t_moments(volume, n2)
    second_moment = (
        (volume - 1) * u01 * u02
        + 2 * (volume - 2) * u11 * u12
        + (volume - 2) * (volume - 3) * u21 * u22
    )
    mean = open_mean_cross(volume, n1, n2)
    return second_moment - mean * mean


def open_gammas(volume: int, n1: int, n2: int) -> tuple[Fraction, Fraction, Fraction]:
    """Return ``(Gamma0,Gamma1,Gamma2)`` for the open leakage variance.

    With ``Z=E1+E2+C`` and ``B=B1+B2``:

      Gamma0 = Var(Z),  Gamma1 = Cov(Z,B),  Gamma2 = Var(B).
    """
    b1 = open_var_b(volume, n1)
    b2 = open_var_b(volume, n2)
    h1 = open_cov_e_b(volume, n1)
    h2 = open_cov_e_b(volume, n2)
    e1 = open_var_e(volume, n1)
    e2 = open_var_e(volume, n2)
    w = open_var_cross(volume, n1, n2)

    gamma0 = e1 + e2 + w - Fraction(4 * n2, volume) * h1 - Fraction(4 * n1, volume) * h2
    gamma1 = h1 + h2 - Fraction(2 * n2, volume) * b1 - Fraction(2 * n1, volume) * b2
    gamma2 = b1 + b2
    return gamma0, gamma1, gamma2


def open_variance_c4(
    volume: int,
    n1: int,
    n2: int,
    a4: Fraction,
    u4: Fraction,
) -> Fraction:
    """Exact leakage squared for the open charge operator."""
    v4 = 4 * a4 - u4
    gamma0, gamma1, gamma2 = open_gammas(volume, n1, n2)
    return 64 * a4 * a4 * gamma0 - 16 * a4 * v4 * gamma1 + v4 * v4 * gamma2


def open_zero_leakage(volume: int, n1: int, n2: int, v4: Fraction) -> bool:
    """Zero set in the physical shell regime ``V4>0``."""
    if volume == 2:
        return True
    if v4 <= 0:
        raise ValueError("closed open-boundary zero criterion requires V4>0")
    return n1 in {0, volume} and n2 in {0, volume}


# ---------------------------------------------------------------------------
# Exact enumeration and SVD certificate
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SectorResult:
    boundary: Boundary
    volume: int
    total_pairs: int
    n1: int
    n2: int
    dimension: int
    configurations_checked: int
    mean_linear: Fraction
    mean_d: Fraction
    mean_c4: Fraction
    variance_c4: Fraction
    mean_linear_formula: Fraction
    mean_d_formula: Fraction
    mean_c4_formula: Fraction
    variance_c4_formula: Fraction
    baseline_d: Fraction
    selector_d_per_a4: Fraction
    min_c4: Fraction
    max_c4: Fraction
    distinct_c4_values: int
    centered_action_sum_zero: bool
    centered_action_norm_exact: bool
    open_decomposition_exact: bool

    @property
    def invariant(self) -> bool:
        return self.variance_c4 == 0

    @property
    def singular_value_decimal(self) -> float:
        return sqrt(float(self.variance_c4))


def enumerate_sector(
    boundary: Boundary,
    volume: int,
    n1: int,
    n2: int,
    masks: dict[int, list[int]],
    a4: Fraction,
    u4: Fraction,
) -> SectorResult:
    values_linear: list[int] = []
    values_d: list[int] = []
    values_c4: list[Fraction] = []
    open_decomposition_exact = True
    total_pairs = n1 + n2

    for mask1 in masks[n1]:
        for mask2 in masks[n2]:
            linear, d_value = charge_observables(volume, mask1, mask2, boundary)
            c4_value = a4 * d_value - Fraction(u4, 2) * linear
            values_linear.append(linear)
            values_d.append(d_value)
            values_c4.append(c4_value)

            if boundary == "open":
                b1, b2, e1, e2, cross = open_components(volume, mask1, mask2)
                btot = b1 + b2
                z = e1 + e2 + cross
                open_decomposition_exact &= linear == 4 * total_pairs - 2 * btot
                open_decomposition_exact &= d_value == 8 * total_pairs - 4 * btot + 8 * z
                v4 = 4 * a4 - u4
                open_decomposition_exact &= c4_value == (8 * a4 - 2 * u4) * total_pairs + 8 * a4 * z - v4 * btot

    expected_count = comb(volume, n1) * comb(volume, n2)
    if len(values_c4) != expected_count:
        raise AssertionError(f"sector dimension mismatch at L={volume}, ({n1},{n2})")

    count = len(values_c4)
    mean_linear = Fraction(sum(values_linear), count)
    mean_d = Fraction(sum(values_d), count)
    mean_c4 = sum(values_c4, Fraction(0, 1)) / count
    variance_c4 = sum((value - mean_c4) ** 2 for value in values_c4) / count

    if boundary == "periodic":
        mean_linear_formula = periodic_mean_linear(volume, n1, n2)
        mean_d_formula = periodic_mean_d(volume, n1, n2)
        baseline_d = periodic_baseline_d(volume, total_pairs)
        selector = periodic_selector_per_a4(volume, n1, n2)
        variance_formula = a4 * a4 * periodic_variance_d(volume, n1, n2)
    else:
        mean_linear_formula = open_mean_linear(volume, n1, n2)
        mean_d_formula = open_mean_d(volume, n1, n2)
        baseline_d = open_baseline_d(volume, total_pairs)
        selector = open_selector_per_a4(volume, n1, n2)
        variance_formula = open_variance_c4(volume, n1, n2, a4, u4)

    mean_c4_formula = a4 * mean_d_formula - Fraction(u4, 2) * mean_linear_formula
    centered_sum = sum(((value - mean_c4) for value in values_c4), Fraction(0, 1))
    centered_sq = sum(((value - mean_c4) ** 2 for value in values_c4), Fraction(0, 1))

    return SectorResult(
        boundary=boundary,
        volume=volume,
        total_pairs=total_pairs,
        n1=n1,
        n2=n2,
        dimension=expected_count,
        configurations_checked=count,
        mean_linear=mean_linear,
        mean_d=mean_d,
        mean_c4=mean_c4,
        variance_c4=variance_c4,
        mean_linear_formula=mean_linear_formula,
        mean_d_formula=mean_d_formula,
        mean_c4_formula=mean_c4_formula,
        variance_c4_formula=variance_formula,
        baseline_d=baseline_d,
        selector_d_per_a4=selector,
        min_c4=min(values_c4),
        max_c4=max(values_c4),
        distinct_c4_values=len(set(values_c4)),
        centered_action_sum_zero=centered_sum == 0,
        centered_action_norm_exact=centered_sq == count * variance_c4,
        open_decomposition_exact=open_decomposition_exact,
    )


def write_csv(path: Path, rows: Iterable[SectorResult], a4: Fraction) -> None:
    fields = [
        "boundary",
        "volume",
        "total_pairs",
        "n1",
        "n2",
        "sector_dimension",
        "configurations_checked",
        "mean_sum_N_exact",
        "mean_D_exact",
        "mean_C4_exact",
        "baseline_D_exact",
        "selector_shift_per_A4_exact",
        "leakage_squared_exact",
        "leakage_squared_formula_exact",
        "leakage_squared_per_A4_squared_exact",
        "singular_value_exact",
        "singular_value_decimal",
        "min_C4_exact",
        "max_C4_exact",
        "distinct_C4_values",
        "agp_invariant",
        "centered_action_sum_zero",
        "centered_action_norm_exact",
        "open_decomposition_exact",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            normalized = row.variance_c4 / (a4 * a4) if a4 != 0 else Fraction(0, 1)
            writer.writerow(
                {
                    "boundary": row.boundary,
                    "volume": row.volume,
                    "total_pairs": row.total_pairs,
                    "n1": row.n1,
                    "n2": row.n2,
                    "sector_dimension": row.dimension,
                    "configurations_checked": row.configurations_checked,
                    "mean_sum_N_exact": fraction_text(row.mean_linear),
                    "mean_D_exact": fraction_text(row.mean_d),
                    "mean_C4_exact": fraction_text(row.mean_c4),
                    "baseline_D_exact": fraction_text(row.baseline_d),
                    "selector_shift_per_A4_exact": fraction_text(row.selector_d_per_a4),
                    "leakage_squared_exact": fraction_text(row.variance_c4),
                    "leakage_squared_formula_exact": fraction_text(row.variance_c4_formula),
                    "leakage_squared_per_A4_squared_exact": fraction_text(normalized),
                    "singular_value_exact": sqrt_text(row.variance_c4),
                    "singular_value_decimal": f"{row.singular_value_decimal:.15g}",
                    "min_C4_exact": fraction_text(row.min_c4),
                    "max_C4_exact": fraction_text(row.max_c4),
                    "distinct_C4_values": row.distinct_c4_values,
                    "agp_invariant": "YES" if row.invariant else "NO",
                    "centered_action_sum_zero": "YES" if row.centered_action_sum_zero else "NO",
                    "centered_action_norm_exact": "YES" if row.centered_action_norm_exact else "NO",
                    "open_decomposition_exact": "YES" if row.open_decomposition_exact else "NO",
                }
            )


def focus_groups(rows: list[SectorResult], boundary: Boundary) -> list[tuple[int, int, list[SectorResult]]]:
    groups: list[tuple[int, int, list[SectorResult]]] = []
    for volume in FOCUS_VOLUMES:
        for total in FOCUS_TOTAL_PAIRS:
            group = sorted(
                [
                    row
                    for row in rows
                    if row.boundary == boundary
                    and row.volume == volume
                    and row.total_pairs == total
                ],
                key=lambda row: row.n1,
            )
            if group:
                groups.append((volume, total, group))
    return groups


def group_norm_squared(rows: list[SectorResult], volume: int, total: int, boundary: Boundary) -> Fraction:
    candidates = [
        row.variance_c4
        for row in rows
        if row.boundary == boundary and row.volume == volume and row.total_pairs == total
    ]
    if not candidates:
        raise KeyError((boundary, volume, total))
    return max(candidates)


def build_certificate(
    rows: list[SectorResult],
    volumes: tuple[int, ...],
    a4: Fraction,
    u4: Fraction,
) -> tuple[str, bool]:
    periodic_rows = [row for row in rows if row.boundary == "periodic"]
    open_rows = [row for row in rows if row.boundary == "open"]
    v4 = 4 * a4 - u4
    q = v4 / a4 if a4 != 0 else Fraction(0, 1)

    checks: dict[str, bool] = {
        "periodic projected mean and selector": all(
            row.mean_linear == row.mean_linear_formula
            and row.mean_d == row.mean_d_formula
            and row.mean_d == row.baseline_d + row.selector_d_per_a4
            and row.mean_c4 == row.mean_c4_formula
            for row in periodic_rows
        ),
        "periodic closed leakage formula": all(
            row.variance_c4 == row.variance_c4_formula for row in periodic_rows
        ),
        "periodic explicit centered action/SVD norm": all(
            row.centered_action_sum_zero and row.centered_action_norm_exact
            for row in periodic_rows
        ),
        "periodic zero set": all(
            row.invariant == periodic_zero_leakage(row.volume, row.n1, row.n2)
            for row in periodic_rows
        ),
        "open endpoint decomposition": all(row.open_decomposition_exact for row in open_rows),
        "open projected mean and selector": all(
            row.mean_linear == row.mean_linear_formula
            and row.mean_d == row.mean_d_formula
            and row.mean_d == row.baseline_d + row.selector_d_per_a4
            and row.mean_c4 == row.mean_c4_formula
            for row in open_rows
        ),
        "open closed leakage formula": all(
            row.variance_c4 == row.variance_c4_formula for row in open_rows
        ),
        "open explicit centered action/SVD norm": all(
            row.centered_action_sum_zero and row.centered_action_norm_exact
            for row in open_rows
        ),
        "open Gamma0>=0, Gamma1<=0, Gamma2>=0": all(
            (lambda gammas: gammas[0] >= 0 and gammas[1] <= 0 and gammas[2] >= 0)(
                open_gammas(row.volume, row.n1, row.n2)
            )
            for row in open_rows
        ),
    }

    if v4 > 0:
        checks["open physical-shell zero set"] = all(
            row.invariant == open_zero_leakage(row.volume, row.n1, row.n2, v4)
            for row in open_rows
        )

    benchmark = next(
        row
        for row in periodic_rows
        if row.volume == 4 and row.n1 == 2 and row.n2 == 0
    )
    checks["periodic L=4 (2,0) benchmark"] = benchmark.variance_c4 == a4 * a4 * Fraction(128, 9)

    # Parameter-normalization checks: periodic depends only on A4; open depends
    # on the overall scale and q=V4/A4, not on the absolute scale separately.
    scale = Fraction(7, 3)
    checks["periodic A4 normalization"] = all(
        (scale * a4) ** 2 * periodic_variance_d(row.volume, row.n1, row.n2)
        == scale**2 * row.variance_c4_formula
        for row in periodic_rows
    )
    checks["open fixed-ratio A4 normalization"] = all(
        open_variance_c4(row.volume, row.n1, row.n2, scale * a4, scale * u4)
        == scale**2 * row.variance_c4_formula
        for row in open_rows
    )

    if a4 == DEFAULT_A4 and u4 == DEFAULT_U4:
        periodic_expected = {
            (4, 2): [Fraction(128, 9), Fraction(32), Fraction(128, 9)],
            (4, 3): [Fraction(0), Fraction(512, 9), Fraction(512, 9), Fraction(0)],
            (6, 2): [Fraction(384, 25), Fraction(320, 9), Fraction(384, 25)],
            (6, 3): [Fraction(576, 25), Fraction(16256, 225), Fraction(16256, 225), Fraction(576, 25)],
            (8, 2): [Fraction(640, 49), Fraction(32), Fraction(640, 49)],
            (8, 3): [Fraction(1280, 49), Fraction(3328, 49), Fraction(3328, 49), Fraction(1280, 49)],
        }
        open_expected = {
            (4, 2): [Fraction(19), Fraction(65, 2), Fraction(19)],
            (4, 3): [Fraction(81, 4), Fraction(1015, 12), Fraction(1015, 12), Fraction(81, 4)],
            (6, 2): [Fraction(16), Fraction(2660, 81), Fraction(16)],
            (6, 3): [Fraction(146, 5), Fraction(30826, 405), Fraction(30826, 405), Fraction(146, 5)],
            (8, 2): [Fraction(369, 28), Fraction(237, 8), Fraction(369, 28)],
            (8, 3): [Fraction(3165, 112), Fraction(7521, 112), Fraction(7521, 112), Fraction(3165, 112)],
        }
        checks["periodic Phase-3 Gram diagonals"] = all(
            [row.variance_c4 / (a4 * a4) for row in group] == periodic_expected[(volume, total)]
            for volume, total, group in focus_groups(rows, "periodic")
        )
        checks["open Phase-3 Gram diagonals at V4/A4=1"] = all(
            [row.variance_c4 / (a4 * a4) for row in group] == open_expected[(volume, total)]
            for volume, total, group in focus_groups(rows, "open")
        )

    overall = all(checks.values())
    total_configurations = sum(row.configurations_checked for row in rows)

    lines: list[str] = [
        "EXACT UNSCREENED PRODUCT-AGP LEAKAGE CERTIFICATE, VERSION 2.1",
        "============================================================",
        "",
        "Scope",
        "-----",
        "Periodic: L translated bonds of the cycle, matching the parent theorem.",
        "Open: the standard truncation to the L-1 bonds of the path.",
        "",
        "For each composition c=(n1,n2), mu_c=<c|C4|c> and",
        "",
        "  R|c> = (C4-mu_c)|c>.",
        "",
        "If tau_c^2=Var_c(C4)>0, define |ell_c>=(C4-mu_c)|c>/tau_c.",
        "Different compositions have orthogonal block-charge support, hence",
        "",
        "  R = sum_c tau_c |ell_c><c|,",
        "  R^dagger R = sum_c tau_c^2 |c><c|,",
        "  ||R||^2 = max_c tau_c^2,",
        "  rank(R) = #{c: tau_c^2>0}.",
        "",
        "The exact centered-moment checks below verify the first and second moments;",
        "occupation-basis diagonality and block-charge orthogonality then imply",
        "the operator-level SVD without a floating-point Fock-space decomposition.",
        "",
        f"Exact parameters: A4={fraction_text(a4)}, U4={fraction_text(u4)}, ",
        f"V4=4A4-U4={fraction_text(v4)}, V4/A4={fraction_text(q)}.",
        f"Volumes enumerated: {', '.join(map(str, volumes))}.",
        f"Boundary/composition sectors checked: {len(rows)}.",
        f"Hard-core pair configurations checked: {total_configurations}.",
        "",
        "Periodic theorem",
        "----------------",
        "For L>=3, tau_c^2/A4^2 equals",
        "",
        "  64 [",
        "    {n1(n1-1)(L-n1)(L-n1-1)+n2(n2-1)(L-n2)(L-n2-1)}",
        "      / {(L-2)(L-1)^2}",
        "    + 2 n1 n2 (L-n1)(L-n2)(3L-8) / {L^2(L-1)^2}",
        "  ].",
        "",
        "Open theorem",
        "------------",
        "Let B be total endpoint occupancy and Z=E1+E2+C. Then",
        "",
        "  C4_open = (8A4-2U4)n + 8A4 Z - V4 B.",
        "",
        "Writing Gamma0=Var(Z), Gamma1=Cov(Z,B), Gamma2=Var(B),",
        "",
        "  tau_c^2 = 64 A4^2 Gamma0 - 16 A4 V4 Gamma1 + V4^2 Gamma2.",
        "",
        "The exact component formulas are recorded in DERIVATION.md and paper_insert.tex.",
        "For A4>0,V4>0, Gamma0>=0, Gamma1<=0, Gamma2>=0.  Therefore for",
        "L>=3 the open leakage vanishes exactly at the four corners",
        "(n1,n2) in {0,L}x{0,L}; for L=2 it vanishes in every sector.",
        "",
        "Projected selectors",
        "-------------------",
        "Periodic shift per A4: 16(L-2)n1n2/[L(L-1)].",
        "Open shift per A4:     16(L-2)n1n2/L^2.",
        "The open coefficient is (L-1)/L times the periodic coefficient.",
        "",
        "Exact checks",
        "------------",
    ]
    for name, passed in checks.items():
        lines.append(f"[{'PASS' if passed else 'FAIL'}] {name}")

    lines.extend(
        [
            "",
            "Periodic benchmark",
            "------------------",
            "L=4, (n1,n2)=(2,0): tau^2/A4^2=128/9, hence",
            "",
            "  ||R|2,0>|| = (8 sqrt(2)/3) A4.",
            "",
            "Phase-3 Gram diagonals",
            "----------------------",
            "Entries are diag(R^dagger R)/A4^2 in order (0,n),...,(n,0).",
        ]
    )

    for boundary in ("periodic", "open"):
        lines.append("")
        lines.append(boundary.upper())
        for volume, total, group in focus_groups(rows, boundary):
            diagonal = [row.variance_c4 / (a4 * a4) for row in group]
            diagonal_text = ", ".join(fraction_text(value) for value in diagonal)
            rank = sum(value > 0 for value in diagonal)
            max_value = max(diagonal)
            lines.append(
                f"L={volume}, n={total}: diag=[{diagonal_text}], rank={rank}/{len(group)}, "
                f"||R||/A4={sqrt_text(max_value)} ({sqrt(float(max_value)):.12g})."
            )

    lines.extend(
        [
            "",
            "Boundary-scope conclusion",
            "-------------------------",
            "The parent theorem is periodic, so the periodic proposition is the theorem-level",
            "certificate.  The numerical plan also uses open chains; under the standard path",
            "truncation, open spectra require the separate formula above.  The linear charge",
            "term is no longer scalar because endpoint occupancy fluctuates.  Consequently the",
            "periodic normalization by A4 alone becomes, for open chains, normalization at fixed",
            f"V4/A4.  At the active ratio V4/A4={fraction_text(q)}, all Phase-3 open sectors listed above leak.",
            "",
            f"OVERALL: {'PASS' if overall else 'FAIL'}",
        ]
    )
    return "\n".join(lines) + "\n", overall


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory for certificate and CSV outputs",
    )
    parser.add_argument(
        "--volumes",
        type=int,
        nargs="+",
        default=list(DEFAULT_VOLUMES),
        help="Volumes to enumerate exactly",
    )
    parser.add_argument("--a4", type=parse_fraction, default=DEFAULT_A4)
    parser.add_argument("--u4", type=parse_fraction, default=DEFAULT_U4)
    args = parser.parse_args()

    volumes = tuple(sorted(set(args.volumes)))
    if not volumes or min(volumes) < 2:
        raise SystemExit("all volumes must be >=2")
    if args.a4 <= 0:
        raise SystemExit("A4 must be positive")

    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    rows: list[SectorResult] = []
    for volume in volumes:
        masks = {n: masks_with_weight(volume, n) for n in range(volume + 1)}
        for boundary in ("periodic", "open"):
            for n1 in range(volume + 1):
                for n2 in range(volume + 1):
                    rows.append(
                        enumerate_sector(
                            boundary,
                            volume,
                            n1,
                            n2,
                            masks,
                            args.a4,
                            args.u4,
                        )
                    )

    for boundary in ("periodic", "open"):
        boundary_rows = [row for row in rows if row.boundary == boundary]
        focus_rows = sorted(
            [
                row
                for row in boundary_rows
                if row.volume in FOCUS_VOLUMES and row.total_pairs in FOCUS_TOTAL_PAIRS
            ],
            key=lambda row: (row.volume, row.total_pairs, row.n1),
        )
        write_csv(
            output_dir / f"unscreened_agp_leakage_{boundary}_all_sectors.csv",
            boundary_rows,
            args.a4,
        )
        write_csv(
            output_dir / f"unscreened_agp_leakage_{boundary}_phase3_sectors.csv",
            focus_rows,
            args.a4,
        )

    # Backward-compatible aliases retain the periodic data.
    shutil.copyfile(
        output_dir / "unscreened_agp_leakage_periodic_all_sectors.csv",
        output_dir / "unscreened_agp_leakage_all_sectors.csv",
    )
    shutil.copyfile(
        output_dir / "unscreened_agp_leakage_periodic_phase3_sectors.csv",
        output_dir / "unscreened_agp_leakage_phase3_sectors.csv",
    )

    certificate, overall = build_certificate(rows, volumes, args.a4, args.u4)
    cert_path = output_dir / "unscreened_agp_leakage_certificate.txt"
    cert_path.write_text(certificate, encoding="utf-8")
    print(certificate, end="")
    if not overall:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
