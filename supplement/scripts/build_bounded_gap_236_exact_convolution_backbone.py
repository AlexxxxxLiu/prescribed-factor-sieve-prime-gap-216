#!/usr/bin/env python3
"""Build the exact truncated-power convolution backbone for the k=48 candidate.

This script does not evaluate I or J and therefore does not prove H_1 <= 236.
It replaces the FFT density recursion by exact FLINT rational polynomials and
records reproducible hashes for the small and big coordinate convolutions that
the final certificate must use.
"""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

from flint import __version__ as flint_version
from flint import fmpq, fmpq_poly

from build_bounded_gap_236_rational_candidate import (
    PROFILE_NUMERATORS,
    combine_basis,
    shifted_chebyshev_basis,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "bounded_gap_236_exact_convolution_backbone.json"

K = 48
RADIUS = Fraction(521, 2000)
DELTA = Fraction(7, 250)
B_SMALL = Fraction(1599, 10000)
B_LARGE = Fraction(1839, 10000)


def q(value: Fraction) -> fmpq:
    return fmpq(value.numerator, value.denominator)


def to_flint_poly(coefficients: list[Fraction]) -> fmpq_poly:
    return fmpq_poly([q(value) for value in coefficients])


def borel(poly: fmpq_poly) -> fmpq_poly:
    """Map sum c_a*x^a to sum c_a*a!*u^(a+1)."""
    result = [fmpq(0)] * (len(poly) + 1)
    for degree in range(len(poly)):
        result[degree + 1] = poly[degree] * math.factorial(degree)
    return fmpq_poly(result)


def inverse_borel(poly: fmpq_poly) -> fmpq_poly:
    """Invert ``borel`` on polynomials with zero constant coefficient."""
    if poly[0] != 0:
        raise ValueError("inverse Borel input must have zero constant term")
    result = [
        poly[degree + 1] / math.factorial(degree)
        for degree in range(len(poly) - 1)
    ]
    return fmpq_poly(result)


def shifted(poly: fmpq_poly, offset: fmpq) -> fmpq_poly:
    return poly(fmpq_poly([offset, fmpq(1)]))


def polynomial_digest(poly: fmpq_poly) -> str:
    digest = hashlib.sha256()
    for coefficient in poly:
        digest.update(str(coefficient.numerator).encode("ascii"))
        digest.update(b"/")
        digest.update(str(coefficient.denominator).encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def encode(value: fmpq | Fraction) -> str:
    if isinstance(value, Fraction):
        return f"{value.numerator}/{value.denominator}"
    return f"{value.numerator}/{value.denominator}"


def evaluate_shifted_pieces(
    pieces: list[tuple[fmpq, fmpq_poly]], value: fmpq
) -> fmpq:
    total = fmpq(0)
    for knot, poly in pieces:
        if value >= knot:
            total += poly(value - knot)
    return total


def direct_two_fold_small(
    density: fmpq_poly, cutoff: fmpq, value: fmpq
) -> fmpq:
    lower = max(fmpq(0), value - cutoff)
    upper = min(cutoff, value)
    if upper <= lower:
        return fmpq(0)
    variable_reflection = density(fmpq_poly([value, fmpq(-1)]))
    antiderivative = (density * variable_reflection).integral()
    return antiderivative(upper) - antiderivative(lower)


def small_pieces(
    density: fmpq_poly,
    cutoff: fmpq,
    order: int,
    maximum_total: fmpq,
) -> list[tuple[fmpq, fmpq_poly]]:
    base = borel(density)
    upper_edge = borel(shifted(density, cutoff))
    maximum_shift = min(order, int(maximum_total / cutoff))
    pieces: list[tuple[fmpq, fmpq_poly]] = []
    for shifted_count in range(maximum_shift + 1):
        transform = (base ** (order - shifted_count)) * (
            upper_edge**shifted_count
        )
        coefficient = (-1) ** shifted_count * math.comb(
            order, shifted_count
        )
        pieces.append(
            (
                shifted_count * cutoff,
                coefficient * inverse_borel(transform),
            )
        )
    return pieces


def big_piece(
    density: fmpq_poly, cutoff: fmpq, order: int
) -> tuple[fmpq, fmpq_poly]:
    if order == 0:
        return fmpq(0), fmpq_poly([1])
    shifted_density = shifted(density, cutoff)
    return order * cutoff, inverse_borel(borel(shifted_density) ** order)


def piece_record(knot: fmpq, poly: fmpq_poly) -> dict[str, object]:
    return {
        "knot": encode(knot),
        "degree": poly.degree(),
        "coefficientCount": len(poly),
        "sha256": polynomial_digest(poly),
    }


def main() -> None:
    profile_power = combine_basis(
        PROFILE_NUMERATORS, shifted_chebyshev_basis(32)
    )
    profile = to_flint_poly(profile_power)
    density = profile * profile
    cutoff = q(DELTA / RADIUS)
    maximum_total = fmpq(1)

    small_2 = small_pieces(density, cutoff, 2, 2 * cutoff)
    sample_points = [cutoff / 2, 3 * cutoff / 2, 5 * cutoff / 2]
    direct_errors = []
    for value in sample_points:
        represented = evaluate_shifted_pieces(small_2, value)
        direct = direct_two_fold_small(density, cutoff, value)
        direct_errors.append(represented - direct)

    small_records: dict[str, object] = {}
    for order in (47, 48):
        pieces = small_pieces(density, cutoff, order, maximum_total)
        small_records[str(order)] = {
            "order": order,
            "maximumTotalRetained": "1/1",
            "maximumActiveShiftCount": len(pieces) - 1,
            "expectedPieceDegree": order * (density.degree() + 1) - 1,
            "pieces": [piece_record(knot, poly) for knot, poly in pieces],
        }

    big_records = []
    for order in range(1, 7):
        knot, poly = big_piece(density, cutoff, order)
        bound = B_SMALL if order <= 2 else B_LARGE
        normalized_bound = q(bound / RADIUS)
        antiderivative = poly.integral()
        retained_mass = (
            antiderivative(normalized_bound - knot) - antiderivative(0)
            if normalized_bound > knot
            else fmpq(0)
        )
        big_records.append(
            {
                "order": order,
                "groupSumBound": encode(bound / RADIUS),
                "retainedMass": encode(retained_mass),
                **piece_record(knot, poly),
            }
        )

    expected_degree_48 = 48 * (density.degree() + 1) - 1
    payload = {
        "name": "Exact FLINT convolution backbone for the k=48 candidate",
        "claimBoundary": (
            "Exact rational convolution data only. I, J, the final integer "
            "margin, and the analytic Type-IIc splice remain unproved."
        ),
        "backend": {"pythonFlintVersion": flint_version},
        "normalizedGeometry": {
            "radius": encode(RADIUS),
            "smallCutoffDeltaOverRadius": encode(DELTA / RADIUS),
            "smallGroupBoundOverRadius": encode(B_SMALL / RADIUS),
            "largeGroupBoundOverRadius": encode(B_LARGE / RADIUS),
        },
        "density": {
            "formula": "f(x)=q(x)^2",
            "degree": density.degree(),
            "coefficientCount": len(density),
            "sha256": polynomial_digest(density),
        },
        "exactChecks": {
            "borelConvolutionDegreeAtOrder48": (
                expected_degree_48 == 3119
            ),
            "twoFoldSmallRepresentationMatchesDirectIntegral": all(
                error == 0 for error in direct_errors
            ),
            "twoFoldSampleErrors": [encode(error) for error in direct_errors],
            "onlyTenSmallShiftsReachNormalizedTotalOne": (
                len(small_records["48"]["pieces"]) == 10
            ),
            "bigCoordinateOrdersStopAtSix": len(big_records) == 6,
        },
        "smallConvolutions": small_records,
        "bigConvolutions": big_records,
        "nextExactStage": (
            "Integrate the m=0..6 I cells and the AA, 2AB, BB J cells by "
            "exact moment recurrences, then verify 48*N_J*D_I-N_I*D_J>0."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(OUTPUT),
                "checks": payload["exactChecks"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
