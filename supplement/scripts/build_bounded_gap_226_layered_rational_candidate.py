#!/usr/bin/env python3
"""Freeze the exact layered-support k=47 candidate.

The radial coefficients are a dyadic quantization of a 90-percent
interpolation toward the generalized eigenvector found on the 1/40000 mesh.
The floating calculation is provenance only; positivity and all subsequent
integrals are checked with exact rational arithmetic.
"""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path

from build_bounded_gap_236_rational_candidate import (
    COMMON_DENOMINATOR,
    PROFILE_NUMERATORS,
    bernstein_positivity_certificate,
    combine_basis,
    encode,
    shifted_chebyshev_basis,
    shifted_legendre_basis,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "bounded_gap_226_layered_rational_candidate.json"

K = 47
RADIUS = Fraction(521, 2000)
DELTA = Fraction(1, 50)
OMEGA = Fraction(57, 10000)
A_ONE = Fraction(1, 4) + OMEGA
SUPPORT_EPSILON = RADIUS - A_ONE
OUTER_RADIUS = A_ONE - SUPPORT_EPSILON
MAX_BIG_COUNT = 10

_ACTIVE_BOUNDS = [157, 161, 179, 179, 188, 199, 203, 208, 208, 211]
B_BY_COUNT = {
    count: Fraction(value, 1000)
    for count, value in enumerate(_ACTIVE_BOUNDS, start=1)
}
for count in range(MAX_BIG_COUNT + 1, int(1 / DELTA) + 1):
    B_BY_COUNT[count] = Fraction(219, 1000)

RADIAL_NUMERATORS = [
    17204829790367296,
    804749628808710,
    7049321079942884,
    -17380069066313472,
    8740170260751930,
    -6919107880673811,
    2084564111565619,
    -1090902544612000,
    -442386526662123,
]


def candidate_binding() -> dict[str, object]:
    return {
        "k": K,
        "denominator": COMMON_DENOMINATOR,
        "profileNumerators": PROFILE_NUMERATORS,
        "radialNumerators": RADIAL_NUMERATORS,
        "support": {
            "radius": encode(RADIUS),
            "outerRadius": encode(OUTER_RADIUS),
            "delta": encode(DELTA),
            "omega": encode(OMEGA),
            "aOne": encode(A_ONE),
            "supportEpsilon": encode(SUPPORT_EPSILON),
            "maximumBigCoordinateCount": MAX_BIG_COUNT,
            "bByCount": {
                str(count): encode(bound)
                for count, bound in sorted(B_BY_COUNT.items())
            },
        },
    }


def candidate_digest() -> str:
    canonical = json.dumps(
        candidate_binding(), sort_keys=True, separators=(",", ":")
    ).encode("ascii")
    return hashlib.sha256(canonical).hexdigest()


def main() -> None:
    profile_power = combine_basis(
        PROFILE_NUMERATORS, shifted_chebyshev_basis(32)
    )
    radial_power = combine_basis(
        RADIAL_NUMERATORS, shifted_legendre_basis(8)
    )
    profile_positivity = bernstein_positivity_certificate(profile_power, 32)
    radial_positivity = bernstein_positivity_certificate(radial_power, 8)

    payload = {
        "name": "Exact layered-support product-radial candidate for k=47",
        "claimBoundary": (
            "This artifact fixes exact coefficients and proves their unit-"
            "interval positivity. The variational inequality and support "
            "partition conditions are certified by companion artifacts."
        ),
        "origin": {
            "sourceMesh": "1/40000",
            "interpolationWeightTowardMeshEigenvector": "9/10",
            "quantizationDenominator": COMMON_DENOMINATOR,
            "floatingScoreAtSourceMesh": 1.0044844789849507,
            "note": (
                "The numerical optimization selected the coefficients but is "
                "not a premise of the exact certificate."
            ),
        },
        "binding": candidate_binding(),
        "sha256": candidate_digest(),
        "profile": {
            "basis": "shifted Chebyshev T_j(2u-1)",
            "degree": 32,
            "numerators": PROFILE_NUMERATORS,
        },
        "radial": {
            "basis": "shifted Legendre P_j(2u-1)",
            "degree": 8,
            "numerators": RADIAL_NUMERATORS,
            "coefficients": [
                encode(Fraction(value, COMMON_DENOMINATOR))
                for value in RADIAL_NUMERATORS
            ],
            "powerBasisCoefficients": [encode(value) for value in radial_power],
        },
        "exactBernsteinPositivityCertificates": {
            "profileOnUnitInterval": profile_positivity,
            "radialOnUnitInterval": radial_positivity,
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(OUTPUT),
                "candidateSha256": candidate_digest(),
                "profilePositive": profile_positivity[
                    "allLeafBernsteinCoefficientsStrictlyPositive"
                ],
                "radialPositive": radial_positivity[
                    "allLeafBernsteinCoefficientsStrictlyPositive"
                ],
                "radialCertifiedLowerBound": radial_positivity[
                    "globalCertifiedLowerBound"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
