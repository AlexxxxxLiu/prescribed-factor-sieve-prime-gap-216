#!/usr/bin/env python3
"""Freeze the 208-dimensional k=46 address-state vector over dyadic rationals."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "outputs"
    / "bounded_gap_216_address_state_order4_joint4_d4_supportfixed_1.25e-5.json"
)
OUTPUT = (
    ROOT
    / "outputs"
    / "bounded_gap_216_address_state_rational_candidate.json"
)
DENOMINATOR = 1 << 50
EXACT_OUTER_RADIUS = "2509/10000"


def canonical_digest(value: object) -> str:
    encoded = json.dumps(
        value, sort_keys=True, separators=(",", ":")
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def main() -> None:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    method = source["method"]
    optimum = source["unconstrainedOptimum"]
    assert isinstance(method, dict)
    assert isinstance(optimum, dict)
    coefficients = [
        float(value) for value in optimum["shiftedLegendreCoefficients"]
    ]
    numerators = [round(value * DENOMINATOR) for value in coefficients]
    quantized = [value / DENOMINATOR for value in numerators]
    errors = [
        abs(original - rounded)
        for original, rounded in zip(coefficients, quantized)
    ]

    binding = {
        "k": int(method["k"]),
        "jDomain": (
            "conservative J truncated at R_o=" + EXACT_OUTER_RADIUS
        ),
        "outerRadius": EXACT_OUTER_RADIUS,
        "supportDigest": method["supportDigest"],
        "radialDegree": 4,
        "addressMassOrder": int(method["addressMassOrder"]),
        "layerDegrees": method["layerDegrees"],
        "massLayerDegrees": method["massLayerDegrees"],
        "dispersionLayerDegrees": method["dispersionLayerDegrees"],
        "coefficientOrder": (
            "layer, then ordinary address mode 0..4, then dispersion mode; "
            "within each active mode, increasing shifted-Legendre degree"
        ),
        "commonDenominator": DENOMINATOR,
        "numerators": numerators,
    }
    payload = {
        "name": "Frozen rational k=46 prime-address moment-state candidate",
        "claimBoundary": (
            "This artifact freezes exact coefficients and metadata only. It "
            "does not certify c^T(J-I)c>0 or any bounded-gap theorem."
        ),
        "origin": {
            "source": str(SOURCE),
            "sourceMesh": method["mesh"],
            "floatingSourceScore46JOverI": optimum["score46JOverI"],
            "floatingFixedVectorCubicDiagnostic": 1.0000266540690037,
            "metadataCorrection": (
                "The source optimizer used outerRadius=0.2509; its legacy "
                "jDomain prose incorrectly named 491/2000. The exact binding "
                "records the value actually used, 2509/10000."
            ),
        },
        "binding": binding,
        "sha256": canonical_digest(binding),
        "quantization": {
            "commonDenominator": DENOMINATOR,
            "coefficientCount": len(numerators),
            "maximumAbsoluteCoefficientError": max(errors),
            "l2CoefficientError": sum(error * error for error in errors)
            ** 0.5,
        },
    }
    if len(numerators) != int(method["matrixDimension"]):
        raise RuntimeError("coefficient count does not match matrix dimension")
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(OUTPUT),
                "sha256": payload["sha256"],
                "coefficientCount": len(numerators),
                "maximumAbsoluteCoefficientError": max(errors),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
