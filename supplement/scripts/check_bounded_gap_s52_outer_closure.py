#!/usr/bin/env python3
"""Exact monomial regression for the S52 outer-inheritance closure."""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "bounded_gap_s52_outer_closure.json"


def add(left: dict[str, Fraction], right: dict[str, Fraction]) -> dict[str, Fraction]:
    keys = set(left) | set(right)
    return {key: left.get(key, Fraction(0)) + right.get(key, Fraction(0)) for key in keys}


def scale(values: dict[str, Fraction], factor: Fraction) -> dict[str, Fraction]:
    return {key: value * factor for key, value in values.items()}


def prune(values: dict[str, Fraction]) -> dict[str, Fraction]:
    return {key: value for key, value in values.items() if value}


def substitute_uv(values: dict[str, Fraction]) -> dict[str, Fraction]:
    result = dict(values)
    power = result.pop("U", Fraction(0))
    if result.pop("V", Fraction(0)) != power:
        raise AssertionError("U and V powers do not form a complete UV factor")
    result["Q"] = result.get("Q", Fraction(0)) + power
    result["q0"] = result.get("q0", Fraction(0)) - power
    return prune(result)


def clean(values: dict[str, Fraction]) -> dict[str, str]:
    return {
        key: f"{value.numerator}/{value.denominator}"
        for key, value in sorted(values.items())
        if value
    }


def main() -> None:
    first_cauchy = {
        "g": Fraction(1),
        "R": Fraction(1),
        "Q": Fraction(1),
        "U": Fraction(1),
        "N": Fraction(1),
        "q0": Fraction(-2),
    }
    duplicated = {
        "g": Fraction(1),
        "R": Fraction(1),
        "Q": Fraction(1),
        "N": Fraction(1),
        "U": Fraction(1),
        "V": Fraction(2),
        "eta": Fraction(-4),
    }
    localized = scale(add(first_cauchy, duplicated), Fraction(1, 2))
    after_factor_interface = substitute_uv(localized)
    expected_inner = {
        "g": Fraction(1),
        "R": Fraction(1),
        "Q": Fraction(2),
        "N": Fraction(1),
        "q0": Fraction(-2),
        "eta": Fraction(-2),
    }
    outer_weight = {
        "q0": Fraction(1),
        "M": Fraction(1),
        "R": Fraction(-1),
        "Q": Fraction(-2),
    }
    before_sums = add(after_factor_interface, outer_weight)
    ell_count = {"N": Fraction(1), "R": Fraction(-1)}
    after_ell = add(before_sums, ell_count)
    expected_after_ell = {
        "gOverQ0Summed": Fraction(1),
        "M": Fraction(1),
        "N": Fraction(2),
        "R": Fraction(-1),
        "eta": Fraction(-2),
    }
    normalized_after_ell = dict(after_ell)
    gcd_power = normalized_after_ell.pop("g")
    q0_power = normalized_after_ell.pop("q0")
    if gcd_power != 1 or q0_power != -1:
        raise AssertionError("outer q0 sum is not gcd(q0,ell)/q0")
    normalized_after_ell["gOverQ0Summed"] = Fraction(1)
    normalized_after_ell = prune(normalized_after_ell)

    checks = {
        "cauchySquareRootCorrect": localized
        == {
            "g": Fraction(1),
            "R": Fraction(1),
            "Q": Fraction(1),
            "U": Fraction(1),
            "N": Fraction(1),
            "q0": Fraction(-1),
            "V": Fraction(1),
            "eta": Fraction(-2),
        },
        "factorInterfaceGivesQ0MinusTwo": after_factor_interface == expected_inner,
        "outerSummandIsGcdOverQ0": gcd_power == 1 and q0_power == -1,
        "ellCountGivesMN2OverR": normalized_after_ell == expected_after_ell,
        "twoHalfEtaLossesLeaveEtaSaving": Fraction(-2) + Fraction(1, 2) + Fraction(1, 2) == -1,
    }
    payload = {
        "name": "S52 outer-closure exact monomial regression",
        "claimBoundary": (
            "This checks the algebra of Proposition 4.1; it does not replace "
            "the analytic Cauchy, completion, or divisor-sum estimates."
        ),
        "localizedCauchyMonomial": clean(localized),
        "afterUVEqualsQOverQ0": clean(after_factor_interface),
        "afterOuterWeightAndEllCount": clean(normalized_after_ell),
        "checks": checks,
        "allChecksPass": all(checks.values()),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if not payload["allChecksPass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
