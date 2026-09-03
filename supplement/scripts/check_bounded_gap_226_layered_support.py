#!/usr/bin/env python3
"""Exact support and partition audit for the layered support geometry.

The inequalities below depend on the support geometry, not on the ambient
Maynard dimension.  The command-line binding keeps that distinction explicit
when the same support is used by the k=47 and k=46 finite candidates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import z3

import check_bounded_gap_236_type_iic as partition_audit
from build_bounded_gap_226_layered_rational_candidate import (
    B_BY_COUNT,
    DELTA as DEFAULT_DELTA,
    K as DEFAULT_K,
    MAX_BIG_COUNT as DEFAULT_MAX_BIG_COUNT,
    OMEGA,
    RADIUS as DEFAULT_RADIUS,
    candidate_digest as support_source_digest,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT / "outputs" / "bounded_gap_226_layered_support_exact.json"
)

ARITHMETIC_EPSILON = Fraction(1, 10_000_000_000)
XI_1 = Fraction(19, 50)
XI_2 = Fraction(2, 5)
XI_3 = Fraction(2, 5)
TAIL_BOUND = Fraction(219, 1000)
DEFAULT_B_BY_COUNT = dict(B_BY_COUNT)


def encode(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def parse_fraction(value: str) -> Fraction:
    return Fraction(value)


def geometry_digest(
    radius: Fraction,
    delta: Fraction,
    omega: Fraction,
    a_one: Fraction,
    support_epsilon: Fraction,
    outer_radius: Fraction,
    maximum_big_count: int,
) -> str:
    if (
        radius == DEFAULT_RADIUS
        and delta == DEFAULT_DELTA
        and omega == OMEGA
        and maximum_big_count == DEFAULT_MAX_BIG_COUNT
        and B_BY_COUNT == DEFAULT_B_BY_COUNT
    ):
        return support_source_digest()
    binding = {
        "radius": encode(radius),
        "delta": encode(delta),
        "omega": encode(omega),
        "aOne": encode(a_one),
        "supportEpsilon": encode(support_epsilon),
        "outerRadius": encode(outer_radius),
        "maximumBigCoordinateCount": maximum_big_count,
        "implicitTailBound": encode(TAIL_BOUND),
        "bByCount": {
            str(count): encode(bound)
            for count, bound in sorted(B_BY_COUNT.items())
        },
    }
    encoded = json.dumps(
        binding, sort_keys=True, separators=(",", ":")
    ).encode("ascii")
    return hashlib.sha256(encoded).hexdigest()


def resolve_binding(
    k: int,
    candidate_json: Path | None,
    geometry_sha: str,
    outer_radius: Fraction,
) -> tuple[str, str]:
    """Return the candidate hash and the geometry-source hash.

    The historical geometry source is the frozen k=47 layered candidate.  A
    later candidate may bind to that exact geometry through supportDigest.
    """

    if candidate_json is None:
        return geometry_sha, geometry_sha

    candidate = json.loads(candidate_json.read_text(encoding="utf-8"))
    binding = candidate.get("binding")
    if not isinstance(binding, dict):
        raise ValueError("candidate JSON has no binding object")
    if int(binding.get("k", -1)) != k:
        raise ValueError("candidate dimension does not match --k")
    if Fraction(str(binding.get("outerRadius"))) != outer_radius:
        raise ValueError("candidate outer radius does not match support geometry")
    if binding.get("supportDigest") != geometry_sha:
        raise ValueError("candidate does not bind to the audited support source")
    candidate_sha = candidate.get("sha256")
    if not isinstance(candidate_sha, str):
        raise ValueError("candidate JSON has no sha256 binding hash")
    return candidate_sha, geometry_sha


def bound_for_count(group_size: int) -> Fraction:
    return B_BY_COUNT.get(group_size, TAIL_BOUND)


def group_bound(group_size: int, _: Fraction) -> Fraction:
    return bound_for_count(group_size)


def audit_fixed_capacities(
    first_size: int,
    second_size: int,
    capacities: list[Fraction],
    delta: Fraction,
    max_iterations: int = 100_000,
) -> dict[str, Any]:
    count = first_size + second_size
    coordinates = z3.Reals(
        " ".join(
            f"fixed_{len(capacities)}_{first_size}_{second_size}_{index}"
            for index in range(count)
        )
    )
    solver = z3.Solver()
    for coordinate in coordinates:
        solver.add(
            coordinate >= partition_audit.z3_rational(delta),
            coordinate <= 1,
        )
    if first_size:
        solver.add(
            z3.Sum(coordinates[:first_size])
            <= partition_audit.z3_rational(bound_for_count(first_size))
        )
    if second_size:
        solver.add(
            z3.Sum(coordinates[first_size:])
            <= partition_audit.z3_rational(bound_for_count(second_size))
        )

    covering_partitions: list[tuple[int, ...]] = []
    for iteration in range(max_iterations):
        status = solver.check()
        if status == z3.unsat:
            return {
                "status": "covered",
                "iterations": iteration,
                "coveringPartitions": [
                    list(item) for item in covering_partitions
                ],
            }
        if status != z3.sat:
            return {"status": "unknown", "iterations": iteration}

        model = solver.model()
        weights = [
            partition_audit.as_fraction(
                model.eval(coordinate, model_completion=True)
            )
            for coordinate in coordinates
        ]
        assignment = partition_audit.exact_pack(weights, capacities)
        if assignment is None:
            return {
                "status": "counterexample",
                "iterations": iteration,
                "weights": [encode(item) for item in weights],
            }

        covering_partitions.append(assignment)
        failures = []
        for bin_index, capacity in enumerate(capacities):
            terms = [
                coordinates[index]
                for index, assigned_bin in enumerate(assignment)
                if assigned_bin == bin_index
            ]
            total = (
                z3.Sum(terms)
                if terms
                else partition_audit.z3_rational(Fraction(0))
            )
            failures.append(total > partition_audit.z3_rational(capacity))
        solver.add(z3.Or(failures))

    return {"status": "iteration-limit", "iterations": max_iterations}


def audit_all_pairs(
    audit,
    maximum_big_count: int,
) -> dict[str, Any]:
    pairs: dict[str, Any] = {}
    for first_size in range(maximum_big_count + 1):
        for second_size in range(
            max(first_size, 1), maximum_big_count + 1
        ):
            key = f"{first_size},{second_size}"
            pairs[key] = audit(first_size, second_size)
    return {
        "allPairsCovered": all(
            result["status"] == "covered" for result in pairs.values()
        ),
        "totalCoveringPartitions": sum(
            len(result.get("coveringPartitions", []))
            for result in pairs.values()
        ),
        "pairs": pairs,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=DEFAULT_K)
    parser.add_argument("--omega", type=parse_fraction, default=OMEGA)
    parser.add_argument("--delta", type=parse_fraction, default=DEFAULT_DELTA)
    parser.add_argument("--radius", type=parse_fraction, default=DEFAULT_RADIUS)
    parser.add_argument("--b1", type=parse_fraction)
    parser.add_argument("--b2", type=parse_fraction)
    parser.add_argument("--b3", type=parse_fraction)
    parser.add_argument("--candidate-json", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.b1 is not None:
        B_BY_COUNT[1] = args.b1
    if args.b2 is not None:
        B_BY_COUNT[2] = args.b2
    if args.b3 is not None:
        B_BY_COUNT[3] = args.b3
    omega = args.omega
    delta = args.delta
    radius = args.radius
    a_one = Fraction(1, 4) + omega
    support_epsilon = radius - a_one
    outer_radius = a_one - support_epsilon
    if delta <= 0 or support_epsilon <= 0 or outer_radius <= 0:
        raise ValueError("omega gives a nonpositive support margin or J radius")
    maximum_geometric_count = int(radius / delta)
    maximum_big_count = max(
        (
            count
            for count in range(1, maximum_geometric_count + 1)
            if count * delta <= bound_for_count(count)
        ),
        default=0,
    )
    geometry_sha = geometry_digest(
        radius,
        delta,
        omega,
        a_one,
        support_epsilon,
        outer_radius,
        maximum_big_count,
    )
    candidate_sha, geometry_sha = resolve_binding(
        args.k, args.candidate_json, geometry_sha, outer_radius
    )

    partition_audit.DELTA = delta
    partition_audit.OMEGA = omega
    partition_audit.group_bound = group_bound

    fixed_conditions = {
        "typeI": [
            XI_1 - 2 * ARITHMETIC_EPSILON,
            Fraction(1, 6) - 4 * omega - 2 * ARITHMETIC_EPSILON,
        ],
        "typeIIa": [
            Fraction(2, 5)
            + Fraction(24, 5) * omega
            + Fraction(7, 5) * delta
            - 2 * ARITHMETIC_EPSILON,
            Fraction(1, 14)
            - Fraction(24, 7) * omega
            - 2 * ARITHMETIC_EPSILON,
        ],
        "typeIIb": [
            Fraction(1, 3)
            + 8 * omega
            + Fraction(7, 3) * delta
            - 4 * ARITHMETIC_EPSILON,
            Fraction(1, 10)
            - Fraction(34, 5) * omega
            - Fraction(7, 5) * delta
            - 4 * ARITHMETIC_EPSILON,
            Fraction(1, 35)
            + Fraction(22, 35) * omega
            + Fraction(21, 35) * delta
            - 4 * ARITHMETIC_EPSILON,
        ],
        "typeIII": [
            1
            - 6 * omega
            - Fraction(3, 2) * XI_3
            - 2 * ARITHMETIC_EPSILON,
            Fraction(5, 2) * omega
            + Fraction(3, 8) * XI_3
            - 2 * ARITHMETIC_EPSILON,
        ],
    }
    fixed_audits = {
        name: {
            "capacities": [encode(value) for value in capacities],
            "audit": audit_all_pairs(
                lambda first, second, caps=capacities: audit_fixed_capacities(
                    first, second, caps, delta
                ),
                maximum_big_count,
            ),
        }
        for name, capacities in fixed_conditions.items()
    }
    type_iic = audit_all_pairs(
        lambda first, second: partition_audit.audit_pair(
            first, second, Fraction(0)
        ),
        maximum_big_count,
    )
    endpoint = audit_all_pairs(
        lambda first, second: partition_audit.audit_endpoint_pair(
            first, second, Fraction(0)
        ),
        maximum_big_count,
    )

    base_margins = {
        "typeIFirst": (
            XI_1
            - 4 * a_one
            + Fraction(2, 3)
            - 2 * ARITHMETIC_EPSILON
            - delta
        ),
        "typeISecond": (
            Fraction(9, 7)
            - Fraction(34, 7) * a_one
            - 2 * ARITHMETIC_EPSILON
            - delta
        ),
        "typeIIZero": (
            Fraction(19, 2)
            - 36 * a_one
            - 13 * delta
            + 100 * ARITHMETIC_EPSILON
        ),
        "typeIIFirst": (
            XI_2 / 10
            - Fraction(32, 10) * a_one
            + Fraction(8, 10)
            - 2 * ARITHMETIC_EPSILON
            - delta
        ),
        "typeIISecond": (
            XI_2 / 4
            + Fraction(11, 16)
            - 3 * a_one
            - 2 * ARITHMETIC_EPSILON
            - delta
        ),
        "typeIII": (
            Fraction(11, 8)
            - Fraction(7, 2) * a_one
            - Fraction(9, 8) * XI_3
            - 2 * ARITHMETIC_EPSILON
            - delta
        ),
    }
    harman_margins = {
        "twoXi1PlusThreeXi2BelowTwo": 2 - 2 * XI_1 - 3 * XI_2,
        "xi3AtLeastXi2": XI_3 - XI_2,
        "xi1PlusNineXi2BelowFour": 4 - XI_1 - 9 * XI_2,
        "twoXi1PlusXi2AboveOne": 2 * XI_1 + XI_2 - 1,
        "seventeenXi2BelowSeven": 7 - 17 * XI_2,
        "roughnessBetaAboveB1": 1 - 2 * XI_2 - B_BY_COUNT[1],
    }
    structural_checks = {
        "radiusIdentity": a_one + support_epsilon == radius,
        "outerRadiusIdentity": a_one - support_epsilon == outer_radius,
        "allBAboveDelta": all(
            delta < bound_for_count(count)
            for count in range(1, int(1 / delta) + 1)
        ),
        "bSequenceNondecreasing": all(
            bound_for_count(count) <= bound_for_count(count + 1)
            for count in range(1, int(1 / delta))
        ),
        "successiveBIncrementsAtMostDelta": all(
            bound_for_count(count + 1) <= bound_for_count(count) + delta
            for count in range(1, int(1 / delta))
        ),
        "maximumBigCoordinateCountIsPossible": (
            maximum_big_count * delta
            <= bound_for_count(maximum_big_count)
        ),
        "nextBigCoordinateCountIsImpossible": (
            (maximum_big_count + 1) * delta
            > bound_for_count(maximum_big_count + 1)
        ),
        "allBaseMarginsPositive": all(
            value > 0 for value in base_margins.values()
        ),
        "allHarmanMarginsNonnegative": all(
            value >= 0 for value in harman_margins.values()
        ),
        "primeMinorantIsPrimeIndicator": XI_2 <= Fraction(2, 5),
    }
    all_partition_audits = [
        item["audit"] for item in fixed_audits.values()
    ] + [type_iic, endpoint]
    structural_checks["allPartitionFamiliesCovered"] = all(
        audit["allPairsCovered"] for audit in all_partition_audits
    )

    payload = {
        "name": f"Exact support audit for the layered k={args.k} candidate",
        "candidateSha256": candidate_sha,
        "supportSourceSha256": geometry_sha,
        "claimBoundary": (
            "This certifies the finite rational support and partition "
            "conditions. The corrected deep Type-IIc equidistribution theorem "
            "S52 is proved separately in the manuscript; this finite audit "
            "does not replace that analytic proof."
        ),
        "parameters": {
            "k": args.k,
            "radius": encode(radius),
            "outerRadius": encode(outer_radius),
            "delta": encode(delta),
            "omega": encode(omega),
            "aOne": encode(a_one),
            "supportEpsilon": encode(support_epsilon),
            "maximumBigCoordinateCount": maximum_big_count,
            "xi1": encode(XI_1),
            "xi2": encode(XI_2),
            "xi3": encode(XI_3),
            "bByCount": {
                str(count): encode(bound)
                for count, bound in sorted(B_BY_COUNT.items())
            },
        },
        "structuralChecks": structural_checks,
        "baseMargins": {key: encode(value) for key, value in base_margins.items()},
        "harmanMargins": {
            key: encode(value) for key, value in harman_margins.items()
        },
        "fixedPartitionConditions": fixed_audits,
        "typeIIcNonnegativeOmega": type_iic,
        "typeIIcRetainedLossEndpoint": {
            "eta": encode(partition_audit.ENDPOINT_ETA),
            "lossConstant": partition_audit.ENDPOINT_LOSS_CONSTANT,
            "audit": endpoint,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "output": str(args.output),
                "candidateSha256": candidate_sha,
                "supportSourceSha256": geometry_sha,
                "checks": structural_checks,
                "partitionCounts": {
                    **{
                        name: item["audit"]["totalCoveringPartitions"]
                        for name, item in fixed_audits.items()
                    },
                    "typeIIc": type_iic["totalCoveringPartitions"],
                    "endpoint": endpoint["totalCoveringPartitions"],
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
