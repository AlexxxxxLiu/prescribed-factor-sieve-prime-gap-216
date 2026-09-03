#!/usr/bin/env python3
"""Exact CEGAR audit of the Type-IIc support condition near H_1 <= 236.

The arithmetic is rational and the symbolic subproblems are quantifier-free
linear real arithmetic.  For each pair of nonempty coordinate groups, the
script alternates between:

1. asking Z3 for a point not covered by any partition found so far; and
2. exactly packing that rational point into the four Type-IIc capacities.

If no packing exists, the point is an exact counterexample.  If Z3 returns
UNSAT, the recorded finite family of partitions covers the whole parameter
polytope.  This audits the nonnegative-omega range used for moduli at or above
the square-root threshold.  The arXiv v1 proposition literally includes a
small negative-omega interval, where its fourth capacity is negative; that
literal statement cannot hold and is reported separately below.
"""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any

try:
    import z3
except ImportError as exc:  # pragma: no cover - environment guidance
    raise SystemExit(
        "z3-solver is required; run with PYTHONPATH=.python-deps or install it"
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "bounded_gap_236_type_iic_exact.json"

DELTA = Fraction(7, 250)
EPSILON = Fraction(1, 10_000_000_000)
XI_1 = Fraction(19, 50)
XI_2 = Fraction(2, 5)
OMEGA = Fraction(3, 1000)
B_SMALL = Fraction(1599, 10_000)
B_CANDIDATE = Fraction(1839, 10_000)
B_FAILURE = Fraction(23, 125)
MAX_GROUP_SIZE = 6
# The Type-IIc lemma statement and the operative factor extraction later in
# its proof use 52*eta.  One earlier display says 100*eta; that isolated
# mismatch must be corrected at source level rather than blended into the
# certificate.  We audit the internally consistent 52*eta formulation.
ENDPOINT_ETA = Fraction(1, 10_000_000_000_000)
ENDPOINT_LOSS_CONSTANT = 52
ENDPOINT_NEGATIVE_WIDTH = Fraction(
    ENDPOINT_LOSS_CONSTANT + 3, 16
) * ENDPOINT_ETA


def encode(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def z3_rational(value: Fraction):
    return z3.Q(value.numerator, value.denominator)


def as_fraction(value) -> Fraction:
    value = z3.simplify(value)
    return Fraction(value.numerator_as_long(), value.denominator_as_long())


def exact_pack(
    weights: list[Fraction], capacities: list[Fraction]
) -> tuple[int, ...] | None:
    """Return one exact four-bin assignment, or None if none exists."""
    order = sorted(range(len(weights)), key=weights.__getitem__, reverse=True)
    ordered_weights = [weights[index] for index in order]
    suffix_mass = [Fraction(0)] * (len(order) + 1)
    for position in range(len(order) - 1, -1, -1):
        suffix_mass[position] = (
            suffix_mass[position + 1] + ordered_weights[position]
        )
    assignment = [-1] * len(weights)
    visited: set[tuple[int, tuple[Fraction, ...]]] = set()

    def search(
        position: int, remaining: tuple[Fraction, ...]
    ) -> tuple[int, ...] | None:
        if position == len(order):
            return tuple(assignment)
        if suffix_mass[position] > sum(remaining, Fraction(0)):
            return None
        weight = ordered_weights[position]
        if weight > max(remaining, default=Fraction(-1)):
            return None

        # Bin labels do not affect feasibility.  Canonicalizing their residual
        # capacities collapses states that differ only by a permutation.
        state = (position, tuple(sorted(remaining, reverse=True)))
        if state in visited:
            return None
        visited.add(state)

        item = order[position]
        equivalent_bins: set[Fraction] = set()
        # Best-fit order usually reaches a witness before creating many loose
        # residual-capacity states; it does not remove any assignment.
        bin_order = sorted(range(len(remaining)), key=remaining.__getitem__)
        for bin_index in bin_order:
            capacity = remaining[bin_index]
            if capacity in equivalent_bins:
                continue
            equivalent_bins.add(capacity)
            if weight > capacity:
                continue
            assignment[item] = bin_index
            updated = list(remaining)
            updated[bin_index] -= weight
            result = search(position + 1, tuple(updated))
            if result is not None:
                return result
        assignment[item] = -1
        return None

    return search(0, tuple(capacities))


def group_bound(group_size: int, large_bound: Fraction) -> Fraction:
    return B_SMALL if group_size <= 2 else large_bound


def audit_pair(
    first_size: int,
    second_size: int,
    large_bound: Fraction,
    max_iterations: int = 100_000,
) -> dict[str, Any]:
    count = first_size + second_size
    coordinates = z3.Reals(" ".join(f"y_{i}" for i in range(count)))
    gamma, omega_zero = z3.Reals("gamma omega_zero")
    solver = z3.Solver()

    for coordinate in coordinates:
        solver.add(coordinate >= z3_rational(DELTA), coordinate <= 1)
    if first_size:
        solver.add(
            z3.Sum(coordinates[:first_size])
            <= z3_rational(group_bound(first_size, large_bound))
        )
    if second_size:
        solver.add(
            z3.Sum(coordinates[first_size:])
            <= z3_rational(group_bound(second_size, large_bound))
        )

    gamma_max = (
        Fraction(1, 3)
        + 8 * OMEGA
        + Fraction(7, 3) * DELTA
        + 3 * EPSILON
    )
    solver.add(
        omega_zero >= 0,
        omega_zero <= z3_rational(OMEGA),
        gamma >= z3_rational(XI_2 - EPSILON),
        gamma <= z3_rational(gamma_max),
    )
    capacities = [
        gamma - z3_rational(2 * DELTA) - 8 * omega_zero - z3_rational(EPSILON),
        z3_rational(Fraction(1, 2)) - gamma - 2 * omega_zero - z3_rational(EPSILON),
        4 * omega_zero + z3_rational(DELTA - EPSILON),
        8 * omega_zero,
    ]

    covering_partitions: list[tuple[int, ...]] = []
    for iteration in range(max_iterations):
        status = solver.check()
        if status == z3.unsat:
            return {
                "status": "covered",
                "iterations": iteration,
                "coveringPartitions": [list(item) for item in covering_partitions],
            }
        if status != z3.sat:
            return {"status": "unknown", "iterations": iteration}

        model = solver.model()
        weights = [
            as_fraction(model.eval(coordinate, model_completion=True))
            for coordinate in coordinates
        ]
        concrete_capacities = [
            as_fraction(model.eval(capacity, model_completion=True))
            for capacity in capacities
        ]
        partition = exact_pack(weights, concrete_capacities)
        if partition is None:
            return {
                "status": "counterexample",
                "iterations": iteration,
                "weights": [encode(item) for item in weights],
                "capacities": [encode(item) for item in concrete_capacities],
                "gamma": encode(model.eval(gamma, model_completion=True).as_fraction()),
                "omegaZero": encode(
                    model.eval(omega_zero, model_completion=True).as_fraction()
                ),
            }

        covering_partitions.append(partition)
        failures = []
        for bin_index in range(4):
            terms = [
                coordinates[index]
                for index, assigned_bin in enumerate(partition)
                if assigned_bin == bin_index
            ]
            total = z3.Sum(terms) if terms else z3_rational(Fraction(0))
            failures.append(total > capacities[bin_index])
        solver.add(z3.Or(failures))

    return {"status": "iteration-limit", "iterations": max_iterations}


def audit_bound(large_bound: Fraction) -> dict[str, Any]:
    pairs: dict[str, Any] = {}
    for first_size in range(0, MAX_GROUP_SIZE + 1):
        for second_size in range(max(first_size, 1), MAX_GROUP_SIZE + 1):
            key = f"{first_size},{second_size}"
            pairs[key] = audit_pair(first_size, second_size, large_bound)
    return {
        "largeGroupBound": encode(large_bound),
        "allPairsCovered": all(
            result["status"] == "covered" for result in pairs.values()
        ),
        "totalCoveringPartitions": sum(
            len(result.get("coveringPartitions", [])) for result in pairs.values()
        ),
        "pairs": pairs,
    }


def audit_endpoint_pair(
    first_size: int,
    second_size: int,
    large_bound: Fraction,
    max_iterations: int = 100_000,
) -> dict[str, Any]:
    """Audit the retained-loss Type-IIc capacities just below sqrt(x)."""
    count = first_size + second_size
    coordinates = z3.Reals(" ".join(f"e_{i}" for i in range(count)))
    gamma, omega_zero, eta = z3.Reals(
        "endpoint_gamma endpoint_omega endpoint_eta"
    )
    solver = z3.Solver()
    for coordinate in coordinates:
        solver.add(coordinate >= z3_rational(DELTA), coordinate <= 1)
    if first_size:
        solver.add(
            z3.Sum(coordinates[:first_size])
            <= z3_rational(group_bound(first_size, large_bound))
        )
    if second_size:
        solver.add(
            z3.Sum(coordinates[first_size:])
            <= z3_rational(group_bound(second_size, large_bound))
        )

    gamma_max = (
        Fraction(1, 3)
        + 8 * OMEGA
        + Fraction(7, 3) * DELTA
        + 3 * EPSILON
    )
    solver.add(
        eta >= 0,
        eta <= z3_rational(ENDPOINT_ETA),
        omega_zero
        >= -z3_rational(Fraction(ENDPOINT_LOSS_CONSTANT + 3, 16))
        * eta,
        omega_zero <= 0,
        gamma >= z3_rational(XI_2 - EPSILON),
        gamma <= z3_rational(gamma_max),
    )
    loss = ENDPOINT_LOSS_CONSTANT
    capacities = [
        gamma
        - z3_rational(2 * DELTA)
        - z3_rational(Fraction(loss + 6)) * eta
        - 8 * omega_zero,
        z3_rational(Fraction(1, 2))
        - gamma
        - 6 * eta
        - 2 * omega_zero,
        z3_rational(DELTA) + 9 * eta + 4 * omega_zero,
        z3_rational(Fraction(loss + 3)) * eta + 8 * omega_zero,
    ]

    covering_partitions: list[tuple[int, ...]] = []
    for iteration in range(max_iterations):
        status = solver.check()
        if status == z3.unsat:
            return {
                "status": "covered",
                "iterations": iteration,
                "coveringPartitions": [list(item) for item in covering_partitions],
            }
        if status != z3.sat:
            return {"status": "unknown", "iterations": iteration}

        model = solver.model()
        weights = [
            as_fraction(model.eval(coordinate, model_completion=True))
            for coordinate in coordinates
        ]
        concrete_capacities = [
            as_fraction(model.eval(capacity, model_completion=True))
            for capacity in capacities
        ]
        partition = exact_pack(weights, concrete_capacities)
        if partition is None:
            return {
                "status": "counterexample",
                "iterations": iteration,
                "weights": [encode(item) for item in weights],
                "capacities": [encode(item) for item in concrete_capacities],
                "gamma": encode(
                    as_fraction(model.eval(gamma, model_completion=True))
                ),
                "omegaZero": encode(
                    as_fraction(model.eval(omega_zero, model_completion=True))
                ),
                "eta": encode(
                    as_fraction(model.eval(eta, model_completion=True))
                ),
            }

        covering_partitions.append(partition)
        failures = []
        for bin_index in range(4):
            terms = [
                coordinates[index]
                for index, assigned_bin in enumerate(partition)
                if assigned_bin == bin_index
            ]
            total = z3.Sum(terms) if terms else z3_rational(Fraction(0))
            failures.append(total > capacities[bin_index])
        solver.add(z3.Or(failures))

    return {"status": "iteration-limit", "iterations": max_iterations}


def audit_endpoint_bound(large_bound: Fraction) -> dict[str, Any]:
    pairs: dict[str, Any] = {}
    for first_size in range(0, MAX_GROUP_SIZE + 1):
        for second_size in range(max(first_size, 1), MAX_GROUP_SIZE + 1):
            key = f"{first_size},{second_size}"
            pairs[key] = audit_endpoint_pair(
                first_size, second_size, large_bound
            )
    return {
        "largeGroupBound": encode(large_bound),
        "allPairsCovered": all(
            result["status"] == "covered" for result in pairs.values()
        ),
        "totalCoveringPartitions": sum(
            len(result.get("coveringPartitions", [])) for result in pairs.values()
        ),
        "pairs": pairs,
    }


def main() -> None:
    candidate = audit_bound(B_CANDIDATE)
    failure_boundary = audit_bound(B_FAILURE)
    endpoint_candidate = audit_endpoint_bound(B_CANDIDATE)
    endpoint_analytic_margins = {
        "typeIIcFirst": Fraction(7, 375) - 6 * EPSILON,
        "typeIIcSecond": Fraction(3, 125) - EPSILON,
        "typeIIcThird": Fraction(1, 125) - 4 * EPSILON,
        "endpointFirstBinLossBudget": (
            Fraction(3, 125)
            - EPSILON
            - (ENDPOINT_LOSS_CONSTANT + 6) * ENDPOINT_ETA
        ),
        "endpointSecondBinLossBudget": (
            Fraction(11, 1500) - 3 * EPSILON - 6 * ENDPOINT_ETA
        ),
    }
    capacity_dominance_margins = {
        # Actual retained-loss capacity minus the simplified capacity used
        # by audit_bound on omega_0 >= 0.
        "first": EPSILON
        - (ENDPOINT_LOSS_CONSTANT + 6) * ENDPOINT_ETA,
        "second": EPSILON - 6 * ENDPOINT_ETA,
        "third": EPSILON + 9 * ENDPOINT_ETA,
        "fourth": (ENDPOINT_LOSS_CONSTANT + 3) * ENDPOINT_ETA,
    }
    omega = OMEGA
    first_bin_bounds = {
        "typeI": XI_1 - 2 * EPSILON,
        "typeIIa": Fraction(2, 5)
        + Fraction(24, 5) * omega
        + Fraction(7, 5) * DELTA
        - 2 * EPSILON,
        "typeIIb": Fraction(1, 3)
        + 8 * omega
        + Fraction(7, 3) * DELTA
        - 4 * EPSILON,
        "typeIII": Fraction(1, 1)
        - 6 * omega
        - Fraction(3, 2) * Fraction(2, 5)
        - 2 * EPSILON,
    }
    structural_checks = {
        "deltaBelowB1": DELTA < B_SMALL,
        "B1AtMostB2": B_SMALL <= B_SMALL,
        "B3AtMostB2PlusDelta": B_CANDIDATE <= B_SMALL + DELTA,
        "BSequenceNondecreasing": B_SMALL <= B_CANDIDATE,
        "typesI_IIa_IIb_IIIAllowAllMassInFirstBin": (
            2 * B_CANDIDATE < min(first_bin_bounds.values())
        ),
    }
    payload = {
        "name": "Exact Type-IIc support audit for the k=48, diameter-236 candidate",
        "arithmetic": "exact rational QF_LRA plus exact rational bin packing",
        "domain": {
            "delta": encode(DELTA),
            "epsilon": encode(EPSILON),
            "xi2": encode(XI_2),
            "omega": encode(OMEGA),
            "omegaZeroRangeAudited": ["0/1", encode(OMEGA)],
            "smallGroupBound": encode(B_SMALL),
            "maximumFeasibleGroupSize": MAX_GROUP_SIZE,
        },
        "structuralChecks": structural_checks,
        "firstBinBounds": {
            key: encode(value) for key, value in first_bin_bounds.items()
        },
        "candidate": candidate,
        "boundaryCheck": failure_boundary,
        "negativeEndpointRepairCandidate": {
            "eta": encode(ENDPOINT_ETA),
            "lossConstant": ENDPOINT_LOSS_CONSTANT,
            "omegaZeroRangeAudited": [
                encode(-ENDPOINT_NEGATIVE_WIDTH),
                "0/1",
            ],
            "etaRangeAudited": ["0/1", encode(ENDPOINT_ETA)],
            "coupledEndpointConstraint": "omega_0 >= -(K+3)*eta/16",
            "retainedCapacities": [
                "gamma-2*delta-(K+6)*eta-8*omega_0",
                "1/2-gamma-6*eta-2*omega_0",
                "delta+9*eta+4*omega_0",
                "(K+3)*eta+8*omega_0",
            ],
            "audit": endpoint_candidate,
            "analyticMargins": {
                key: encode(value)
                for key, value in endpoint_analytic_margins.items()
            },
            "allDisplayedAnalyticMarginsPositive": all(
                value > 0 for value in endpoint_analytic_margins.values()
            ),
            "nonnegativeCapacityDominanceMargins": {
                key: encode(value)
                for key, value in capacity_dominance_margins.items()
            },
            "allNonnegativeRetainedCapacitiesDominateSimplifiedAudit": all(
                value > 0 for value in capacity_dominance_margins.values()
            ),
            "bombieriVinogradovSplice": (
                "For omega_0 <= -(K+3)*eta/16, q is at most a constant "
                "times x^(1/2-(K+3)*eta/8), hence eventually at most "
                "x^(1/2)/(log x)^B for every fixed B.  The retained-loss "
                "audit covers the complementary band up to omega_0=0."
            ),
            "remainingAnalyticStep": (
                "Independently validate the deep Type-IIc equidistribution "
                "lemma after correcting its isolated 100-versus-52 display. "
                "The elementary endpoint/Bombieri-Vinogradov coverage is "
                "closed by the recorded inequalities."
            ),
        },
        "sourceIssueAndRepair": (
            "arXiv:2608.31126v1 states omega_0 in [-epsilon, omega]. "
            "For omega_0 < 0 its fourth capacity 8*omega_0 is negative, so "
            "even an empty fourth bin violates the simplified displayed "
            "inequality.  The certificate retains the omitted 55*eta loss "
            "capacity near zero and uses Bombieri-Vinogradov below that fixed "
            "band.  It also consistently uses the loss 52 appearing in the "
            "lemma statement and operative factor extraction; the isolated "
            "100 display still requires correction in the source theorem."
        ),
        "claimBoundary": (
            "A covered result certifies the finite partition property on the "
            "stated rational polytope. It does not certify the variational "
            "inequality 48*J(F)>I(F), DHL[48,2], or H_1<=236."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "candidateBound": candidate["largeGroupBound"],
                "candidateCovered": candidate["allPairsCovered"],
                "candidatePartitions": candidate["totalCoveringPartitions"],
                "boundaryBound": failure_boundary["largeGroupBound"],
                "boundaryCovered": failure_boundary["allPairsCovered"],
                "boundaryCounterexamplePairs": [
                    key
                    for key, value in failure_boundary["pairs"].items()
                    if value["status"] == "counterexample"
                ],
                "negativeEndpointCovered": endpoint_candidate["allPairsCovered"],
                "negativeEndpointPartitions": endpoint_candidate[
                    "totalCoveringPartitions"
                ],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
