#!/usr/bin/env python3
"""Finite regression audit for the generalized gcd-average repair.

The mathematical proof is Lemma 4.7 in the companion note.  This script
exhaustively checks its explicit finite inequality over a deliberately broad
set of small parameters and records the q_0 exponent budget used afterwards.
It is a regression check, not a substitute for the divisor-expansion proof.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def divisor_count(n: int) -> int:
    count = 1
    p = 2
    while p * p <= n:
        exponent = 0
        while n % p == 0:
            n //= p
            exponent += 1
        if exponent:
            count *= exponent + 1
        p += 1
    if n > 1:
        count *= 2
    return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="outputs/bounded_gap_s52_gcd_repair.json",
    )
    args = parser.parse_args()

    cases = 0
    maximum_ratio = 0.0
    maximizer: dict[str, int] | None = None
    for m in range(1, 121):
        tau_m = divisor_count(m)
        for a in range(-10, 11):
            if a == 0:
                continue
            gcd_am = math.gcd(a, m)
            for b in range(-10, 11):
                for k_bound in range(0, 13):
                    gcd_values = sorted(
                        math.gcd(a * k + b, m)
                        for k in range(-k_bound, k_bound + 1)
                    )
                    cursor = 0
                    lhs = 0
                    for threshold in range(1, m + 2):
                        while (
                            cursor < len(gcd_values)
                            and gcd_values[cursor] <= threshold
                        ):
                            lhs += gcd_values[cursor]
                            cursor += 1
                        rhs = tau_m * (
                            2 * gcd_am * k_bound + threshold
                        )
                        cases += 1
                        if lhs > rhs:
                            raise SystemExit(
                                "gcd repair failed at "
                                f"m={m}, A={a}, B={b}, K={k_bound}, "
                                f"T={threshold}: {lhs}>{rhs}"
                            )
                        ratio = lhs / rhs if rhs else 0.0
                        if ratio > maximum_ratio:
                            maximum_ratio = ratio
                            maximizer = {
                                "m": m,
                                "A": a,
                                "B": b,
                                "K": k_bound,
                                "T": threshold,
                                "lhs": lhs,
                                "rhs": rhs,
                            }

    original_q0 = [-3, -6, -5, -8, -7, -10, -6, -9]
    repaired_q0 = [value + 1 for value in original_q0]
    report = {
        "name": "finite regression audit for the S52 gcd-average repair",
        "claimBoundary": (
            "The exhaustive check guards the implementation; Lemma 4.7's "
            "divisor expansion is the proof."
        ),
        "casesChecked": cases,
        "maximumLhsOverRhs": maximum_ratio,
        "maximizer": maximizer,
        "explicitBound": "tau(m) * (2*gcd(A,m)*K + T)",
        "application": {
            "oldNormalization": "gcd(w2,m)*T",
            "newLoss": "1 + K/T up to an absolute factor",
            "originalResidualQ0Exponents": original_q0,
            "repairedResidualQ0Exponents": repaired_q0,
            "allRepairedQ0ExponentsNonpositive": all(
                value <= 0 for value in repaired_q0
            ),
            "etaGainInSecondBranch": -100,
        },
        "allChecksPass": True,
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "allChecksPass": True,
        "casesChecked": cases,
        "maximumLhsOverRhs": maximum_ratio,
    }, sort_keys=True))
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
