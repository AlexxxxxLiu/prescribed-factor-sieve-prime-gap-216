#!/usr/bin/env python3
"""Finite model audit for the prime-support transport lemma.

This script is a mechanical regression check, not the proof.  The proof is
the support-set argument in Lemma 4.3 of the S52 factor-interface note.  A
prime satisfying the lemma's hypotheses occupies exactly one of eight
disjoint support cells (or no cell), so a finite universe checks every local
incidence pattern.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path


PRIMES = (2, 3, 5, 7, 11)
CELLS = (
    "none",
    "d1",
    "r1",
    "q0",
    "u1",
    "g",
    "v1_only",
    "v2_only",
    "q2",
)


def product(values: list[int]) -> int:
    answer = 1
    for value in values:
        answer *= value
    return answer


def squarefree(value: int) -> bool:
    return all(value % (prime * prime) for prime in PRIMES)


def divisors(value: int) -> list[int]:
    factors = [prime for prime in PRIMES if value % prime == 0]
    return [
        product([prime for prime, keep in zip(factors, mask, strict=True) if keep])
        for mask in itertools.product((False, True), repeat=len(factors))
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="outputs/bounded_gap_s52_support_transport.json",
    )
    args = parser.parse_args()

    cases = 0
    q3_cases = 0
    for allocation in itertools.product(CELLS, repeat=len(PRIMES)):
        buckets = {
            cell: product(
                [prime for prime, assigned in zip(PRIMES, allocation, strict=True)
                 if assigned == cell]
            )
            for cell in CELLS
        }
        d1 = buckets["d1"]
        r1 = buckets["r1"]
        r = d1 * r1
        q0 = buckets["q0"]
        u1 = buckets["u1"]
        g = buckets["g"]
        v1_star = buckets["v1_only"]
        v2_star = buckets["v2_only"]
        v1 = g * v1_star
        v2 = g * v2_star
        q2 = buckets["q2"]
        lcm_v = math.lcm(v1, v2)
        m = r1 * q0 * u1 * lcm_v * q2

        assert squarefree(q0 * u1 * v1 * r)
        assert squarefree(q0 * u1 * v2 * r)
        assert squarefree(q0 * q2 * r)
        assert math.gcd(u1 * v1 * v2, q0 * q2) == 1
        assert squarefree(m)
        assert math.gcd(d1, m) == 1
        assert math.gcd(r1 * q0 * u1 * q2, lcm_v) == 1

        for q3 in divisors(m // q0):
            q3_cases += 1
            assert m % (q0 * q3) == 0
            assert squarefree(q0 * q3)
        cases += 1

    status = {
        "allChecksPass": True,
        "primeUniverse": list(PRIMES),
        "supportCells": list(CELLS),
        "supportAllocationsChecked": cases,
        "terminalQ3DivisorsChecked": q3_cases,
    }
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
    print(json.dumps(status, sort_keys=True))
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
