#!/usr/bin/env python3
"""Finite regression audit for the reversible small-prime preprocessing.

The mathematical statement is Lemma 3.1 and equations (8.2)--(8.3) in the
S52 factor-interface note.  This script checks every squarefree incidence
pattern on a five-prime universe.  It is not a substitute for the proof.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path


PRIMES = (2, 3, 5, 7, 11)
SMALL_PRIMES = (2, 3)


def product(values: list[int]) -> int:
    answer = 1
    for value in values:
        answer *= value
    return answer


def divisors(value: int) -> list[int]:
    factors = [prime for prime in PRIMES if value % prime == 0]
    return [
        product([prime for prime, keep in zip(factors, mask, strict=True) if keep])
        for mask in itertools.product((False, True), repeat=len(factors))
    ]


def squarefree_numbers() -> list[int]:
    return [
        product([prime for prime, keep in zip(PRIMES, mask, strict=True) if keep])
        for mask in itertools.product((False, True), repeat=len(PRIMES))
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="outputs/bounded_gap_s52_reversible_preprocessing.json",
    )
    args = parser.parse_args()

    cases = 0
    moved_prime_counts: dict[int, int] = {}
    small_prime_product = product(list(SMALL_PRIMES))

    for d in squarefree_numbers():
        small_part = math.gcd(d, small_prime_product)
        for r0 in divisors(d):
            q_raw = d // r0
            for u0 in divisors(q_raw):
                for d10 in divisors(r0):
                    t = math.gcd(q_raw, small_part)
                    b = math.gcd(u0, t)
                    q = q_raw // t
                    r = r0 * t
                    u_prime = u0 // b

                    assert q * r == d
                    assert math.gcd(q, small_prime_product) == 1
                    assert q % u_prime == 0
                    assert r % d10 == 0

                    # Retaining (t,b) makes the normalization invertible.
                    assert q * t == q_raw
                    assert r // t == r0
                    assert u_prime * b == u0

                    # These are the two exact interface changes.
                    assert Fraction(u_prime, u0) == Fraction(1, b)
                    assert Fraction(r0 * r0, d**4) == Fraction(r * r, t * t * d**4)

                    moved = sum(t % prime == 0 for prime in SMALL_PRIMES)
                    moved_prime_counts[moved] = moved_prime_counts.get(moved, 0) + 1
                    cases += 1

    status = {
        "allChecksPass": True,
        "primeUniverse": list(PRIMES),
        "smallPrimeUniverse": list(SMALL_PRIMES),
        "tuplesChecked": cases,
        "movedPrimeHistogram": moved_prime_counts,
        "verifiedIdentities": [
            "qr=d",
            "gcd(q,small-prime-product)=1",
            "uPrime|q",
            "d10|r",
            "qRaw=qt",
            "r0=r/t",
            "u0=uPrime*b",
            "uPrime/u0=1/b",
            "r0^2/d^4=r^2/(t^2*d^4)",
        ],
    }
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
    print(json.dumps(status, sort_keys=True))
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
