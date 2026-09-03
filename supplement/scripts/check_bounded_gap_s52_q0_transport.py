#!/usr/bin/env python3
"""Audit the q_0 transport in the corrected S52 q-van der Corput step.

This is a symbolic power ledger, not a numerical experiment.  It records the
only places where the enlarged V-window can enter the proof that removes the
large-gcd and diagonal tuples.  Powers are relative to fixed H; the final
column also substitutes H <= x^(4 omega + delta + 7 eta) / q_0 when needed.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Power:
    q0: int = 0
    h: int = 0

    def __add__(self, other: "Power") -> "Power":
        return Power(self.q0 + other.q0, self.h + other.h)

    def residual_q0(self) -> int:
        """Substitute one q_0^{-1} for every positive power of H."""

        return self.q0 - self.h


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="outputs/bounded_gap_s52_q0_transport.json",
    )
    args = parser.parse_args()

    # Corrected interfaces, with H left as an independent scale.
    delta = Power(q0=-2)
    delta1 = Power(q0=-2)
    y = Power(q0=1, h=2)  # Y <= q0 x^(delta+5 eta) H^2 / g.
    inverse_large_gcd = delta1 + Power(h=-2)
    y_over_large_gcd = y + inverse_large_gcd

    branches = {
        # Delta H^2 N.
        "exact_diagonal_error": delta + Power(h=2),
        # max(H^2/t,H) N^2 (Y/q); the max/min/H^2 product is 1.
        "equal_y_unequal_n_error": y_over_large_gcd,
        # Delta Y^2/t is the longer of the two n=ntilde terms.
        "equal_n_unequal_y_long_term": delta + y + y,
        # The possible C-values cost at most q0; Delta supplies q0^-2.
        "fully_off_diagonal_error": delta + Power(q0=1, h=2),
        # This is the one main-term loss that must be retained.
        "diagonal_pair_boundary": y + y + Power(h=-1),
    }

    records: dict[str, object] = {}
    errors_nonpositive = True
    for name, power in branches.items():
        residual = power.residual_q0()
        is_main = name == "diagonal_pair_boundary"
        if not is_main:
            errors_nonpositive &= residual <= 0
        records[name] = {
            "q0PowerAtFixedH": power.q0,
            "hPower": power.h,
            "q0PowerAfterHUpper": residual,
            "role": "retained main-term loss" if is_main else "discarded error",
            "passes": residual <= 2 if is_main else residual <= 0,
        }

    status = {
        "YOverLargeGcdHasQ0Saving": y_over_large_gcd.q0 == -1,
        "allDiscardedErrorsHaveNonpositiveResidualQ0": errors_nonpositive,
        "onlyRetainedPositiveCostIsAtMostQ0Squared": (
            records["diagonal_pair_boundary"]["q0PowerAtFixedH"] <= 2
        ),
    }
    status["allChecksPass"] = all(status.values())

    output = {
        "name": "S52 intermediate q0 transport ledger",
        "interfaces": {
            "VUpper": "q0 x^(delta+5 eta) H",
            "YUpper": "q0 x^(delta+5 eta) H^2 / g",
            "DeltaUpper": "N / (q0^2 x^(50 eta) H^2)",
            "Delta1Upper": "N / (q0^2 x^(55 eta) H^2)",
            "largeGcdThreshold": (
                "max(1/t,1/H) x^(delta+100 eta) H^2 N / (g Delta1)"
            ),
            "YOverLargeGcd": (
                "q0^-1 x^(-150 eta) min(t,H) / H^2"
            ),
        },
        "branches": records,
        "status": status,
        "conclusion": (
            "The non-diagonal q0 losses are absorbed by Delta or Delta1. "
            "Only the diagonal pair boundary retains the q0^2 cost."
        ),
    }

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    if not status["allChecksPass"]:
        raise SystemExit("S52 q0 transport audit failed")
    print(json.dumps(status, sort_keys=True))
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
