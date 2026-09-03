#!/usr/bin/env python3
"""Exact exponent ledger for the corrected S52 factor-interface proof.

The script does not numerically test an exponential-sum estimate.  It checks
the deterministic algebra after the q-van der Corput reduction: the three
remaining error families are dominated by the three strict hypotheses in
the S. Type I lemma, even when Delta* is kept as a genuine minimum.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path


@dataclass(frozen=True)
class LinearForm:
    constant: Fraction = Fraction(0)
    omega: Fraction = Fraction(0)
    delta: Fraction = Fraction(0)
    gamma: Fraction = Fraction(0)

    def __sub__(self, other: "LinearForm") -> "LinearForm":
        return LinearForm(
            self.constant - other.constant,
            self.omega - other.omega,
            self.delta - other.delta,
            self.gamma - other.gamma,
        )

    def evaluate(self, *, omega: Fraction, delta: Fraction, gamma: Fraction) -> Fraction:
        return (
            self.constant
            + self.omega * omega
            + self.delta * delta
            + self.gamma * gamma
        )

    def as_dict(self) -> dict[str, str]:
        return {
            "constant": frac(self.constant),
            "omega": frac(self.omega),
            "delta": frac(self.delta),
            "gamma": frac(self.gamma),
        }

    def display(self) -> str:
        terms = [frac(self.constant)]
        for coefficient, name in (
            (self.omega, "omega"),
            (self.delta, "delta"),
            (self.gamma, "gamma"),
        ):
            if coefficient:
                terms.append(f"({frac(coefficient)})*{name}")
        return " + ".join(terms)


@dataclass(frozen=True)
class RawExponent:
    """Exponent before substituting H <= x^(4w+d+7e)/q0."""

    constant: Fraction = Fraction(0)
    omega: Fraction = Fraction(0)
    delta: Fraction = Fraction(0)
    gamma: Fraction = Fraction(0)
    eta: Fraction = Fraction(0)
    q0: Fraction = Fraction(0)
    h: Fraction = Fraction(0)

    def __add__(self, other: "RawExponent") -> "RawExponent":
        return RawExponent(*(
            left + right
            for left, right in zip(self.values(), other.values(), strict=True)
        ))

    def __sub__(self, other: "RawExponent") -> "RawExponent":
        return RawExponent(*(
            left - right
            for left, right in zip(self.values(), other.values(), strict=True)
        ))

    def values(self) -> tuple[Fraction, ...]:
        return (
            self.constant,
            self.omega,
            self.delta,
            self.gamma,
            self.eta,
            self.q0,
            self.h,
        )

    def substitute_h_upper(self) -> tuple[LinearForm, int, int]:
        """Return (eta-free form, eta coefficient, residual q0 exponent)."""

        return (
            LinearForm(
                self.constant,
                self.omega + 4 * self.h,
                self.delta + self.h,
                self.gamma,
            ),
            int(self.eta + 7 * self.h),
            int(self.q0 - self.h),
        )


def raw(
    *,
    constant: int = 0,
    omega: int = 0,
    delta: int = 0,
    gamma: int = 0,
    eta: int = 0,
    q0: int = 0,
    h: int = 0,
) -> RawExponent:
    return RawExponent(
        Fraction(constant),
        Fraction(omega),
        Fraction(delta),
        Fraction(gamma),
        Fraction(eta),
        Fraction(q0),
        Fraction(h),
    )


def frac(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def nonnegative_coefficient_difference(form: LinearForm) -> bool:
    """Certify nonnegativity for omega, delta >= 0 and free gamma.

    Target-minus-derived differences used below have zero constant and gamma
    coefficients, so coefficientwise checking is exact.
    """

    return (
        form.constant == 0
        and form.gamma == 0
        and form.omega >= 0
        and form.delta >= 0
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="outputs/bounded_gap_s52_exponent_ledger.json",
        help="JSON output path",
    )
    args = parser.parse_args()

    targets = {
        "WTS1": LinearForm(-1, 8, 4, 2),
        "WTS2": LinearForm(0, 32, 10, -1),
        "WTS3": LinearForm(1, 48, 16, -4),
    }

    # The eta coefficients come from retaining every q_0 factor and using
    # H <= x^(4 omega + delta + 7 eta)/q_0.  The two WTS3 pairs correspond
    # to 1/Delta* = max(x^(5 eta)|Lambda|/N, 1/Delta_1).
    branches = {
        "WTS1_H": ("WTS1", LinearForm(-1, 4, 3, 2), 94),
        "WTS1_H2": ("WTS1", LinearForm(-1, 8, 3, 2), 91),
        "WTS2_H7": ("WTS2", LinearForm(0, 28, 10, -1), 245),
        "WTS2_H8": ("WTS2", LinearForm(0, 32, 10, -1), 242),
        "WTS3_Lambda_H11": ("WTS3", LinearForm(1, 44, 16, -4), 282),
        "WTS3_Lambda_H12": ("WTS3", LinearForm(1, 48, 16, -4), 279),
        "WTS3_Delta1_H11": ("WTS3", LinearForm(1, 44, 16, -4), 327),
        "WTS3_Delta1_H12": ("WTS3", LinearForm(1, 48, 16, -4), 324),
    }

    # Mechanical reconstruction from the primitive inequalities.  Variables
    # h and q0 denote log_x H and log_x q_0 before the H upper bound is used.
    common = raw(delta=1, eta=131, q0=-2)
    oscillatory_h5 = raw(delta=1, eta=10, q0=2, h=5)
    oscillatory_h6 = raw(h=6)
    m_lower = raw(constant=1, gamma=-2, eta=54, q0=2, h=4)
    delta_star_reciprocal = raw(delta=1, gamma=-1, eta=55, q0=2, h=2)

    # m <= q0 x^(delta-eta) M H^2 / Delta1, with M=x/N.
    m_upper_without_delta1 = raw(
        constant=1, delta=1, gamma=-1, eta=-1, q0=1, h=2
    )
    n_reciprocal = raw(gamma=-1)
    delta_product_lambda = raw(delta=2, gamma=-2, eta=65, q0=3, h=4)
    delta_product_delta1 = raw(delta=2, gamma=-2, eta=110, q0=4, h=4)

    mechanically_derived = {
        "WTS1_H": common + oscillatory_h5 - m_lower,
        "WTS1_H2": common + oscillatory_h6 - m_lower,
        "WTS2_H7": common + oscillatory_h5 + delta_star_reciprocal,
        "WTS2_H8": common + oscillatory_h6 + delta_star_reciprocal,
        "WTS3_Lambda_H11": (
            common
            + n_reciprocal
            + m_upper_without_delta1
            + oscillatory_h5
            + delta_product_lambda
        ),
        "WTS3_Lambda_H12": (
            common
            + n_reciprocal
            + m_upper_without_delta1
            + oscillatory_h6
            + delta_product_lambda
        ),
        "WTS3_Delta1_H11": (
            common
            + n_reciprocal
            + m_upper_without_delta1
            + oscillatory_h5
            + delta_product_delta1
        ),
        "WTS3_Delta1_H12": (
            common
            + n_reciprocal
            + m_upper_without_delta1
            + oscillatory_h6
            + delta_product_delta1
        ),
    }

    branch_records: dict[str, object] = {}
    all_dominated = True
    all_mechanically_matched = True
    all_q0_powers_nonpositive = True
    all_gcd_repair_branches_pass = True
    max_eta_coefficient = 0
    for name, (target_name, form, eta_coefficient) in branches.items():
        derived_form, derived_eta, residual_q0 = mechanically_derived[
            name
        ].substitute_h_upper()
        mechanically_matched = derived_form == form and derived_eta == eta_coefficient
        q0_nonpositive = residual_q0 <= 0
        gcd_repair_q0 = residual_q0 + 1
        gcd_repair_eta = eta_coefficient - 100
        gcd_repair_pass = gcd_repair_q0 <= 0 and gcd_repair_eta <= eta_coefficient
        all_mechanically_matched &= mechanically_matched
        all_q0_powers_nonpositive &= q0_nonpositive
        all_gcd_repair_branches_pass &= gcd_repair_pass
        difference = targets[target_name] - form
        dominated = nonnegative_coefficient_difference(difference)
        all_dominated &= dominated
        max_eta_coefficient = max(max_eta_coefficient, eta_coefficient)
        branch_records[name] = {
            "target": target_name,
            "baseExponent": form.as_dict(),
            "baseExponentDisplay": form.display(),
            "etaCoefficient": eta_coefficient,
            "residualQ0Exponent": residual_q0,
            "residualQ0Nonpositive": q0_nonpositive,
            "gcdAverageRepair": {
                "multiplier": "q0*x^(-100 eta)",
                "etaCoefficient": gcd_repair_eta,
                "residualQ0Exponent": gcd_repair_q0,
                "passes": gcd_repair_pass,
            },
            "mechanicallyDerived": mechanically_matched,
            "targetMinusBase": difference.as_dict(),
            "coefficientwiseDominated": dominated,
        }

    # Moving the product of primes <= D_0 between the two initial factors
    # costs x^o(1).  For fixed eta and sufficiently large x it is harmless;
    # eight eta units are reserved here so the ledger has an explicit slot.
    subpower_eta_reserve = 8
    safe_eta_coefficient = 400
    eta_budget_ok = max_eta_coefficient + subpower_eta_reserve <= safe_eta_coefficient

    omega = Fraction(57, 10_000)
    delta = Fraction(1, 50)
    gamma_lower_2 = 32 * omega + 10 * delta
    gamma_lower_3 = (1 + 48 * omega + 16 * delta) / 4
    gamma_lower = max(gamma_lower_2, gamma_lower_3)
    gamma_upper = (1 - 8 * omega - 4 * delta) / 2
    gamma_midpoint = (gamma_lower + gamma_upper) / 2

    slacks = {
        "WTS1": -targets["WTS1"].evaluate(
            omega=omega, delta=delta, gamma=gamma_midpoint
        ),
        "WTS2": -targets["WTS2"].evaluate(
            omega=omega, delta=delta, gamma=gamma_midpoint
        ),
        "WTS3": -targets["WTS3"].evaluate(
            omega=omega, delta=delta, gamma=gamma_midpoint
        ),
    }
    eta_upper = min(delta, *slacks.values()) / (2 * safe_eta_coefficient)

    output = {
        "name": "corrected S52 factor-interface exponent ledger",
        "status": {
            "allBaseExponentsDominated": all_dominated,
            "allBranchesMechanicallyDerived": all_mechanically_matched,
            "allResidualQ0PowersNonpositive": all_q0_powers_nonpositive,
            "allGcdAverageRepairBranchesPass": all_gcd_repair_branches_pass,
            "etaBudgetFits": eta_budget_ok,
            "allChecksPass": (
                all_dominated
                and all_mechanically_matched
                and all_q0_powers_nonpositive
                and all_gcd_repair_branches_pass
                and eta_budget_ok
                and gamma_lower < gamma_upper
            ),
        },
        "interfaceIdentities": {
            "H": "x^eta R Q^2 / (q0 M)",
            "extractedDelta": "N / (q0^2 x^(50 eta) H^2)",
            "Delta1Range": [
                "N / (q0^2 x^(delta+55 eta) H^2)",
                "N / (q0^2 x^(55 eta) H^2)",
            ],
            "LambdaUpper": "q0 x^(delta+5 eta) H^2 / (w1 gcd(v1,v2))",
            "DeltaStar": "min(N/(|Lambda| x^(5 eta)), Delta1)",
            "DeltaStarLower": "N / (q0^2 x^(delta+55 eta) H^2)",
            "generalGcdAverage": "sum gcd <= x^eta ((A,m) K + T)",
            "gcdAverageRepairFactor": "1 + S_x^O(1) q0 x^(-100 eta)",
        },
        "targets": {
            name: {
                "exponent": form.as_dict(),
                "display": form.display(),
                "requiredSign": "< 0",
            }
            for name, form in targets.items()
        },
        "branches": branch_records,
        "etaAccounting": {
            "maximumDerivedCoefficient": max_eta_coefficient,
            "subpowerTransportReserve": subpower_eta_reserve,
            "safeCoefficient": safe_eta_coefficient,
            "fits": eta_budget_ok,
            "sufficientRule": (
                "eta < min(delta, -WTS1, -WTS2, -WTS3) / "
                f"{2 * safe_eta_coefficient}"
            ),
        },
        "diameter226ParameterSlice": {
            "omega": frac(omega),
            "delta": frac(delta),
            "gammaLowerFromWTS2": frac(gamma_lower_2),
            "gammaLowerFromWTS3": frac(gamma_lower_3),
            "gammaLower": frac(gamma_lower),
            "gammaUpper": frac(gamma_upper),
            "intervalNonempty": gamma_lower < gamma_upper,
            "sampleGamma": frac(gamma_midpoint),
            "sampleSlacks": {name: frac(value) for name, value in slacks.items()},
            "sampleSufficientEtaUpper": frac(eta_upper),
        },
    }

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    if not output["status"]["allChecksPass"]:
        raise SystemExit("S52 exponent ledger failed")
    print(json.dumps(output["status"], sort_keys=True))
    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
