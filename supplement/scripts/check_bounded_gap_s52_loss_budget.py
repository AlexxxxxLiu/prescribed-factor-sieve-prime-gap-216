#!/usr/bin/env python3
"""Parameter audit for the S52 exponent losses.

This script separates principal (omega, delta, gamma) exponents from the
technical eta budget.  It does not prove that the source argument remains
valid after changing the divisor-window constant.  It answers the narrower
algebraic question: after the two factor interfaces have been obtained, which
terminal exponents would change if the constant 52 or the smoothing losses
were changed?
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path


def f(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def decimal(value: Fraction) -> float:
    return float(value.numerator / value.denominator)


def parse_fraction(value: str) -> Fraction:
    return Fraction(value)


def eta_coefficients(
    *,
    window_constant: int,
    smoothing_loss: int,
    terminal_loss: int,
    h_loss: int,
    diagonal_loss: int,
) -> dict[str, int]:
    # d1 has central scale x^{-(C-2) eta}; shortening the d interval by
    # x^{smoothing_loss eta} gives Delta_1 exponent kappa=C-2+smoothing.
    kappa = window_constant - 2 + smoothing_loss
    return {
        "WTS1_H": terminal_loss - kappa + diagonal_loss + h_loss + 1,
        "WTS1_H2": terminal_loss - kappa + 2 * h_loss + 1,
        "WTS2_H7": (
            terminal_loss + kappa + diagonal_loss + 7 * h_loss
        ),
        "WTS2_H8": terminal_loss + kappa + 8 * h_loss,
        "WTS3_Lambda_H11": (
            terminal_loss
            + kappa
            + 2 * diagonal_loss
            + 11 * h_loss
            - 1
        ),
        "WTS3_Lambda_H12": (
            terminal_loss
            + kappa
            + diagonal_loss
            + 12 * h_loss
            - 1
        ),
        "WTS3_Delta1_H11": (
            terminal_loss
            + 2 * kappa
            + diagonal_loss
            + 11 * h_loss
            - 1
        ),
        "WTS3_Delta1_H12": (
            terminal_loss + 2 * kappa + 12 * h_loss - 1
        ),
    }


def parameter_slice(*, omega: Fraction, delta: Fraction) -> dict[str, object]:
    lower_wts2 = 32 * omega + 10 * delta
    lower_wts3 = Fraction(1, 4) + 12 * omega + 4 * delta
    lower = max(lower_wts2, lower_wts3)
    upper = Fraction(1, 2) - 4 * omega - 2 * delta
    midpoint = (lower + upper) / 2
    return {
        "omega": f(omega),
        "delta": f(delta),
        "gammaLowerWTS2": f(lower_wts2),
        "gammaLowerWTS3": f(lower_wts3),
        "activeLowerBranch": "WTS2" if lower_wts2 >= lower_wts3 else "WTS3",
        "gammaLower": f(lower),
        "gammaUpper": f(upper),
        "feasible": lower < upper,
        "sampleGamma": f(midpoint),
        "width": f(upper - lower),
        "widthDecimal": decimal(upper - lower),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--window-constant", type=int, default=52)
    parser.add_argument("--smoothing-loss", type=int, default=5)
    parser.add_argument("--terminal-loss", type=int, default=131)
    parser.add_argument("--h-loss", type=int, default=7)
    parser.add_argument("--diagonal-loss", type=int, default=10)
    parser.add_argument("--gcd-repair-gain", type=int, default=100)
    parser.add_argument("--omega", type=parse_fraction, default=Fraction(57, 10_000))
    parser.add_argument("--delta", type=parse_fraction, default=Fraction(1, 50))
    parser.add_argument(
        "--output",
        default="outputs/bounded_gap_s52_loss_budget.json",
    )
    args = parser.parse_args()

    kappa = args.window_constant - 2 + args.smoothing_loss
    coeffs = eta_coefficients(
        window_constant=args.window_constant,
        smoothing_loss=args.smoothing_loss,
        terminal_loss=args.terminal_loss,
        h_loss=args.h_loss,
        diagonal_loss=args.diagonal_loss,
    )

    expected_default = {
        "WTS1_H": 94,
        "WTS1_H2": 91,
        "WTS2_H7": 245,
        "WTS2_H8": 242,
        "WTS3_Lambda_H11": 282,
        "WTS3_Lambda_H12": 279,
        "WTS3_Delta1_H11": 327,
        "WTS3_Delta1_H12": 324,
    }
    is_default = (
        args.window_constant,
        args.smoothing_loss,
        args.terminal_loss,
        args.h_loss,
        args.diagonal_loss,
    ) == (52, 5, 131, 7, 10)
    default_regression_pass = (not is_default) or coeffs == expected_default

    principal_conditions = {
        "upperGamma": "gamma < 1/2 - 4 omega - 2 delta",
        "lowerGammaWTS2": "gamma > 32 omega + 10 delta",
        "lowerGammaWTS3": "gamma > 1/4 + 12 omega + 4 delta",
        "combinedFeasibility": "72 omega + 24 delta < 1",
        "secondaryFeasibility": "64 omega + 24 delta < 1",
    }

    max_omega_wts2 = (1 - 24 * args.delta) / 72
    max_omega_wts3 = (1 - 24 * args.delta) / 64
    max_omega = min(max_omega_wts2, max_omega_wts3)

    slices = {
        "requested": parameter_slice(omega=args.omega, delta=args.delta),
        "diameter226": parameter_slice(
            omega=Fraction(57, 10_000), delta=Fraction(1, 50)
        ),
        "diameter216Legacy": parameter_slice(
            omega=Fraction(3, 1000), delta=Fraction(7, 250)
        ),
    }

    output = {
        "name": "S52 principal-versus-eta loss budget",
        "status": {
            "defaultRegressionPass": default_regression_pass,
            "requestedPrincipalRegionFeasible": slices["requested"]["feasible"],
            "algebraicAuditPass": default_regression_pass and kappa > 0,
            "claimBoundary": (
                "Changing C is only a post-interface sensitivity calculation; "
                "it is not a proof of an S_C theorem."
            ),
        },
        "parameters": {
            "windowConstantC": args.window_constant,
            "smoothingLoss": args.smoothing_loss,
            "delta1EtaExponentKappa": kappa,
            "terminalLoss": args.terminal_loss,
            "hUpperEtaLoss": args.h_loss,
            "diagonalEtaLoss": args.diagonal_loss,
            "gcdRepairGain": args.gcd_repair_gain,
        },
        "interfaceFormulae": {
            "extractedDelta": "N/(q0^2 x^((C-2) eta) H^2)",
            "delta1": "N/(q0^2 x^((delta+kappa eta)) H^2)",
            "kappa": "C-2+smoothingLoss",
        },
        "principalConditions": principal_conditions,
        "principalIndependence": {
            "windowConstantAffectsPrincipalRegion": False,
            "reason": (
                "C multiplies eta only; all C-dependence vanishes in the "
                "eta-to-zero principal exponent polytope."
            ),
            "maximumOmegaAtRequestedDelta": f(max_omega),
            "maximumOmegaAtRequestedDeltaDecimal": decimal(max_omega),
            "maximumOmegaFromWTS2": f(max_omega_wts2),
            "maximumOmegaFromWTS3": f(max_omega_wts3),
        },
        "etaCoefficients": coeffs,
        "maximumEtaCoefficient": max(coeffs.values()),
        "gcdRepairEtaCoefficients": {
            name: value - args.gcd_repair_gain for name, value in coeffs.items()
        },
        "slices": slices,
        "nextTrueImprovementTarget": {
            "activeGlobalBottleneck": "WTS2_H8 at the edge of the region",
            "requiredTypeOfGain": (
                "reduce a principal coefficient in 32 omega + 10 delta - gamma, "
                "or improve the finite k=45/k=46 Maynard functional"
            ),
            "doesLowering52AloneSuffice": False,
        },
    }

    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    if not output["status"]["algebraicAuditPass"]:
        raise SystemExit("loss-budget audit failed")
    print(json.dumps(output["status"], sort_keys=True))
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
