#!/usr/bin/env python3
"""Build the frozen rational product-radial k=48 candidate certificate.

The rational candidate and its Bernstein positivity certificates are exact.
The FFT convolution scores are explicitly exploratory floating-point data and
do not certify the variational inequality, DHL[48,2], or H_1 <= 236.
"""

from __future__ import annotations

import hashlib
import json
import math
from fractions import Fraction
from math import comb
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs" / "bounded_gap_236_rational_candidate.json"

K = 48
RADIUS = Fraction(521, 2000)
OUTER_RADIUS = Fraction(491, 2000)
DELTA = Fraction(7, 250)
B_SMALL = Fraction(1599, 10000)
B_LARGE = Fraction(1839, 10000)
COMMON_DENOMINATOR = 2**60

# q(u) = sum_j PROFILE_NUMERATORS[j] / 2^60 * T_j(2u-1).
PROFILE_NUMERATORS = [
    84265234121726416,
    -142893298421964544,
    122629789604345376,
    -105729321642428720,
    91391889518002576,
    -79128886273929664,
    68590101744362928,
    -59504938689358944,
    51655742830090568,
    -44863267008961344,
    38977683148286936,
    -33872504633309072,
    29440204798230304,
    -25588908921300000,
    22239811571273840,
    -19325110051517680,
    16786320281221538,
    -14572885215433662,
    12641012637711650,
    -10952696276384644,
    9474885620540164,
    -8178777720987322,
    7039209913546395,
    -6034136553410870,
    5144175970185577,
    -4352216245054130,
    3643070273313602,
    -3003172041936671,
    2420307216208698,
    -1883372057830567,
    1382155436599828,
    -907139283406681,
    449313289073579,
]

# r(u) = sum_j RADIAL_NUMERATORS[j] / 2^60 * P_j(2u-1).
RADIAL_NUMERATORS = [
    75971692445695152,
    -159107193767709472,
    223644392570526624,
    -242290931779389792,
    197096508760712256,
    -138174241580346720,
    75812797481343664,
    -32349295341155232,
    8474450879968022,
]

STRESS_MESHES = [1.0e-4, 5.0e-5, 2.5e-5, 1.25e-5]
QUANTIZED_UNION_STEP = 5.0e-5


def encode(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def poly_add(left: list[Fraction], right: list[Fraction]) -> list[Fraction]:
    size = max(len(left), len(right))
    return [
        (left[index] if index < len(left) else Fraction(0))
        + (right[index] if index < len(right) else Fraction(0))
        for index in range(size)
    ]


def poly_scale(poly: list[Fraction], scalar: Fraction) -> list[Fraction]:
    return [scalar * coefficient for coefficient in poly]


def poly_multiply(
    left: list[Fraction], right: list[Fraction]
) -> list[Fraction]:
    result = [Fraction(0)] * (len(left) + len(right) - 1)
    for first_index, first in enumerate(left):
        for second_index, second in enumerate(right):
            result[first_index + second_index] += first * second
    return result


def shifted_chebyshev_basis(degree: int) -> list[list[Fraction]]:
    affine = [Fraction(-1), Fraction(2)]
    basis = [[Fraction(1)]]
    if degree == 0:
        return basis
    basis.append(affine)
    for _ in range(2, degree + 1):
        basis.append(
            poly_add(
                poly_scale(poly_multiply(affine, basis[-1]), Fraction(2)),
                poly_scale(basis[-2], Fraction(-1)),
            )
        )
    return basis


def shifted_legendre_basis(degree: int) -> list[list[Fraction]]:
    affine = [Fraction(-1), Fraction(2)]
    basis = [[Fraction(1)]]
    if degree == 0:
        return basis
    basis.append(affine)
    for index in range(1, degree):
        numerator = poly_add(
            poly_scale(
                poly_multiply(affine, basis[-1]), Fraction(2 * index + 1)
            ),
            poly_scale(basis[-2], Fraction(-index)),
        )
        basis.append(poly_scale(numerator, Fraction(1, index + 1)))
    return basis


def combine_basis(
    numerators: list[int], basis: list[list[Fraction]]
) -> list[Fraction]:
    result = [Fraction(0)]
    for numerator, basis_polynomial in zip(numerators, basis):
        result = poly_add(
            result,
            poly_scale(
                basis_polynomial, Fraction(numerator, COMMON_DENOMINATOR)
            ),
        )
    return result


def power_to_bernstein(
    power_coefficients: list[Fraction], degree: int
) -> list[Fraction]:
    return [
        sum(
            power_coefficients[power]
            * Fraction(comb(index, power), comb(degree, power))
            for power in range(index + 1)
        )
        for index in range(degree + 1)
    ]


def split_bernstein_midpoint(
    coefficients: list[Fraction],
) -> tuple[list[Fraction], list[Fraction]]:
    rows = [coefficients]
    while len(rows[-1]) > 1:
        rows.append(
            [
                (left + right) / 2
                for left, right in zip(rows[-1], rows[-1][1:])
            ]
        )
    left = [row[0] for row in rows]
    right = [row[-1] for row in reversed(rows)]
    return left, right


def bernstein_positivity_certificate(
    power_coefficients: list[Fraction], degree: int
) -> dict[str, object]:
    root = power_to_bernstein(power_coefficients, degree)
    pending = [(root, 0, Fraction(0), Fraction(1))]
    leaves: list[dict[str, object]] = []
    while pending:
        coefficients, depth, lower, upper = pending.pop()
        minimum = min(coefficients)
        if minimum > 0:
            leaves.append(
                {
                    "interval": [encode(lower), encode(upper)],
                    "depth": depth,
                    "bernsteinCoefficients": [
                        encode(value) for value in coefficients
                    ],
                    "minimumCoefficient": encode(minimum),
                    "minimumCoefficientDecimal": float(minimum),
                }
            )
            continue
        if depth >= 20:
            raise RuntimeError("Bernstein positivity subdivision did not close")
        left, right = split_bernstein_midpoint(coefficients)
        midpoint = (lower + upper) / 2
        pending.append((right, depth + 1, midpoint, upper))
        pending.append((left, depth + 1, lower, midpoint))
    global_minimum = min(
        Fraction(leaf["minimumCoefficient"])
        for leaf in leaves
    )
    return {
        "method": "exact rational midpoint de Casteljau subdivision",
        "allLeafBernsteinCoefficientsStrictlyPositive": True,
        "leafCount": len(leaves),
        "maximumDepth": max(int(leaf["depth"]) for leaf in leaves),
        "globalCertifiedLowerBound": encode(global_minimum),
        "globalCertifiedLowerBoundDecimal": float(global_minimum),
        "leaves": leaves,
    }


def evaluate_profile(normalized_x: np.ndarray) -> np.ndarray:
    coefficients = np.asarray(PROFILE_NUMERATORS, dtype=float)
    coefficients /= float(COMMON_DENOMINATOR)
    return np.polynomial.chebyshev.chebval(2 * normalized_x - 1, coefficients)


def causal_convolution(
    left: np.ndarray, right: np.ndarray, dx: float, nfft: int
) -> np.ndarray:
    return (
        np.fft.irfft(
            np.fft.rfft(left, nfft) * np.fft.rfft(right, nfft), nfft
        )[: len(left)]
        * dx
    )


def supported_score(
    dx: float,
    radial_power: np.ndarray,
) -> dict[str, object]:
    """Return the exploratory 48 J/I score with its nonnegative channels."""
    radius = float(RADIUS)
    outer_radius = float(OUTER_RADIUS)
    delta = float(DELTA)
    b_small = float(B_SMALL)
    b_large = float(B_LARGE)
    intervals = round(radius / dx)
    if abs(intervals * dx - radius) > 1.0e-13:
        raise ValueError("the stress mesh must divide 521/2000")

    x = np.arange(intervals + 1) * dx
    normalized_x = x / radius
    profile = evaluate_profile(normalized_x)
    squared_profile = profile**2
    profile_mass = float(np.sum(squared_profile) * dx)
    small = squared_profile * (x <= delta + 1.0e-14)
    big = squared_profile * (x > delta + 1.0e-14)
    small_mass = float(np.sum(small) * dx)
    big_mass = float(np.sum(big) * dx)
    small /= small_mass
    big /= big_mass
    probability_small = small_mass / profile_mass
    probability_big = big_mass / profile_mass
    nfft = 1 << (2 * (intervals + 1) - 1).bit_length()

    small_convolutions: list[np.ndarray | None] = [None] * (K + 1)
    small_convolutions[1] = small
    for order in range(2, K + 1):
        previous = small_convolutions[order - 1]
        assert previous is not None
        small_convolutions[order] = causal_convolution(
            previous, small, dx, nfft
        )

    big_convolutions: list[np.ndarray | None] = [None] * 7
    big_convolutions[1] = big
    for order in range(2, 7):
        previous = big_convolutions[order - 1]
        assert previous is not None
        big_convolutions[order] = causal_convolution(previous, big, dx, nfft)

    allowed_density_by_count = np.zeros((7, intervals + 1))
    for big_count in range(7):
        small_count = K - big_count
        component_weight = (
            comb(K, big_count)
            * probability_big**big_count
            * probability_small**small_count
        )
        small_density = small_convolutions[small_count]
        assert small_density is not None
        if big_count == 0:
            component_density = small_density
        else:
            big_density = big_convolutions[big_count]
            assert big_density is not None
            big_density = big_density.copy()
            boundary = b_small if big_count <= 2 else b_large
            big_density[x > boundary + 1.0e-14] = 0.0
            component_density = causal_convolution(
                big_density, small_density, dx, nfft
            )
        allowed_density_by_count[big_count] = (
            component_weight * component_density
        )

    allowed_density = np.sum(allowed_density_by_count, axis=0)

    radial_factor = np.polynomial.polynomial.polyval(
        normalized_x, radial_power
    )
    expectation_i = float(np.sum(allowed_density * radial_factor**2) * dx)
    expectation_i_by_count = np.sum(
        allowed_density_by_count * radial_factor[np.newaxis, :] ** 2,
        axis=1,
    ) * dx

    radial_degree = len(radial_power) - 1
    moments = np.zeros((radial_degree + 1, intervals + 1))
    for degree in range(radial_degree + 1):
        moments[degree] = np.cumsum(normalized_x**degree * profile) * dx

    binomial_derivatives = np.zeros_like(moments)
    for derivative in range(radial_degree + 1):
        coefficients = [
            radial_power[degree] * comb(degree, derivative)
            for degree in range(derivative, len(radial_power))
        ]
        binomial_derivatives[derivative] = np.polynomial.polynomial.polyval(
            normalized_x, coefficients
        )

    delta_index = round(delta / dx)
    outer_index = round(outer_radius / dx)

    def interval_integral(
        sum_indices: np.ndarray, upper_index: int, lower_index: int = 0
    ) -> np.ndarray:
        interval_moments = moments[:, upper_index].copy()
        if lower_index:
            interval_moments -= moments[:, lower_index - 1]
        return interval_moments @ binomial_derivatives[:, sum_indices]

    all_indices = np.arange(intervals + 1)
    small_channel = np.zeros(intervals + 1)
    constant_small = all_indices <= intervals - delta_index
    small_channel[constant_small] = interval_integral(
        all_indices[constant_small], delta_index
    )
    for sum_index in all_indices[~constant_small]:
        small_channel[sum_index] = interval_integral(
            np.asarray([sum_index]), intervals - sum_index
        )[0]

    big_boundary_channel = np.zeros(intervals + 1)
    for sum_index in range(intervals - delta_index + 1):
        big_boundary_channel[sum_index] = interval_integral(
            np.asarray([sum_index]),
            intervals - sum_index,
            delta_index + 1,
        )[0]

    # The last-coordinate amplitude is A+B, where A is its small-coordinate
    # channel and B is its big-coordinate channel.  All three pieces of
    # (A+B)^2 must be retained: the eventual exact certificate cannot infer
    # the variational inequality from A^2 alone.
    expectation_j_by_count = np.zeros((7, 3))
    quantized_union_j_by_count = np.zeros(7)
    for big_count in range(7):
        small_count = K - 1 - big_count
        component_weight = (
            comb(K - 1, big_count)
            * probability_big**big_count
            * probability_small**small_count
        )
        small_density = small_convolutions[small_count]
        assert small_density is not None
        if big_count == 0:
            big_indices = np.asarray([0])
            big_weights = np.asarray([1.0])
        else:
            big_density = big_convolutions[big_count]
            assert big_density is not None
            big_indices = np.flatnonzero(big_density > 1.0e-300)
            big_weights = big_density[big_indices] * dx

        current_boundary = b_small if 0 < big_count <= 2 else b_large
        next_boundary = b_small if big_count + 1 <= 2 else b_large
        for big_index, big_weight in zip(big_indices, big_weights):
            max_small_index = min(
                outer_index - big_index, intervals - big_index
            )
            if max_small_index < 0:
                continue
            small_indices = np.arange(max_small_index + 1)
            sum_indices = big_index + small_indices
            small_last_channel = np.zeros(max_small_index + 1)
            big_last_channel = np.zeros(max_small_index + 1)

            if big_count == 0 or x[big_index] <= current_boundary + 1.0e-14:
                small_last_channel += small_channel[sum_indices]

            upper = next_boundary - x[big_index]
            if upper > delta:
                upper_index = min(
                    intervals, int(math.floor((upper + 1.0e-14) / dx))
                )
                constant_big = (
                    x[small_indices] < radius - next_boundary - 1.0e-14
                )
                if np.any(constant_big):
                    big_last_channel[constant_big] += interval_integral(
                        sum_indices[constant_big],
                        upper_index,
                        delta_index + 1,
                    )
                if np.any(~constant_big):
                    big_last_channel[~constant_big] += big_boundary_channel[
                        sum_indices[~constant_big]
                    ]

            common = (
                component_weight
                * big_weight
                * small_density[small_indices]
                * dx
            )
            expectation_j_by_count[big_count, 0] += float(
                np.sum(common * small_last_channel**2)
            )
            expectation_j_by_count[big_count, 1] += float(
                np.sum(
                    common
                    * (2.0 * small_last_channel * big_last_channel)
                )
            )
            expectation_j_by_count[big_count, 2] += float(
                np.sum(common * big_last_channel**2)
            )

            # Rigorous-certificate candidate: on the subregion where the
            # complete small interval [0,delta] is available, round the
            # guaranteed union endpoint down to a fixed rational grid.  The
            # resulting [0,a_cell] amplitude is pointwise no larger than the
            # true A+B amplitude because q and r are positive.
            union_lower_channel = np.zeros(max_small_index + 1)
            total_values = x[sum_indices]
            current_big_allowed = (
                big_count == 0
                or x[big_index] <= current_boundary + 1.0e-14
            )
            if current_big_allowed:
                eligible = total_values < radius - 1.0e-14
                total_endpoint = radius - total_values
                big_endpoint = next_boundary - x[big_index]
                guaranteed_endpoint = np.where(
                    total_endpoint < delta,
                    total_endpoint,
                    np.maximum(
                        delta,
                        np.minimum(total_endpoint, big_endpoint),
                    ),
                )
                below_delta = guaranteed_endpoint < delta
                quantized_endpoint = np.empty_like(guaranteed_endpoint)
                quantized_endpoint[below_delta] = (
                    np.floor(
                        (guaranteed_endpoint[below_delta] + 1.0e-14)
                        / QUANTIZED_UNION_STEP
                    )
                    * QUANTIZED_UNION_STEP
                )
                extra_steps = np.floor(
                    (
                        guaranteed_endpoint[~below_delta]
                        - delta
                        + 1.0e-14
                    )
                    / QUANTIZED_UNION_STEP
                ).astype(int)
                quantized_endpoint[~below_delta] = (
                    delta + extra_steps * QUANTIZED_UNION_STEP
                )
                upper_indices = np.minimum(
                    intervals,
                    np.floor(
                        (quantized_endpoint + 1.0e-14) / dx
                    ).astype(int),
                )
                for upper_index in np.unique(upper_indices[eligible]):
                    mask = eligible & (upper_indices == upper_index)
                    union_lower_channel[mask] = interval_integral(
                        sum_indices[mask], int(upper_index)
                    )
            quantized_union_j_by_count[big_count] += float(
                np.sum(common * union_lower_channel**2)
            )

    normalization = K / (profile_mass * expectation_i)
    normalized = expectation_j_by_count * normalization
    normalized_quantized_union = quantized_union_j_by_count * normalization
    channel_totals = np.sum(normalized, axis=0)
    total = float(np.sum(channel_totals))
    return {
        "mesh": dx,
        "score48JOverI": total,
        "integralDiagnostics": {
            "profileMassPhysical": profile_mass,
            "allowedNormalizedExpectationI": expectation_i,
            "physicalI": profile_mass**K * expectation_i,
            "physicalIByBigCoordinateCount": [
                float(profile_mass**K * value)
                for value in expectation_i_by_count
            ],
        },
        "channels": {
            "smallSquare": float(channel_totals[0]),
            "twiceSmallBigCross": float(channel_totals[1]),
            "bigSquare": float(channel_totals[2]),
            "sum": total,
        },
        "byOuterBigCoordinateCount": [
            {
                "bigCount": big_count,
                "smallSquare": float(normalized[big_count, 0]),
                "twiceSmallBigCross": float(normalized[big_count, 1]),
                "bigSquare": float(normalized[big_count, 2]),
                "sum": float(np.sum(normalized[big_count])),
            }
            for big_count in range(7)
        ],
        "quantizedUnionLowerCandidate": {
            "physicalEndpointStep": QUANTIZED_UNION_STEP,
            "score48JLowerOverI": float(
                np.sum(normalized_quantized_union)
            ),
            "byOuterBigCoordinateCount": [
                float(value) for value in normalized_quantized_union
            ],
            "warning": (
                "This remains a floating outer quadrature diagnostic.  Its "
                "inner amplitude is a pointwise analytic lower channel, but "
                "the displayed score is not yet an interval enclosure."
            ),
        },
    }


def candidate_digest() -> str:
    canonical = json.dumps(
        {
            "denominator": COMMON_DENOMINATOR,
            "profileNumerators": PROFILE_NUMERATORS,
            "radialNumerators": RADIAL_NUMERATORS,
            "support": {
                "radius": encode(RADIUS),
                "outerRadius": encode(OUTER_RADIUS),
                "delta": encode(DELTA),
                "bSmall": encode(B_SMALL),
                "bLarge": encode(B_LARGE),
            },
        },
        sort_keys=True,
        separators=(",", ":"),
    ).encode("ascii")
    return hashlib.sha256(canonical).hexdigest()


def main() -> None:
    profile_power = combine_basis(
        PROFILE_NUMERATORS, shifted_chebyshev_basis(32)
    )
    radial_power_exact = combine_basis(
        RADIAL_NUMERATORS, shifted_legendre_basis(8)
    )
    profile_positivity = bernstein_positivity_certificate(profile_power, 32)
    radial_positivity = bernstein_positivity_certificate(radial_power_exact, 8)
    radial_power = np.asarray([float(value) for value in radial_power_exact])

    scores = [supported_score(mesh, radial_power) for mesh in STRESS_MESHES]
    values = [item["score48JOverI"] for item in scores]
    first_order = [2 * values[index + 1] - values[index] for index in range(3)]
    second_order = (4 * first_order[-1] - first_order[-2]) / 3
    small_channel_scores = [
        {
            "mesh": item["mesh"],
            "score48JSmallSquareOverI": item["channels"]["smallSquare"],
        }
        for item in scores
    ]
    small_values = [
        item["score48JSmallSquareOverI"] for item in small_channel_scores
    ]
    small_first_order = [
        2 * small_values[index + 1] - small_values[index]
        for index in range(3)
    ]
    small_second_order = (
        4 * small_first_order[-1] - small_first_order[-2]
    ) / 3

    sample = np.linspace(0.0, 1.0, 1_000_001)
    profile_sample = evaluate_profile(sample)
    radial_sample = np.polynomial.polynomial.polyval(sample, radial_power)

    payload = {
        "name": "Frozen rational product-radial candidate for k=48",
        "claimBoundary": {
            "status": "candidate-only; not a proof",
            "statement": (
                "The rational coefficients and Bernstein positivity certificates "
                "are exact. The convolution scores use IEEE-754 double precision, "
                "rectangle quadrature, FFT convolution, and Richardson diagnostics. "
                "They are not interval enclosures and do not prove 48*J(F)>I(F), "
                "DHL[48,2], or H_1<=236."
            ),
        },
        "definition": {
            "k": K,
            "formula": (
                "F(t)=1_Omega(t)*r((sum_i t_i)/(521/2000))*"
                "product_i q(t_i/(521/2000))"
            ),
            "profileFormula": "q(u)=sum_{j=0}^{32} a_j*T_j(2u-1)",
            "radialFormula": "r(u)=sum_{j=0}^{8} b_j*P_j(2u-1)",
            "support": {
                "radius": encode(RADIUS),
                "restrictedOuterRadiusUsedInJStressTest": encode(OUTER_RADIUS),
                "largeCoordinateThreshold": encode(DELTA),
                "largeCoordinateSumBoundForOneOrTwo": encode(B_SMALL),
                "largeCoordinateSumBoundForThreeThroughSix": encode(B_LARGE),
                "sevenLargeCoordinatesImpossible": bool(7 * DELTA > B_LARGE),
            },
        },
        "rationalCandidate": {
            "commonDenominator": COMMON_DENOMINATOR,
            "sha256": candidate_digest(),
            "profile": {
                "basis": "shifted Chebyshev T_j(2u-1)",
                "degree": 32,
                "numerators": PROFILE_NUMERATORS,
                "coefficients": [
                    encode(Fraction(value, COMMON_DENOMINATOR))
                    for value in PROFILE_NUMERATORS
                ],
                "powerBasisCoefficients": [
                    encode(value) for value in profile_power
                ],
            },
            "radial": {
                "basis": "shifted Legendre P_j(2u-1)",
                "degree": 8,
                "numerators": RADIAL_NUMERATORS,
                "coefficients": [
                    encode(Fraction(value, COMMON_DENOMINATOR))
                    for value in RADIAL_NUMERATORS
                ],
                "powerBasisCoefficients": [
                    encode(value) for value in radial_power_exact
                ],
            },
        },
        "exactBernsteinPositivityCertificates": {
            "profileOnUnitInterval": profile_positivity,
            "radialOnUnitInterval": radial_positivity,
        },
        "floatingPointStressTest": {
            "arithmetic": "IEEE-754 float64 with NumPy FFT",
            "quadrature": "aligned left rectangle rule used consistently",
            "meshScores": scores,
            "firstOrderRichardsonSequence": first_order,
            "secondOrderRichardsonEstimate": second_order,
            "sampledProfileMinimum": float(np.min(profile_sample)),
            "sampledProfileMinimumLocation": float(
                sample[int(np.argmin(profile_sample))]
            ),
            "sampledRadialMinimum": float(np.min(radial_sample)),
            "sampledRadialMinimumLocation": float(
                sample[int(np.argmin(radial_sample))]
            ),
            "candidateAppearsAboveThreshold": bool(second_order > 1.0),
            "smallLastCoordinateLowerChannel": {
                "definition": (
                    "Retain only the A^2 term from the t_k<=delta amplitude. "
                    "It is nonnegative but is numerically insufficient for "
                    "48*J/I>1.  The exact certificate must also retain 2*A*B "
                    "and B^2 for every outer big-coordinate count 0 through 6."
                ),
                "meshScores": small_channel_scores,
                "firstOrderRichardsonSequence": small_first_order,
                "secondOrderRichardsonEstimate": small_second_order,
                "candidateAppearsAboveThreshold": bool(
                    small_second_order > 1.0
                ),
            },
            "nonProofWarning": (
                "Mesh refinement and Richardson extrapolation are diagnostics only. "
                "A proof still requires exact or outward-rounded evaluation of I and J."
            ),
        },
        "remainingExactProofTarget": (
            "Compute I=N_I/D_I and the full nonnegative truncated integral "
            "J_{<=R_o}=N_J/D_J for this exact rational candidate, retaining "
            "A^2, 2*A*B, and B^2 for outer big-coordinate counts 0 through 6, "
            "then verify 48*N_J*D_I-N_I*D_J>0 as an integer inequality."
        ),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), "sha256": candidate_digest()}, indent=2))


if __name__ == "__main__":
    main()
