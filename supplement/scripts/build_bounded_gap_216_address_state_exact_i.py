#!/usr/bin/env python3
"""Exactly evaluate I_46 for the frozen 208-dimensional address state.

The calculation uses FLINT rational polynomials.  It certifies only the
denominator integral I.  The companion J calculation and the analytic S52
input remain separate proof obligations.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

from flint import __version__ as flint_version
from flint import fmpq, fmpq_poly

from build_bounded_gap_226_layered_rational_candidate import (
    B_BY_COUNT,
    DELTA,
    MAX_BIG_COUNT,
    PROFILE_NUMERATORS,
    RADIUS,
    candidate_digest as support_digest,
)
from build_bounded_gap_236_exact_convolution_backbone import (
    borel,
    encode,
    inverse_borel,
    q,
    shifted,
    to_flint_poly,
)
from build_bounded_gap_236_rational_candidate import (
    combine_basis,
    shifted_chebyshev_basis,
    shifted_legendre_basis,
)


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = (
    ROOT / "outputs" / "bounded_gap_216_address_state_rational_candidate.json"
)
OUTPUT = ROOT / "outputs" / "bounded_gap_216_address_state_exact_i.json"

K = 46
Bivar = dict[tuple[int, int], fmpq]
TransformPieces = dict[int, fmpq_poly]
FACTORIALS = [1]


def rational_digest(value: fmpq) -> str:
    digest = hashlib.sha256()
    digest.update(str(value.numerator).encode("ascii"))
    digest.update(b"/")
    digest.update(str(value.denominator).encode("ascii"))
    return digest.hexdigest()


def monomial(power: int) -> fmpq_poly:
    return fmpq_poly([fmpq(0)] * power + [fmpq(1)])


def ensure_factorials(maximum: int) -> None:
    while len(FACTORIALS) <= maximum:
        FACTORIALS.append(FACTORIALS[-1] * len(FACTORIALS))


def cached_borel(poly: fmpq_poly) -> fmpq_poly:
    ensure_factorials(len(poly) - 1)
    return fmpq_poly(
        [fmpq(0)]
        + [poly[degree] * FACTORIALS[degree] for degree in range(len(poly))]
    )


def cached_inverse_borel(poly: fmpq_poly) -> fmpq_poly:
    if poly[0] != 0:
        raise ValueError("inverse Borel input must have zero constant term")
    ensure_factorials(len(poly) - 2)
    return fmpq_poly(
        [
            poly[degree + 1] / FACTORIALS[degree]
            for degree in range(len(poly) - 1)
        ]
    )


def transform_add(
    left: TransformPieces,
    right: TransformPieces,
    scale: fmpq = fmpq(1),
) -> TransformPieces:
    result: dict[int, fmpq_poly] = dict(left)
    for shift_index, poly in right.items():
        value = result.get(shift_index, fmpq_poly([0])) + scale * poly
        if value == 0:
            result.pop(shift_index, None)
        else:
            result[shift_index] = value
    return result


def transform_multiply(
    left: TransformPieces,
    right: TransformPieces,
    maximum_shift: int,
) -> TransformPieces:
    result: dict[int, fmpq_poly] = {}
    for left_shift, left_poly in left.items():
        for right_shift, right_poly in right.items():
            shift_index = left_shift + right_shift
            if shift_index > maximum_shift:
                continue
            value = result.get(shift_index, fmpq_poly([0]))
            value += left_poly * right_poly
            result[shift_index] = value
    return {key: value for key, value in result.items() if value != 0}


def two_piece_power(
    base: TransformPieces, order: int, maximum_shift: int
) -> TransformPieces:
    if order == 0:
        return {0: fmpq_poly([1])}
    if set(base) != {0, 1}:
        raise ValueError("two_piece_power requires shifts zero and one")
    return {
        shifted_count: (
            math.comb(order, shifted_count)
            * base[0] ** (order - shifted_count)
            * base[1] ** shifted_count
        )
        for shifted_count in range(min(order, maximum_shift) + 1)
    }


def compact_coordinate_transform(
    density: fmpq_poly, cutoff: fmpq
) -> TransformPieces:
    return {
        0: borel(density),
        1: -borel(shifted(density, cutoff)),
    }


def compact_group_moments(
    densities: tuple[fmpq_poly, fmpq_poly, fmpq_poly],
    cutoff: fmpq,
    order: int,
    maximum_shift: int,
) -> list[dict[int, fmpq_poly]]:
    transforms = [
        compact_coordinate_transform(density, cutoff)
        for density in densities
    ]
    base_n = two_piece_power(transforms[0], order, maximum_shift)
    base_n_minus_1 = two_piece_power(
        transforms[0], order - 1, maximum_shift
    )
    first = transform_multiply(
        transforms[1], base_n_minus_1, maximum_shift
    )
    first = {key: order * value for key, value in first.items()}

    second = transform_multiply(
        transforms[2], base_n_minus_1, maximum_shift
    )
    second = {key: order * value for key, value in second.items()}
    if order >= 2:
        marked_pair = transform_multiply(
            transforms[1], transforms[1], maximum_shift
        )
        marked_pair = transform_multiply(
            marked_pair,
            two_piece_power(transforms[0], order - 2, maximum_shift),
            maximum_shift,
        )
        second = transform_add(
            second, marked_pair, fmpq(order * (order - 1))
        )

    return [
        {
            key: inverse_borel(value)
            for key, value in transform.items()
        }
        for transform in (base_n, first, second)
    ]


def big_group_moments(
    densities: tuple[fmpq_poly, fmpq_poly, fmpq_poly],
    cutoff: fmpq,
    order: int,
) -> list[fmpq_poly]:
    if order < 1:
        raise ValueError("a big group must be nonempty")
    transforms = [borel(shifted(density, cutoff)) for density in densities]
    base = transforms[0]
    zeroth = base**order
    first = order * transforms[1] * base ** (order - 1)
    second = order * transforms[2] * base ** (order - 1)
    if order >= 2:
        second += (
            order
            * (order - 1)
            * transforms[1]
            * transforms[1]
            * base ** (order - 2)
        )
    return [inverse_borel(value) for value in (zeroth, first, second)]


def bivar_add_term(
    target: Bivar, u_power: int, b_power: int, value: fmpq
) -> None:
    if value == 0:
        return
    key = (u_power, b_power)
    target[key] = target.get(key, fmpq(0)) + value
    if target[key] == 0:
        del target[key]


def bivar_multiply(left: Bivar, right: Bivar) -> Bivar:
    result: defaultdict[tuple[int, int], fmpq] = defaultdict(fmpq)
    for (left_u, left_b), left_value in left.items():
        for (right_u, right_b), right_value in right.items():
            result[(left_u + right_u, left_b + right_b)] += (
                left_value * right_value
            )
    return {key: value for key, value in result.items() if value != 0}


def bivar_scale(poly: Bivar, scalar: fmpq) -> Bivar:
    return {key: scalar * value for key, value in poly.items() if value != 0}


def combine_flint_basis(
    numerators: list[int], denominator: int, basis: list[fmpq_poly]
) -> fmpq_poly:
    result = fmpq_poly([0])
    for numerator, basis_poly in zip(numerators, basis):
        result += fmpq(numerator, denominator) * basis_poly
    return result


def decode_amplitudes(
    binding: dict[str, object], cutoff: fmpq
) -> list[tuple[Bivar, fmpq_poly]]:
    numerators = [int(value) for value in binding["numerators"]]
    denominator = int(binding["commonDenominator"])
    layer_degrees = [int(value) for value in binding["layerDegrees"]]
    mass_degrees = [int(value) for value in binding["massLayerDegrees"]]
    dispersion_degrees = [
        int(value) for value in binding["dispersionLayerDegrees"]
    ]
    radial_basis = [
        to_flint_poly(poly) for poly in shifted_legendre_basis(4)
    ]
    address_basis = radial_basis
    cursor = 0
    amplitudes: list[tuple[Bivar, fmpq_poly]] = []

    for layer in range(len(layer_degrees)):
        ordinary: Bivar = {}
        degree = layer_degrees[layer]
        radial = combine_flint_basis(
            numerators[cursor : cursor + degree + 1],
            denominator,
            radial_basis,
        )
        cursor += degree + 1
        for u_power, value in enumerate(radial):
            bivar_add_term(ordinary, u_power, 0, value)

        mass_degree = mass_degrees[layer]
        if mass_degree >= 0:
            boundary = q(B_BY_COUNT[layer] / RADIUS)
            headroom = boundary - layer * cutoff
            if layer == 0 or headroom <= 0:
                raise ArithmeticError("invalid positive-layer address interval")
            address_coordinate = fmpq_poly(
                [-layer * cutoff / headroom, fmpq(1) / headroom]
            )
            for mode in range(1, 5):
                radial = combine_flint_basis(
                    numerators[cursor : cursor + mass_degree + 1],
                    denominator,
                    radial_basis,
                )
                cursor += mass_degree + 1
                mass_poly = address_basis[mode](address_coordinate)
                for u_power, u_value in enumerate(radial):
                    for b_power, b_value in enumerate(mass_poly):
                        bivar_add_term(
                            ordinary,
                            u_power,
                            b_power,
                            u_value * b_value,
                        )

        dispersion_degree = dispersion_degrees[layer]
        if dispersion_degree >= 0:
            dispersion = combine_flint_basis(
                numerators[cursor : cursor + dispersion_degree + 1],
                denominator,
                radial_basis,
            )
            cursor += dispersion_degree + 1
        else:
            dispersion = fmpq_poly([0])
        amplitudes.append((ordinary, dispersion))

    if cursor != len(numerators):
        raise ArithmeticError(
            f"decoded {cursor} coefficients, expected {len(numerators)}"
        )
    return amplitudes


def amplitude_powers(ordinary: Bivar, dispersion: fmpq_poly) -> list[Bivar]:
    dispersion_bivar = {
        (degree, 0): value
        for degree, value in enumerate(dispersion)
        if value != 0
    }
    return [
        bivar_multiply(ordinary, ordinary),
        bivar_scale(bivar_multiply(ordinary, dispersion_bivar), fmpq(2)),
        bivar_multiply(dispersion_bivar, dispersion_bivar),
    ]


def integral_interval(poly: fmpq_poly, upper: fmpq) -> fmpq:
    if upper <= 0:
        return fmpq(0)
    antiderivative = poly.integral()
    return antiderivative(upper) - antiderivative(0)


def convolve_from_left_borel(
    left_borel: fmpq_poly, right: fmpq_poly
) -> fmpq_poly:
    return cached_inverse_borel(left_borel * cached_borel(right))


def integrate_sum_amplitude(
    density: fmpq_poly,
    amplitude: fmpq_poly,
    offset: fmpq,
    upper: fmpq,
) -> fmpq:
    if upper <= 0 or density == 0 or amplitude == 0:
        return fmpq(0)
    return integral_interval(density * shifted(amplitude, offset), upper)


def integrate_delta_big(
    small: fmpq_poly,
    amplitude: Bivar,
    total_offset: fmpq,
    total_cap: fmpq,
) -> fmpq:
    if total_cap <= 0:
        return fmpq(0)
    univariate = fmpq_poly([0])
    affine_total = fmpq_poly([total_offset, 1])
    for (u_power, b_power), coefficient in amplitude.items():
        if b_power == 0:
            univariate += coefficient * affine_total**u_power
    return integral_interval(small * univariate, total_cap)


def integrate_triangle_bivar(
    small: fmpq_poly,
    big: fmpq_poly,
    amplitude: Bivar,
    total_offset: fmpq,
    big_offset: fmpq,
    total_cap: fmpq,
    big_cap: fmpq,
) -> fmpq:
    if total_cap <= 0 or big_cap <= 0 or not amplitude:
        return fmpq(0)
    maximum_u_power = max(key[0] for key in amplitude)
    big_affine = fmpq_poly([big_offset, 1])
    maximum_b_power = max(key[1] for key in amplitude)
    big_powers = [big_affine**power for power in range(maximum_b_power + 1)]
    grouped_by_b = [fmpq_poly([0]) for _ in range(maximum_b_power + 1)]
    for (u_power, b_power), coefficient in amplitude.items():
        grouped_by_b[b_power] += coefficient * monomial(u_power)

    # Integrate by total sum z=x+y.  The Borel transform turns convolution
    # into polynomial multiplication.  Capping y at Y is inclusion-exclusion:
    # subtract the translated tail y=Y+y'.  This avoids composing a degree
    # roughly 3000 antiderivative with T-y.
    left_borel = cached_borel(small)
    result = fmpq(0)
    for b_power, u_amplitude in enumerate(grouped_by_b):
        if u_amplitude == 0:
            continue
        weighted_big = big * big_powers[b_power]
        total_density = convolve_from_left_borel(left_borel, weighted_big)
        result += integrate_sum_amplitude(
            total_density, u_amplitude, total_offset, total_cap
        )
        if big_cap < total_cap:
            tail_big = shifted(big, big_cap)
            tail_affine = fmpq_poly([big_offset + big_cap, 1])
            tail_big *= tail_affine**b_power
            tail_density = convolve_from_left_borel(left_borel, tail_big)
            result -= integrate_sum_amplitude(
                tail_density,
                u_amplitude,
                total_offset + big_cap,
                total_cap - big_cap,
            )
    return result


def exact_layer_integral(
    layer: int,
    density_moments: tuple[fmpq_poly, fmpq_poly, fmpq_poly],
    cutoff: fmpq,
    amplitudes: list[Bivar],
) -> fmpq:
    small_count = K - layer
    maximum_shift = int((Fraction(1) - layer * (DELTA / RADIUS)) / (DELTA / RADIUS))
    small_moments = compact_group_moments(
        density_moments, cutoff, small_count, maximum_shift
    )

    if layer == 0:
        subtotal = fmpq(0)
        for moment_power, amplitude in enumerate(amplitudes):
            for shift_index, small_poly in small_moments[moment_power].items():
                total_offset = shift_index * cutoff
                subtotal += integrate_delta_big(
                    small_poly,
                    amplitude,
                    total_offset,
                    fmpq(1) - total_offset,
                )
        return subtotal

    big_moments = big_group_moments(density_moments, cutoff, layer)
    big_offset = layer * cutoff
    big_cap = q(B_BY_COUNT[layer] / RADIUS) - big_offset
    subtotal = fmpq(0)
    for moment_power, amplitude in enumerate(amplitudes):
        for small_power in range(moment_power + 1):
            big_power = moment_power - small_power
            multiplicity = math.comb(moment_power, small_power)
            for shift_index, small_poly in small_moments[small_power].items():
                total_offset = big_offset + shift_index * cutoff
                subtotal += multiplicity * integrate_triangle_bivar(
                    small_poly,
                    big_moments[big_power],
                    amplitude,
                    total_offset,
                    big_offset,
                    fmpq(1) - total_offset,
                    big_cap,
                )
    return math.comb(K, layer) * subtotal


def main() -> None:
    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    binding = candidate["binding"]
    if candidate["sha256"] != hashlib.sha256(
        json.dumps(binding, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest():
        raise ArithmeticError("address-state binding hash mismatch")
    if int(binding["k"]) != K:
        raise ArithmeticError("wrong k in address-state candidate")
    if binding["supportDigest"] != support_digest():
        raise ArithmeticError("address-state support digest mismatch")
    if binding["outerRadius"] != "2509/10000":
        raise ArithmeticError("address-state outer radius metadata mismatch")

    profile = to_flint_poly(
        combine_basis(PROFILE_NUMERATORS, shifted_chebyshev_basis(32))
    )
    density = profile * profile
    variable = monomial(1)
    density_moments = (
        density,
        variable**2 * density,
        variable**4 * density,
    )
    cutoff = q(DELTA / RADIUS)
    decoded = decode_amplitudes(binding, cutoff)

    layer_values: list[fmpq] = []
    for layer, (ordinary, dispersion) in enumerate(decoded):
        value = exact_layer_integral(
            layer,
            density_moments,
            cutoff,
            amplitude_powers(ordinary, dispersion),
        )
        layer_values.append(value)
        print(
            json.dumps(
                {
                    "layer": layer,
                    "decimal": float(value),
                    "numeratorDigits": len(str(value.numerator)),
                    "denominatorDigits": len(str(value.denominator)),
                }
            ),
            flush=True,
        )

    normalized_i = sum(layer_values, fmpq(0))
    physical_i = q(RADIUS) ** K * normalized_i
    payload = {
        "name": "Exact rational I_46 for the prime-address moment state",
        "candidateSha256": candidate["sha256"],
        "supportSha256": support_digest(),
        "claimBoundary": (
            "This artifact certifies I_46 only. Exact J_46 positivity and the "
            "source-level S52 proof are separate obligations."
        ),
        "backend": {"pythonFlintVersion": flint_version},
        "formula": (
            "I=R^46 sum_m binom(46,m) integral "
            "|A_m(U,B)+D C_m(U)|^2 product q(x_i)^2 dx"
        ),
        "geometry": {
            "radius": encode(RADIUS),
            "delta": encode(DELTA),
            "normalizedCutoff": encode(DELTA / RADIUS),
            "maximumBigCoordinateCount": MAX_BIG_COUNT,
        },
        "normalizedContributionsByBigCoordinateCount": [
            {
                "bigCount": layer,
                "exact": encode(value),
                "decimal": float(value),
                "sha256": rational_digest(value),
            }
            for layer, value in enumerate(layer_values)
        ],
        "normalizedI": {
            "exact": encode(normalized_i),
            "decimal": float(normalized_i),
            "sha256": rational_digest(normalized_i),
            "numeratorDigits": len(str(normalized_i.numerator)),
            "denominatorDigits": len(str(normalized_i.denominator)),
        },
        "physicalI": {
            "exact": encode(physical_i),
            "decimal": float(physical_i),
            "sha256": rational_digest(physical_i),
        },
        "checks": {
            "allLayerContributionsPositive": all(value > 0 for value in layer_values),
            "sumMatchesRecordedTotal": sum(layer_values, fmpq(0)) == normalized_i,
            "physicalScalingIsRTo46": physical_i == q(RADIUS) ** K * normalized_i,
            "coefficientCountIs208": len(binding["numerators"]) == 208,
        },
        "nextExactStage": (
            "Evaluate the conservative AA, 2AB, and BB pieces of J_46 with "
            "the same address and dispersion moments, then prove "
            "46*R*J-I>0 as an exact rational inequality."
        ),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(OUTPUT),
                "normalizedIDecimal": float(normalized_i),
                "checks": payload["checks"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
