#!/usr/bin/env python3
"""Rigorous Arb enclosure of I_46 for the frozen address-state candidate."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

from flint import arb, arb_poly, ctx, fmpq_poly

import build_bounded_gap_216_address_state_exact_i as exact
from build_bounded_gap_226_layered_rational_candidate import (
    B_BY_COUNT,
    DELTA,
    MAX_BIG_COUNT,
    PROFILE_NUMERATORS,
    RADIUS,
    candidate_digest as support_digest,
)
from build_bounded_gap_236_exact_convolution_backbone import q, to_flint_poly
from build_bounded_gap_236_rational_candidate import (
    combine_basis,
    shifted_chebyshev_basis,
)


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = exact.CANDIDATE
OUTPUT = ROOT / "outputs" / "bounded_gap_216_address_state_i_arb.json"
K = 46

BallBivar = dict[tuple[int, int], arb]
BallPieces = dict[int, arb_poly]
FACTORIALS = [1]


def a(value: object) -> arb:
    return arb(value)


def to_ball_poly(poly: fmpq_poly) -> arb_poly:
    return arb_poly([a(value) for value in poly])


def ensure_factorials(maximum: int) -> None:
    while len(FACTORIALS) <= maximum:
        FACTORIALS.append(FACTORIALS[-1] * len(FACTORIALS))


def borel(poly: arb_poly) -> arb_poly:
    ensure_factorials(len(poly) - 1)
    return arb_poly(
        [a(0)]
        + [poly[degree] * FACTORIALS[degree] for degree in range(len(poly))]
    )


def inverse_borel(poly: arb_poly) -> arb_poly:
    ensure_factorials(len(poly) - 2)
    return arb_poly(
        [
            poly[degree + 1] / FACTORIALS[degree]
            for degree in range(len(poly) - 1)
        ]
    )


def shifted(poly: arb_poly, offset: arb) -> arb_poly:
    return poly(arb_poly([offset, a(1)]))


def monomial(power: int) -> arb_poly:
    return arb_poly([a(0)] * power + [a(1)])


def transform_add(
    left: BallPieces, right: BallPieces, scale: int = 1
) -> BallPieces:
    result = dict(left)
    for shift_index, poly in right.items():
        result[shift_index] = result.get(shift_index, arb_poly([0])) + scale * poly
    return result


def transform_multiply(
    left: BallPieces, right: BallPieces, maximum_shift: int
) -> BallPieces:
    result: dict[int, arb_poly] = {}
    for left_shift, left_poly in left.items():
        for right_shift, right_poly in right.items():
            shift_index = left_shift + right_shift
            if shift_index <= maximum_shift:
                result[shift_index] = (
                    result.get(shift_index, arb_poly([0]))
                    + left_poly * right_poly
                )
    return result


def two_piece_power(
    base: BallPieces, order: int, maximum_shift: int
) -> BallPieces:
    if order == 0:
        return {0: arb_poly([1])}
    return {
        shifted_count: (
            math.comb(order, shifted_count)
            * base[0] ** (order - shifted_count)
            * base[1] ** shifted_count
        )
        for shifted_count in range(min(order, maximum_shift) + 1)
    }


def compact_group_moments(
    densities: tuple[arb_poly, arb_poly, arb_poly],
    cutoff: arb,
    order: int,
    maximum_shift: int,
) -> list[BallPieces]:
    coordinate = [
        {0: borel(density), 1: -borel(shifted(density, cutoff))}
        for density in densities
    ]
    base_n = two_piece_power(coordinate[0], order, maximum_shift)
    base_n1 = two_piece_power(coordinate[0], order - 1, maximum_shift)
    first = transform_multiply(coordinate[1], base_n1, maximum_shift)
    first = {key: order * value for key, value in first.items()}
    second = transform_multiply(coordinate[2], base_n1, maximum_shift)
    second = {key: order * value for key, value in second.items()}
    if order >= 2:
        pair = transform_multiply(coordinate[1], coordinate[1], maximum_shift)
        pair = transform_multiply(
            pair,
            two_piece_power(coordinate[0], order - 2, maximum_shift),
            maximum_shift,
        )
        second = transform_add(second, pair, order * (order - 1))
    return [
        {key: inverse_borel(poly) for key, poly in pieces.items()}
        for pieces in (base_n, first, second)
    ]


def big_group_moments(
    densities: tuple[arb_poly, arb_poly, arb_poly],
    cutoff: arb,
    order: int,
) -> list[arb_poly]:
    transforms = [borel(shifted(density, cutoff)) for density in densities]
    base = transforms[0]
    zeroth = base**order
    first = order * transforms[1] * base ** (order - 1)
    second = order * transforms[2] * base ** (order - 1)
    if order >= 2:
        second += (
            order
            * (order - 1)
            * transforms[1] ** 2
            * base ** (order - 2)
        )
    return [inverse_borel(poly) for poly in (zeroth, first, second)]


def ball_bivar(poly: exact.Bivar) -> BallBivar:
    return {key: a(value) for key, value in poly.items()}


def ball_amplitude_powers(
    ordinary: exact.Bivar, dispersion: fmpq_poly
) -> list[BallBivar]:
    return [ball_bivar(poly) for poly in exact.amplitude_powers(ordinary, dispersion)]


def integral_interval(poly: arb_poly, upper: arb) -> arb:
    antiderivative = poly.integral()
    return antiderivative(upper) - antiderivative(a(0))


def integrate_delta_big(
    small: arb_poly,
    amplitude: BallBivar,
    total_offset: arb,
    total_cap: arb,
) -> arb:
    affine = arb_poly([total_offset, a(1)])
    univariate = arb_poly([0])
    for (u_power, b_power), coefficient in amplitude.items():
        if b_power == 0:
            univariate += coefficient * affine**u_power
    return integral_interval(small * univariate, total_cap)


def convolve_from_left_borel(left_borel: arb_poly, right: arb_poly) -> arb_poly:
    return inverse_borel(left_borel * borel(right))


def integrate_sum_amplitude(
    density: arb_poly,
    amplitude: arb_poly,
    offset: arb,
    upper: arb,
) -> arb:
    return integral_interval(density * shifted(amplitude, offset), upper)


def integrate_triangle_bivar(
    small: arb_poly,
    big: arb_poly,
    amplitude: BallBivar,
    total_offset: arb,
    big_offset: arb,
    total_cap: arb,
    big_cap: arb,
    has_tail: bool,
) -> arb:
    if not amplitude:
        return a(0)
    maximum_u = max(key[0] for key in amplitude)
    maximum_b = max(key[1] for key in amplitude)
    grouped = [arb_poly([0]) for _ in range(maximum_b + 1)]
    for (u_power, b_power), coefficient in amplitude.items():
        grouped[b_power] += coefficient * monomial(u_power)

    big_affine = arb_poly([big_offset, a(1)])
    left_borel = borel(small)
    result = a(0)
    for b_power, u_amplitude in enumerate(grouped):
        weighted_big = big * big_affine**b_power
        density = convolve_from_left_borel(left_borel, weighted_big)
        result += integrate_sum_amplitude(
            density, u_amplitude, total_offset, total_cap
        )
        if has_tail:
            translated_big = shifted(big, big_cap)
            tail_affine = arb_poly([big_offset + big_cap, a(1)])
            translated_big *= tail_affine**b_power
            tail_density = convolve_from_left_borel(left_borel, translated_big)
            result -= integrate_sum_amplitude(
                tail_density,
                u_amplitude,
                total_offset + big_cap,
                total_cap - big_cap,
            )
    return result


def layer_integral(
    layer: int,
    density_moments: tuple[arb_poly, arb_poly, arb_poly],
    cutoff: arb,
    amplitudes: list[BallBivar],
) -> arb:
    cutoff_fraction = DELTA / RADIUS
    maximum_shift = int(
        (1 - layer * cutoff_fraction) / cutoff_fraction
    )
    small_moments = compact_group_moments(
        density_moments, cutoff, K - layer, maximum_shift
    )
    if layer == 0:
        result = a(0)
        for power, amplitude in enumerate(amplitudes):
            for shift_index, small in small_moments[power].items():
                offset = shift_index * cutoff
                result += integrate_delta_big(
                    small, amplitude, offset, a(1) - offset
                )
        return result

    big_moments = big_group_moments(density_moments, cutoff, layer)
    big_offset_fraction = layer * cutoff_fraction
    big_cap_fraction = B_BY_COUNT[layer] / RADIUS - big_offset_fraction
    big_offset = a(q(big_offset_fraction))
    big_cap = a(q(big_cap_fraction))
    result = a(0)
    for power, amplitude in enumerate(amplitudes):
        for small_power in range(power + 1):
            big_power = power - small_power
            multiplicity = math.comb(power, small_power)
            for shift_index, small in small_moments[small_power].items():
                total_offset_fraction = (
                    big_offset_fraction + shift_index * cutoff_fraction
                )
                total_cap_fraction = 1 - total_offset_fraction
                total_offset = a(q(total_offset_fraction))
                total_cap = a(q(total_cap_fraction))
                result += multiplicity * integrate_triangle_bivar(
                    small,
                    big_moments[big_power],
                    amplitude,
                    total_offset,
                    big_offset,
                    total_cap,
                    big_cap,
                    big_cap_fraction < total_cap_fraction,
                )
    return math.comb(K, layer) * result


def record(value: arb) -> dict[str, object]:
    midpoint, radius, exponent = value.mid_rad_10exp()
    return {
        "ball": str(value),
        "lower": str(value.lower()),
        "upper": str(value.upper()),
        "midRad10Exp": [str(midpoint), str(radius), int(exponent)],
        "strictlyPositive": value.lower() > 0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--precision", type=int, default=384)
    args = parser.parse_args()
    ctx.prec = args.precision

    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    binding = candidate["binding"]
    binding_hash = hashlib.sha256(
        json.dumps(binding, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if binding_hash != candidate["sha256"]:
        raise ArithmeticError("candidate binding hash mismatch")
    if binding["supportDigest"] != support_digest():
        raise ArithmeticError("support hash mismatch")

    exact_profile = to_flint_poly(
        combine_basis(PROFILE_NUMERATORS, shifted_chebyshev_basis(32))
    )
    profile = to_ball_poly(exact_profile)
    density = profile * profile
    variable = monomial(1)
    density_moments = (density, variable**2 * density, variable**4 * density)
    cutoff = a(q(DELTA / RADIUS))
    decoded = exact.decode_amplitudes(binding, q(DELTA / RADIUS))

    layers: list[arb] = []
    for layer, (ordinary, dispersion) in enumerate(decoded):
        value = layer_integral(
            layer,
            density_moments,
            cutoff,
            ball_amplitude_powers(ordinary, dispersion),
        )
        layers.append(value)
        print(json.dumps({"layer": layer, **record(value)}), flush=True)

    normalized_i = sum(layers, a(0))
    physical_i = a(q(RADIUS)) ** K * normalized_i
    payload = {
        "name": "Rigorous Arb I_46 certificate for the address-state candidate",
        "candidateSha256": candidate["sha256"],
        "supportSha256": support_digest(),
        "precisionBits": args.precision,
        "claimBoundary": (
            "This is a rigorous outward-rounded enclosure of I_46 only. "
            "The J margin and the source-level S52 proof remain separate."
        ),
        "layers": [
            {"bigCount": index, **record(value)}
            for index, value in enumerate(layers)
        ],
        "normalizedI": record(normalized_i),
        "physicalI": record(physical_i),
        "checks": {
            "allLayerLowerBoundsPositive": all(
                value.lower() > 0 for value in layers
            ),
            "normalizedILowerBoundPositive": normalized_i.lower() > 0,
            "coefficientCountIs208": len(binding["numerators"]) == 208,
            "layerCountIs11": len(layers) == MAX_BIG_COUNT + 1,
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), **payload["checks"]}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"certificate failed: {error}", file=sys.stderr)
        raise
