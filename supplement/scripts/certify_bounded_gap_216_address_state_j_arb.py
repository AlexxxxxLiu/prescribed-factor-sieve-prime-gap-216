#!/usr/bin/env python3
"""Rigorous Arb certificate for the conservative k=46 address-state J form."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

from flint import arb, arb_poly, ctx, fmpq, fmpq_poly

import build_bounded_gap_216_address_state_exact_i as exact
import certify_bounded_gap_216_address_state_i_arb as ai
from build_bounded_gap_226_layered_rational_candidate import (
    B_BY_COUNT,
    DELTA,
    MAX_BIG_COUNT,
    OUTER_RADIUS,
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
I_CERTIFICATE = ROOT / "outputs" / "bounded_gap_216_address_state_i_arb.json"
OUTPUT = ROOT / "outputs" / "bounded_gap_216_address_state_j_arb.json"
K = 46

Trivar = dict[tuple[int, int, int], fmpq]
RatBivar = dict[tuple[int, int], fmpq]
BallBivar = dict[tuple[int, int], arb]
MOMENT_CACHE: dict[tuple[tuple[str, ...], str, int], list[arb]] = {}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def trivar_add(
    target: Trivar,
    s_power: int,
    b_power: int,
    z_power: int,
    value: fmpq,
) -> None:
    if value == 0:
        return
    key = (s_power, b_power, z_power)
    target[key] = target.get(key, fmpq(0)) + value
    if target[key] == 0:
        del target[key]


def primitive_pair(
    ordinary: exact.Bivar,
    dispersion: fmpq_poly,
    profile: fmpq_poly,
    shift_big_address: bool,
) -> tuple[Trivar, Trivar]:
    """Return primitives for the constant and outer-D channel parts."""
    constant: Trivar = {}
    dlinear: Trivar = {}
    for (u_degree, b_degree), amplitude_value in ordinary.items():
        for u_mark in range(u_degree + 1):
            u_factor = math.comb(u_degree, u_mark)
            b_marks = range(b_degree + 1) if shift_big_address else (0,)
            for b_mark in b_marks:
                b_factor = (
                    math.comb(b_degree, b_mark)
                    if shift_big_address
                    else 1
                )
                residual_b_degree = (
                    b_degree - b_mark
                    if shift_big_address
                    else b_degree
                )
                for profile_degree, profile_value in enumerate(profile):
                    t_degree = u_mark + b_mark + profile_degree
                    trivar_add(
                        constant,
                        u_degree - u_mark,
                        residual_b_degree,
                        t_degree + 1,
                        amplitude_value
                        * u_factor
                        * b_factor
                        * profile_value
                        / (t_degree + 1),
                    )

    for u_degree, amplitude_value in enumerate(dispersion):
        if amplitude_value == 0:
            continue
        for u_mark in range(u_degree + 1):
            u_factor = math.comb(u_degree, u_mark)
            for profile_degree, profile_value in enumerate(profile):
                base_degree = u_mark + profile_degree
                coefficient = amplitude_value * u_factor * profile_value
                trivar_add(
                    dlinear,
                    u_degree - u_mark,
                    0,
                    base_degree + 1,
                    coefficient / (base_degree + 1),
                )
                trivar_add(
                    constant,
                    u_degree - u_mark,
                    0,
                    base_degree + 3,
                    coefficient / (base_degree + 3),
                )
    return constant, dlinear


def bivar_add_term(
    target: RatBivar, s_power: int, b_power: int, value: fmpq
) -> None:
    if value == 0:
        return
    key = (s_power, b_power)
    target[key] = target.get(key, fmpq(0)) + value
    if target[key] == 0:
        del target[key]


def evaluate_primitive(
    primitive: Trivar, endpoint: str, value: fmpq | None = None
) -> RatBivar:
    result: RatBivar = {}
    for (s_power, b_power, z_power), coefficient in primitive.items():
        if endpoint == "constant":
            assert value is not None
            bivar_add_term(
                result, s_power, b_power, coefficient * value**z_power
            )
        elif endpoint == "radial":
            for mark in range(z_power + 1):
                bivar_add_term(
                    result,
                    s_power + mark,
                    b_power,
                    coefficient * math.comb(z_power, mark) * (-1) ** mark,
                )
        elif endpoint == "address":
            assert value is not None
            for mark in range(z_power + 1):
                bivar_add_term(
                    result,
                    s_power,
                    b_power + mark,
                    coefficient
                    * math.comb(z_power, mark)
                    * value ** (z_power - mark)
                    * (-1) ** mark,
                )
        else:
            raise ValueError(f"unknown endpoint {endpoint}")
    return result


def bivar_add(
    left: RatBivar, right: RatBivar, scale: fmpq = fmpq(1)
) -> RatBivar:
    result = dict(left)
    for key, value in right.items():
        result[key] = result.get(key, fmpq(0)) + scale * value
        if result[key] == 0:
            del result[key]
    return result


def bivar_multiply(left: RatBivar, right: RatBivar) -> RatBivar:
    result: defaultdict[tuple[int, int], fmpq] = defaultdict(fmpq)
    for (left_s, left_b), left_value in left.items():
        for (right_s, right_b), right_value in right.items():
            result[(left_s + right_s, left_b + right_b)] += (
                left_value * right_value
            )
    return {key: value for key, value in result.items() if value != 0}


def channel_product(
    left: tuple[RatBivar, RatBivar],
    right: tuple[RatBivar, RatBivar],
    scale: int = 1,
) -> list[RatBivar]:
    zeroth = bivar_multiply(left[0], right[0])
    first = bivar_add(
        bivar_multiply(left[0], right[1]),
        bivar_multiply(left[1], right[0]),
    )
    second = bivar_multiply(left[1], right[1])
    if scale != 1:
        for poly in (zeroth, first, second):
            for key in list(poly):
                poly[key] *= scale
    return [zeroth, first, second]


def endpoint_channels(
    decoded: list[tuple[exact.Bivar, fmpq_poly]],
    profile: fmpq_poly,
    cutoff: fmpq,
) -> tuple[
    list[tuple[RatBivar, RatBivar]],
    list[tuple[RatBivar, RatBivar]],
    list[tuple[RatBivar, RatBivar] | None],
    list[tuple[RatBivar, RatBivar] | None],
]:
    small_constant = []
    small_radial = []
    big_address: list[tuple[RatBivar, RatBivar] | None] = [None] * len(decoded)
    big_radial: list[tuple[RatBivar, RatBivar] | None] = [None] * len(decoded)

    for ordinary, dispersion in decoded:
        primitive = primitive_pair(ordinary, dispersion, profile, False)
        small_constant.append(
            tuple(
                evaluate_primitive(part, "constant", cutoff)
                for part in primitive
            )
        )
        small_radial.append(
            tuple(evaluate_primitive(part, "radial") for part in primitive)
        )

    for layer in range(1, len(decoded)):
        ordinary, dispersion = decoded[layer]
        primitive = primitive_pair(ordinary, dispersion, profile, True)
        lower = tuple(
            evaluate_primitive(part, "constant", cutoff)
            for part in primitive
        )
        boundary = q(B_BY_COUNT[layer] / RADIUS)
        big_address[layer] = tuple(
            bivar_add(
                evaluate_primitive(part, "address", boundary),
                lower_part,
                fmpq(-1),
            )
            for part, lower_part in zip(primitive, lower)
        )
        big_radial[layer] = tuple(
            bivar_add(
                evaluate_primitive(part, "radial"),
                lower_part,
                fmpq(-1),
            )
            for part, lower_part in zip(primitive, lower)
        )
    return small_constant, small_radial, big_address, big_radial


def to_ball_bivars(polys: list[RatBivar]) -> list[BallBivar]:
    return [{key: arb(value) for key, value in poly.items()} for poly in polys]


def interval_delta(
    density: arb_poly,
    amplitude: BallBivar,
    offset: arb,
    lower: arb,
    upper: arb,
) -> arb:
    if not amplitude:
        return arb(0)
    affine = arb_poly([offset, arb(1)])
    poly = arb_poly([0])
    for (s_power, b_power), coefficient in amplitude.items():
        if b_power == 0:
            poly += coefficient * affine**s_power
    antiderivative = (density * poly).integral()
    return antiderivative(upper) - antiderivative(lower)


def rectangle_bivar(
    small: arb_poly,
    big: arb_poly,
    amplitude: BallBivar,
    total_offset: arb,
    big_offset: arb,
    x_cap: arb,
    y_cap: arb,
) -> arb:
    if not amplitude:
        return arb(0)
    maximum_s = max(key[0] for key in amplitude)
    maximum_b = max(key[1] for key in amplitude)

    def moments(poly: arb_poly, cap: arb, maximum: int) -> list[arb]:
        # Object addresses can be reused after an arb_poly is destroyed.  A
        # coefficient signature makes cache reuse depend on mathematical
        # content instead of Python allocator history.
        signature = tuple(str(poly[index]) for index in range(len(poly)))
        key = (signature, str(cap), maximum)
        cached = MOMENT_CACHE.get(key)
        if cached is not None:
            return cached
        values = []
        for power in range(maximum + 1):
            antiderivative = (poly * ai.monomial(power)).integral()
            values.append(antiderivative(cap) - antiderivative(arb(0)))
        MOMENT_CACHE[key] = values
        return values

    x_moments = moments(small, x_cap, maximum_s)
    y_moments = moments(big, y_cap, maximum_s + maximum_b)
    total_powers = [total_offset**power for power in range(maximum_s + 1)]
    big_powers = [big_offset**power for power in range(maximum_b + 1)]

    # E[y^j(B_0+y)^b] on the big-coordinate rectangle.
    shifted_y_moments: list[list[arb]] = []
    for b_power in range(maximum_b + 1):
        shifted_y_moments.append(
            [
                sum(
                    (
                        math.comb(b_power, mark)
                        * big_powers[b_power - mark]
                        * y_moments[y_power + mark]
                    )
                    for mark in range(b_power + 1)
                )
                for y_power in range(maximum_s + 1)
            ]
        )

    result = arb(0)
    for (s_power, b_power), coefficient in amplitude.items():
        for x_power in range(s_power + 1):
            for y_power in range(s_power - x_power + 1):
                constant_power = s_power - x_power - y_power
                multinomial = (
                    math.factorial(s_power)
                    // (
                        math.factorial(x_power)
                        * math.factorial(y_power)
                        * math.factorial(constant_power)
                    )
                )
                result += (
                    coefficient
                    * multinomial
                    * total_powers[constant_power]
                    * x_moments[x_power]
                    * shifted_y_moments[b_power][y_power]
                )
    return result


def moment_sum(
    power_amplitudes: list[BallBivar],
    small_moments: list[dict[int, arb_poly]],
    big_moments: list[arb_poly] | None,
    shift_index: int,
    region,
) -> arb:
    total = arb(0)
    for d_power, amplitude in enumerate(power_amplitudes):
        if not amplitude:
            continue
        if big_moments is None:
            small = small_moments[d_power].get(shift_index)
            if small is not None:
                total += region(small, None, amplitude)
            continue
        for small_power in range(d_power + 1):
            small = small_moments[small_power].get(shift_index)
            if small is None:
                continue
            total += math.comb(d_power, small_power) * region(
                small, big_moments[d_power - small_power], amplitude
            )
    return total


def layer_channels(
    layer: int,
    density_moments: tuple[arb_poly, arb_poly, arb_poly],
    cutoff: arb,
    channel_data,
) -> tuple[arb, arb, arb]:
    small_constant, small_radial, big_address, big_radial = channel_data
    cutoff_fraction = DELTA / RADIUS
    outer_fraction = OUTER_RADIUS / RADIUS
    small_count = K - 1 - layer
    maximum_shift = int(
        (outer_fraction - layer * cutoff_fraction) / cutoff_fraction
    )
    small_moments = ai.compact_group_moments(
        density_moments, cutoff, small_count, maximum_shift
    )
    big_moments = (
        None
        if layer == 0
        else ai.big_group_moments(density_moments, cutoff, layer)
    )
    big_offset_fraction = layer * cutoff_fraction
    big_offset = arb(q(big_offset_fraction))
    current_big_cap_fraction = (
        None
        if layer == 0
        else B_BY_COUNT[layer] / RADIUS - big_offset_fraction
    )
    current_big_cap = (
        None
        if current_big_cap_fraction is None
        else arb(q(current_big_cap_fraction))
    )

    aa_low = to_ball_bivars(
        channel_product(small_constant[layer], small_constant[layer])
    )
    aa_high = to_ball_bivars(
        channel_product(small_radial[layer], small_radial[layer])
    )
    target = layer + 1
    has_big_channel = (
        target <= MAX_BIG_COUNT and big_address[target] is not None
    )
    if has_big_channel:
        assert big_address[target] is not None
        assert big_radial[target] is not None
        cap_cross = to_ball_bivars(
            channel_product(
                small_constant[layer], big_address[target], scale=2
            )
        )
        cap_square = to_ball_bivars(
            channel_product(big_address[target], big_address[target])
        )
        radial_cross = to_ball_bivars(
            channel_product(
                small_constant[layer], big_radial[target], scale=2
            )
        )
        radial_square = to_ball_bivars(
            channel_product(big_radial[target], big_radial[target])
        )
    aa = arb(0)
    ab = arb(0)
    bb = arb(0)
    combination = math.comb(K - 1, layer)

    for shift_index in sorted(small_moments[0]):
        total_offset_fraction = big_offset_fraction + shift_index * cutoff_fraction
        total_offset = arb(q(total_offset_fraction))
        low_cap_fraction = 1 - cutoff_fraction - total_offset_fraction
        outer_cap_fraction = outer_fraction - total_offset_fraction
        if outer_cap_fraction <= 0:
            continue
        low_cap = arb(q(low_cap_fraction))
        outer_cap = arb(q(outer_cap_fraction))

        def triangle_region(
            cap: arb, cap_fraction, big_cap: arb | None
        ):
            def integrate(small, big, amplitude):
                if big is None:
                    return interval_delta(
                        small, amplitude, total_offset, arb(0), cap
                    )
                assert big_cap is not None
                return ai.integrate_triangle_bivar(
                    small,
                    big,
                    amplitude,
                    total_offset,
                    big_offset,
                    cap,
                    big_cap,
                    current_big_cap_fraction < cap_fraction,
                )
            return integrate

        if low_cap_fraction > 0:
            aa += moment_sum(
                aa_low,
                small_moments,
                big_moments,
                shift_index,
                triangle_region(low_cap, low_cap_fraction, current_big_cap),
            )
            aa -= moment_sum(
                aa_high,
                small_moments,
                big_moments,
                shift_index,
                triangle_region(low_cap, low_cap_fraction, current_big_cap),
            )
        aa += moment_sum(
            aa_high,
            small_moments,
            big_moments,
            shift_index,
            triangle_region(outer_cap, outer_cap_fraction, current_big_cap),
        )

        if not has_big_channel:
            continue
        next_cap_fraction = B_BY_COUNT[target] / RADIUS - target * cutoff_fraction
        if next_cap_fraction <= 0 or low_cap_fraction <= 0:
            continue
        x_boundary_fraction = low_cap_fraction - next_cap_fraction
        x_boundary = arb(q(max(x_boundary_fraction, 0)))
        next_cap = arb(q(next_cap_fraction))

        if x_boundary_fraction > 0:
            def rectangle_region(small, big, amplitude):
                if big is None:
                    return interval_delta(
                        small,
                        amplitude,
                        total_offset,
                        arb(0),
                        x_boundary,
                    )
                return rectangle_bivar(
                    small,
                    big,
                    amplitude,
                    total_offset,
                    big_offset,
                    x_boundary,
                    next_cap,
                )

            ab += moment_sum(
                cap_cross,
                small_moments,
                big_moments,
                shift_index,
                rectangle_region,
            )
            bb += moment_sum(
                cap_square,
                small_moments,
                big_moments,
                shift_index,
                rectangle_region,
            )

        lower_fraction = max(x_boundary_fraction, 0)
        radial_cap_fraction = low_cap_fraction - lower_fraction
        if radial_cap_fraction > 0:
            lower = arb(q(lower_fraction))
            radial_cap = arb(q(radial_cap_fraction))

            def lower_triangle_region(small, big, amplitude):
                shifted_small = ai.shifted(small, lower)
                if big is None:
                    return interval_delta(
                        shifted_small,
                        amplitude,
                        total_offset + lower,
                        arb(0),
                        radial_cap,
                    )
                return ai.integrate_triangle_bivar(
                    shifted_small,
                    big,
                    amplitude,
                    total_offset + lower,
                    big_offset,
                    radial_cap,
                    radial_cap,
                    False,
                )

            ab += moment_sum(
                radial_cross,
                small_moments,
                big_moments,
                shift_index,
                lower_triangle_region,
            )
            bb += moment_sum(
                radial_square,
                small_moments,
                big_moments,
                shift_index,
                lower_triangle_region,
            )
    return combination * aa, combination * ab, combination * bb


def rec(value: arb) -> dict[str, object]:
    midpoint, radius, exponent = value.mid_rad_10exp()
    return {
        "ball": str(value),
        "lower": str(value.lower()),
        "upper": str(value.upper()),
        "midRad10Exp": [str(midpoint), str(radius), int(exponent)],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--precision", type=int, default=1024)
    args = parser.parse_args()
    ctx.prec = args.precision

    candidate = json.loads(CANDIDATE.read_text(encoding="utf-8"))
    binding = candidate["binding"]
    if binding["supportDigest"] != support_digest():
        raise ArithmeticError("support hash mismatch")
    expected_hash = hashlib.sha256(
        json.dumps(binding, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if expected_hash != candidate["sha256"]:
        raise ArithmeticError("candidate binding hash mismatch")

    i_payload = json.loads(I_CERTIFICATE.read_text(encoding="utf-8"))
    if i_payload["candidateSha256"] != candidate["sha256"]:
        raise ArithmeticError("I and J candidate hashes differ")
    if i_payload["supportSha256"] != support_digest():
        raise ArithmeticError("I certificate support hash mismatch")
    if int(i_payload["precisionBits"]) < args.precision:
        raise ArithmeticError("I certificate precision is below requested J precision")
    normalized_i = arb(i_payload["normalizedI"]["ball"])

    profile_exact = to_flint_poly(
        combine_basis(PROFILE_NUMERATORS, shifted_chebyshev_basis(32))
    )
    decoded = exact.decode_amplitudes(binding, q(DELTA / RADIUS))
    channel_data = endpoint_channels(
        decoded, profile_exact, q(DELTA / RADIUS)
    )
    profile = ai.to_ball_poly(profile_exact)
    density = profile * profile
    variable = ai.monomial(1)
    density_moments = (density, variable**2 * density, variable**4 * density)
    cutoff = arb(q(DELTA / RADIUS))

    layer_values = []
    totals = [arb(0), arb(0), arb(0)]
    for layer in range(MAX_BIG_COUNT + 1):
        channels = layer_channels(
            layer, density_moments, cutoff, channel_data
        )
        layer_values.append(channels)
        for index, value in enumerate(channels):
            totals[index] += value
        print(
            json.dumps(
                {
                    "layer": layer,
                    "aa": str(channels[0]),
                    "twiceAb": str(channels[1]),
                    "bb": str(channels[2]),
                }
            ),
            flush=True,
        )

    normalized_j = sum(totals, arb(0))
    ratio = K * arb(q(RADIUS)) * normalized_j / normalized_i
    margin = K * arb(q(RADIUS)) * normalized_j - normalized_i
    payload = {
        "name": "Rigorous Arb J_46 and variational-margin certificate",
        "candidateSha256": candidate["sha256"],
        "supportSha256": support_digest(),
        "iCertificateFileSha256": sha256_file(I_CERTIFICATE),
        "precisionBits": args.precision,
        "claimBoundary": (
            "This certifies the finite conservative variational inequality. "
            "The manuscript proves the analytic S52 theorem separately and "
            "combines the two components to claim unconditional H_1<=216. "
            "This numerical certificate does not replace that analytic proof."
        ),
        "geometry": {
            "radius": str(RADIUS),
            "outerRadius": str(OUTER_RADIUS),
            "delta": str(DELTA),
        },
        "channelTotals": {
            "aa": rec(totals[0]),
            "twiceAb": rec(totals[1]),
            "bb": rec(totals[2]),
        },
        "normalizedJ": rec(normalized_j),
        "normalizedI": rec(normalized_i),
        "score46JOverI": rec(ratio),
        "margin46RJMinusI": rec(margin),
        "layers": [
            {
                "bigCount": layer,
                "aa": rec(values[0]),
                "twiceAb": rec(values[1]),
                "bb": rec(values[2]),
            }
            for layer, values in enumerate(layer_values)
        ],
        "checks": {
            "normalizedJLowerBoundPositive": normalized_j.lower() > 0,
            "scoreLowerBoundExceedsOne": ratio.lower() > 1,
            "marginLowerBoundPositive": margin.lower() > 0,
            "candidateHashesMatch": True,
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(OUTPUT), **payload["checks"]}, indent=2))


if __name__ == "__main__":
    main()
