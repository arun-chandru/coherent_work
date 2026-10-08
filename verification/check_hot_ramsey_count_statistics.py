"""Certified finite-grid shot budgets and optional actual-count intervals.

All claimed inequalities use FLINT Arb. No observed counts are generated.
The deterministic continuum certificate is an input, not reproved here.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path

from flint import arb, ctx  # noqa: E402

ctx.prec = 192


def rational(value):
    f = Fraction(str(value))
    return arb(f.numerator) / f.denominator


def upper_integer(value):
    # Float locates a candidate only; Arb certifies the final strict bound.
    candidate = math.ceil(float(value.upper())) + 1
    if not arb(candidate) > value:
        raise ArithmeticError("Integer planning bound not certified")
    return candidate


def cdf(n, k, p):
    """P(Bin(n,p)<=k), positive recurrence suited to rare successes."""
    if k < 0:
        return arb(0)
    if k >= n:
        return arb(1)
    if p == 0:
        return arb(1)
    if p == 1:
        return arb(0)
    term = (1 - p) ** n
    total = term
    ratio = p / (1 - p)
    for j in range(1, k + 1):
        term *= rational(Fraction(n - j + 1, j)) * ratio
        total += term
    return total


def cdf_root_bracket(n, k, target, iterations=96):
    """Enclose root of decreasing binomial CDF by exact dyadic brackets."""
    lo, hi = Fraction(0), Fraction(1)
    for _ in range(iterations):
        mid = (lo + hi) / 2
        gap = cdf(n, k, rational(mid)) - target
        if gap > 0:
            lo = mid
        elif gap < 0:
            hi = mid
        else:
            # Preserve the certified bracket rather than guess a sign.
            break
    return rational(lo), rational(hi)


def count_interval(n, k, alpha):
    """Outer enclosure of the equal-tailed exact binomial interval."""
    if not (n >= 1 and 0 <= k <= n):
        raise ValueError("Require n>=1 and an integer 0<=k<=n")
    if k > n // 2:
        lo, hi = count_interval(n, n - k, alpha)
        return 1 - hi, 1 - lo
    if k > 100000:
        raise ValueError("This positive-sum checker is intended for rare counts")
    lo = arb(0) if k == 0 else cdf_root_bracket(n, k - 1, 1 - alpha / 2)[0]
    hi = cdf_root_bracket(n, k, alpha / 2)[1]
    return lo, hi


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--certificate", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--alpha", default="1/20")
    parser.add_argument("--epsilon", default="1/20")
    parser.add_argument("--count", type=int)
    parser.add_argument("--trials", type=int)
    args = parser.parse_args()
    if (args.count is None) != (args.trials is None):
        parser.error("Provide both --count and --trials, or neither")
    alpha, epsilon = rational(args.alpha), rational(args.epsilon)
    if not (alpha > 0 and alpha < 1 and epsilon > 0):
        parser.error("Require 0<alpha<1 and epsilon>0")
    raw = args.certificate.read_bytes()
    cert = json.loads(raw)
    params, budgets = cert["parameters"], cert["budgets"]
    lam = rational(Fraction(params["s"]) / Fraction(params["r"]))
    target = arb(budgets["target_cap_source"])
    source_mean = arb(budgets["finite_device_weighted_mean"])
    direct_bias = arb(budgets["total_absolute_systematic_bias_upper"]).upper()
    conservative_bias = (
        arb(budgets["conservative_absolute_systematic_bias_upper"]).upper()
        if "conservative_absolute_systematic_bias_upper" in budgets else
        (arb(budgets["source_quadrature_minus_target"]).abs_upper()
         + arb(budgets["weighted_absolute_sample_bias_upper"]).upper()).upper()
    )

    # Actual design: exact positive rational weights on the predeclared grid.
    denominator = 2 ** 80
    numerators, implemented, upper_probs, times, time_fractions = [], [], [], [], []
    weight_error = arb(0)
    for item in cert["samples"]:
        exact_weight = arb(item["weight"])
        integer = (exact_weight * denominator).floor().unique_fmpz()
        if integer is None or integer <= 0:
            raise ArithmeticError("Weight quantization not uniquely certified")
        numerator = int(integer)
        wi = rational(Fraction(numerator, denominator))
        ui = arb(item["scaled_probability"]).upper()
        if not ui >= 0:
            raise ArithmeticError("Uncertified probability upper bound")
        weight_error += (wi - exact_weight).abs_upper() * ui
        numerators.append(numerator)
        implemented.append(wi)
        upper_probs.append(ui)
        time_fraction = Fraction(item["time_exact"])
        time_fractions.append(time_fraction)
        times.append(rational(time_fraction))
    total_numerator = sum(numerators)
    weight_sum = rational(Fraction(total_numerator, denominator))
    weight_error = weight_error.upper()
    scaled_mean_upper = (source_mean.upper() + weight_error).upper()
    direct_bias = (direct_bias + weight_error).upper()
    conservative_bias = (conservative_bias + weight_error).upper()
    ell = (2 / alpha).log()
    required = 2 * lam * weight_sum * ell * (scaled_mean_upper + epsilon / 3) / epsilon ** 2
    preparations = upper_integer(required)
    exponent = preparations * epsilon ** 2 / (
        2 * lam * weight_sum * (scaled_mean_upper + epsilon / 3))
    failure_bound = 2 * (-exponent).exp()
    if not failure_bound < alpha:
        raise ArithmeticError("Statistical failure probability not certified")
    mean_time = sum((w * t for w, t in zip(implemented, times)), arb(0)) / weight_sum
    max_time = rational(max(time_fractions))
    r = rational(params["r"])
    direct_total = (direct_bias + epsilon).upper()
    conservative_total = (conservative_bias + epsilon).upper()

    # Comparison with the same aggregate Bernstein bound at equal shots.
    # This is a bound comparison, not an oracle lower bound or measured gain.
    m = len(implemented)
    equal_A = m * sum((w * w * u for w, u in zip(implemented, upper_probs)), arb(0))
    equal_B = rational(Fraction(m * max(numerators), denominator))
    equal_n0 = upper_integer(2 * lam * ell * (equal_A + epsilon * equal_B / 3) / epsilon ** 2)
    equal_total = m * ((equal_n0 + m - 1) // m)
    result = {
        "status": "Certified shot budget; no experimental counts supplied" if args.count is None
                  else "Certified shot budget and interval for user-supplied integer count",
        "arithmetic": "FLINT Arb 192 bits; exact rational implemented weights",
        "certificate_path": str(args.certificate),
        "certificate_sha256": hashlib.sha256(raw).hexdigest(),
        "statistics_checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "parameters": params,
        "lambda": str(lam),
        "alpha_exact": args.alpha,
        "epsilon_stat_exact": args.epsilon,
        "source_target_ball": str(target),
        "implemented_weight_sum": str(weight_sum),
        "weight_rounding_bias_upper": str(weight_error),
        "scaled_actual_weighted_mean_upper_D": str(scaled_mean_upper),
        "required_preparation_expression_ball": str(required),
        "preparations": preparations,
        "fresh_pure_ancilla_uses": preparations,
        "ramsey_sequences": preparations,
        "ancilla_readouts": preparations,
        "time_label_draws": preparations,
        "failure_probability_upper": str(failure_bound.upper()),
        "point_specific_direct_systematic_upper": str(direct_bias),
        "conservative_systematic_upper_from_serialized_enclosures": str(conservative_bias),
        "point_specific_direct_total_error_upper": str(direct_total),
        "conservative_total_error_upper": str(conservative_total),
        "direct_total_below_quarter_target": bool(direct_total < target / 4),
        "conservative_total_below_quarter_target": bool(conservative_total < target / 4),
        "pooled_acceptance_probability_enclosure": str(
            (source_mean + arb(0, weight_error)) / (lam * weight_sum)),
        "expected_success_count_enclosure": str(
            preparations * (source_mean + arb(0, weight_error)) / (lam * weight_sum)),
        "mean_dimensionless_time_per_trial": str(mean_time),
        "max_dimensionless_time_per_trial": str(max_time),
        "max_total_control_duration_in_hbar_over_Ephys0": str(preparations * r * r * max_time),
        "expected_total_control_duration_in_hbar_over_Ephys0": str(preparations * r * r * mean_time),
        "equal_shot_aggregate_bernstein_preparation_bound": equal_total,
        "sampling_design": {
            "weight_denominator": denominator,
            "probability_denominator": total_numerator,
            "nodes": [{"time_exact": item["time_exact"], "probability_numerator": n}
                      for item, n in zip(cert["samples"], numerators)],
        },
        "scope": "IID ideal sharp-wall full-Gibbs Ramsey trials; deterministic certificate is a prerequisite. Physical drift/noise, preparation and readout overheads are not newly bounded here.",
    }
    if args.count is not None:
        lo, hi = count_interval(args.trials, args.count, alpha)
        scale = lam * weight_sum
        result["user_supplied_observation"] = {
            "trials": args.trials, "successes": args.count,
            "pooled_probability_confidence_lower": str(lo),
            "pooled_probability_confidence_upper": str(hi),
            "raw_point_estimate": str(scale * args.count / args.trials),
            "source_interval_direct_lower": str((scale * lo - direct_bias).lower()),
            "source_interval_direct_upper": str((scale * hi + direct_bias).upper()),
            "source_interval_conservative_lower": str((scale * lo - conservative_bias).lower()),
            "source_interval_conservative_upper": str((scale * hi + conservative_bias).upper()),
        }
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in result.items() if k not in ("sampling_design", "parameters")}, indent=2))


if __name__ == "__main__":
    main()
