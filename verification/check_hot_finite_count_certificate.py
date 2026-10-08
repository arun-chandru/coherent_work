"""Actual sharp-wall Ramsey probabilities and finite positive reconstruction.

All mathematical claims use FLINT Arb/Acb. Ordinary floats only locate roots,
time execution, and print noncertifying resource/central-value diagnostics.
The complete Gibbs and final spectral tails are paid; there is no mesh.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parent
PREVIOUS = ROOT
from check_hot_ramsey_finite_certificate import (  # noqa: E402
    acb, arb, ball, complete_prefix, ctx, free_gibbs_tail,
    serialize_mode, signed_full_overlap, strict_positive,
)

ctx.prec = 192


def interval(lo, hi):
    return lo.lower().union(hi.upper())


def integrate(fn, lo, hi):
    value = acb.integral(fn, lo, hi, abs_tol=ball("1e-11"),
                         rel_tol=ball("1e-11"), eval_limit=200000,
                         depth_limit=40)
    if not value.is_finite():
        raise ArithmeticError("Nonfinite rigorous quadrature")
    if not value.imag.contains(0):
        raise ArithmeticError("Real integral has nonzero imaginary enclosure")
    if not value.real.rad() < ball("1e-8"):
        raise ArithmeticError(f"Quadrature accuracy not certified: {value}")
    return value.real


def high_density_q(q, analytic):
    k = (q*q + 4).sqrt(analytic=analytic)
    return 18*q*(1+k*q)/(arb.pi()*(q*q+4)**2*(k+q)*(q*q+1)**2)


def tail_density_u(u, analytic):
    h = (1+4*u*u).sqrt(analytic=analytic)
    return 18*u**4*(u*u+h)/(arb.pi()*h**4*(h+1)*(1+u*u)**2)


def source_tail(qmax):
    """Exact source mass above E=4+qmax^2, by q=1/u."""
    return integrate(tail_density_u, 0, 1/qmax)


def source_value(t, qmax, tail):
    def low(theta, analytic):
        sn, cs = theta.sin(), theta.cos()
        z = 2*t*sn*sn
        return 9*t*t/arb.pi()*sn*sn*cs*z.sinc()**2/(1+2*cs)**2
    def high(q, analytic):
        return high_density_q(q, analytic)*(t*(q*q+4)/2).sin()**2
    retained = integrate(low, 0, arb.pi()/2) + integrate(high, 0, qmax)
    return interval(retained.lower(), retained.upper()+tail.upper())


def capped_source(qmax, tail):
    """Cap is exactly one; E=1 occurs at theta=pi/6."""
    def below(theta, analytic):
        cs = theta.cos()
        return 9*cs/(arb.pi()*(1+2*cs)**2)
    def above(theta, analytic):
        sn, cs = theta.sin(), theta.cos()
        return 9*cs/(4*arb.pi()*sn*sn*(1+2*cs)**2)
    return (integrate(below, 0, arb.pi()/6)
            + integrate(above, arb.pi()/6, arb.pi()/2)
            + integrate(high_density_q, 0, qmax) + tail)


def primitive(t):
    if t == 0:
        return ball(0)
    return 4/arb.pi()*(t.si()-(1-t.cos())/t)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=256)
    parser.add_argument("--horizon", default="32")
    parser.add_argument("--s", default="12")
    parser.add_argument("--inverse-r", type=int, default=50)
    parser.add_argument("--initial-count", type=int, default=40)
    parser.add_argument("--final-scaled-cut", default="100")
    parser.add_argument("--source-qmax", default="4")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    start = time.perf_counter()
    s, horizon, qmax = map(ball, (args.s, args.horizon, args.source_qmax))
    if args.inverse_r < 1:
        raise ValueError("Positive integer inverse-r required")
    initial_barrier = args.inverse_r**2
    final_barrier = 4*initial_barrier
    gamma = args.inverse_r*s
    theta = s/initial_barrier
    count = args.samples
    if count < 1:
        raise ValueError("Positive sample count required")
    strict_positive(horizon, "positive horizon")
    strict_positive(qmax, "positive source cutoff")
    initial_all = complete_prefix(initial_barrier,
        args.inverse_r+(args.initial_count+3)*math.pi/2)
    if len(initial_all) <= args.initial_count:
        raise ArithmeticError("Missing initial omitted neighbor")
    initial = initial_all[:args.initial_count]
    final_cut = initial_barrier*ball(args.final_scaled_cut)
    final_all = complete_prefix(final_barrier, math.sqrt(float(final_cut))+4)
    final = []
    for mode in final_all:
        if mode["energy"] < final_cut:
            final.append(mode)
        elif mode["energy"] > final_cut:
            final_neighbor = mode
            break
        else:
            raise ArithmeticError("Unresolved final cutoff")
    if final[-1]["n"]+1 != final_neighbor["n"]:
        raise ArithmeticError("Missing final omitted neighbor")
    square_overlap = [[signed_full_overlap(i, f)**2 for f in final] for i in initial]
    overlap_tail, independent_spectral_tail = [], []
    for mode, row in zip(initial, square_overlap):
        residual = (1-sum(row, ball(0))).upper()
        strict_positive(residual, "nonnegative final projection residual")
        gap = final_neighbor["energy"].lower()-mode["energy"]
        strict_positive(gap, "omitted spectrum above every retained input")
        bound = ((final_barrier-initial_barrier)**2*mode["source"]/gap**2).upper()
        overlap_tail.append(min(residual, bound))
        independent_spectral_tail.append(bound)
    thermal_weights = [(-theta*i["energy"]).exp() for i in initial]
    partial_z = sum(thermal_weights, ball(0))
    omitted_gibbs = free_gibbs_tail(theta, len(initial))
    zlo, zhi = partial_z.lower(), (partial_z+omitted_gibbs).upper()
    omitted_final_trace = sum((w*t for w,t in zip(thermal_weights, overlap_tail)), ball(0))
    omitted_probability = (omitted_gibbs/zlo).upper()
    final_probability_tail = (omitted_final_trace/zlo).upper()
    activated_weight = sum((w for w,i in zip(thermal_weights,initial)
                            if i["energy"] > initial_barrier), ball(0))
    for i in initial:
        if not (i["energy"] > initial_barrier or i["energy"] < initial_barrier):
            raise ArithmeticError("Initial barrier threshold is unresolved")
    activated_probability = ((activated_weight+omitted_gibbs)/zlo).upper()
    energies = [[(f["energy"]-i["energy"])/initial_barrier for f in final] for i in initial]
    weighted_overlaps = [[w*x for x in row] for w,row in zip(thermal_weights,square_overlap)]
    source_tail_ball = source_tail(qmax)
    phi = capped_source(qmax, source_tail_ball)
    strict_positive(phi, "positive capped source")
    samples = []
    measured_sum = ball(0)
    source_sum = ball(0)
    aggregate_difference = ball(0)
    weighted_abs_bias = ball(0)
    max_sample_bias = ball(0)
    capped_finite_sum = ball(0)
    for row_e,row_w in zip(energies,weighted_overlaps):
        for gap,weight in zip(row_e,row_w):
            absolute = abs(gap)
            clipped = interval(min(ball(1), absolute.lower()),
                               min(ball(1), absolute.upper()))
            capped_finite_sum += weight*clipped
    capped_finite = gamma*interval((capped_finite_sum/zhi).lower(),
        ((capped_finite_sum+omitted_final_trace+omitted_gibbs)/zlo).upper())
    print(json.dumps({"progress":"spectral setup complete", "initial_count":len(initial),
        "final_count":len(final), "source_cap":str(phi),
        "full_time_finite_device_capped_work":str(capped_finite),
        "source_high_energy_mass_tail":str(source_tail_ball),
        "elapsed_seconds_diagnostic":time.perf_counter()-start}), flush=True)
    for index in range(count):
        time_exact = Fraction(2*index+1, 2*count)*Fraction(args.horizon)
        t = ball(time_exact)
        finite_trace = sum((weight*(t*gap/2).sin()**2
            for row_e,row_w in zip(energies,weighted_overlaps)
            for gap,weight in zip(row_e,row_w)), ball(0))
        probability = interval((finite_trace/zhi).lower(),
            ((finite_trace+omitted_final_trace+omitted_gibbs)/zlo).upper())
        source = source_value(t,qmax,source_tail_ball)
        wt = primitive(horizon*(index+1)/count)-primitive(horizon*index/count)
        strict_positive(wt,"positive exact quadrature weight")
        scaled_probability = gamma*probability
        bias = scaled_probability-source
        measured_sum += wt*scaled_probability
        source_sum += wt*source
        aggregate_difference += wt*bias
        weighted_abs_bias += wt*abs(bias)
        max_sample_bias = max(max_sample_bias, abs(bias).upper())
        samples.append({"time_exact":str(time_exact), "weight":str(wt),
            "probability":str(probability), "scaled_probability":str(scaled_probability),
            "source":str(source), "signed_finite_device_bias":str(bias)})
        if (index+1)%32 == 0 or index+1 == count:
            print(json.dumps({"progress":"samples", "completed":index+1,
                "total":count,"elapsed_seconds_diagnostic":time.perf_counter()-start}),flush=True)
    total_weight = primitive(horizon)
    total_bias = measured_sum-phi
    conservative_bias = abs(source_sum-phi)+weighted_abs_bias
    rational_claims = {}
    if (args.inverse_r == 50 and args.samples == 256 and Fraction(args.horizon) == 32
            and Fraction(args.s) == 12 and args.initial_count == 40
            and Fraction(args.final_scaled_cut) == 100):
        rational_claims = {
            "target_lower": "4040654186/10000000000",
            "target_upper": "4040654188/10000000000",
            "finite_weighted_mean_lower": "51836/100000",
            "finite_weighted_mean_upper": "51838/100000",
            "raw_absolute_systematic_bias_upper": "114314/1000000",
            "raw_bias_exceeds_one_quarter_of_target": True,
            "activated_weighted_contribution_upper": "8087/1000000",
        }
        strict_positive(phi-ball(rational_claims["target_lower"]), "cap lower")
        strict_positive(ball(rational_claims["target_upper"])-phi, "cap upper")
        strict_positive(measured_sum-ball(rational_claims["finite_weighted_mean_lower"]),
                        "raw mean lower")
        strict_positive(ball(rational_claims["finite_weighted_mean_upper"])-measured_sum,
                        "raw mean upper")
        strict_positive(ball(rational_claims["raw_absolute_systematic_bias_upper"])
                        -abs(total_bias), "raw bias upper")
        strict_positive(total_bias-phi/4, "raw bias alone exceeds 25 percent target")
        strict_positive(ball(rational_claims["activated_weighted_contribution_upper"])
                        -gamma*total_weight*activated_probability, "activated upper")
    if (args.inverse_r == 200 and args.samples == 256 and Fraction(args.horizon) == 32
            and Fraction(args.s) == 128 and args.initial_count == 64
            and Fraction(args.final_scaled_cut) == 100):
        rational_claims = {
            "target_lower": "4040654186/10000000000",
            "target_upper": "4040654188/10000000000",
            "finite_weighted_mean_lower": "40358/100000",
            "finite_weighted_mean_upper": "40360/100000",
            "raw_absolute_systematic_bias_upper": "48/100000",
            "conservative_absolute_systematic_bias_upper": "95/1000",
            "conservative_bias_plus_005_below_one_quarter_of_target": True,
        }
        strict_positive(phi-ball(rational_claims["target_lower"]), "follow-up cap lower")
        strict_positive(ball(rational_claims["target_upper"])-phi, "follow-up cap upper")
        strict_positive(measured_sum-ball(rational_claims["finite_weighted_mean_lower"]),
                        "follow-up raw mean lower")
        strict_positive(ball(rational_claims["finite_weighted_mean_upper"])-measured_sum,
                        "follow-up raw mean upper")
        strict_positive(ball(rational_claims["raw_absolute_systematic_bias_upper"])
                        -abs(total_bias), "follow-up raw bias upper")
        strict_positive(ball(rational_claims["conservative_absolute_systematic_bias_upper"])
                        -conservative_bias, "follow-up conservative bias upper")
        strict_positive(phi/4-conservative_bias-ball(".005"),
                        "conservative bias plus .005 is below one quarter target")
    budget = {
        "target_cap_source":str(phi), "finite_device_weighted_mean":str(measured_sum),
        "source_weighted_quadrature":str(source_sum),
        "source_quadrature_minus_target":str(source_sum-phi),
        "source_quadrature_absolute_bias_upper":str(abs(source_sum-phi).upper()),
        "aggregate_device_minus_source":str(aggregate_difference),
        "weighted_absolute_sample_bias_upper":str(weighted_abs_bias.upper()),
        "maximum_absolute_sample_bias_upper":str(max_sample_bias),
        "total_signed_systematic_bias":str(total_bias),
        "total_absolute_systematic_bias_upper":str(abs(total_bias).upper()),
        "conservative_absolute_systematic_bias_upper":str(
            conservative_bias.upper()),
        "sum_weights":str(total_weight),
        "sum_computed_weights":str(sum((arb(x["weight"]) for x in samples),ball(0))),
        "raw_scaled_omitted_gibbs_probability":str(gamma*omitted_probability),
        "raw_scaled_final_spectral_probability_tail":str(gamma*final_probability_tail),
        "weighted_scaled_gibbs_tail_upper":str((gamma*total_weight*omitted_probability).upper()),
        "weighted_scaled_final_tail_upper":str((gamma*total_weight*final_probability_tail).upper()),
        "weighted_source_high_energy_tail_upper":str((total_weight*source_tail_ball).upper()),
        "activated_input_probability_upper":str(activated_probability),
        "weighted_activated_contribution_upper":str((gamma*total_weight*activated_probability).upper()),
        "full_time_finite_device_capped_work":str(capped_finite),
        "full_time_capped_work_minus_source":str(capped_finite-phi),
    }
    source_paths = [Path(__file__), PREVIOUS/"check_hot_ramsey_finite_certificate.py",
        ROOT/"check_phy_info_finite_purity_certificate.py",
        ROOT/"check_phy_info_continuum_certificate.py"]
    result = {"status":"validated exact-continuum finite-sample Ramsey reconstruction, all omitted modes paid",
        "arithmetic":"FLINT Arb/Acb 192 bits; certified rational Sturm root brackets; rigorous analytic quadrature",
        "parameters":{"r":str(Fraction(1,args.inverse_r)),"s":args.s,"cap":"1","horizon":args.horizon,
            "samples":count,"mesh":str(Fraction(args.horizon)/count),
            "source_qmax":args.source_qmax,"source_energy_cut":str(4+Fraction(args.source_qmax)**2),
            "initial_count":len(initial),"final_count":len(final),
            "initial_barrier":initial_barrier,"final_barrier":final_barrier},
        "budgets":budget,"samples":samples,"partition":str(interval(zlo,zhi)),
        "certified_rational_claims":rational_claims,
        "omitted_initial_gibbs_unnormalized":str(omitted_gibbs),
        "omitted_initial_probability_upper":str(omitted_probability),
        "final_overlap_tail_bounds":[str(x) for x in overlap_tail],
        "independent_spectral_tail_bounds":[str(x) for x in independent_spectral_tail],
        "source_high_energy_mass":str(source_tail_ball),
        "initial_modes":[serialize_mode(x) for x in initial],
        "final_modes":[serialize_mode(x) for x in final],
        "first_omitted_initial_mode":serialize_mode(initial_all[len(initial)]),
        "first_omitted_final_mode":serialize_mode(final_neighbor),
        "source_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
        "elapsed_seconds_diagnostic":time.perf_counter()-start,
        "scope":"Ideal one-round actual sharp-wall Ramsey counts; complete initial Gibbs preparation; no output projector; deterministic model bias certificate, not observed counts or apparatus calibration."}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["status"],"parameters":result["parameters"],"budgets":budget,
        "elapsed_seconds_diagnostic":result["elapsed_seconds_diagnostic"]},indent=2),flush=True)


if __name__ == "__main__":
    main()
