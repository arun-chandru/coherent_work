"""Finite-height continuum/full-Gibbs band-purity certificate.

Only floating root proposals use ordinary arithmetic. Accepted roots, their
Sturm indices, full band enumeration, signed continuum overlaps, complete
Gibbs tails, and all claimed inequalities use exact inputs and Arb balls.
No spatial mesh, numerical quadrature, or propagated truncated state is used.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from check_phy_info_continuum_certificate import (
    arb, ball, certified_mode, ctx, hull, secular, strict_positive,
)

ctx.prec = 256


def super_secular(k, barrier):
    q = (k*k-ball(barrier)).sqrt()
    return k*k.cos()*q.sin()+q*k.sin()*q.cos()


def certified_floor_pi(x, name):
    candidate = math.floor(float(x)/math.pi)
    strict_positive(x-candidate*arb.pi(), f"{name}: lower integer bracket")
    strict_positive((candidate+1)*arb.pi()-x, f"{name}: upper integer bracket")
    return candidate


def certified_super_modes(barrier, end):
    """Propose roots, then certify their exact nodal labels independently."""
    def fn(k):
        q = math.sqrt(k*k-barrier)
        return k*math.cos(k)*math.sin(q)+q*math.sin(k)*math.cos(q)

    candidates = []
    x = math.sqrt(barrier)+1e-6
    fx = fn(x)
    while x < end:
        y = min(x+.002, end)
        fy = fn(y)
        if fx*fy < 0:
            lo, hi, fl = x, y, fx
            for _ in range(55):
                mid = (lo+hi)/2
                fm = fn(mid)
                if fl*fm > 0:
                    lo, fl = mid, fm
                else:
                    hi = mid
            center = round((lo+hi)/2*10**8)
            candidates.append((Fraction(center-2, 10**8), Fraction(center+2, 10**8)))
        x, fx = y, fy
    result = []
    for lo, hi in candidates:
        fl, fh = super_secular(ball(lo), barrier), super_secular(ball(hi), barrier)
        strict_positive(-fl*fh, "superbarrier candidate sign change")
        sign_left = 1 if fl > 0 else -1
        for _ in range(80):
            mid = (lo+hi)/2
            fm = super_secular(ball(mid), barrier)
            if fm > 0:
                sign_mid = 1
            elif fm < 0:
                sign_mid = -1
            else:
                raise ArithmeticError("Indeterminate bisection sign")
            if sign_mid == sign_left:
                lo = mid
            else:
                hi = mid
        strict_positive(-super_secular(ball(lo), barrier)*super_secular(ball(hi), barrier),
                        "superbarrier final sign change")
        k = hull(lo, hi)
        strict_positive(k*k-barrier, "superbarrier condition")
        q = (k*k-barrier).sqrt()
        strict_positive(abs(k.sin()), "nonzero interface amplitude")
        strict_positive(abs(q.sin()), "nonzero outer sine denominator")
        # No zero at the interface: these are ALL interior zeros of the
        # piecewise sine eigenfunction. Sturm oscillation fixes its label.
        inside_nodes = certified_floor_pi(k, "interior nodal count")
        outside_nodes = certified_floor_pi(q, "exterior nodal count")
        n = inside_nodes+outside_nodes+1
        ext_shape = (ball("1/2")-(2*q).sin()/(4*q))/q.sin()**2
        exterior = k.sin()**2*ext_shape
        norm = ball("1/2")-(2*k).sin()/(4*k)+exterior
        strict_positive(norm, "positive superbarrier normalization")
        result.append(dict(n=n, kind="superbarrier", barrier=barrier,
                           lo=lo, hi=hi, k=k, q=q, energy=k*k,
                           exterior=exterior, norm=norm, source=exterior/norm))
    return result


def signed_overlap(i, f):
    """Analytic Green-identity overlap; both eigenfunctions have sine inside."""
    a = i["kap"]
    if f["kind"] == "subbarrier":
        b = f["kap"]
        strict_positive(abs(a-b), "nonzero hyperbolic cross denominator")
        integral = ((a+b).sinh()/(2*(a+b))-(a-b).sinh()/(2*(a-b))) \
                   / (a.sinh()*b.sinh())
    else:
        b = f["q"]
        integral = (a*a.coth()-b*b.cot())/(a*a+b*b)
    gap = f["energy"]-i["energy"]
    strict_positive(gap, "strict endpoint energy gap")
    return ((f["barrier"]-i["barrier"])*i["k"].sin()*f["k"].sin()*integral
            / (gap*(i["norm"]*f["norm"]).sqrt()))


def serialize_mode(m):
    fn = secular if m["kind"] == "subbarrier" else super_secular
    return dict(n=m["n"], kind=m["kind"], barrier=m["barrier"],
                k_lower_exact=str(m["lo"]), k_upper_exact=str(m["hi"]),
                energy_ball=str(m["energy"]),
                secular_lower_ball=str(fn(ball(m["lo"]), m["barrier"])),
                secular_upper_ball=str(fn(ball(m["hi"]), m["barrier"])),
                normalization_integral_ball=str(m["norm"]),
                exterior_probability_ball=str(m["source"]))


def gibbs_tail(theta, barrier, retained):
    """Neumann interface split. First retained levels must all be left branch."""
    pi2 = arb.pi()**2
    last_left = pi2*ball(Fraction(2*retained-1, 2))**2
    first_right = barrier+pi2/4
    strict_positive(first_right-last_left, "first retained split levels all left")
    # Left tail starts k=retained+1; right tail starts k=1. For k=q+1+j,
    # (k-1/2)^2 >= (q+1/2)^2 + 2(q+1)j.
    left = (-theta*pi2*ball(Fraction(2*retained+1, 2))**2).exp() \
           / (1-(-2*theta*pi2*(retained+1)).exp())
    right = (-theta*(barrier+pi2/4)).exp()/(1-(-2*theta*pi2).exp())
    return left+right


def action_bound(action, weights, projected, residual_sq, source, tail, zlo, zhi):
    static_q_lo = (sum((w*q for w,q in zip(weights, projected)), ball(0))/zhi).lower()
    static_q_hi = ((sum((w*q for w,q in zip(weights, projected)), ball(0))+tail)/zlo).upper()
    source_hi = ((sum((w*s for w,s in zip(weights, source)), ball(0))+tail)/zlo).upper()
    amplitude = static_q_lo.sqrt()-action*source_hi.sqrt()
    strict_positive(amplitude, "positive action-uniform population amplitude")
    yield_lo = (amplitude**2).lower()
    yield_hi = ((static_q_hi.sqrt()+action*source_hi.sqrt())**2).upper()
    # Each omitted state's perpendicular band mass is at most one, even
    # after the drive, so add its Gibbs weight once, with no Duhamel penalty.
    off_numerator = sum((w*(t.sqrt()+action*s.sqrt())**2
                         for w,t,s in zip(weights, residual_sq, source)), ball(0))+tail
    off_hi = (off_numerator/zlo).upper()
    fidelity_lo = (1-off_hi/yield_lo).lower()
    strict_positive(fidelity_lo, "positive common-mode fidelity")
    purity_lo = (fidelity_lo**2).lower()
    return dict(action_exact=str(action), static_band_probability_lower=str(static_q_lo),
                static_band_probability_upper=str(static_q_hi), source_probability_upper=str(source_hi),
                yield_lower=str(yield_lo), yield_upper=str(yield_hi),
                orthogonal_band_mass_upper=str(off_hi),
                common_mode_fidelity_lower=str(fidelity_lo), purity_lower=str(purity_lo)), purity_lo, yield_lo


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    initial_barrier, final_barrier, retained = 2500, 10000, 16
    theta = ball("1/100")
    initial = [certified_mode(n, initial_barrier) for n in range(1, retained+1)]
    # All subbarrier states plus superbarrier candidates through beyond the
    # upper band edge. Contiguous exact Sturm labels certify completeness.
    final_modes = [certified_mode(n, final_barrier) for n in range(1, 33)]
    final_modes += certified_super_modes(final_barrier, math.sqrt(8*initial_barrier)+4)
    labels = [m["n"] for m in final_modes]
    if labels != list(range(1, len(final_modes)+1)):
        raise ArithmeticError(f"Final Sturm labels have a gap or duplicate: {labels}")
    for previous, current in zip(final_modes, final_modes[1:]):
        if not previous["hi"] < current["lo"]:
            raise ArithmeticError("Final spectral brackets are not ordered/disjoint")
    band_lo, band_hi = 2*initial_barrier, 8*initial_barrier
    final = []
    below, above = [], []
    for m in final_modes:
        if m["energy"] < band_lo:
            below.append(m)
        elif m["energy"] > band_hi:
            above.append(m)
        elif m["energy"] > band_lo and m["energy"] < band_hi:
            final.append(m)
        else:
            raise ArithmeticError("Band edge does not separate from root enclosure")
    if not below or not above or not final:
        raise ArithmeticError("Missing bracketing states outside the complete band")
    if not (below[-1]["n"]+1 == final[0]["n"] and final[-1]["n"]+1 == above[0]["n"]):
        raise ArithmeticError("Full band and neighboring outside indices not contiguous")
    overlap = [[signed_overlap(i, f) for f in final] for i in initial]
    projected = [sum((x*x for x in row), ball(0)) for row in overlap]
    strict_positive(projected[0], "reference projected mode nonzero")
    residual_sq = [ball(0)]
    for n, row in enumerate(overlap[1:], 2):
        cross = sum((a*b for a,b in zip(overlap[0], row)), ball(0))
        t2 = projected[n-1]-cross**2/projected[0]
        strict_positive(t2, f"strict static perpendicular residual, label {n}")
        residual_sq.append(t2)
    weights = [(-theta*m["energy"]).exp() for m in initial]
    tail = gibbs_tail(theta, initial_barrier, retained)
    tail2 = gibbs_tail(2*theta, initial_barrier, retained)
    partial_z = sum(weights, ball(0))
    zlo, zhi = partial_z.lower(), (partial_z+tail).upper()
    input_purity_lo = (sum((w*w for w in weights), ball(0))/zhi**2).lower()
    input_purity_hi = ((sum((w*w for w in weights), ball(0))+tail2)/zlo**2).upper()
    strict_positive(ball("3/10")-input_purity_hi, "complete input purity < 0.3")
    reference_weights = [x*x/projected[0] for x in overlap[0]]
    mu_hi = max(x.upper() for x in reference_weights)
    strict_positive(ball("78652/1000000")-mu_hi, "reference maximum final-energy weight < 0.078652")
    static_dephased_numerator = sum((wm*sum((wn*row[m]**2
                                     for wn,row in zip(weights,overlap)), ball(0))
                                     for m,wm in enumerate(reference_weights)), ball(0))
    static_dephased_lo = (static_dephased_numerator/zhi).lower()
    static_dephased_hi = ((static_dephased_numerator+mu_hi*tail)/zlo).upper()
    source = [m["source"] for m in initial]
    actions = []
    principal = None
    for rational in ("0", "1/80", "1/40", "1/20", "1/8"):
        report, purity, probability = action_bound(ball(rational), weights, projected,
                                                  residual_sq, source, tail, zlo, zhi)
        report["action_rational"] = rational
        actions.append(report)
        if rational == "1/80":
            strict_positive(purity-ball("9972/10000"), "principal output purity > 0.9972")
            strict_positive(probability-ball("16608/100000000"), "principal yield > 0.00016608")
            principal = report
    # Entire initial-energy tail for matched monotone-classical envelope.
    # Standard Chernoff + length-two free-box Gaussian partition bound,
    # evaluated at a=1/(2u), u=1250: theta-a>0.
    u = ball(1250)
    strict_positive(theta-1/(2*u), "Chernoff inverse-temperature admissibility")
    classical_hi = (2*u.sqrt()/(2*arb.pi()).sqrt()/zlo
                    *(ball("1/2")-theta*u).exp()).upper()
    finite_tail_numerator = tail
    for w,m in zip(weights,initial):
        if m["energy"] > u:
            finite_tail_numerator += w
        elif not m["energy"] < u:
            raise ArithmeticError("Initial mode interval straddles classical threshold")
    classical_finite_hi = (finite_tail_numerator/zlo).upper()
    principal_yield_lo = arb(principal["yield_lower"])
    principal_yield_hi = arb(principal["yield_upper"])
    principal_off_hi = arb(principal["orthogonal_band_mass_upper"])
    bright_rate_lo = (principal_yield_lo-principal_off_hi).lower()
    diagonal_rate_hi = (mu_hi*principal_yield_hi).upper()
    readout_gap_lo = (bright_rate_lo-diagonal_rate_hi).lower()
    strict_positive(arb(principal["common_mode_fidelity_lower"])-ball("9986/10000"),
                    "principal common-mode fidelity > 0.9986")
    strict_positive(ball("18066/100000000")-principal_yield_hi, "principal yield < 0.00018066")
    strict_positive(bright_rate_lo-ball("16585/100000000"), "driven bright rate > 0.00016585")
    strict_positive(ball("1421/100000000")-diagonal_rate_hi, "diagonal bright rate < 0.00001421")
    strict_positive(readout_gap_lo-ball("15164/100000000"), "bright rate gap > 0.00015164")
    strict_positive(ball("58/100000000")-classical_finite_hi, "classical initial tail < 0.00000058")
    result = dict(
        status="validated continuum full-Gibbs coherent-band certificate: all assertions passed",
        arithmetic="FLINT Arb, 256 bits; exact rational inputs; certified signed root brackets",
        parameters=dict(r="1/50", theta="1/100", initial_barrier=initial_barrier,
                        final_barrier=final_barrier, scaled_band=[2,8],
                        principal_action="1/80", linear_ramp_scaled_duration="1/120",
                        linear_ramp_duration="1/300000"),
        complete_band=dict(rank=len(final), first_index=final[0]["n"], last_index=final[-1]["n"],
                           lower_neighbor=serialize_mode(below[-1]), upper_neighbor=serialize_mode(above[0])),
        partition_lower=str(zlo), partition_upper=str(zhi),
        complete_omitted_unnormalized_gibbs_bound=str(tail),
        complete_omitted_probability_upper=str((tail/zlo).upper()),
        input_purity_lower=str(input_purity_lo), input_purity_upper=str(input_purity_hi),
        principal_certificate=principal, action_tradeoff=actions,
        matched_monotone_classical_tail_upper=str(classical_hi),
        sharper_matched_monotone_classical_tail_upper=str(classical_finite_hi),
        reference_readout=dict(reference_projected_norm_squared=str(projected[0]),
                              max_energy_weight_upper=str(mu_hi),
                              energy_weight_balls=[str(x) for x in reference_weights],
                              signed_reference_coefficients=[str(x/projected[0].sqrt()) for x in overlap[0]],
                              static_dephased_bright_rate_lower=str(static_dephased_lo),
                              static_dephased_bright_rate_upper=str(static_dephased_hi),
                              driven_bright_rate_lower=str(bright_rate_lo),
                              any_same_mass_diagonal_bright_rate_upper=str(diagonal_rate_hi),
                              unconditional_bright_rate_gap_lower=str(readout_gap_lo)),
        initial_modes=[serialize_mode(m) for m in initial],
        all_final_modes_through_upper_neighbor=[serialize_mode(m) for m in final_modes
                                               if m["n"] <= above[0]["n"]],
        static_projected_probability_balls=[str(x) for x in projected],
        static_orthogonal_residual_squared_balls=[str(x) for x in residual_sq],
        reference_mode="w=P_[5000,20000](H_f) phi_1^i / ||P_[5000,20000](H_f) phi_1^i||",
        certified_rational_claims=dict(input_purity_upper="3/10", output_purity_lower="9972/10000",
                                       band_yield_lower="16608/100000000", band_yield_upper="18066/100000000",
                                       common_mode_fidelity_lower="9986/10000",
                                       max_reference_energy_weight_upper="78652/1000000",
                                       driven_bright_rate_lower="16585/100000000",
                                       same_mass_diagonal_bright_rate_upper="1421/100000000",
                                       unconditional_bright_rate_gap_lower="15164/100000000",
                                       matched_monotone_classical_tail_upper="58/100000000"),
        scope="Single particle, exact step, full Gibbs; coherent band block; every measurable height control with endpoints 1 and 4 and integrated action <=1/80. Monotonicity only for classical envelope.",
    )
    output = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(output, encoding="utf-8")
    print(json.dumps({k:result[k] for k in ("status", "parameters", "complete_band",
                     "complete_omitted_probability_upper", "input_purity_lower", "input_purity_upper",
                     "principal_certificate", "action_tradeoff", "sharper_matched_monotone_classical_tail_upper")}, indent=2))
    print(json.dumps({k:v for k,v in result["reference_readout"].items()
                      if k not in ("energy_weight_balls", "signed_reference_coefficients")}, indent=2))


if __name__ == "__main__":
    main()
