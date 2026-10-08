"""Validated continuum endpoint certificate, using FLINT/Arb ball arithmetic.

Floating arithmetic is used only to propose root brackets. Every root sign,
spectral membership, normalization, overlap, thermal bound, and final inequality
is subsequently checked with Arb enclosures of exact rational inputs.

Requires the externally installed python-flint distribution.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from flint import arb, ctx, fmpq

ctx.prec = 256


def ball(q):
    """Enclose an exact Fraction, integer, or rational string."""
    if isinstance(q, arb):
        return q
    if isinstance(q, float):
        raise TypeError("Use an exact rational input, never a Python float")
    q = Fraction(q)
    return arb(fmpq(q.numerator, q.denominator))


def hull(lo, hi):
    if lo > hi:
        raise ArithmeticError("Reversed rational bracket")
    # Arb's union encloses both rounded exact endpoints and everything between.
    return ball(lo).union(ball(hi))


def strict_positive(x, message):
    if not x > 0:
        raise ArithmeticError(f"Unproved positive inequality: {message}: {x}")


def secular(k, barrier):
    kap = (ball(barrier) - k*k).sqrt()
    return k*k.cos() + kap*kap.coth()*k.sin()


def propose_root(n, barrier):
    """Noncertifying localization only; returned bracket still needs signs."""
    lo = (n-.5)*math.pi + 1e-8
    hi = min(n*math.pi - 1e-8, math.sqrt(barrier)-1e-8)
    def fn(k):
        kap = math.sqrt(barrier-k*k)
        return k*math.cos(k)+kap/math.tanh(kap)*math.sin(k)
    flo = fn(lo)
    if not flo*fn(hi) < 0:
        raise ArithmeticError("Floating localization did not find a candidate")
    for _ in range(60):
        mid = (lo+hi)/2
        if flo*fn(mid) > 0:
            lo, flo = mid, fn(mid)
        else:
            hi = mid
    center = round((lo+hi)/2 * 10**8)
    return Fraction(center-2, 10**8), Fraction(center+2, 10**8)


def certified_mode(n, barrier):
    lo, hi = propose_root(n, barrier)
    f_lo, f_hi = secular(ball(lo), barrier), secular(ball(hi), barrier)
    strict_positive(-f_lo*f_hi, "opposite root bracket signs")
    sign_left = 1 if f_lo > 0 else -1
    # Exact-rational bisection; every decision excludes zero in the ball.
    for _ in range(80):
        mid = (lo+hi)/2
        f_mid = secular(ball(mid), barrier)
        if f_mid > 0:
            sign_mid = 1
        elif f_mid < 0:
            sign_mid = -1
        else:
            raise ArithmeticError("Precision insufficient for a root sign")
        if sign_mid == sign_left:
            lo = mid
        else:
            hi = mid
    strict_positive(-secular(ball(lo),barrier)*secular(ball(hi),barrier),
                    "final exact bracket signs")
    k = hull(lo, hi)
    strict_positive(k-(ball(n)-ball("1/2"))*arb.pi(), "left nodal interval")
    strict_positive(ball(n)*arb.pi()-k, "right nodal interval")
    strict_positive(ball(barrier)-k*k, "subbarrier condition")
    # Existence: opposite signs and continuity. Sturm oscillation: each root
    # has n-1 interior nodes and no exterior nodes; hence it is eigenvalue n.
    kap = (ball(barrier)-k*k).sqrt()
    hsquare = kap.coth()/(2*kap)-1/(2*kap.sinh()**2)
    exterior = k.sin()**2*hsquare
    norm = ball("1/2")-(2*k).sin()/(4*k)+exterior
    strict_positive(norm, "positive eigenfunction norm")
    return dict(n=n, kind="subbarrier", barrier=barrier, lo=lo, hi=hi, k=k, kap=kap,
                energy=k*k, exterior=exterior, norm=norm,
                source=exterior/norm)


def super_secular(k):
    q = (k*k-2500).sqrt()
    return k*k.cos()*q.sin()+q*k.sin()*q.cos()


def certified_super_modes():
    """Find a subset of roots above Bf; completeness is not needed."""
    def fn(k):
        q = math.sqrt(k*k-2500)
        return k*math.cos(k)*math.sin(q)+q*math.sin(k)*math.cos(q)
    candidates = []
    step = .001
    x = 50.000001
    fx = fn(x)
    end = math.sqrt(625*7.9)
    while x < end:
        y = min(x+step,end)
        fy = fn(y)
        if fx*fy < 0:
            left, right, fleft = x,y,fx
            for _ in range(50):
                mid = (left+right)/2
                fm = fn(mid)
                if fleft*fm > 0:
                    left, fleft = mid,fm
                else:
                    right = mid
            center = round((left+right)/2*10**8)
            candidates.append((Fraction(center-2,10**8),Fraction(center+2,10**8)))
        x,fx = y,fy
    modes = []
    for index,(lo,hi) in enumerate(candidates,1):
        flo,fhi = super_secular(ball(lo)),super_secular(ball(hi))
        strict_positive(-flo*fhi,"suprabarrier initial sign change")
        sign_left = 1 if flo > 0 else -1
        for _ in range(80):
            mid = (lo+hi)/2
            fm = super_secular(ball(mid))
            if fm > 0:
                sign_mid = 1
            elif fm < 0:
                sign_mid = -1
            else:
                raise ArithmeticError("Unresolved suprabarrier root sign")
            if sign_mid == sign_left:
                lo = mid
            else:
                hi = mid
        strict_positive(-super_secular(ball(lo))*super_secular(ball(hi)),
                        "suprabarrier final sign change")
        if modes:
            if not modes[-1]["hi"] < lo:
                raise ArithmeticError("Selected suprabarrier brackets overlap")
        k = hull(lo,hi)
        strict_positive(k*k-2500,"strictly above final barrier")
        q = (k*k-2500).sqrt()
        strict_positive(abs(q.sin()),"nonzero exterior sine denominator")
        hsquare = (ball("1/2")-(2*q).sin()/(4*q))/q.sin()**2
        exterior = k.sin()**2*hsquare
        norm = ball("1/2")-(2*k).sin()/(4*k)+exterior
        strict_positive(norm,"suprabarrier positive normalization")
        modes.append(dict(n=None,subset_index=index,kind="superbarrier",barrier=2500,
                          lo=lo,hi=hi,k=k,q=q,energy=k*k,exterior=exterior,
                          norm=norm,source=exterior/norm))
    return modes


def overlap_square(i, f):
    a = i["kap"]
    if f["kind"] == "subbarrier":
        b = f["kap"]
        strict_positive(abs(a-b), "nonzero sinh product denominator")
        integral = ((a+b).sinh()/(2*(a+b))
                    -(a-b).sinh()/(2*(a-b)))/(a.sinh()*b.sinh())
    else:
        b = f["q"]
        integral = (a*a.coth()-b*b.cot())/(a*a+b*b)
    cross = i["k"].sin()*f["k"].sin()*integral
    # Green's identity: (Ef-Ei)<f,i>=(Bf-Bi)<f,Q i>.
    difference = f["energy"]-i["energy"]
    strict_positive(difference, "endpoint energy gap")
    return ((ball(f["barrier"]-i["barrier"])*cross/difference)**2
            /(i["norm"]*f["norm"]))


def serial_mode(m):
    secular_fn = (lambda k: secular(k,m["barrier"])) if m["kind"] == "subbarrier" else super_secular
    return {"n": m["n"], "kind":m["kind"], "barrier": m["barrier"],
            "k_lower_exact": str(m["lo"]), "k_upper_exact": str(m["hi"]),
            "energy_ball": str(m["energy"]),
            "scaled_energy_ball": str(m["energy"]/625),
            "secular_lower_ball": str(secular_fn(ball(m["lo"]))),
            "secular_upper_ball": str(secular_fn(ball(m["hi"]))),
            "normalization_integral_ball": str(m["norm"]),
            "exterior_probability_ball": str(m["source"])}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    initial = [certified_mode(n,625) for n in range(1,9)]
    final_sub = [certified_mode(n,2500) for n in range(12,17)]
    final = final_sub+certified_super_modes()
    for m in final:
        # This also places each selected mode in O^{-h}=[2.1,7.9], for
        # O=[2,8] and readout guard h=1/10.
        strict_positive(m["energy"]/625-ball("21/10"), "lower guarded membership")
        strict_positive(ball("79/10")-m["energy"]/625, "upper guarded membership")
    theta = ball("1/10")
    weights = [(-theta*m["energy"]).exp() for m in initial]
    partial_z = sum(weights, ball(0))
    c = theta*arb.pi()**2/4
    remainder = (-c*81).exp()/(1-(-c*19).exp())
    # Point endpoints are rigorous outward bounds and remain Arb objects.
    zlo, zhi = partial_z.lower(), (partial_z+remainder).upper()
    projected = [sum((overlap_square(i,f) for f in final),ball(0))
                 for i in initial]
    q_lower = (sum((w*a for w,a in zip(weights,projected)),ball(0))/zhi).lower()
    source_upper = ((sum((w*m["source"] for w,m in zip(weights,initial)),ball(0))
                     +remainder)/zlo).upper()
    action = ball("3/8")
    amplitude = q_lower.sqrt()-action*source_upper.sqrt()
    strict_positive(amplitude, "action floor amplitude")
    floor = amplitude**2
    # FE51: theta=1/10, L=2, threshold u=2/(4r^2)=625/2.
    u = ball("625/2")
    classical_upper = (2*u.sqrt()/(2*arb.pi()).sqrt()/zlo
                       *(ball("1/2")-theta*u).exp())
    gap = floor-classical_upper
    strict_positive(gap, "strict continuum/classical separation")
    # Coarse exact rational inequalities are the reported certificate.
    strict_positive(q_lower-ball("22980717/100000000000"), "Q > 0.00022980717")
    strict_positive(ball("674811/1000000000")-source_upper, "q2 < 0.000674811")
    strict_positive(floor-ball("293544/10000000000"), "floor > 2.93544e-5")
    strict_positive(ball("1457/1000000000000000")-classical_upper, "tail < 1.457e-12")
    strict_positive(gap-ball("293544/10000000000"), "gap > 2.93544e-5")
    guarded_u = ball("19/10")*625/4
    guarded_tail = (2*guarded_u.sqrt()/(2*arb.pi()).sqrt()/zlo
                    *(ball("1/2")-theta*guarded_u).exp())
    epsilon = ball("1/100000")
    delta_drive = delta_null = ball("1/1000000")
    guarded_gap = floor-guarded_tail-2*epsilon-delta_drive-delta_null
    strict_positive(ball("7/1000000000000")-guarded_tail, "guarded tail < 7e-12")
    strict_positive(guarded_gap-ball("73544/10000000000"), "guarded gap > 7.3544e-6")
    sub_projected = [sum((overlap_square(i,f) for f in final_sub),ball(0)) for i in initial]
    sub_q_lower = (sum((w*a for w,a in zip(weights,sub_projected)),ball(0))/zhi).lower()
    sub_amplitude = sub_q_lower.sqrt()-action*source_upper.sqrt()
    strict_positive(sub_amplitude,"five-state base amplitude")
    sub_floor = sub_amplitude**2
    strict_positive(sub_floor-ball("5695/10000000000"),"five-state base floor > 5.695e-7")
    result = {
        "status": "validated continuum certificate: all assertions passed",
        "arithmetic": "FLINT Arb rigorous balls; precision 256 bits; exact rational inputs",
        "parameters": {"r": "1/25", "theta": "1/10", "initial_barrier":625,
                       "final_barrier":2500,"scaled_band":[2,8],"action":"3/8"},
        "initial_modes": [serial_mode(m) for m in initial],
        "selected_final_modes": [serial_mode(m) for m in final],
        "unnormalized_omitted_gibbs_bound_ball": str(remainder),
        "partition_lower": str(zlo), "partition_upper": str(zhi),
        "band_probability_lower": str(q_lower), "source_probability_upper": str(source_upper),
        "uniform_band_floor_ball": str(floor),
        "full_initial_tail_upper_ball": str(classical_upper),
        "separation_gap_ball": str(gap),
        "five_subbarrier_state_base_certificate": {"scaled_band":[2,4],
            "band_probability_lower":str(sub_q_lower),"band_floor_ball":str(sub_floor)},
        "certified_rational_bounds": {
            "Q_lower_strict": "22980717/100000000000",
            "q2_upper_strict": "674811/1000000000",
            "band_floor_lower_strict": "293544/10000000000",
            "initial_tail_upper_strict": "1457/1000000000000000",
            "gap_lower_strict": "293544/10000000000"},
        "conditional_detector_budget": {
            "observed_band": [2,8], "guard_h": "1/10",
            "failure_probability": "1/100000",
            "delta_drive": "1/1000000", "delta_null": "1/1000000",
            "guarded_initial_tail_upper_ball": str(guarded_tail),
            "predicted_observed_gap_ball": str(guarded_gap),
            "predicted_observed_gap_lower_strict": "73544/10000000000",
            "scope": "Sufficient hypothetical calibrated budgets, not apparatus validation"},
        "initial_projected_probability_balls": [str(v) for v in projected],
    }
    output = json.dumps(result,indent=2)+"\n"
    if args.output:
        args.output.write_text(output,encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
