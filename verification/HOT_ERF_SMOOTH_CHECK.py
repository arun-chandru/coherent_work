"""Validated actual analytic erf wall, full Gibbs preparation and all tails.

Exact Taylor ODE enclosures on the resolved collar; analytic Cauchy remainders,
Sturm node labels, polynomial integral overlaps, and full-space Gram bounds.
Floats only propose root brackets. No physical spectral or spatial projection.
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
from check_hot_ramsey_finite_certificate import ball, free_gibbs_tail, lower_nonnegative, strict_positive
from flint import arb, acb, arb_poly, ctx
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

ctx.prec = 256
R = Fraction(1, 20)
WIDTH = R / 2
LEFT = 1 - 5 * WIDTH
CELL = WIDTH / 5
NCELL = 50
ORDER = 48
CAUCHY_RADIUS = 4
INITIAL_BARRIER = int(1 / (R * R))
FINAL_BARRIER = 4 * INITIAL_BARRIER
S = 8
INPUT_COUNT = 30
E_CUT = 100 * INITIAL_BARRIER
ERF5 = ball(5).erf()


def sign(x, name):
    if x > 0:
        return 1
    if x < 0:
        return -1
    raise ArithmeticError(f"Unresolved sign {name}: {x}")


def inflate(x, error):
    return x + arb(0, error.upper())


def floor_pi(x, name):
    proposed = math.floor(float(x) / math.pi)
    strict_positive(x - proposed * arb.pi(), name + " lower")
    strict_positive((proposed + 1) * arb.pi() - x, name + " upper")
    return proposed


def erf_coefficients(cell):
    """Coefficients of [erf(-5+(cell+x)/5)+erf(5)]/[2erf(5)]."""
    z, dz = ball(Fraction(-25 + cell, 5)), ball(Fraction(1, 5))
    g = [(-z*z).exp()]
    # g(x)=exp[-(z+dz*x)^2]; g'=-2dz(z+dz*x)g.
    for n in range(ORDER):
        g.append((-2*dz*z*g[n] - (2*dz*dz*g[n-1] if n else 0))/(n+1))
    f = [(z.erf()+ERF5)/(2*ERF5)]
    for n in range(ORDER):
        f.append(dz*g[n]/(arb.pi().sqrt()*ERF5*(n+1)))
    return f


PROFILE_COEFF = [erf_coefficients(j) for j in range(NCELL)]
# For |x|<=R_c, |Im(z+x/5)|<=R_c/5. Integrate erf vertically
# from its real part: |erf(zeta)|<=1+2 b exp(b^2)/sqrt(pi).
B_IM = ball(Fraction(CAUCHY_RADIUS, 5))
F_COMPLEX_BOUND = (1 + (1 + 2*B_IM*(B_IM*B_IM).exp()/arb.pi().sqrt())/ERF5)/2


def ode_cell(u, v, energy, barrier, cell, keep=False):
    """Analytic solution coefficients and rigorous complex-Cauchy remainder.

    x is the local unit coordinate, v is physical y derivative. Y=(u,h*v)
    solves Y'=[[0,1],[h^2(aF-E),0]]Y. Entire F and Gronwall on radial
    complex paths bound ||Y||_infty by exp(R_c*M)||Y(0)||_infty.
    """
    h = ball(CELL)
    a = [u, h*v]
    f = PROFILE_COEFF[cell]
    h2 = h*h
    for n in range(ORDER):
        convolution = sum((f[j]*a[n-j] for j in range(n+1)), ball(0))
        a.append(h2*(barrier*convolution-energy*a[n])/((n+1)*(n+2)))
    matrix_bound = max(ball(1), (h2*(abs(energy)+barrier*F_COMPLEX_BOUND)).upper())
    initial_bound = max(abs(u).upper(), abs(h*v).upper())
    analytic_bound = (initial_bound*(CAUCHY_RADIUS*matrix_bound).exp()).upper()
    remainder = (analytic_bound / ((CAUCHY_RADIUS-1)*ball(CAUCHY_RADIUS)**ORDER)).upper()
    # Both components are truncated through x^ORDER. u needs a[0:ORDER+1]
    # and h*v needs (j+1)a[j+1], for which a[ORDER+1] was computed.
    up = sum(a[:ORDER+1], ball(0))
    vp = sum(((j+1)*a[j+1] for j in range(ORDER+1)), ball(0))/h
    ue, ve = inflate(up, remainder), inflate(vp, remainder/h)
    if not keep:
        return ue, ve
    coefficients = a[:ORDER+1]
    polynomial = arb_poly(coefficients)
    polynomial_bound = sum((abs(c).upper() for c in coefficients), ball(0)).upper()
    return ue, ve, {"poly": polynomial, "bound": polynomial_bound, "error": remainder}


def collar(k, barrier, keep=False):
    energy, length = k*k, ball(LEFT)
    u, v = (k*length).sin(), k*(k*length).cos()
    # Center the numerical balls after every cell and carry one Euclidean
    # error in (u,v/kappa). The real logarithmic norm is small even when
    # the solution oscillates; this prevents artificial box wrapping.
    kappa = ball(math.ceil(math.sqrt(float(energy)+barrier)))
    logarithmic_norm = ((abs(kappa*kappa-energy)+barrier)/(2*kappa)).upper()
    growth = (ball(CELL)*logarithmic_norm).exp().upper()
    error = (u.rad()**2+(v.rad()/kappa)**2).sqrt().upper()
    u, v = u.mid(), v.mid()
    states, cells = [(inflate(u,error),inflate(v,kappa*error))], []
    for j in range(NCELL):
        ue, ve, data = ode_cell(u, v, energy, barrier, j, True)
        propagated_error = (growth*error).upper()
        data["error"] = (data["error"]+propagated_error).upper()
        local_error = (ue.rad()**2+(ve.rad()/kappa)**2).sqrt().upper()
        error = (propagated_error+local_error).upper()
        u, v = ue.mid(), ve.mid()
        if keep:
            cells.append(data)
        states.append((inflate(u,error),inflate(v,kappa*error)))
    return states, cells


def secular(k, barrier):
    u, v = collar(k, barrier)[0][-1]
    d, length = k*k-barrier, ball(LEFT)
    if d > 0:
        q = d.sqrt()
        return u*(q*length).cos()+v*length*(q*length).sinc()
    if d < 0:
        q = (-d).sqrt()
        return u+v*(q*length).tanh()/q
    raise ArithmeticError("Threshold unresolved in secular function")


def float_secular(k, barrier):
    length, width = float(LEFT), float(WIDTH)
    normalization = math.erf(5)
    def rhs(y, state):
        f = (math.erf((y-1)/width)+normalization)/(2*normalization)
        return state[1], (barrier*f-k*k)*state[0]
    sol = solve_ivp(rhs, (length, 2-length),
                    (math.sin(k*length), k*math.cos(k*length)),
                    method="DOP853", rtol=2e-13, atol=2e-14)
    if not sol.success:
        raise ArithmeticError("Float root proposal integration failed")
    u, v = sol.y[:, -1]
    d = k*k-barrier
    if d > 0:
        q = math.sqrt(d)
        return u*math.cos(q*length)+v*math.sin(q*length)/q
    if d < 0:
        q = math.sqrt(-d)
        return u+v*math.tanh(q*length)/q
    return u+v*length


def product_integral(ca, cb):
    integrated = (ca["poly"]*cb["poly"]).integral()(ball(1))*ball(CELL)
    error = ball(CELL)*(ca["error"]*cb["bound"]+cb["error"]*ca["bound"]
                        +ca["error"]*cb["error"])
    return inflate(integrated, error)


def certify_mode(lo, hi, barrier):
    fl, fh = secular(ball(lo), barrier), secular(ball(hi), barrier)
    strict_positive(-fl*fh, "signed rational root bracket")
    # Secants only propose much narrower rational brackets. Every candidate
    # endpoint is checked again with the rigorous ODE before being accepted.
    proposed = (ball(lo)*fh-ball(hi)*fl)/(fh-fl)
    center = Fraction(proposed.mid().str(70, radius=False))
    for denominator in (10**16,10**14,10**12):
        newlo,newhi = center-Fraction(1,denominator),center+Fraction(1,denominator)
        if not lo < newlo < newhi < hi:
            continue
        nfl,nfh = secular(ball(newlo),barrier),secular(ball(newhi),barrier)
        if -nfl*nfh > 0:
            lo,hi,fl,fh = newlo,newhi,nfl,nfh
            break
    if hi-lo > Fraction(1,10**13):
        raise ArithmeticError("Could not certify a narrow root bracket")
    k = ball(lo).union(ball(hi))
    states, cells = collar(k, barrier, True)
    nodes = floor_pi(k*ball(LEFT), "left Sturm floor")
    strict_positive(arb.pi()-k*ball(CELL), "at most one zero in each cell")
    for u, _ in states:
        strict_positive(abs(u), "nonzero cell endpoint")
    node_counts = [int(sign(a[0], "cell entrance") != sign(b[0], "cell exit"))
                   for a, b in zip(states, states[1:])]
    nodes += sum(node_counts)
    u, _ = states[-1]
    d, length = k*k-barrier, ball(LEFT)
    if d > 0:
        q = d.sqrt()
        strict_positive(abs((q*length).sin()), "right sine denominator")
        right_nodes = floor_pi(q*length, "right Sturm floor")
        right_norm = u*u*(length/2-(2*q*length).sin()/(4*q))/(q*length).sin()**2
        kind = "superbarrier"
    elif d < 0:
        q = (-d).sqrt()
        right_nodes = 0
        right_norm = u*u*((q*length).coth()/(2*q)-length/(2*(q*length).sinh()**2))
        kind = "subbarrier"
    else:
        raise ArithmeticError("Right threshold unresolved")
    left_norm = length/2-(2*k*length).sin()/(4*k)
    cell_norms = [product_integral(c, c) for c in cells]
    for val in cell_norms:
        strict_positive(val, "positive cell norm")
    norm = left_norm+sum(cell_norms, ball(0))+right_norm
    strict_positive(norm, "positive full norm")
    # Since 0<=F<=1 and F=0 on the left, this safely bounds ||F phi||^2.
    source_upper = ((sum(cell_norms, ball(0))+right_norm)/norm).upper()
    return {"n": nodes+right_nodes+1, "kind": kind, "lo": lo, "hi": hi,
            "k": k, "energy": k*k, "barrier": barrier, "states": states,
            "cells": cells, "norm": norm, "source_upper": source_upper,
            "left_norm": left_norm, "right_norm": right_norm,
            "cell_norms": cell_norms, "node_counts": node_counts,
            "right_nodes": right_nodes, "secular_lo": fl, "secular_hi": fh}


def complete_prefix(barrier, end, cache_path):
    if cache_path.exists():
        raw = json.loads(cache_path.read_text(encoding="utf-8"))
        expected = {"barrier": barrier, "end": end, "width": str(WIDTH)}
        if raw["parameters"] != expected:
            raise ArithmeticError("Float proposal cache belongs to another profile")
        proposals = [(Fraction(a), Fraction(b)) for a,b in raw["brackets"]]
    else:
        proposals = []
        x, fx = .01, float_secular(.01, barrier)
        while x < end:
            y = min(x+.2, end)
            fy = float_secular(y, barrier)
            if fx*fy < 0:
                root = brentq(lambda k: float_secular(k, barrier), x, y, xtol=2e-13)
                center = round(root*10**10)
                proposals.append((Fraction(center-2, 10**10), Fraction(center+2, 10**10)))
            x, fx = y, fy
        cache_path.write_text(json.dumps({"parameters": {"barrier": barrier, "end": end, "width": str(WIDTH)},
            "warning": "Untrusted floating proposals only; every use rechecks all root signs and Sturm labels.",
            "brackets": [[str(a),str(b)] for a,b in proposals]}, indent=2)+"\n", encoding="utf-8")
    modes = []
    for i, (lo,hi) in enumerate(proposals):
        mode = certify_mode(lo,hi,barrier)
        if mode["n"] != i+1:
            raise ArithmeticError(f"Missing or repeated exact Sturm label: expected {i+1}, got {mode['n']}")
        if modes and not modes[-1]["hi"] < lo:
            raise ArithmeticError("Overlapping or unordered exact root brackets")
        modes.append(mode)
        if len(modes) % 20 == 0:
            print(f"Certified barrier {barrier}: {len(modes)} contiguous Sturm modes", flush=True)
    return modes


def signed_overlap(a, b):
    length = ball(LEFT)
    left = length*((length*(a["k"]-b["k"])).sinc()
                   -(length*(a["k"]+b["k"])).sinc())/2
    qa, qb = acb(a["energy"]-a["barrier"]).sqrt(), acb(b["energy"]-b["barrier"]).sqrt()
    denom = 2*(qa*length).sin()*(qb*length).sin()
    strict_positive(abs(denom), "right overlap sine denominator")
    right_shape = (length*((length*(qa-qb)).sinc()-(length*(qa+qb)).sinc())/denom).real
    right = a["states"][-1][0]*b["states"][-1][0]*right_shape
    middle = sum((product_integral(ca,cb) for ca,cb in zip(a["cells"],b["cells"])),ball(0))
    return (left+middle+right)/(a["norm"]*b["norm"]).sqrt()


def serialize_mode(m):
    return {"n": m["n"], "kind": m["kind"], "barrier": m["barrier"],
        "k_lower_exact": str(m["lo"]), "k_upper_exact": str(m["hi"]),
        "energy_ball": str(m["energy"]), "normalization_ball": str(m["norm"]),
        "left_norm_ball": str(m["left_norm"]), "right_norm_ball": str(m["right_norm"]),
        "cell_norm_balls": [str(v) for v in m["cell_norms"]],
        "F_squared_mass_upper": str(m["source_upper"]),
        "cell_node_counts": m["node_counts"], "right_nodes": m["right_nodes"],
        "largest_cell_Taylor_remainder": str(max(c["error"] for c in m["cells"])),
        "secular_lo_ball": str(m["secular_lo"]), "secular_hi_ball": str(m["secular_hi"])}


def evaluate(initial, final, neighbor, s):
    overlap = [[signed_overlap(a,b) for b in final] for a in initial]
    tails = []
    for a, row in zip(initial,overlap):
        gap = neighbor["energy"].lower()-a["energy"]
        strict_positive(gap,"final tail above retained initial mode")
        spectral = (9*INITIAL_BARRIER**2*a["source_upper"]/gap**2).upper()
        residual = 1-sum((v*v for v in row),ball(0))
        # The exact residual is nonnegative; an enclosure can cross zero.
        strict_positive(residual.upper(),"nonnegative overlap residual enclosure")
        tails.append(min(spectral,residual.upper()))
    theta, duration = ball(s)/INITIAL_BARRIER, 2*arb.pi()/(3*INITIAL_BARRIER)
    weights = [(-theta*a["energy"]).exp() for a in initial]
    tail, tail2 = free_gibbs_tail(theta,len(initial)),free_gibbs_tail(2*theta,len(initial))
    zp = sum(weights,ball(0))
    zlo,zhi = zp.lower(),(zp+tail).upper()
    pin = ((sum((w*w for w in weights),ball(0))+tail2)/zlo**2).upper()
    phase_b = [acb(0,-duration*b["energy"]).exp() for b in final]
    phase_a = [acb(0,duration*a["energy"]).exp() for a in initial]
    trlo,trhi,hslo = ball(0),ball(0),ball(0)
    for m in range(len(initial)):
        for q in range(m,len(initial)):
            v = sum((acb(overlap[m][j]*overlap[q][j])*phase_b[j] for j in range(len(final))),acb(0))
            g = (acb(ball("1/2")) if m == q else acb(0))-(phase_a[m]*v+(phase_a[q]*v).conjugate())/4
            err = ((tails[m]*tails[q]).sqrt()/2).upper()
            absolute = lower_nonnegative(abs(g)-err)
            hslo += weights[m]*weights[q]*(1 if m == q else 2)*absolute**2
            if m == q:
                trlo += weights[m]*lower_nonnegative(g.real-err)
                trhi += weights[m]*min(ball(1),(g.real+err).upper())
    trlo,trhi = trlo.lower(),(trhi+tail).upper()
    strict_positive(trlo,"positive comparator yield")
    return {"input_purity_upper": str(pin),"ideal_yield_lower": str((trlo/zhi).lower()),
            "ideal_yield_upper": str((trhi/zlo).upper()),
            "ideal_output_purity_lower": str((hslo/trhi**2).lower()),
            "partition_lower": str(zlo),"partition_upper": str(zhi),
            "omitted_initial_probability_upper": str((tail/zlo).upper()),
            "weighted_final_Gram_error_upper": str((sum((w*t for w,t in zip(weights,tails)),ball(0))/(2*zlo)).upper()),
            "overlap_tail_bounds": [str(x) for x in tails],
            "unnormalized_trace_lower": str(trlo),"unnormalized_trace_upper": str(trhi),
            "unnormalized_HS_squared_lower": str(hslo.lower())}


def extend_to_actual_erf(comparator):
    epsilon = ((1-ERF5)/2).upper()
    tau, dt, nu = 2*arb.pi()/3,ball("1/10000"),ball("9999/10000")
    # Bounded uniform perturbation of Gibbs square roots by Schatten Holder.
    preparation = (S*epsilon*(S*epsilon/2).exp()).upper()
    propagation = (ball("5/2")*tau*epsilon).upper()
    h = (preparation+propagation).upper()
    p0lo,p0hi = arb(comparator["ideal_yield_lower"]).lower(),arb(comparator["ideal_yield_upper"]).upper()
    pure0 = arb(comparator["ideal_output_purity_lower"]).lower()
    strict_positive(p0lo.sqrt()-h,"actual analytic wall yield")
    plo,phi = ((p0lo.sqrt()-h)**2).lower(),((p0hi.sqrt()+h)**2).upper()
    pure = (pure0-2*h/p0lo.sqrt()).lower()
    input_purity = (arb(comparator["input_purity_upper"])+2*preparation).upper()
    # Timing is paid after passing to the actual analytic F, in its own frame.
    timing_error = 3*dt/2
    timing_plo = ((plo.sqrt()-timing_error)**2).lower()
    timing_phi = ((phi.sqrt()+timing_error)**2).upper()
    timing_pure = (pure-2*timing_error/plo.sqrt()).lower()
    noisy_lo = (nu*timing_plo+(1-nu)/2).lower()
    noisy_hi = (nu*timing_phi+(1-nu)/2).upper()
    noisy_pure = ((nu*timing_plo/(nu*timing_plo+(1-nu)/2))**2*timing_pure).lower()
    strict_positive(ball("1/2")-input_purity,"predeclared actual input purity below .5")
    strict_positive(pure-ball("9/10"),"predeclared analytic wall ideal purity above .9")
    strict_positive(noisy_pure-input_purity,"predeclared actual noisy purity improvement")
    thresholds = {}
    for target in ("0.8","0.9","0.95"):
        if timing_pure > ball(target):
            alpha = (ball(target)/timing_pure).sqrt()
            required = (alpha/(alpha+2*timing_plo*(1-alpha))).upper()
            thresholds[target] = str(required)
    alpha = (input_purity/timing_pure).sqrt()
    advantage_visibility = (alpha/(alpha+2*timing_plo*(1-alpha))).upper()
    # Actual 10--90 width follows from erf(z)=.8. Bracket and rational checks.
    zlo,zhi = ball("0.9061938"),ball("0.9061939")
    strict_positive(ball("4/5")-zlo.erf(),"10--90 erf root lower")
    strict_positive(zhi.erf()-ball("4/5"),"10--90 erf root upper")
    return {"uniform_profile_error_upper": str(epsilon),
        "Gibbs_square_root_HS_distance_upper": str(preparation),
        "profile_propagation_operator_error_upper": str(propagation),
        "actual_erf_input_purity_upper": str(input_purity),
        "actual_erf_ideal_yield_lower": str(plo),"actual_erf_ideal_yield_upper": str(phi),
        "actual_erf_ideal_purity_lower": str(pure),
        "timing_ideal_yield_lower": str(timing_plo),"timing_ideal_yield_upper": str(timing_phi),
        "timing_ideal_purity_lower": str(timing_pure),
        "noisy_yield_lower": str(noisy_lo),"noisy_yield_upper": str(noisy_hi),
        "noisy_purity_lower": str(noisy_pure),"mean_preparations_upper": str((1/noisy_lo).upper()),
        "sufficient_visibility_for_purity": thresholds,
        "sufficient_visibility_for_strict_purity_improvement": str(advantage_visibility),
        "10_to_90_width_over_r_lower": "0.9061938", "10_to_90_width_over_r_upper": "0.9061939",
        "timing_error_dimensionless": "1/10000", "visibility": "9999/10000"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path(__file__).with_name("HOT_ERF_SMOOTH_CERTIFICATE.json"))
    args = parser.parse_args()
    start = time.monotonic()
    base = Path(__file__).parent
    print("Fixed analytic erf target: r=.05, s=8, erf scale=.025; beginning rigorous certificate.",flush=True)
    initial_all = complete_prefix(INITIAL_BARRIER,80,base/"HOT_ERF_FLOAT_PROPOSALS_A.json")
    final_all = complete_prefix(FINAL_BARRIER,204,base/"HOT_ERF_FLOAT_PROPOSALS_B.json")
    initial = initial_all[:INPUT_COUNT]
    if len(initial_all) <= INPUT_COUNT:
        raise ArithmeticError("Missing omitted initial neighbor")
    final=[]
    for mode in final_all:
        if mode["energy"] < E_CUT:
            final.append(mode)
        elif mode["energy"] > E_CUT:
            neighbor=mode
            break
        else:
            raise ArithmeticError("Final cutoff crosses mode enclosure")
    else:
        raise ArithmeticError("Missing omitted final neighbor")
    print(f"Certified contiguous prefixes A={len(initial)}, B={len(final)}; forming signed full overlaps.",flush=True)
    comparator = evaluate(initial,final,neighbor,S)
    actual = extend_to_actual_erf(comparator)
    sources = [Path(__file__),
        ROOT/"check_hot_ramsey_finite_certificate.py",
        ROOT/"check_phy_info_finite_purity_certificate.py",
        ROOT/"check_phy_info_continuum_certificate.py"]
    result={"status":"validated analytic erf wall; full Gibbs, continuum tails, visibility, deterministic timing",
        "parameters":{"r":str(R),"s":str(S),"erf_scale":str(WIDTH),"left_collar_boundary":str(LEFT),
                      "collar_cell_width":str(CELL),"cells":NCELL,"Taylor_order":ORDER,
                      "Cauchy_radius":CAUCHY_RADIUS,"initial_modes":len(initial),"final_modes":len(final)},
        "Cauchy_profile_bound":str(F_COMPLEX_BOUND),"comparator":comparator,"actual_analytic_erf":actual,
        "initial_modes":[serialize_mode(m) for m in initial],
        "final_modes":[serialize_mode(m) for m in final],
        "first_omitted_initial_mode":serialize_mode(initial_all[INPUT_COUNT]),
        "first_omitted_final_mode":serialize_mode(neighbor),
        "source_sha256":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        "elapsed_seconds":time.monotonic()-start,
        "scope":"Actual analytic erf profile on the full Dirichlet box, its own full Gibbs state, and actual unprojected output. Timing claim is each deterministic time, not unrecorded random jitter. Visibility is the declared real positive channel; no apparatus realization or generic noise assertion."}
    args.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:result[k] for k in ("status","parameters","actual_analytic_erf","elapsed_seconds")},indent=2))


if __name__ == "__main__":
    main()
