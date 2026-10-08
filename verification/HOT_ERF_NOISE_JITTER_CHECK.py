"""Scalar full-Gibbs jitter/readout extension of frozen analytic erf certificate.

No new spectral solve or operating-point search. Independent duration law on
tau +/- .001, visibility .9999, false positive 1e-5, false negative .01.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
from check_phy_info_continuum_certificate import ball,strict_positive
from flint import arb,ctx
ctx.prec=256


def free_tail(s,r,n):
    c=s*r*r*arb.pi()**2/4
    return (-c*(n+1)**2).exp()/(1-(-c*(2*n+3)).exp())


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",type=Path,default=Path(__file__).with_name("HOT_ERF_NOISE_JITTER_CERTIFICATE.json"))
    args=parser.parse_args()
    baseline=Path(__file__).with_name("HOT_ERF_SMOOTH_CERTIFICATE.json")
    data=json.loads(baseline.read_text(encoding="utf-8"))
    for path,digest in data["source_sha256"].items():
        if hashlib.sha256((ROOT/path).read_bytes()).hexdigest()!=digest:
            raise ArithmeticError("Frozen smooth source hash mismatch: "+path)
    if any(data["parameters"][key]!=value for key,value in {"r":"1/20","s":"8","erf_scale":"1/40"}.items()):
        raise ArithmeticError("Unexpected smooth baseline parameters")
    r,s=ball("1/20"),ball(8)
    actual=data["actual_analytic_erf"]
    epsilon=arb(actual["uniform_profile_error_upper"]).upper()
    preparation=arb(actual["Gibbs_square_root_HS_distance_upper"]).upper()
    zlo=arb(data["comparator"]["partition_lower"]).lower()
    actual_zlo=((-s*epsilon).exp()*zlo).lower()
    erf5=ball(5).erf()
    energy2,qmass=ball(0),ball(0)
    for mode in data["initial_modes"]:
        e=arb(mode["energy_ball"])*r*r
        # Actual erf eigenvalues differ from comparator by at most epsilon.
        energy2 += (e+epsilon)**2*(-s*(e-epsilon)).exp()
        norms=[arb(x) for x in mode["cell_norm_balls"]]
        weighted=arb(mode["right_norm_ball"])
        for j,norm in enumerate(norms):
            z=ball(Fraction(-24+j,5))
            fmax=(z.erf()+erf5)/(2*erf5)
            weighted += fmax*fmax*norm
        qmass += (-s*e).exp()*weighted/arb(mode["normalization_ball"])
    n=len(data["initial_modes"])
    tail=free_tail(s,r,n)
    c=s*r*r*arb.pi()**2/4
    x0=r*r*arb.pi()**2*(n+1)**2/4
    strict_positive(x0-2/s,"monotone second-moment free tail")
    q0=(-c*(2*n+3)).exp()
    moment_tail=x0*x0*(-c*(n+1)**2).exp()*(1+11*q0+11*q0*q0+q0**3)/(1-q0)**5
    energy2=((energy2+moment_tail)/actual_zlo).upper()
    # rho actual vs comparator costs <=2 preparation in trace norm;
    # ||F^2-F0^2||<=2epsilon because both multipliers lie in [0,1].
    qmass=((qmass+tail)/zlo+2*preparation+2*epsilon).upper()
    lipschitz=(energy2.sqrt()+ball("3/2")*qmass.sqrt()).upper()
    delta,nu,fp,fn=ball("1/1000"),ball("9999/10000"),ball("1/100000"),ball("1/100")
    p0lo,p0hi=arb(actual["actual_erf_ideal_yield_lower"]).lower(),arb(actual["actual_erf_ideal_yield_upper"]).upper()
    purity0=arb(actual["actual_erf_ideal_purity_lower"]).lower()
    d=(lipschitz*delta).upper()
    strict_positive(p0lo.sqrt()-d,"nonzero probability for every timing law")
    plo,phi=((p0lo.sqrt()-d)**2).lower(),((p0hi.sqrt()+d)**2).upper()
    purity=(purity0-2*d/p0lo.sqrt()).lower()
    coefficient=(1-fn-fp)*nu
    background=(1-fn-fp)*(1-nu)/2+fp
    noisy_lo,noisy_hi=(coefficient*plo+background).lower(),(coefficient*phi+background).upper()
    noisy_purity=((coefficient*plo/(coefficient*plo+background))**2*purity).lower()
    preparations=(1/noisy_lo).upper()
    pnu_lo=(nu*plo+(1-nu)/2).lower()
    false_fraction=(fp*(1-pnu_lo)/noisy_lo).upper()
    sufficient={}
    for target in ("0.8","0.9","0.95"):
        if purity > ball(target):
            sufficient[target]=str(((1+2*fp/(1-fn-fp))/(1+2*plo*((purity/ball(target)).sqrt()-1))).upper())
    input_purity=arb(actual["actual_erf_input_purity_upper"]).upper()
    # A model-independent necessary ceiling for this same device at exact
    # nominal time, nu=.999 and ideal readout: ||D||_2<=sqrt(P_input).
    poor_nu=ball("999/1000")
    alpha_upper=(poor_nu*p0hi/(poor_nu*p0hi+(1-poor_nu)/2)).upper()
    noisy_ceiling=(alpha_upper+(1-alpha_upper)*input_purity.sqrt())**2
    strict_positive(ball("881/1000")-noisy_ceiling,"nu=.999 cannot reach .90 purity at nominal time")
    # Strict purification is a less demanding target than .90 purity.
    moderate_nu=ball("998/1000")
    ma=(1-fn-fp)*moderate_nu
    mb=(1-fn-fp)*(1-moderate_nu)/2+fp
    moderate_purity=((ma*plo/(ma*plo+mb))**2*purity).lower()
    strict_positive(moderate_purity-input_purity,"strict smooth-wall purity gain at nu=.998")
    # Fixed threshold before the scalar run, not adjusted after seeing output.
    strict_positive(noisy_purity-ball("9/10"),"predeclared analytic-wall purity >.90 with jitter and readout")
    strict_positive(ball(400)-preparations,"predeclared fewer than 400 preparations")
    result={"status":"validated actual analytic erf wall with independent arbitrary bounded duration law and specified detector/dephasing channel",
        "parameters":{"r":"1/20","s":"8","erf_scale":"1/40","timing_halfwidth":"1/1000",
                      "visibility":"9999/10000","false_positive":"1/100000","false_negative":"1/100"},
        "actual_second_energy_moment_upper":str(energy2),"actual_F_squared_moment_upper":str(qmass),
        "laboratory_frame_HS_Lipschitz_upper":str(lipschitz),"full_second_energy_tail_unnormalized":str(moment_tail),
        "jittered_ideal_yield_lower":str(plo),"jittered_ideal_yield_upper":str(phi),
        "jittered_ideal_purity_lower":str(purity),"recorded_yield_lower":str(noisy_lo),
        "recorded_yield_upper":str(noisy_hi),"recorded_purity_lower":str(noisy_purity),
        "mean_preparations_upper":str(preparations),"detector_false_accept_fraction_upper":str(false_fraction),
        "sufficient_visibility_for_purity_with_readout":sufficient,
        "nominal_time_purity_upper_at_visibility_0999_ideal_readout":str(noisy_ceiling.upper()),
        "purity_lower_at_visibility_0998_with_same_jitter_and_readout":str(moderate_purity),
        "input_purity_upper":actual["actual_erf_input_purity_upper"],
        "baseline_sha256":hashlib.sha256(baseline.read_bytes()).hexdigest(),
        "baseline_source_sha256":data["source_sha256"],
        "extension_source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "scope":"All initial/final modes retained mathematically. Arbitrary duration law supported in tau +/- .001, independent of prepared state. Visibility .9999 and independent detector false-positive 1e-5/false-negative .01. No apparatus or generic correlated error assertion."}
    args.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
