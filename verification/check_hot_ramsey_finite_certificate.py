"""Full-space continuum Ramsey certificate, with certified spectral tails.

Inherited exact Sturm root certification is imported without modification.
All overlap integrals and certificate arithmetic use Arb/Acb at 256 bits.
The --diagnostic option uses ordinary arithmetic AFTER certified root/overlap
construction and produces no validated claim.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
import json
import math
from pathlib import Path

from check_phy_info_finite_purity_certificate import (  # noqa: E402
    ball, certified_mode, certified_super_modes, ctx, serialize_mode,
    strict_positive,
)
from flint import arb, acb  # noqa: E402

ctx.prec = 256


def sinhc(value):
    return acb(0, value).sinc().real


def signed_full_overlap(initial, final):
    """Direct integral, admitting either input kind and every energy ordering."""
    left = ((initial["k"] - final["k"]).sinc()
            - (initial["k"] + final["k"]).sinc()) / 2
    ikind, fkind = initial["kind"], final["kind"]
    if ikind == fkind == "subbarrier":
        a, b = initial["kap"], final["kap"]
        right = (sinhc(a + b) - sinhc(a - b)) / (2 * a.sinh() * b.sinh())
    elif ikind == fkind == "superbarrier":
        a, b = initial["q"], final["q"]
        right = ((a - b).sinc() - (a + b).sinc()) / (2 * a.sin() * b.sin())
    else:
        sub, sup = (initial, final) if ikind == "subbarrier" else (final, initial)
        a, b = sub["kap"], sup["q"]
        right = (a * a.coth() - b * b.cot()) / (a * a + b * b)
    return (left + initial["k"].sin() * final["k"].sin() * right) / (
        initial["norm"] * final["norm"]
    ).sqrt()


def complete_prefix(barrier, end):
    """Exact, disjoint brackets and contiguous nodal labels certify completeness."""
    last_sub = math.floor(math.sqrt(barrier) / math.pi + 0.5)
    modes = [certified_mode(n, barrier) for n in range(1, last_sub + 1)]
    modes += certified_super_modes(barrier, end)
    labels = [mode["n"] for mode in modes]
    if labels != list(range(1, len(modes) + 1)):
        raise ArithmeticError(f"Missing/duplicate Sturm index: {labels}")
    for previous, current in zip(modes, modes[1:]):
        if not previous["hi"] < current["lo"]:
            raise ArithmeticError("Unordered/overlapping spectral brackets")
    return modes


def free_gibbs_tail(theta, retained):
    """A >= free Dirichlet length-two box, including every omitted input."""
    coefficient = theta * arb.pi() ** 2 / 4
    return (-coefficient * (retained + 1) ** 2).exp() / (
        1 - (-coefficient * (2 * retained + 3)).exp()
    )


def lower_nonnegative(value):
    lower = value.lower()
    return lower if lower > 0 else ball(0)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--s", default="12", help="exact rational scaled inverse temperature")
    parser.add_argument("--initial-count", type=int, default=40)
    parser.add_argument("--final-scaled-cut", default="100")
    parser.add_argument("--diagnostic", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    initial_barrier, final_barrier = 2500, 10000
    retained = args.initial_count
    s = ball(args.s)
    theta = s / initial_barrier
    duration = 2 * arb.pi() / (3 * initial_barrier)
    final_cut = ball(args.final_scaled_cut) * initial_barrier

    all_initial = complete_prefix(initial_barrier, 50 + (retained + 3) * math.pi / 2)
    if len(all_initial) <= retained:
        raise ArithmeticError("Initial prefix does not include omitted neighbor")
    initial, input_neighbor = all_initial[:retained], all_initial[retained]
    all_final = complete_prefix(final_barrier, math.sqrt(float(final_cut)) + 4)
    final = []
    for mode in all_final:
        if mode["energy"] < final_cut:
            final.append(mode)
        elif mode["energy"] > final_cut:
            final_neighbor = mode
            break
        else:
            raise ArithmeticError("Requested final cutoff straddles a root")
    else:
        raise ArithmeticError("Missing first final omitted neighbor")
    if not final or final[-1]["n"] + 1 != final_neighbor["n"]:
        raise ArithmeticError("Final prefix/neighbor not contiguous")

    overlaps = [[signed_full_overlap(i, f) for f in final] for i in initial]
    omitted_lambda = final_neighbor["energy"].lower()
    tails = []
    spectral_tails = []
    for mode, row in zip(initial, overlaps):
        gap = omitted_lambda - mode["energy"]
        strict_positive(gap, "Final tail above every retained input energy")
        spectral = ((final_barrier - initial_barrier) ** 2 * mode["source"] / gap**2).upper()
        residual = (1 - sum((x * x for x in row), ball(0))).upper()
        if not residual >= 0:
            raise ArithmeticError("Overlap prefix norm exceeds one")
        tails.append(min(spectral, residual))
        spectral_tails.append(spectral)

    weights = [(-theta * mode["energy"]).exp() for mode in initial]
    tail = free_gibbs_tail(theta, retained)
    tail2 = free_gibbs_tail(2 * theta, retained)
    zpartial = sum(weights, ball(0))
    zlo, zhi = zpartial.lower(), (zpartial + tail).upper()
    input_purity_lo = (sum((w * w for w in weights), ball(0)) / zhi**2).lower()
    input_purity_hi = ((sum((w * w for w in weights), ball(0)) + tail2) / zlo**2).upper()

    if args.diagnostic:
        import numpy as np
        matrix = np.array([[float(x) for x in row] for row in overlaps])
        ein = np.array([float(m["energy"]) for m in initial])
        efinal = np.array([float(m["energy"]) for m in final])
        unitaries = (matrix * np.exp(-1j * float(duration) * efinal)) @ matrix.T
        phase = np.exp(1j * float(duration) * ein)
        gram = 0.5 * np.eye(retained) - 0.25 * (
            phase[:, None] * unitaries + np.conj(phase[None, :] * unitaries)
        )
        w = np.array([float(x) for x in weights])
        trace = np.dot(w, np.diag(gram).real)
        purity = np.sum(w[:, None] * w[None, :] * np.abs(gram) ** 2) / trace**2
        print(json.dumps({
            "status": "NONCERTIFIED diagnostic; finite Gram central values only",
            "s": args.s, "initial_count": retained, "final_count": len(final),
            "yield_central": trace / float(zpartial), "purity_central": purity,
            "input_purity_upper_certified": str(input_purity_hi),
            "omitted_input_probability_upper": str((tail / zlo).upper()),
            "final_gram_trace_error_upper": str((sum((w * t for w, t in zip(weights, tails)), ball(0)) / 2).upper()),
            "final_overlap_tail_max": str(max(tails)),
        }, indent=2))
        return

    final_phase = [acb(0, -duration * f["energy"]).exp() for f in final]
    initial_phase = [acb(0, duration * i["energy"]).exp() for i in initial]
    trace_lo = ball(0)
    trace_hi = ball(0)
    numerator_lo = ball(0)
    gram_diagonal = []
    max_gram_error = ball(0)
    for m in range(retained):
        for n in range(m, retained):
            unitary = sum((acb(overlaps[m][j] * overlaps[n][j]) * final_phase[j]
                           for j in range(len(final))), acb(0))
            finite = (acb(ball("1/2")) if m == n else acb(0)) - (
                initial_phase[m] * unitary + (initial_phase[n] * unitary).conjugate()
            ) / 4
            error = ((tails[m] * tails[n]).sqrt() / 2).upper()
            max_gram_error = max(max_gram_error, error)
            absolute_lower = lower_nonnegative(abs(finite) - error)
            coefficient = weights[m] * weights[n] * (1 if m == n else 2)
            numerator_lo += coefficient * absolute_lower**2
            if m == n:
                g_lo = lower_nonnegative(finite.real - error)
                g_hi = min(ball(1), (finite.real + error).upper())
                strict_positive(g_hi - g_lo, "ordered diagonal Gram enclosure")
                trace_lo += weights[m] * g_lo
                trace_hi += weights[m] * g_hi
                gram_diagonal.append({"lower": str(g_lo), "upper": str(g_hi),
                                      "final_tail_error_upper": str(error)})

    trace_lo = trace_lo.lower()
    trace_hi = (trace_hi + tail).upper()
    strict_positive(trace_lo, "positive complete dark-port yield")
    purity_lo = (numerator_lo / trace_hi**2).lower()
    yield_lo, yield_hi = (trace_lo / zhi).lower(), (trace_hi / zlo).upper()
    strict_positive(ball(1) - purity_lo, "nonvacuous purity lower bound below one")
    strict_positive(purity_lo - input_purity_hi, "strict full-state purity improvement")

    rational_claims = {}
    if Fraction(args.s) == 12 and retained == 40 and Fraction(args.final_scaled_cut) == 100:
        rational_claims = {
            "input_purity_upper": "182445/1000000",
            "output_purity_lower": "98819/100000",
            "dark_port_yield_lower": "45025/100000000",
            "dark_port_yield_upper": "45027/100000000",
            "omitted_initial_probability_upper": "1/1000000000",
        }
        strict_positive(ball(rational_claims["input_purity_upper"]) - input_purity_hi,
                        "complete input purity < 0.182445")
        strict_positive(purity_lo - ball(rational_claims["output_purity_lower"]),
                        "complete output purity > 0.98819")
        strict_positive(yield_lo - ball(rational_claims["dark_port_yield_lower"]),
                        "complete dark-port yield > 0.00045025")
        strict_positive(ball(rational_claims["dark_port_yield_upper"]) - yield_hi,
                        "complete dark-port yield < 0.00045027")
        strict_positive(ball(rational_claims["omitted_initial_probability_upper"]) - tail / zlo,
                        "complete omitted input probability < 1e-9")
    if Fraction(args.s) == 8 and retained == 48 and Fraction(args.final_scaled_cut) == 100:
        rational_claims = {
            "input_purity_upper": "146928/1000000",
            "output_purity_lower": "96852/100000",
            "dark_port_yield_lower": "72122/100000000",
            "dark_port_yield_upper": "72125/100000000",
            "omitted_initial_probability_upper": "24/10000000000",
        }
        strict_positive(ball(rational_claims["input_purity_upper"]) - input_purity_hi,
                        "complete input purity < 0.146928")
        strict_positive(purity_lo - ball(rational_claims["output_purity_lower"]),
                        "complete output purity > 0.96852")
        strict_positive(yield_lo - ball(rational_claims["dark_port_yield_lower"]),
                        "complete dark-port yield > 0.00072122")
        strict_positive(ball(rational_claims["dark_port_yield_upper"]) - yield_hi,
                        "complete dark-port yield < 0.00072125")
        strict_positive(ball(rational_claims["omitted_initial_probability_upper"]) - tail / zlo,
                        "complete omitted input probability < 2.4e-9")

    result = {
        "status": "validated full-space continuum Gibbs Ramsey purity improvement",
        "arithmetic": "FLINT Arb/Acb, 256 bits; exact rational inputs and certified signed Sturm brackets",
        "parameters": {"r": "1/50", "s": args.s, "theta": str(Fraction(args.s) / 2500),
                       "scaled_duration": "2*pi/3", "physical_dimensionless_duration": "2*pi/7500",
                       "initial_barrier": initial_barrier, "final_barrier": final_barrier,
                       "initial_retained_count": retained, "final_retained_count": len(final),
                       "requested_final_scaled_cut": args.final_scaled_cut},
        "input_purity_lower": str(input_purity_lo), "input_purity_upper": str(input_purity_hi),
        "dark_port_yield_lower": str(yield_lo), "dark_port_yield_upper": str(yield_hi),
        "output_purity_lower": str(purity_lo),
        "unnormalized_output_trace_lower": str(trace_lo), "unnormalized_output_trace_upper": str(trace_hi),
        "unnormalized_output_hilbert_schmidt_squared_lower": str(numerator_lo.lower()),
        "partition_lower": str(zlo), "partition_upper": str(zhi),
        "omitted_initial_unnormalized_gibbs_bound": str(tail),
        "omitted_initial_probability_upper": str((tail / zlo).upper()),
        "first_omitted_initial_mode": serialize_mode(input_neighbor),
        "first_omitted_final_mode": serialize_mode(final_neighbor),
        "final_tail_lambda_lower": str(omitted_lambda),
        "max_gram_tail_error_upper": str(max_gram_error),
        "certified_rational_claims": rational_claims,
        "retained_diagonal_gram_enclosures": gram_diagonal,
        "final_overlap_tail_bounds": [str(x) for x in tails],
        "separate_spectral_tail_bounds": [str(x) for x in spectral_tails],
        "initial_modes": [serialize_mode(m) for m in initial],
        "final_modes": [serialize_mode(m) for m in final],
        "scope": "One actual K=(exp(-i tau B)-exp(-i tau A))/2 on full Hilbert space; full initial Gibbs state; no energy-band projector or spatial discretization; ideal coherent ancilla/relative phase only.",
    }
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "parameters", "input_purity_lower", "input_purity_upper",
        "dark_port_yield_lower", "dark_port_yield_upper", "output_purity_lower",
        "omitted_initial_probability_upper", "max_gram_tail_error_upper",
    )}, indent=2))


if __name__ == "__main__":
    main()
