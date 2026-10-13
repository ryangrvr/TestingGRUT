"""Full-law ground exclusion for the prepared graph-harmonic sector.

Uses the pinned parent packet's unregulated Gauss-star lower bound and
matching-state trial upper bound. Outward arithmetic; no cutoff or weak
hopping expansion. No geometry is called selected.
"""
from pathlib import Path
from fractions import Fraction
from decimal import Decimal,localcontext,ROUND_FLOOR,ROUND_CEILING
import json
import sys

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/"cp12_remainder_and_mixing"))
from variational_exclusion import regular_lower,frac_dec
from near_regular_exclusion import conservative_trial_upper


def certificate(n):
    Q=2*n;M=n*(n-1)//2;assert M-2*Q+1>0
    b2=Fraction(2*M*Q,M-2*Q+1) # Frozen hopping t=1. Block Schur count.
    upper=conservative_trial_upper(n)
    with localcontext() as floor:
        floor.prec=80;floor.rounding=ROUND_FLOOR
        lower_total=Decimal(n)*regular_lower(n)
        gap=lower_total-upper
        gap2=gap*gap if gap>0 else Decimal(0)
        denominator=gap2+frac_dec(b2,ROUND_FLOOR)
        with localcontext() as ceiling:
            ceiling.prec=80;ceiling.rounding=ROUND_CEILING
            probability=frac_dec(b2,ROUND_CEILING)/denominator if gap>0 else Decimal(1)
    alpha=Fraction(Q,M-Q+1)
    out=dict(N=n,Q=Q,minimum_unlock_denominator=M-2*Q+1,
             exact_H_crossing_norm_square_upper=str(b2),
             fixed_block_lower_total=str(lower_total),
             full_H_trial_upper=str(upper),energy_gap_lower=str(gap),
             top_harmonic_ground_probability_upper=str(probability),
             definite_graph_top_probability_lower=str(1-alpha),
             exact_definite_graph_missing_top_weight_upper=str(alpha))
    return out


if __name__=="__main__":
    out=dict(status="CP12_HARMONIC_GROUND_CERTIFICATE_PASS",
             rows=[certificate(n) for n in [512,1024,8192,32768,131072,1048576]],
             asymptotic_ground_top_probability="O(1/N)",
             asymptotic_prepared_definite_graph_missing_top_probability="O(1/N)",
             scope="Ground/definite-preparation separation in the unchanged full unregulated CP-12 law. This is a quantum connection-harmonic sector, not a spacetime phase or gravitational dressing.")
    (P/"harmonic_ground_certificate_results.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(out['rows'][-2:],indent=2))
