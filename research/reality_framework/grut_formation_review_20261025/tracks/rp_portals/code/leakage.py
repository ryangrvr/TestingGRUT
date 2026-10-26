"""rp_portals / leakage.py (Claude, own code; rates reused from our xcharge and hw_gates tracks, stated).

Cold-object X^6 leakage at the heavy spot point with the full heavy-Y denominator:
alpha_6 = kappa^2 mu h_Sigma / (sqrt2 (M_Y^4 - 2 mu^2 h_Sigma^2)), c6 = alpha_6 M_red^2, Gamma t_0 = c6^2 x (unit-c6 rate) x t_0.
The unit-coefficient rates for the Q = 200000 ground object are our own xcharge numbers (coherent channel from
forced.py, two-quantum anti-X pair channel from pairs.py, free-Born and the two LDA variants) [Q own tracks];
hw_gates gate5 is the second own evaluation [Q own].  Codex's PAIR_RESULTS 'included_channels_rate' is quoted [Q].

Usage: python3 -B leakage.py <out.json>
"""
import json, math, sys, pathlib

HI = 1.3334961939833566e11
MRED = 2.435e18
T0 = 13.8e9*365.25*86400

def main(out):
    M, h, mu, kap = 30*HI, 40*HI, 7.159456159513794*HI, 0.0034586847211042523
    r = math.sqrt(2)*mu*h/M**2
    alpha_full = kap**2*mu*h/(math.sqrt(2)*(M**4 - 2*mu*mu*h*h))
    alpha_lead = kap**2*mu*h/(math.sqrt(2)*M**4)
    c6 = alpha_full*MRED**2
    root = pathlib.Path(__file__).resolve().parents[2]
    xc = json.load(open(root/'xcharge/results/pairs_Q200000.json'))['6']
    fo = [row for row in json.load(open(root/'xcharge/results/forced_Q200000.json'))['rows'] if row['d'] == 6][0]
    hw = json.load(open(root/'hw_gates/results/gate5_light.json'))['survival']
    coh = fo['frac_rate_per_s']
    pairs = {k: xc[k] for k in ('frac_core_rate_free_born', 'frac_core_rate_lda_vacuum', 'frac_core_rate_lda_local_mass')}
    res = dict(r=r, exact_over_leading=alpha_full/alpha_lead, one_over_1_minus_r2=1/(1-r*r),
               alpha6_full_GeVm2=alpha_full, alpha6_leading_GeVm2=alpha_lead, c6_full=c6, c6_leading=alpha_lead*MRED**2,
               codex_alpha6=2.1088597931217505e-31, codex_c6=1250390.420687231, codex_Gamma_t0=7.08367857848652e-17,
               hw_gates_alpha6=hw['alpha6_GeVm2'], hw_gates_c6=hw['c6'], hw_gates_Gamma_t0=hw['rate_times_age'],
               xcharge_unit_rates_per_s=dict(coherent=coh, **pairs),
               codex_included_unit_rate=1.0403620880845258e-46,
               Gamma_t0=dict(coherent_plus_pair_free_born=c6*c6*(coh+pairs['frac_core_rate_free_born'])*T0,
                             coherent_plus_pair_lda_vacuum=c6*c6*(coh+pairs['frac_core_rate_lda_vacuum'])*T0,
                             coherent_plus_pair_lda_local_mass=c6*c6*(coh+pairs['frac_core_rate_lda_local_mass'])*T0,
                             leading_denominator_free_born=(alpha_lead*MRED**2)**2*(coh+pairs['frac_core_rate_free_born'])*T0),
               pair_to_coherent=pairs['frac_core_rate_free_born']/coh,
               c6_ceiling_1pct=math.sqrt(0.01/((coh+pairs['frac_core_rate_free_born'])*T0)),
               note='alpha_6 convention: V_eff contains alpha_6 X^6 + h.c. after integrating out both real Y components at Sigma = h/sqrt2; '
                    'the pair channel is a free-Born/LDA estimate (xcharge: O(1) medium corrections not computed).')
    json.dump(res, open(out, 'w'), indent=1)
    print(json.dumps(res, indent=1))

if __name__ == '__main__':
    main(sys.argv[1])
