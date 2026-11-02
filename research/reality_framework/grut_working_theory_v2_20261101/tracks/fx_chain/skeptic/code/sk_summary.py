"""fx_chain_skeptic: arithmetic that composes the skeptic numbers (frozen test with the recalibrated coefficients,
eta = 0 case, nucleation work estimate for late re-formation, p versus s/n_Q exponent).  Reads ../results/*.json."""
import json, math, os
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
kin = json.load(open(R + '/sk_kin.json')); se = json.load(open(R + '/start_entropy.json'))
out = {}
# frozen test coefficients for an object family member with E/Q = x m_X: A = 1/x - 1, B = 2/x
def coef(x): return 1 / x - 1, 2 / x
for lab, x, p in (('Q2e5_supplied', 0.93739999, kin['m1_gas_relaxed_R497']['p']),
                  ('Q~1e12_slow_passage', 0.7770, kin['m1_gas_relaxed_R497_YQv6.128e-09']['p'])):
    A, B = coef(x)
    out['frozen_' + lab] = dict(E_over_Q=x, A=A, B=B, p_gas=p, eta_min_gas=1 - (0.01 - B * p) / A, eta_min_capture=1 - 0.01 / A,
                                test_at_eta0=A + B * p, fail_factor_eta0=(A + B * p) / 0.01)
# p versus entropy per charge
S0 = se['S_pre_hartree']; p0 = kin['m1_gas_relaxed_R497']['p']
ex = []
for n, k, f in (('sk_gas_mult1.0_R497_V0x0', 1.0, 0.0), ('sk_gas_mult1.5_R497', 1.5, 1.0), ('sk_gas_mult1.2_R497', 1.2, 1.0)):
    r = [x for x in se['rows'] if abs(x['k'] - k) < 1e-9 and abs(x['f'] - f) < 1e-9][0]
    ex.append(dict(case=n, s_ratio=r['s_per_net'] / S0, p_ratio=kin[n]['p'] / p0, exponent=math.log(kin[n]['p'] / p0) / math.log(r['s_per_net'] / S0)))
ar = se['a_ref_start']; ex.append(dict(case='m1_gas_aref_R497', s_ratio=ar['ratio'], p_ratio=kin['m1_gas_aref_R497']['p'] / p0,
                                       exponent=math.log(kin['m1_gas_aref_R497']['p'] / p0) / math.log(ar['ratio'])))
out['p_vs_entropy'] = ex
out['entropy_for_eta1_limit'] = dict(p_limit=0.01 / 2.13356, s_ratio_needed=(0.01 / 2.13356 / p0) ** (1 / 9.7))
# late re-formation from the gas (all charge in gas, Boltzmann): mu = m - T (s/n - 5/2); work to assemble the smallest
# bound objects W = Q (E/Q - mu) with E/Q(7e4) = 1.00106, Q_min = 5.67e4 (E/Q >~ 1)
sn = S0
rows = []
for T in (0.09, 0.05, 0.03, 0.02, 0.01, 0.005):
    mu = 1 - T * (sn - 2.5)
    W = 7e4 * (1.00106 - mu)
    rows.append(dict(T=T, mu_gas=mu, fugacity=math.exp(-(sn - 2.5)), W_7e4=W, W_over_T=W / T))
out['late_nucleation_estimate'] = rows
json.dump(out, open(R + '/skeptic_summary.json', 'w'), indent=1)
print(json.dumps(out, indent=1))
