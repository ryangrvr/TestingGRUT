#!/usr/bin/env python3
"""Collect track-xcharge numbers into results/summary.json (no new physics; arithmetic on saved results)."""
import json, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parents[1]
RES = HERE / 'results'
T_REF = 13.8e9 * 365.25 * 86400
MRED = 2.435e18
HBAR = 6.582119569e-25
MX = 0.0924089637375362
LAMX = 0.001

forced = json.load(open(RES / 'forced_Q200000.json'))
pairs = json.load(open(RES / 'pairs_Q200000.json'))
scan = json.load(open(RES / 'scan_ground_family.json'))
td = {d: json.load(open(RES / f'timedomain_d{d}.json')) for d in (5, 6)}

out = dict(codex=dict(frac5=1.5359968620859373e-12, frac6=3.91958726382564e-50, cmax5=1.222681409729816e-4,
                      M5=1.99e22, M6=8.8e10, age_frac5=668918.769, age_frac6=1.7069601851e-32))
best = {r['d']: r for r in forced['rows'] if r.get('R') == 55.0 and not r.get('free')}
coh = {}
for d in (5, 6):
    r = best[d]
    W = np.array([(d + 1), (d - 1), d], float)
    coh[d] = dict(frac=r['frac_rate_per_s'], age=r['frac_times_age'], cmax=r['c_max_1pct'], M=r['M_unit_coeff_1pct'],
                  flux=r['flux'], channel_core_share=(W * np.array(r['flux']) / np.dot(W, r['flux'])).tolist(),
                  drift_coeff_dimless=8 * np.pi * float(np.dot(W, r['flux'])))
out['coherent_single_quantum'] = coh
out['ratio_d5_d6_coherent'] = coh[5]['frac'] / coh[6]['frac']
X0 = 0.924089637375362 * 2.619923134249244 / np.sqrt(2)
out['central_|X0|_GeV'] = X0
out['ratio_over_(M/|X0|)^2_coherent'] = out['ratio_d5_d6_coherent'] / (MRED / X0)**2

# time-domain check
tdc = []
for d in (5, 6):
    epsfac = d * 2**((1 - d) / 2)
    for r in td[d]['runs']:
        eps = epsfac * r['eta']
        pred = -coh[d]['drift_coeff_dimless'] * eps**2
        tdc.append(dict(d=d, eta=r['eta'], eps=eps, even_drift=r['even_drift'], predicted_drift=pred,
                        rel_dev=r['even_drift'] / pred - 1, osc=r['odd_osc_amp'], osc_pred=r['predicted_first_order_amp'],
                        osc_rel_dev=r['odd_osc_amp'] / r['predicted_first_order_amp'] - 1, odd_drift=r['odd_drift']))
out['time_domain'] = tdc

# pairs and totals
tot = {}
for d in (5, 6):
    p = pairs[str(d)]
    variants = {k: p['frac_core_rate_' + k] for k in ('free_born', 'lda_vacuum', 'lda_local_mass')}
    t = coh[d]['frac'] + variants['free_born']
    tot[d] = dict(pair_frac=variants, pair_over_coherent=variants['free_born'] / coh[d]['frac'],
                  total_frac=t, total_age=t * T_REF,
                  cmax_total=float(np.sqrt(0.01 / (t * T_REF))),
                  M_unit_total=float(MRED * (t * T_REF / 0.01)**(1 / (2 * (d - 4)))),
                  cmax_total_poulin=float(np.sqrt(0.087 / (t * T_REF))))
    tot[d]['S_needed'] = float(-np.log(tot[d]['cmax_total']))
out['with_pairs'] = tot
out['S_needed_codex_coherent_only'] = float(-np.log(coh[5]['cmax']))
out['ratio_d5_d6_pairs'] = tot[5]['pair_frac']['free_born'] / tot[6]['pair_frac']['free_born']
out['ratio_over_(M/|X0|)^2_pairs'] = out['ratio_d5_d6_pairs'] / (MRED / X0)**2

# scan
rows = scan['rows']
out['scan'] = [dict(Omega=r['Omega'], Q=r['Q'], F0=r['F0'], E_over_QmX=r['E_over_QmX'],
                    d5_single=r['d5_single_frac'], d5_pair=r['d5_pair_frac'],
                    d6_single=r['d6_single_frac'], d6_pair=r['d6_pair_frac']) for r in rows]

# formation: homogeneous supplied condensate at T = 20 MeV
YQ = 4.8294009076988175e-9
T = 0.020; gs = 10.75; gstar = 10.75; MPL = 1.220890e19
s = 2 * np.pi**2 / 45 * gs * T**3
n = YQ * s
X2 = n / (2 * MX)
H = 1.66 * np.sqrt(gstar) * T**2 / MPL
form = {}
from math import comb
for d in (5, 6):
    lam = 1 / MRED**(d - 4)
    beta = np.sqrt(1 - 4 / (d - 2)**2)
    g2 = comb(d, 2)**2 * lam**2 * X2**(d - 2)
    rate_per_vol = g2 * beta / (4 * np.pi)
    frac = (d - 2) * rate_per_vol / n / HBAR
    form[d] = dict(frac_rate_per_s=frac, Gamma_over_H=frac * HBAR / H)
out['formation_homogeneous_20MeV'] = dict(n_GeV3=n, X_amp_GeV=float(np.sqrt(X2)), H_per_s=H / HBAR, per_d=form,
                                          note='coherent single-quantum emission forbidden for a homogeneous source (k = 0)')

json.dump(out, open(RES / 'summary.json', 'w'), indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
