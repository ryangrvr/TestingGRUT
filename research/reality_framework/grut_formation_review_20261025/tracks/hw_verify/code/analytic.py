#!/usr/bin/env python3
"""hw_verify (Claude, own code): closed-form cross-checks.

1. Reproduction [C of Codex's quoted formula] of Codex 2908e91's conditional envelope bound at the FIXED stated
   couplings (r_i = 0.1 M_Y and M_Y, c_H = 1 and 10) with our own implementation, and the same bound evaluated at the
   ls_scaling heavy spot point and ladder (the heavy points lie inside Codex's "controlled matching" diagnostics).
2. Why the two packets agree: at fixed couplings the quiet-regime yield needs r_i/M_Y ~ 12, i.e. M_Y/H ~ 6
   (resonant, not quiet); the heavy window instead raises M_Y, h, mu together with H (the scaling symmetry) and
   kappa, r.  The yield at the quiet spot point is 24% of Codex's envelope bound.
3. Ladder closed forms from our ladder.json: R per Y quantum, R per vector, the window ends.
Usage: python3 -B analytic.py ../results   (reads ladder.json, writes analytic.json)
"""
import json, sys
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hw

R = Path(sys.argv[1])
out = {}

# 1. Codex's bound at the fixed couplings, own implementation of the quoted formula
cp0 = dict(hw.STATED)
rows = []
for cH in (1.0, 10.0):
    for ratio in (0.1, 1.0):
        r_i = ratio * cp0['M_GeV']
        Hi = np.sqrt(2 * cp0['lx'] / cH) * r_i            # Codex's leading radial relation (x0 = sqrt(cH/(2 lx)) H)
        b = hw.codex_envelope_bound(cp0, Hi, cH, r_i / Hi)
        rows.append(dict(cH=cH, r_i_over_M=ratio, H_i_GeV=Hi, M_over_H=cp0['M_GeV'] / Hi, yield_bound_over_YQ=b['yield_bound_over_YQ'],
                         codex_quoted={(1.0, 0.1): 1.2201e-8, (10.0, 0.1): 5.4890e-8, (1.0, 1.0): 3.8583e-5, (10.0, 1.0): 1.7358e-4}[(cH, ratio)],
                         diagnostics={k: b[k] for k in ('kappa_X2_over_M2', 'X_freq_over_Y_gap', 'X_freq_over_Sigma_mass',
                                                        'Y_quartic_relative_force', 'Sigma_relative_forced_shift')}))
out['codex_bound_fixed_couplings'] = rows
# the r_i that would give Y_Q at fixed couplings, c_H = 10, and the corresponding M_Y/H
from scipy.optimize import brentq
f = lambda lr: np.log(hw.codex_envelope_bound(cp0, np.sqrt(2 * cp0['lx'] / 10) * np.exp(lr) * 1e12, 10.0,
                                             1 / np.sqrt(2 * cp0['lx'] / 10))['yield_bound_over_YQ'])
lr = brentq(f, np.log(0.1), np.log(1e3))
ri = np.exp(lr); Hi = np.sqrt(2 * cp0['lx'] / 10) * ri * 1e12
out['codex_bound_reaches_YQ_at'] = dict(cH=10.0, r_i_over_M=ri, H_i_GeV=Hi, M_over_H=1e12 / Hi, codex_quoted_r_i_over_M=11.87,
                                        note='M_Y/H ~ 6 there: the resonant regime (ls_scaling m = 10 already fails on debris)')

# 2. the bound at the heavy ladder (own ladder.json)
L = json.load(open(R / 'ladder.json'))['rows']
lad = []
for r in sorted([r for r in L if r['name'].startswith('m') and abs(r['delta'] - np.pi / 2) < 1e-6], key=lambda r: r['couplings_over_H']['M']):
    b = r['codex_envelope']
    lad.append(dict(m=r['couplings_over_H']['M'], eps=r['eps'], NX=r['NX'], NX_envelope_bound=b['NX_bound'],
                    NX_over_bound=abs(r['NX']) / b['NX_bound'], yield_over_YQ=r['yield_over_YQ'],
                    yield_bound_over_YQ=b['yield_bound_over_YQ'],
                    diagnostics={k: b[k] for k in ('kappa_X2_over_M2', 'X_freq_over_Y_gap', 'X_freq_over_Sigma_mass',
                                                   'Y_quartic_relative_force', 'Sigma_relative_forced_shift')},
                    R_central=r['R_central'], R_L=r['R_L'], R_per_Y_quantum=r['deposit_per_Y'] / r['EX'],
                    R_per_vector=r['deposit_per_A'] / r['EX'], NY=r['NYH'] + r['NYL'], c6=r['c6'], kappa=r['kappa'],
                    kappa_over_max=r['kappa'] / r['kappa_max']))
for a, b in zip(lad, lad[1:]):
    b['slope_NY'] = -np.log(b['NY'] / a['NY']) / np.log(b['m'] / a['m'])
    b['slope_R'] = -np.log(b['R_central'] / a['R_central']) / np.log(b['m'] / a['m'])
out['heavy_ladder'] = lad
spot = next(r for r in lad if r['m'] == 30)
out['window'] = dict(m_min_debris_between=[10, 15], m_max_stability=hw.kappa_max(hw.heavy()) / (1e-3 * np.sqrt(0.5 * (1 - 0.45**2) / 30)),
                     vector_allowance_per_net_charge=(spot['R_L'] - spot['R_central']) / (spot['R_per_vector'] * abs(spot['NX'])),
                     Y_allowance_per_net_charge=(spot['R_L'] - spot['R_central']) / (spot['R_per_Y_quantum'] * abs(spot['NX'])))
(R / 'analytic.json').write_text(json.dumps(out, indent=1, default=float))
print(json.dumps(out, indent=1, default=float))
