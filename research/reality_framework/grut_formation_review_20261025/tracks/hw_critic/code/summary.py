#!/usr/bin/env python3
"""hw_critic (Claude): print every number quoted in NOTES.md / CONSOLIDATED.md from results/*.json.
Usage: python3 -B summary.py ../results"""
import json, math, sys
from pathlib import Path
R = Path(sys.argv[1])
ck = json.load(open(R / 'checks.json')); th = json.load(open(R / 'thresholds.json'))
ho = json.load(open(R / 'holder.json')); we = json.load(open(R / 'window_edges.json'))
rows = {r['tag']: r for r in ck['rows']}
s = rows['spot_relaxed']
print('SPOT (own, Y-relaxed): eps %.5f  NX %.2f  EX %.1f  NYH %.5f  NYL %.4f  Nsig %.1e  R_Y %.3f  R_casc %.3f  R %.4f  yield %.4f  R_L(table) %.2f'
      % (s['eps'], s['NX'], s['EX'], s['NYH'], s['NYL'], s['Nsig'], s['R_Y'], s['R_casc'], s['R_central'], s['yield_over_YQ'], s['R_L_table']))
c = rows['spot_codex']
print('SPOT (own, Codex prep): eps %.5f  NYH %.2f  R %.0f' % (c['eps'], c['NYH'], c['R_central']))
a = rows['spot_a100_None']; print('window control a_end 100: eps %.5f R %.4f' % (a['eps'], a['R_central']))
print('\nHELD relaxation from Codex state (N, eps, NY, R):')
for N in (0.5, 1.0, 2.0):
    r = rows['codex_held_%s' % N]; print('  N %.1f  eps %+.4f  NY %7.3f  R %7.2f' % (N, r['eps'], r['NY'], r['R_central']))
print('\nFREE (unheld) de Sitter stage from the relaxed state, then release:')
print('  N     delta_phys(end)  nX_phys(end)  eps(with velocities)  eps(velocities zeroed)  R(with)  R(zeroed)')
for N in (0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0):
    v = rows['free_vel_%s' % N]; z = rows['free_novel_%s' % N]
    d = (v['stage_end']['delta'] - math.pi) % (2 * math.pi); d = d - 2 * math.pi if d > math.pi else d   # S = -s e^{i delta}: delta_phys = arg S - 6 arg X - pi
    print('  %-5g %10.3f %14.0f %14.4f %18.4f %12.3f %8.3f' % (N, d, v['stage_end']['nX_phys'], v['eps'], z['eps'], v['R_central'], z['R_central']))
print('\nFINITE RELEASE (Gamma_phi/H_i, eps, NY, R_central, R_reheat, diluted yield/Y_Q):')
for g in (10.0, 1.0, 0.1):
    r = rows['gamma_%s' % g]
    print('  %5g  %+.4f  %.4g  %.4f  %.4f  %.3f' % (g, r['eps'], r['NY'], r['R_central'], r['R_central_reheat'], r['yield_diluted_over_YQ']))
print('\nD1 HOLDER switch-off (own): ' + '; '.join('%s/%s eps %+.4f R %.3g' % (r['mode'], r['Dl'], r['eps'], r['R_central'])
                                               for r in sorted(ho['rows'], key=lambda r: r['tag'])))
print('\nTHRESHOLDS: ' + ', '.join('%s %.2f' % (k, v) for k, v in th['thresholds'].items()))
print('GMR formula N_A = %.4f H_i^3 (R_A %.1f); headroom per net charge: %s' % (
    th['gmr']['N_A_per_Hi3'], th['gmr']['R_A_from_formula'], ', '.join('%s %.2e' % (k, v) for k, v in th['vector_headroom_per_net_charge'].items())))
print('FULL DEBRIS / threshold (vector from own modes | GMR formula):')
for e in (1e-7, 3e-7, 1e-6, 3e-6, 1e-5):
    m = [r for r in th['rows'] if r['e'] == e and r['vec'] == 'modes'][0]; g = [r for r in th['rows'] if r['e'] == e and r['vec'] == 'GMR'][0]
    print('  e %.0e  R_tot %5.1f|%5.1f  proxy(direct 81.5) %.2f|%.2f  ke(8.26) %.2f|%.2f  ke(6.97) %.2f|%.2f  ke+s(19.6) %.2f|%.2f' % (
        e, m['R_total'], g['R_total'], m['ratio']['direct_proxy_hw_verify'], g['ratio']['direct_proxy_hw_verify'],
        m['ratio']['ke_muc0.861'], g['ratio']['ke_muc0.861'], m['ratio']['ke_muc0.923'], g['ratio']['ke_muc0.923'],
        m['ratio']['ke_plus_s'], g['ratio']['ke_plus_s']))
print('\nWINDOW LOWER EDGE m (mediator only | with vector e=1e-7, 1e-6, 3e-6):')
for k, v in we.items():
    f = lambda x: ('%.1f' % x) if isinstance(x, float) else str(x)
    print('  %-24s %s | %s, %s, %s' % (k, f(v['mediator_only']), f(v['with_vector_e=1e-07']), f(v['with_vector_e=1e-06']), f(v['with_vector_e=3e-06'])))
