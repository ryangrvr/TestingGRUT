"""R8 retention sensitivity, object-only energy drift, Cartesian R36 symmetry residuals.
Usage: python -I misc_checks.py <raw_dir> <protocol.json> <TRANSVERSE_PROTOCOL.json>"""
import sys, json, numpy as np
R = sys.argv[1]; P = json.load(open(sys.argv[2])); T = json.load(open(sys.argv[3]))
L = lambda n: np.genfromtxt(f'{R}/{n}.csv', delimiter=',', names=True)
objs = [c for c in P['cases'] if c['case'] in ('single', 'pair')] + T['cases']
e = max((np.max(abs(L(c['id'])['E']/L(c['id'])['E'][0]-1)), c['id']) for c in objs)
print('object runs only: max |dE/E0| = %.3e (%s)' % e)
r8 = []
for c in objs:
    a = L(c['id']); m = a['t'] >= 120; Q0 = a['Q'][0]
    v = (a['c0_Q8'][m]+a['c1_Q8'][m])/Q0 if c['case'] == 'pair' else np.maximum(a['c0_Q8'][m], a['c1_Q8'][m])/Q0
    r8.append((v.min(), c['id']))
print('R8 retention: min over 26 runs = %.4f%% (%s); max = %.4f%%' % (100*min(r8)[0], min(r8)[1], 100*max(r8)[0]))
print('Cartesian R36 symmetry residual |res|/|gen| (translation, rotation) at t=80,120,160:')
for n in ['pair_reference_coarse', 'pair_reference_medium', 'pair_reference_fine', 'pair_reference_box', 'single_generic_coarse', 'single_generic_medium', 'single_generic_fine', 'single_generic_box']:
    a = L(n)
    print('  %-24s' % n, ' '.join('t%d:(%.4f,%.4f)' % (t, a['res0_norm'][a['t'] == t][0]/a['gen0_norm'][a['t'] == t][0], a['res1_norm'][a['t'] == t][0]/a['gen1_norm'][a['t'] == t][0]) for t in (80, 120, 160)))
