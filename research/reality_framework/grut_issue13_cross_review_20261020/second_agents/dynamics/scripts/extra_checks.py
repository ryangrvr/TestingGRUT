"""Extra checks: wall proximity, injected perturbation energy, endpoint components, merger phase.
Usage: python -I extra_checks.py <raw_dir>"""
import sys, numpy as np
R = sys.argv[1]
L = lambda n: np.genfromtxt(f'{R}/{n}.csv', delimiter=',', names=True, ndmin=1)
print('-- closest approach of an R12 ball to a fixed-vacuum face, late window (pairs with transverse drift)')
for n, Lx, Lz in [('pair_negative_fine', 24, 48), ('pair_generic_fine', 24, 48), ('pair_generic_medium', 24, 48), ('pair_generic_coarse', 24, 48), ('pair_rotated_fine', 24, 48)]:
    a = L(n); m = a['t'] >= 120
    for j in (0, 1):
        x, y, z = a[f'c{j}_x'][m], a[f'c{j}_y'][m], a[f'c{j}_z'][m]
        g = np.stack([Lx-abs(x), Lx-abs(y), Lz-abs(z)], 1)-12
        k = np.unravel_index(np.argmin(g), g.shape)
        print('  %-20s core%d end centre=(%.2f,%.2f,%.2f) min gap=%.2f along %s at t=%d; late centre velocity=(%.4f,%.4f,%.4f)' % (
            n, j, x[-1], y[-1], z[-1], g.min(), 'xyz'[k[1]], a['t'][m][k[0]], *(np.polyfit(a['t'][m], v, 1)[0] for v in (x, y, z))))
print('-- injected perturbation: initial E/Q of perturbed vs reference')
for p, r in [('single_generic_fine', 'single_reference_fine'), ('single_negative_fine', 'single_reference_fine'), ('pair_generic_fine', 'pair_reference_fine'), ('pair_negative_fine', 'pair_reference_fine')]:
    a, b = L(p), L(r)
    print('  %-22s E0=%.6f Q0=%.6f  dE/E_ref=%.4f%%  dQ/Q_ref=%.5f%%  E/Q=%.6f vs %.6f' % (p, a['E'][0], a['Q'][0], 100*(a['E'][0]/b['E'][0]-1), 100*(a['Q'][0]/b['Q'][0]-1), a['E'][0]/a['Q'][0], b['E'][0]/b['Q'][0]))
print('-- endpoint full-grid components (single runs), 2% threshold: count, top-2 Q / Q0')
for n in ['single_reference_coarse', 'single_generic_coarse', 'single_reference_medium', 'single_generic_medium', 'single_reference_fine', 'single_generic_fine', 'single_negative_fine', 'single_generic_box', 'single_generic_dt']:
    a = L(n); c = np.genfromtxt(f'{R}/{n}.csv.components.csv', delimiter=',', names=True, ndmin=1)
    s = c[np.isclose(c['threshold'], .02)]; q = np.sort(s['Q'])[::-1]
    print('  %-24s n=%d top=%s' % (n, len(s), np.round(q[:2]/a['Q'][0], 6)))
print('-- merger phase in pair_reference_fine: centre separation and Q12 per core')
a = L('pair_reference_fine')
for r in a:
    if 25 <= r['t'] <= 80:
        sep = np.sqrt((r['c0_x']-r['c1_x'])**2+(r['c0_y']-r['c1_y'])**2+(r['c0_z']-r['c1_z'])**2)
        print('  t=%3d sep=%6.2f Q12=(%.2f, %.2f) Q8=(%.2f, %.2f) maxF2=%.3f minS=%.3f' % (r['t'], sep, r['c0_Q12'], r['c1_Q12'], r['c0_Q8'], r['c1_Q8'], r['maxF2'], r['minS']))
