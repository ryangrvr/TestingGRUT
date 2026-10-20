"""Cross-code check: old cylindrical m=1 relative-shift response (large-domain replays) vs Cartesian direct-impact x/delta.
Usage: python -I old_m1_vs_cartesian.py <raw_dir>"""
import json, sys, numpy as np
R = sys.argv[1]
H = [0.5, 1/3, 0.25]
def rich(v, h):
    d1, d2 = v[0]-v[1], v[1]-v[2]
    f = lambda p: (h[0]**p-h[1]**p)/(h[1]**p-h[2]**p) - d1/d2
    lo, hi = 0.05, 8
    for _ in range(200):
        m = .5*(lo+hi)
        if f(lo)*f(m) <= 0: hi = m
        else: lo = m
    p = .5*(lo+hi); C = d2/(h[1]**p-h[2]**p)
    return p, v[2]-C*h[2]**p
tr = {h: json.load(open(f'{R}/large_modes_h{h:.6f}.json'))['trace'] for h in H}
amp = [tr[h][-1]['norm'][2]/tr[h][0]['norm'][2] for h in H]
res = [tr[h][-1]['relative_symmetry_residual'] for h in H]
disp = [[c['fitted_displacement'] for c in tr[h][-1]['late_core_shape']] for h in H]
print('old relative-shift amplification t160:', np.round(amp, 6), ' Richardson p=%.2f continuum=%.3f' % rich(amp, H))
for j, nm in enumerate(['translation', 'rotation']):
    r = [x[j] for x in res]
    print('old %s residual t160:' % nm, np.round(r, 5), ' observed orders: %.2f, %.2f' % (np.log(r[0]/r[1])/np.log(1.5), np.log(r[1]/r[2])/np.log(4/3)))
for j, nm in enumerate(['smaller core (z<0)', 'larger core (z>0)']):
    v = [d[j] for d in disp]
    p, f0 = rich(v, H)
    print('old fitted displacement %s:' % nm, np.round(v, 4), ' Richardson p=%.2f continuum=%.3f' % (p, f0))
# time-history of old fitted displacement at finest grid vs Cartesian impact_fine x/delta
cart = np.genfromtxt(f'{R}/impact_fine.csv', delimiter=',', names=True)
print('t   old(h=.25) core0, core1    Cartesian impact_fine x/delta core0, core1')
for r in tr[0.25]:
    t = r['time']
    if t % 20 == 0 and t >= 60:
        k = np.argmin(abs(cart['t']-t))
        print('%3d  %8.3f %8.3f      %8.3f %8.3f' % (t, r['late_core_shape'][0]['fitted_displacement'], r['late_core_shape'][1]['fitted_displacement'], cart['c0_x'][k]/0.02, cart['c1_x'][k]/0.02))
