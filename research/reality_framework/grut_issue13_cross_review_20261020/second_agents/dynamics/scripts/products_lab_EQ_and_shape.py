"""Late lab-frame E12/Q12 of the outgoing products and R8 anisotropy change (impact vs reference).
Usage: python -I products_lab_EQ_and_shape.py <raw_dir>"""
import sys, numpy as np
R = sys.argv[1]
L = lambda n: np.genfromtxt(f'{R}/{n}.csv', delimiter=',', names=True)
def aniso(a, i):
    cov = np.empty((len(a), 3, 3)); q = a[f'c{i}_Q8']
    for x, y, z in [(0, 0, 'xx'), (0, 1, 'xy'), (0, 2, 'xz'), (1, 1, 'yy'), (1, 2, 'yz'), (2, 2, 'zz')]:
        cov[:, x, y] = cov[:, y, x] = a[f'c{i}_{z}']/q
    ev = np.linalg.eigvalsh(cov); return np.sqrt(ev[:, -1]/ev[:, 0])
for n in ['pair_reference_fine', 'pair_generic_fine', 'pair_negative_fine', 'pair_rotated_fine']:
    a = L(n); m = a['t'] >= 120
    out = []
    for j in (0, 1):
        v = np.sqrt(sum(np.polyfit(a['t'][m], a[f'c{j}_{x}'][m], 1)[0]**2 for x in 'xyz'))
        eq = a[f'c{j}_E12'][m]/a[f'c{j}_Q12'][m]
        out.append('core%d Q12=%.1f lab E/Q=[%.4f,%.4f] speed=%.4f rest-frame E/Q~%.4f' % (j, a[f'c{j}_Q12'][m].mean(), eq.min(), eq.max(), v, eq.mean()*np.sqrt(1-v*v)))
    print(n, ' | '.join(out))
for imp, ref in [('impact_coarse', 'pair_reference_coarse'), ('impact_medium', 'pair_reference_medium'), ('impact_fine', 'pair_reference_fine'), ('impact_half_fine', 'pair_reference_fine'), ('impact_box', 'pair_reference_box'), ('impact_dt', 'pair_reference_dt')]:
    a, b = L(imp), L(ref); m = a['t'] >= 120
    print('%-17s max |d anisotropy| vs ref: %s' % (imp, [float('%.6g' % np.max(abs(aniso(a[m], i)-aniso(b[m], i)))) for i in (0, 1)]))
