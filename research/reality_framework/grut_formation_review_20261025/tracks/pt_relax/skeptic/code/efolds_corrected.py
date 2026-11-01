"""pt_relax skeptic (Claude, own code): pt_relax's 2<->4 e-fold table recomputed with the Monte Carlo c_eff
(mc24.py: one real scalar, Bose, mu = 0, Hartree mass on every line) in place of the parametric c = (10-1e3)/(4 pi)^5.
Same V^4 T formula and s_onset = 2.63e12 as pt_relax (results/scales.json item2_rates).  Writes results/efolds_corrected.json.
"""
import json, math, os
HERE = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(HERE, '..', 'results')
mc = {(r['m_over_T'], r['stat']): r['c_eff'] for r in json.load(open(os.path.join(R, 'mc24.json')))['rows']}
T, s_on = 17.556, 2.629e12
LX, G2, LP = 1e-3, 1e-2, 3e-2
c_lo, c_hi = 10 / (4 * math.pi) ** 5, 1e3 / (4 * math.pi) ** 5
cX, cP = mc[(0.16, 'bose')], mc[(0.36, 'bose')]          # M_X/T = 2.84/17.56 = 0.16, M_phi/T = 6.40/17.56 = 0.36
rows = []
for name, V4, cs in [('2<->4 X only, V = 4 lambda_X', (4 * LX) ** 4, (cX, cX)),
                     ('2<->4 two g^2 vertices, V = 2 g^2', (2 * G2) ** 4, (cP, cX)),
                     ('2<->4 mixed, V^4 = (2g^2)^2 (6 lambda_phi)^2', (2 * G2) ** 2 * (6 * LP) ** 2, (cP, cX)),
                     ('2<->4 phi self, V = 6 lambda_phi', (6 * LP) ** 4, (cP, cX))]:
    rows.append(dict(channel=name, pt_relax_efolds=[c_lo * V4 * T * s_on, c_hi * V4 * T * s_on],
                     mc_c_eff_range=list(cs), corrected_efolds=[cs[0] * V4 * T * s_on, cs[1] * V4 * T * s_on]))
    print('%-45s pt_relax %.1e - %.1e   corrected %.1e - %.1e' % (name, *rows[-1]['pt_relax_efolds'], *rows[-1]['corrected_efolds']))
json.dump(dict(T=T, s_onset=s_on, c_pt_relax=[c_lo, c_hi], c_eff_X=cX, c_eff_phi=cP, rows=rows,
               caveat='c_eff from a single real scalar with lambda phi^4/4!; complex-field and mixed-vertex combinatorics change it by O(1-10), not 1e2-1e5'),
          open(os.path.join(R, 'efolds_corrected.json'), 'w'), indent=1)
