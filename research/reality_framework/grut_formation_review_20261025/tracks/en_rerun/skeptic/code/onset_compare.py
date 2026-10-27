"""en_rerun_skeptic -- compare the own-engine N48 history tau 20-85 (own_engine.py from the tau-20 Codex initial field)
with Codex's PUBLISHED FIELD_CURVES (N48, full resolution): onset thresholds and curves. Own code.
usage: python3 -I onset_compare.py <own_20_85.json> <FIELD_CURVES.npz> <out_json>"""
import sys, json
import numpy as np
LX, G2, LP = .001, .01, .03
rows = json.load(open(sys.argv[1]))['rows']; C = np.load(sys.argv[2]); outp = sys.argv[3]
tau = np.array([r['tau'] for r in rows])
g = lambda k: np.array([r[k] for r in rows])
mr2 = g('q0')**2 + g('q1')**2
contrast = np.sqrt(np.maximum(g('charge_RMS')**2 - g('nQ')**2, 0))/abs(g('nQ'))
fluct = np.maximum(g('X2') - mr2, 0)/g('X2')
EX0 = .5*(g('p0')**2 + g('p1')**2) + LX/4*mr2**2 + G2/2*g('qphi')**2*mr2
EP0 = .5*g('pphi')**2 + LP/4*g('qphi')**4
outside = 1 - (EX0 + EP0)/g('rho')
phimix = (g('Kphi') + g('Gphi') + g('Uphi') + g('Umix'))/g('rho')
pt = C['N48_tau']; idx = np.searchsorted(pt, tau - 1e-7); assert np.max(abs(pt[idx] - tau)) < 1e-6
cmp = {}
for name, mine in [('charge_contrast_RMS', contrast), ('X_field_fluctuation_variance_fraction', fluct),
                   ('energy_outside_spatial_mean_fraction', outside), ('phi_plus_interaction_energy_fraction', phimix)]:
    pub = C['N48_' + name][idx]
    cmp[name] = float(np.max(abs(mine - pub))/np.max(abs(pub)))
first = lambda v, lev: float(tau[np.argmax(v >= lev)])
# sub-sample interpolated crossings (linear) -- the published thresholds are first SAMPLE at/above the level
def interp_cross(v, lev):
    j = int(np.argmax(v >= lev)); return float(tau[j-1] + (lev - v[j-1])*(tau[j] - tau[j-1])/(v[j] - v[j-1]))
res = dict(own_vs_published_curve_maxrel_tau20_85=cmp,
           energy_drift_own=float(np.max(abs(g('rho')/g('rho')[0] - 1))), charge_drift_own=float(np.max(abs(g('nQ')/g('nQ')[0] - 1))),
           first_contrast_one_sample=first(contrast, 1.), first_contrast_one_interp=interp_cross(contrast, 1.),
           first_X_variance_half_sample=first(fluct, .5), first_X_variance_half_interp=interp_cross(fluct, .5),
           crossings_contrast_one=[float(tau[j]) for j in np.flatnonzero((contrast[1:] >= 1) & (contrast[:-1] < 1)) + 1][:8],
           crossings_X_variance_half=[float(tau[j]) for j in np.flatnonzero((fluct[1:] >= .5) & (fluct[:-1] < .5)) + 1][:8])
json.dump(res, open(outp, 'w'), indent=1); print(json.dumps(res, indent=1))
