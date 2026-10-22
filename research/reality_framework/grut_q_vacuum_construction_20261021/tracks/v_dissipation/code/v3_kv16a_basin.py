"""v_dissipation, script 3: KV16a basin of the Minkowski end state (track NOTES sec 4.2), recomputed independently.

Method (not the track's): the drain per e-fold does not depend on H (|H| dt = dN for H > 0), so in E_* units
    d r_V/dN = -(gamma/4) (r_d,i e^{-3N} + 4 r_V)^2                         (v1 check D_Estar_dust_drain)
holds for radiation + dust + vacuum whatever the radiation content, as long as H > 0.  The basin boundary is the
trajectory with r_V(N -> inf) = 0+, obtained by integrating this Riccati equation BACKWARD from N = 40 (where
r_V = (gamma/24) r_d^2 e^{-6N}) to N = 0; backward integration is stable because forward neighbours diverge.
It is therefore independent of lambda AND of r_r,i.  A forward check in tau (state h, r_V, r_d; radiation from the
constraint) at r_r,i = 10 and r_r,i = 1000 brackets the same boundary.
Grid fractions: lambda in [-1, 1] step 0.05 (41) x f_i in [0.20, 3.00] step 0.05 (57), r_V = lambda - f^2 + 3/f^2.
Track values [NOTES sec 4.2 table]: fractions 0.4018 / 0.4014 / 0.3188 / 0 (gamma_phys; r_d = 0, 10, 100, 1000),
0.4018 / 0.3992 / 0.1904 / 0 (gamma = 0.01); boundaries 0, 0.01329, 1.419, above box; 0, 0.04261, 5.368, above box.
Output: ../results/v3_kv16a_basin.json.  Run: python3 -I v3_kv16a_basin.py
"""
import json, os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

GPH = 1 / (32 * np.pi**2)
lams = np.round(np.arange(-1.0, 1.0 + 1e-9, 0.05), 10)
fis = np.round(np.arange(0.20, 3.0 + 1e-9, 0.05), 10)
LL, FF = np.meshgrid(lams, fis, indexing='ij')
RV = LL - FF**2 + 3 / FF**2
ntot = RV.size


def boundary_N(rd, gam, Nbig=40.0):
    if rd == 0:
        return 0.0
    def rhs(N, r):
        return [-(gam / 4) * (rd * np.exp(-3 * N) + 4 * r[0])**2]
    def blow(N, r):          # backward blow-up: the separatrix does not reach N = 0 (no finite boundary)
        return r[0] - 1e8
    blow.terminal = True
    r_end = gam / 24 * rd**2 * np.exp(-6 * Nbig)
    s = solve_ivp(rhs, (Nbig, 0.0), [r_end], method='LSODA', rtol=1e-12, atol=1e-30, events=[blow])
    if s.status == 1:
        return float('inf')
    return float(s.y[0, -1])


def forward_class(rV1, rr1, rd1, gam, Nmax=60.0):
    """forward integration in N of (r_V, r_r, r_d); classify by sign of r_V once dust is negligible."""
    def rhs(N, y):
        rV, rr, rdd = y
        k = 2 * rV + rdd / 2                    # hdot + 2h^2 in E_* units
        drain = gam * k * k                     # (gamma/4)(r_d + 4 r_V)^2 = gamma k^2  per e-fold
        return [-drain, -4 * rr + drain, -3 * rdd]
    def ev(N, y):                 # r_V crosses zero from above (monotone drain: no return)
        return y[0]
    ev.terminal = True; ev.direction = -1

    def ev_h(N, y):               # h^2 = r_V + r_r + r_d reaches zero: turnaround, crunch
        return y[0] + y[1] + y[2]
    ev_h.terminal = True; ev_h.direction = -1
    s = solve_ivp(rhs, (0, Nmax), [rV1, rr1, rd1], method='LSODA', rtol=1e-11,
                  atol=[1e-14 * max(abs(rV1), 1e-3), 1e-30, 1e-30], events=[ev, ev_h])
    if rV1 < 0 and s.status == 0:
        return 'undecided'
    if s.status == 1:
        return 'crunch'
    return 'mink' if s.y[0, -1] > 0 else 'crunch'


out = {'grid': {'lambda': [float(lams[0]), float(lams[-1]), 0.05, len(lams)], 'f_i': [0.2, 3.0, 0.05, len(fis)], 'points': int(ntot)}}
track_frac = {('phys', 0): 0.4018, ('phys', 10): 0.4014, ('phys', 100): 0.3188, ('phys', 1000): 0.0,
              ('0.01', 0): 0.4018, ('0.01', 10): 0.3992, ('0.01', 100): 0.1904, ('0.01', 1000): 0.0}
track_bnd = {('phys', 10): 0.01329, ('phys', 100): 1.419, ('0.01', 10): 0.04261, ('0.01', 100): 5.368}
rows = []
for gname, gam in (('phys', GPH), ('0.01', 0.01)):
    for rd in (0, 10, 100, 1000):
        b = boundary_N(rd, gam)
        frac = float(np.mean(RV > b)) if rd > 0 else float(np.mean(RV >= 0))
        row = {'gamma': gname, 'r_d,i': rd, 'boundary_rV_i': b, 'leading_gamma_rd^2/24': gam * rd**2 / 24,
               'fraction': frac, 'track_fraction': track_frac[(gname, rd)],
               'track_boundary': track_bnd.get((gname, rd)), 'max_rV_in_box': float(RV.max())}
        # forward cross-check at two radiation contents: bracket the boundary at +-1e-6 relative
        if rd > 0 and np.isfinite(b):
            chk = {}
            for rr1 in (10.0, 1000.0):
                lo = forward_class(b * (1 - 1e-6), rr1, rd, gam)
                hi = forward_class(b * (1 + 1e-6), rr1, rd, gam)
                chk['r_r=%g' % rr1] = [lo, hi]
            row['forward_bracket_(-1e-6,+1e-6)'] = chk
        rows.append(row)
out['rows'] = rows
# radiation-only: analytic boundary f*(lambda) = sqrt((lambda + sqrt(lambda^2 + 12))/2)
out['f_star(lambda)'] = {str(l): float(np.sqrt((l + np.sqrt(l * l + 12)) / 2)) for l in (-1, -0.5, 0, 0.5, 1)}
# spot-check radiation-only classification by forward integration on 300 random grid points
rng = np.random.default_rng(1)
idx = rng.choice(ntot, 300, replace=False)
mism = 0
for i in idx:
    rv = RV.flat[i]
    c = forward_class(rv, 10.0, 0.0, GPH, Nmax=200.0) if rv != 0 else 'mink'
    if (c == 'mink') != (rv >= 0):
        mism += 1
out['radiation_only_spotcheck_300_points_mismatches'] = mism
with open(os.path.join(os.path.dirname(__file__), '..', 'results', 'v3_kv16a_basin.json'), 'w') as fh:
    json.dump(out, fh, indent=1)
for r in rows:
    print(r['gamma'], r['r_d,i'], 'bnd=%.5g (lead %.4g, track %s)' % (r['boundary_rV_i'], r['leading_gamma_rd^2/24'], r['track_boundary']),
          'frac=%.4f (track %.4f)' % (r['fraction'], r['track_fraction']), r.get('forward_bracket_(-1e-6,+1e-6)', ''))
print('spotcheck mismatches', mism, 'max rV in box', RV.max())
