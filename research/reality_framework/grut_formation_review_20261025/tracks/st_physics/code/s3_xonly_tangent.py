"""st_physics step 3 (own code): time-domain test that the late full-system tangent is the X-only quartic
linearisation.

Take Codex's saved 12-component tangent at k = 2.13 and 5.8 (and 0.02) at tau0 = 200 (selected_vectors in
quiet*.npz: q = chi-components of (x, y, s), p = conformal momenta), keep only the X part, and propagate
    d^2 dchi/d eta^2 = -[k^2 + lam(|chi|^2 I + 2 chi chi^T)] dchi
on (i) a cubic spline of the published X background in eta = a/C and (ii) our own pure-quartic orbit started from
the background at tau0.  Compare with Codex's X components at later samples (relative vector difference).
Also reports the drift of the pure-quartic orbit from the published background (phase error), which measures
the non-quartic (Y/Sigma) influence on the background over this interval.

usage: python3 -I s3_xonly_tangent.py <codex timing packet dir> <out.json>
"""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
sys.path.insert(0, __import__('os').path.dirname(__file__))
import orbit_tools as O

pkt, out = sys.argv[1:3]
rr = {r['tag']: r for r in json.load(open(f'{pkt}/RESPONSE_RESULTS.json'))['rows']}
LAM = O.LAM
res = {}
for tag in ('quiet3', 'quiet5'):
    d = np.load(f'{pkt}/{tag}.npz')
    times, bg, sv, svk = d['times'], d['bg'], d['selected_vectors'], d['selected_vector_ks']
    a, c, cp, E, L = O.conformal_X(bg)
    C = rr[tag]['late_C']
    eta = a / C
    m = times >= 150
    spl = CubicSpline(eta[m], c[m])
    i0 = int(np.argmin(abs(times - 200)))
    e0 = eta[i0]
    # (ii) own pure-quartic background from tau0
    def bgrhs(t, y):
        q = y[:2]
        return np.r_[y[2:], -LAM * (q @ q) * q]
    own = solve_ivp(bgrhs, (e0, eta[-1]), np.r_[c[i0], cp[i0]], method='DOP853', rtol=1e-12, atol=1e-10,
                    dense_output=True)
    later = times >= 200
    dif = np.array([np.linalg.norm(own.sol(e)[:2] - c[j]) for e, j in zip(eta[later], np.flatnonzero(later))])
    drift = dict(max_abs_dchi=float(dif.max()), rmax=float(np.max(np.hypot(*c[later].T))),
                 dchi_at_800=float(dif[-1]))
    row = dict(background_drift_pure_quartic_vs_published=drift, vectors=[])
    for j, k in enumerate(svk):
        v0 = sv[i0, j]
        y0 = np.r_[v0[0:2], v0[6:8]]
        for label, chifun in (('spline_bg', lambda e: spl(e)), ('own_quartic_bg', lambda e: own.sol(e)[:2])):
            def rhs(t, y):
                q = chifun(t)
                M = k * k * np.eye(2) + LAM * ((q @ q) * np.eye(2) + 2 * np.outer(q, q))
                return np.r_[y[2:], -M @ y[:2]]
            s = solve_ivp(rhs, (e0, eta[-1]), y0, method='DOP853', rtol=1e-11, atol=1e-14, t_eval=eta[later],
                          dense_output=False)
            ours = s.y.T
            theirs = np.c_[sv[later, j, 0:2], sv[later, j, 6:8]]
            nrm = np.linalg.norm(theirs, axis=1)
            rel = np.linalg.norm(ours - theirs, axis=1) / nrm
            # amplitude response from X part (dchi . chi / |chi|^2), compared at the samples
            row['vectors'].append(dict(k=float(k), background=label, rel_vector_diff_max=float(rel.max()),
                                       rel_vector_diff_at_800=float(rel[-1]),
                                       growth_factor_codex_200_to_800=float(nrm[-1] / nrm[0])))
            print(tag, 'k', k, label, 'max rel diff %.3e at800 %.3e growth %.3e' % (rel.max(), rel[-1], nrm[-1] / nrm[0]),
                  flush=True)
    print(tag, 'background drift', drift)
    res[tag] = row
json.dump(res, open(out, 'w'), indent=1)
