"""st_physics step 1 (own code): conformal X invariants along Codex's quiet backgrounds and Floquet spectra of the
matched isolated quartic orbits.

usage: python3 -I s1_orbits_floquet.py <codex timing packet dir> <grut_reheating_modes QUARTIC_ORBIT.json> <out.json>
"""
import sys, json, time
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
import orbit_tools as O

pkt, qorb, out = sys.argv[1:4]
KS = np.array([.02, .05, .1, .2, .4, .8, 1.2, 1.6, 2., 2.13, 2.4, 2.8, 3.2, 4., 5.6, 5.8, 6., 8.])
res = dict(scope='Isolated complex quartic orbit (lam=1e-3, conformal units) matched to the X field of each quiet '
                 'background; Floquet exponents per conformal time. Auxiliary diagnostic, own code.')
t0 = time.time()
for tag in ('quiet3', 'quiet5'):
    d = np.load(f'{pkt}/{tag}.npz')
    times, bg = d['times'], d['bg']
    a, c, cp, E, L = O.conformal_X(bg)
    late = times >= 20
    inv = dict(E_tau10=float(E[np.argmin(abs(times - 10))]), L_tau10=float(L[np.argmin(abs(times - 10))]),
               E_late_mean=float(E[late].mean()), E_late_rel_range=float(np.ptp(E[late]) / E[late].mean()),
               L_late_mean=float(L[late].mean()), L_late_rel_range=float(np.ptp(L[late]) / L[late].mean()),
               E_800=float(E[-1]), L_800=float(L[-1]))
    row = dict(invariants=inv, orbits={})
    for label, (Em, Lm) in dict(tau10=(inv['E_tau10'], inv['L_tau10']),
                                late_mean=(inv['E_late_mean'], inv['L_late_mean'])).items():
        orb = O.Orbit(Em, Lm)
        s = orb.summary()
        # Floquet on Codex's 18 response momenta, a fine scan, the IR and both band edges
        s['KS'] = [dict(k=float(k), mu=orb.mu(k)[0]) for k in KS]
        fine = np.r_[np.linspace(0.05, 7.0, 140)]
        s['fine'] = [dict(k=float(k), mu=orb.mu(k)[0]) for k in fine]
        # refine peak by golden search on the fine grid around its max
        from scipy.optimize import minimize_scalar
        kf = fine[int(np.argmax([r['mu'] for r in s['fine']]))]
        pk = minimize_scalar(lambda k: -orb.mu(k)[0], bounds=(kf - 0.1, kf + 0.1), method='bounded',
                             options=dict(xatol=1e-5))
        s['peak'] = dict(k=float(pk.x), mu=float(-pk.fun), multiplier=float(np.exp(-pk.fun * orb.T)))
        s['ir'] = [dict(k=float(k), mu=orb.mu(k)[0], mu_over_k=orb.mu(k)[0] / k) for k in (0.002, 0.005, 0.01, 0.02, 0.05)]
        mu0, m0, sy0 = orb.mu(0.0)
        s['k0'] = dict(mu=mu0, multiplier_minus_1=m0 - 1, symplectic=sy0)
        # band edges by bisection on mu > threshold
        def edge(lo, hi, thr=1e-7):
            flo = orb.mu(lo)[0] > thr
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if (orb.mu(mid)[0] > thr) == flo:
                    lo = mid
                else:
                    hi = mid
            return 0.5 * (lo + hi)
        mus = np.array([r['mu'] for r in s['fine']])
        on = mus > 1e-7
        edges = [float(edge(fine[i], fine[i + 1])) for i in range(len(fine) - 1) if on[i] != on[i + 1]]
        s['band_edges'] = edges
        s['second_band_peak'] = None
        if len(edges) >= 3:
            sel = (fine > edges[1]) & (fine < (edges[2] if len(edges) > 2 else 7))
            if sel.any():
                j = int(np.argmax(np.where(sel, mus, -1)))
                s['second_band_peak'] = dict(k=float(fine[j]), mu=float(mus[j]))
        s['radau_control_peak'] = orb.mu(s['peak']['k'], method='Radau')[0]
        s['max_symplectic'] = max(orb.mu(k)[2] for k in (0.5, 2.0, 4.0, 6.0))
        row['orbits'][label] = s
        print(tag, label, {k: s[k] for k in ('E', 'L', 'rmin', 'rmax', 'radial_period_eta', 'E_circ_over_E',
                                             'eps_charge_fraction', 'peak', 'band_edges')}, 'IR', s['ir'][2],
              'k0', s['k0'], f'{time.time() - t0:.0f}s', flush=True)
    res[tag] = row

# reproduction check of Codex's tau10 orbit for (3,0.5): compare with QUARTIC_ORBIT.json nodes
pub = json.load(open(qorb))
o3 = O.Orbit(res['quiet3']['invariants']['E_tau10'], res['quiet3']['invariants']['L_tau10'])
cmp = [(r['k'], r['mu_conformal'], o3.mu(r['k'])[0]) for r in pub['rows'][::10]]
res['codex_quartic_orbit_check'] = dict(period_codex=pub['radial_period'], period_ours=o3.T,
                                        rotation_codex=pub['rotation_angle'], rotation_ours=o3.dth,
                                        max_rel_mu_diff=max(abs(b - c) / max(b, 1e-3) for _, b, c in cmp), nodes=cmp)
# pure-quartic scaling: chi -> s chi maps k -> s k, mu -> s mu, eps invariant
s35 = (res['quiet5']['orbits']['tau10']['E'] / res['quiet3']['orbits']['tau10']['E']) ** 0.25
res['scaling_check'] = dict(s_5_over_3=s35,
                            eps3=res['quiet3']['orbits']['tau10']['eps_charge_fraction'],
                            eps5=res['quiet5']['orbits']['tau10']['eps_charge_fraction'],
                            mu_peak3_times_s=res['quiet3']['orbits']['tau10']['peak']['mu'] * s35,
                            k_peak3_times_s=res['quiet3']['orbits']['tau10']['peak']['k'] * s35)
json.dump(res, open(out, 'w'), indent=1)
print(json.dumps(res['codex_quartic_orbit_check'], indent=1)[:600], res['scaling_check'])
