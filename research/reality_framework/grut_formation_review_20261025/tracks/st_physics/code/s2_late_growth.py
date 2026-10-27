"""st_physics step 2 (own analysis of Codex's published response arrays).

(a) Effective growth exponent per conformal time between the published 0.03 rad and 1e-6 rad unit-amplitude samples:
    a^2 = 2C(tau - tau_*) => eta = a/C + const, so mu_eff = ln(0.03/1e-6) / ((a_1e-6 - a_0.03)/C).
(b) Late envelope growth rate per momentum: max|response| in windows of one radial period (in eta) of the matched
    quartic orbit, least-squares slope of ln(envelope) against eta over tau in [200, 790].  Compared with the
    Floquet exponent of the matched quartic orbit (results/orbits_floquet.json) at the same k.
(c) Where on the X orbit the published unit crossings fall (|chi_X| relative to r_min, r_max).

usage: python3 -I s2_late_growth.py <codex timing packet dir> <orbits_floquet.json> <out.json>
"""
import sys, json
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
import orbit_tools as O

pkt, orbf, out = sys.argv[1:4]
orb = json.load(open(orbf))
rr = json.load(open(f'{pkt}/RESPONSE_RESULTS.json'))
rows = {r['tag']: r for r in rr['rows']}
res = {}
for tag in ('quiet3', 'quiet5'):
    d = np.load(f'{pkt}/{tag}.npz')
    times, bg, ks = d['times'], d['bg'], d['ks']
    amp, chg = d['amplitude_per_radian'], d['charge_fraction_per_radian']
    a, c, cp, E, L = O.conformal_X(bg)
    row = rows[tag]
    C = row['late_C']
    tau_star = row['matched_radiation_tau_offset']
    # radiation check of a(tau) after tau = 80
    m80 = times >= 80
    rad_dev = float(np.max(abs(a[m80] / np.sqrt(2 * C * (times[m80] - tau_star)) - 1)))
    eta = a / C
    cr = {c_['phase_seed_radians']: c_ for c_ in row['crossings']}
    a03 = cr[0.03]['amplitude_unity']['a_bracket'][1]
    a06 = cr[1e-6]['amplitude_unity']['a_bracket'][1]
    mu_eff = np.log(0.03 / 1e-6) / ((a06 - a03) / C)
    # same with bracket ends to show the sampling spread
    mu_eff_lo = np.log(0.03 / 1e-6) / ((cr[1e-6]['amplitude_unity']['a_bracket'][1] - cr[0.03]['amplitude_unity']['a_bracket'][0]) / C)
    mu_eff_hi = np.log(0.03 / 1e-6) / ((cr[1e-6]['amplitude_unity']['a_bracket'][0] - cr[0.03]['amplitude_unity']['a_bracket'][1]) / C)
    a03c = cr[0.03]['local_charge_unity']['a_bracket'][1]
    a06c = cr[1e-6]['local_charge_unity']['a_bracket'][1]
    mu_eff_charge = np.log(0.03 / 1e-6) / ((a06c - a03c) / C)
    lm = orb[tag]['orbits']['late_mean']
    T = lm['radial_period_eta']
    flo = {round(r['k'], 4): r['mu'] for r in lm['KS']}
    flo10 = {round(r['k'], 4): r['mu'] for r in orb[tag]['orbits']['tau10']['KS']}
    sel = (times >= 200) & (times <= 790)
    edges = np.arange(eta[sel][0], eta[sel][-1], T)
    fits = []
    for j, k in enumerate(ks):
        out_k = dict(k=float(k), floquet_mu_late_orbit=flo[round(float(k), 4)], floquet_mu_tau10_orbit=flo10[round(float(k), 4)])
        for name, arr in (('amplitude', amp[:, j]), ('charge', chg[:, j])):
            env, ec = [], []
            for lo in edges[:-1]:
                w = sel & (eta >= lo) & (eta < lo + T)
                if w.sum() < 3:
                    continue
                env.append(np.max(abs(arr[w])))
                ec.append(lo + T / 2)
            env, ec = np.array(env), np.array(ec)
            p = np.polyfit(ec, np.log(env), 1)
            resid = np.log(env) - np.polyval(p, ec)
            out_k[f'{name}_mu_fit'] = float(p[0])
            out_k[f'{name}_fit_rms_resid'] = float(np.sqrt(np.mean(resid ** 2)))
        fits.append(out_k)
    # vector norm growth of the X part at the selected k (selected_vectors: q(6), p(6))
    sv, svk = d['selected_vectors'], d['selected_vector_ks']
    vec_fits = []
    for j, k in enumerate(svk):
        nX = np.sqrt(sv[:, j, 0] ** 2 + sv[:, j, 1] ** 2 + (sv[:, j, 6] ** 2 + sv[:, j, 7] ** 2) / 4.0)
        env, ec = [], []
        for lo in edges[:-1]:
            w = sel & (eta >= lo) & (eta < lo + T)
            if w.sum() < 3:
                continue
            env.append(np.max(nX[w])); ec.append(lo + T / 2)
        p = np.polyfit(np.array(ec), np.log(env), 1)
        vec_fits.append(dict(k=float(k), X_norm_mu_fit=float(p[0])))
    # where the published crossings sit on the orbit
    def at(tau):
        i = int(np.argmin(abs(times - tau)))
        r = float(np.hypot(*c[i]))
        return dict(tau=float(times[i]), a=float(a[i]), r_chi=r, r_over_rmin=r / lm['rmin'], r_over_rmax=r / lm['rmax'],
                    phys_absx=float(r / (np.sqrt(2) * a[i])))
    cross_pos = {f'{s}_amplitude_unity': at(cr[s]['amplitude_unity']['tau_bracket'][1]) for s in (0.03, 1e-6)}
    cross_pos.update({f'{s}_charge_unity': at(cr[s]['local_charge_unity']['tau_bracket'][1]) for s in (0.03, 1e-6)})
    # alternative, pericentre-free measure: |delta chi_X| / r_max at the peak k vector (k = 2.13 saved)
    jk = int(np.argmin(abs(svk - 2.13)))
    dchi = np.hypot(sv[:, jk, 0], sv[:, jk, 1])
    first = {}
    for s in (0.03, 1e-6):
        for name, val in (('dchi_over_rmax', s * dchi / lm['rmax']), ('dchi_over_rcirc', s * dchi / lm['r_circ'])):
            w = np.flatnonzero(val >= 1)
            first[f'{s}_{name}'] = None if not len(w) else dict(tau=float(times[w[0]]), a=float(a[w[0]]))
    peak_fit = max(fits, key=lambda f: f['amplitude_mu_fit'])
    res[tag] = dict(C=C, tau_star=tau_star, max_rel_dev_from_radiation_a_after_tau80=rad_dev,
                    crossing_a=dict(a_003=a03, a_1e6=a06, a_003_charge=a03c, a_1e6_charge=a06c),
                    mu_eff_amplitude=float(mu_eff), mu_eff_amplitude_bracket=[float(mu_eff_lo), float(mu_eff_hi)],
                    mu_eff_charge=float(mu_eff_charge),
                    floquet_peak_late_orbit=lm['peak'], floquet_peak_tau10_orbit=orb[tag]['orbits']['tau10']['peak'],
                    late_fits=fits, peak_k_by_fit=peak_fit, vector_fits=vec_fits,
                    crossing_positions_on_orbit=cross_pos, pericentre_free_first_unit=first,
                    late_mean_orbit=dict(rmin=lm['rmin'], rmax=lm['rmax'], r_circ=lm['r_circ'], T_eta=T))
    print(tag, 'C', C, 'rad_dev', rad_dev, 'mu_eff', mu_eff, [mu_eff_lo, mu_eff_hi], 'mu_eff_charge', mu_eff_charge)
    for f in fits:
        print('  k %.2f  fitA %.4f fitQ %.4f  Floquet(late) %.4f (tau10) %.4f  rmsA %.3f' % (
            f['k'], f['amplitude_mu_fit'], f['charge_mu_fit'], f['floquet_mu_late_orbit'], f['floquet_mu_tau10_orbit'],
            f['amplitude_fit_rms_resid']))
    print('  vector fits', vec_fits)
    print('  crossings on orbit', json.dumps(cross_pos))
    print('  pericentre-free', first)
json.dump(res, open(out, 'w'), indent=1)
