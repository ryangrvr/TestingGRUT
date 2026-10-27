"""st_rerun: own arithmetic for the decay clocks and a growth-rate diagnostic, on Codex's background arrays.

Independent of Codex's Clock class (no PCHIP, no Gauss panels, no closed-form tail): early energy integral by
trapezoid on the dense a(tau) samples (dtau=0.01, from our dense rerun with Codex's integrator), tail by direct
numerical quadrature of gamma*sqrt(2C(s-tau*))*exp(-gamma*s) in the variable u=gamma(s-tau*).
Also: ln(response envelope) vs conformal time eta=int dtau/a at k=2.13 (2.0), compared with the auxiliary quartic
mu_conformal; a labelled [E] extrapolation of when a 1e-9 rad seed would reach unity if that growth continued.
Usage: python3 -I clocks_check.py <published_dir> <dense_dir> <out_json>
"""
import sys, json, math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
pub, dense, outp = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
SC = json.loads((pub / 'SOURCE_CLOCKS.json').read_text())
RR = {r['tag']: r for r in json.loads((pub / 'RESPONSE_RESULTS.json').read_text())['rows']}
QS = json.loads((pub / 'QUARTIC_SCAN.json').read_text())
mu_peak = max(r['mu_conformal'] for r in QS['rows'])
out = dict(quartic_mu_peak=mu_peak)
for h in SC['histories']:
    tag = h['tag']; b = 'b3' if tag == 'quiet3' else 'b5'
    z = np.load(pub / (tag + '.npz')); bg = z['bg']; a_pub = np.exp(bg[:, 15]); H = bg[:, 16]; t_pub = z['times']
    C_arr = a_pub ** 2 * H
    C = float(C_arr[-1]); off = float(t_pub[-1] - a_pub[-1] ** 2 / (2 * C))
    d = np.load(dense / f'time_{b}_r10.npz'); t = d['times']; a = d['a']
    row = dict(C_own=C, C_codex=h['late_C'], offset_own=off, offset_codex=h['matched_radiation_tau_offset'],
               C_range_last10pct=float(np.ptp(C_arr[t_pub >= 720]) / C),
               dense_a800_vs_published=float(a[-1] / a_pub[-1] - 1), sources={})
    for name, s in h['sources'].items():
        g = s['gamma_over_Hi']
        def early(tau):
            m = t <= tau; tt = np.r_[t[m], tau] if t[m][-1] < tau else t[m]
            aa = np.interp(tt, t, a)
            return float(np.trapezoid(g * aa * np.exp(-g * tt), tt))
        u0 = g * (800. - off)
        tail = quad(lambda u: math.sqrt(2 * C * u / g) * math.exp(-u - g * off), u0, np.inf, epsabs=0, epsrel=1e-12, limit=500)[0]
        total = early(800.) + tail
        med_tau = math.log(2) / g; med_a = math.sqrt(2 * C * (med_tau - off))
        cross = []
        for c in s['at_response_crossings']:
            fN = -math.expm1(-g * c['tau']); fE = early(c['tau']) / total
            cross.append(dict(seed=c['phase_seed'], kind=c['kind'], tau=c['tau'], fN_own=fN, fN_codex=c['decayed_number_fraction_rest_clock'],
                              fE_own=fE, fE_codex=c['cumulative_conformal_energy_fraction_prescribed_source'],
                              fE_rel_diff=fE / c['cumulative_conformal_energy_fraction_prescribed_source'] - 1))
        row['sources'][name] = dict(gamma=g, median_tau_own=med_tau, median_tau_codex=s['number_quantiles'][2]['tau'],
                                    median_a_own=med_a, median_a_codex=s['number_quantiles'][2]['a'],
                                    mean_a_own=total, mean_a_codex=s['mean_a_per_decay'], mean_a_rel_diff=total / s['mean_a_per_decay'] - 1,
                                    max_fE_rel_diff=max(abs(x['fE_rel_diff']) for x in cross), crossings=cross)
    # growth-rate diagnostic on dense responses
    amp = d['amp']; ks = d['ks']; eta = np.r_[0, np.cumsum(np.diff(t) * 0.5 * (1 / a[1:] + 1 / a[:-1]))]
    j = int(np.argmin(abs(ks - (2.13 if tag == 'quiet3' else 2.0))))
    env = np.abs(amp).max(axis=1)
    fits = {}
    for lo, hi in [(150, 400), (400, 800), (150, 800)]:
        # envelope: running maximum in windows of 20 tau
        sel = (t >= lo) & (t <= hi); tt = t[sel]; ee = env[sel]; et = eta[sel]
        edges = np.arange(lo, hi + 1e-9, 20.); px = []; py = []
        for e0, e1 in zip(edges[:-1], edges[1:]):
            m = (tt >= e0) & (tt < e1)
            if m.any():
                i = np.argmax(ee[m]); px.append(et[m][i]); py.append(np.log(ee[m][i]))
        slope = float(np.polyfit(px, py, 1)[0]); fits[f'{lo}-{hi}'] = slope
    # [E] extrapolation for 1e-9: last envelope max, continue with slope of 400-800
    s = fits['400-800']; late = t >= 700; il = int(np.argmax(env[late])); lnA800 = math.log(env[late][il])
    eta_ref = float(eta[late][il]); a_ref = float(a[late][il])
    deta = (math.log(1e9) - lnA800) / s; eta_c = eta_ref + deta
    # matched radiation tail: d eta = d tau / a and a^2 = 2C(tau - tau*)  =>  a(eta) = a_ref + C (eta - eta_ref)
    a_c = a_ref + C * deta; tau_c = off + a_c ** 2 / (2 * C)
    row['growth'] = dict(k=float(ks[j]), eta_800=float(eta[-1]), slope_ln_envelope_per_eta=fits, quartic_mu_peak=mu_peak,
                         ln_env_per_radian_max_700_800=lnA800, tau_of_that_max=float(t[late][il]),
                         E_extrapolated_1e9_unity=dict(eta=eta_c, a=a_c, tau=tau_c,
                                                       label='[E] assumes the 400-800 envelope slope persists; not computed'))
    out[tag] = row
    print(tag, json.dumps({k: v for k, v in row.items() if k != 'sources'}, indent=None, default=float))
    for name, s_ in row['sources'].items():
        print('  ', name, 'median', s_['median_tau_own'], s_['median_tau_codex'], s_['median_a_own'], s_['median_a_codex'],
              'mean_a rel', s_['mean_a_rel_diff'], 'max fE rel', s_['max_fE_rel_diff'])
outp.write_text(json.dumps(out, indent=1, default=float) + '\n')
