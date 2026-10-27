"""st_physics skeptic, check k6: arithmetic re-derivations (own code; [D] on published numbers, [Q] widths).

P6/P7 extrapolations with our quadrature eta; P11/P12 circular orbit and Xbar budget; P13 wavelength-cube charge;
Y-sextic ratio at the actual a = 1 state; Y_L clocks at the published unit samples by our own trapezoid quadrature of
the published a(tau) plus the analytic radiation tail.

usage: python3 -I k6_arith.py <codex timing packet dir> <k2_strobe_growth.json> <k1_corot_floquet.json> <out.json>
"""
import sys, json
import numpy as np
from scipy.special import gammaincc, gamma as G

pkt, k2f, k1f, out = sys.argv[1:5]
k2, k1 = json.load(open(k2f)), json.load(open(k1f))
rr = {r['tag']: r for r in json.load(open(f'{pkt}/RESPONSE_RESULTS.json'))['rows']}
sc = {h['tag']: h for h in json.load(open(f'{pkt}/SOURCE_CLOCKS.json'))['histories']}
LAM, KAP, ML2 = 1e-3, 0.0034586847211042523, 900 - 2 * 7.159456159513794 * 28.284271247461902
res = {}
for tag in ('quiet3', 'quiet5'):
    d = np.load(f'{pkt}/{tag}.npz'); t = d['times']; a = np.exp(d['bg'][:, 15])
    C, ts = k2[tag]['C_end'], k2[tag]['tau_star']
    mu = k2[tag]['mu_eff_quadrature_eta']['amplitude_unity']['mu']
    cr = {c['phase_seed_radians']: c for c in rr[tag]['crossings']}
    a06 = cr[1e-6]['amplitude_unity']['a_bracket'][1]
    a9 = a06 + C / mu * np.log(1e3)
    gam = sc[tag]['sources']['YL']['gamma_over_Hi']
    tauYL = np.log(2) / gam
    aYL = np.sqrt(2 * C * (tauYL - ts))
    log10s = -6 - mu * (aYL - a06) / C / np.log(10)
    # Y_L clocks
    def FE(tau):
        w = t <= tau
        num = np.trapezoid(gam * a[w] * np.exp(-gam * t[w]), t[w])
        den = np.trapezoid(gam * a * np.exp(-gam * t), t) + np.sqrt(2 * C / gam) * np.exp(-gam * ts) * G(1.5) * gammaincc(1.5, gam * (t[-1] - ts))
        return float(num / den)
    tc = cr[0.03]['amplitude_unity']['tau_bracket'][1]
    # circular orbit and Xbar (tau10 orbit)
    o = k1[tag + '_tau10']; E, L = o['E'], o['L']
    Ec = 0.75 * LAM ** (1 / 3) * L ** (4 / 3); rc = (L / np.sqrt(LAM)) ** (1 / 3); kp = k1[tag + '_tau10']['peak']['k']
    w_lo, w_hi = kp, np.sqrt(kp ** 2 + 3 * LAM * rc ** 2)
    Tth = (15 * (E - Ec) / np.pi ** 2) ** 0.25
    norm = abs(rr[tag]['late_NX_normalization_through80'])
    res[tag] = dict(
        P6_unit_1e9=dict(a=a9, tau=ts + a9 ** 2 / (2 * C)), P7=dict(tau_YL=tauYL, a_YL=aYL, log10_seed=log10s),
        YL_at_003_unit=dict(tau=tc, F_N=float(1 - np.exp(-gam * tc)), F_E=FE(tc),
                            codex=[x for x in sc[tag]['sources']['YL']['at_response_crossings'] if x['kind'] == 'amplitude_unity' and x['phase_seed'] == 0.03]),
        P11=dict(E_circ=Ec, E_circ_over_E=Ec / E, release=1 - Ec / E),
        P12=dict(pairs=[(1 - Ec / E) * E / (2 * w_hi * L), (1 - Ec / E) * E / (2 * w_lo * L)],
                 gas=[(E / (w_hi * L) - 1) / 2, (E / (w_lo * L) - 1) / 2],
                 equilibrium_mu0=float(1.2020569 * Tth ** 3 / np.pi ** 2 / L), T=Tth),
        P13=dict(Q_cube_2p13=norm * (2 * np.pi / 2.13) ** 3, Q_cube_own_peak=norm * (2 * np.pi / kp) ** 3, k_peak=kp,
                 Q_cube_2p13_late_L=k1[tag + '_late']['L'] * (2 * np.pi / 2.13) ** 3),
        Y_sextic_over_quartic_actual_a1=float(KAP ** 2 * 76.39978266 ** 2 / (ML2 * LAM)))
    print(tag, json.dumps(res[tag], default=float)[:1800])
json.dump(res, open(out, 'w'), indent=1, default=float)
