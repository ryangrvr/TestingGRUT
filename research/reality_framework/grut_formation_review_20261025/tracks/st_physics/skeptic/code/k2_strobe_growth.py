"""st_physics skeptic, check k2 (own analysis of Codex's published arrays, different method from the reviewer).

(a) Conformal time by direct quadrature eta = int dtau / a over the published samples (no radiation-law mapping).
(b) Stroboscopic late growth: locate every pericentre of the published X orbit (local minima of |chi_X| refined by a
    parabola), take the published response at the sample nearest each pericentre (amplitude and charge arrays, all
    18 k) and the X part of the saved 12-vectors (k = 0.02, 2.13, 5.8) in the CO-ROTATING frame at the pericentre,
    and regress ln|.| on eta over tau in [200, 800].  Compared with our own co-rotating Floquet exponents (k1).
(c) Effective growth per conformal time from the published unit samples, using our quadrature eta.
(d) L(t) and E(t) oscillation and drift of the conformal X invariants; C = a^2 H and the radiation tail.

usage: python3 -I k2_strobe_growth.py <codex timing packet dir> <k1_corot_floquet.json> <out.json>
"""
import sys, json
import numpy as np
from scipy.integrate import cumulative_trapezoid

LAM = 1e-3
pkt, k1f, out = sys.argv[1:4]
k1 = json.load(open(k1f))
rr = {r['tag']: r for r in json.load(open(f'{pkt}/RESPONSE_RESULTS.json'))['rows']}
res = {}
for tag in ('quiet3', 'quiet5'):
    d = np.load(f'{pkt}/{tag}.npz')
    t, bg, ks = d['times'], d['bg'], d['ks']
    amp, chg, sv, svk = d['amplitude_per_radian'], d['charge_fraction_per_radian'], d['selected_vectors'], d['selected_vector_ks']
    a = np.exp(bg[:, 15]); H = bg[:, 16]
    x = bg[:, 2] + 1j * bg[:, 3]; xd = bg[:, 8] + 1j * bg[:, 9]
    chi = np.sqrt(2) * a * x; chip = np.sqrt(2) * a * a * (xd + H * x)
    r = abs(chi)
    E = 0.5 * abs(chip) ** 2 + LAM / 4 * r ** 4
    L = np.imag(np.conj(chi) * chip)
    # (a) conformal time: trapezoid of 1/a in tau; accuracy check against Simpson-like refinement on halves
    eta = np.r_[0, cumulative_trapezoid(1 / a, t)]
    # radiation-law mapping for comparison
    C_end = float(a[-1] ** 2 * H[-1])
    late = t >= 80
    eta_rad = a / C_end
    off = eta[late] - eta_rad[late]
    # (d) C and tail
    Cs = a * a * H
    tau_star = t[-1] - a[-1] ** 2 / (2 * C_end)
    tail_dev = float(np.max(abs(a[late] / np.sqrt(2 * C_end * (t[late] - tau_star)) - 1)))
    # pericentres: local minima of r over tau >= 20
    idx = [i for i in range(1, len(t) - 1) if t[i] >= 20 and r[i] < r[i - 1] and r[i] <= r[i + 1]]
    peri = []
    for i in idx:
        # parabola in eta through (i-1, i, i+1) for r^2
        e3, r3 = eta[i - 1:i + 2], r[i - 1:i + 2] ** 2
        p = np.polyfit(e3 - e3[1], r3, 2)
        e_min = -p[1] / (2 * p[0]) if p[0] > 0 else 0.0
        peri.append(dict(i=i, tau=float(t[i]), eta=float(eta[i]), eta_min=float(e3[1] + e_min),
                         r_sample=float(r[i]), r2_min_parab=float(np.polyval(p, e_min))))
    # stroboscopic regression over tau in [200, 800]
    sel = [p for p in peri if 200 <= p['tau'] <= 800]
    ii = np.array([p['i'] for p in sel]); ee = eta[ii]
    rows = []
    flo = {round(q['k'], 4): q['mu'] for q in k1[tag + '_late']['KS']}
    flo10 = {round(q['k'], 4): q['mu'] for q in k1[tag + '_tau10']['KS']}
    for j, k in enumerate(ks):
        row = dict(k=float(k), floquet_late=flo[round(float(k), 4)], floquet_tau10=flo10[round(float(k), 4)])
        for name, arr in (('amp', amp[:, j]), ('chg', chg[:, j])):
            y = np.log(abs(arr[ii]) + 1e-300)
            pfit = np.polyfit(ee, y, 1)
            row[name + '_strobe_mu'] = float(pfit[0])
            row[name + '_strobe_rms'] = float(np.sqrt(np.mean((y - np.polyval(pfit, ee)) ** 2)))
        rows.append(row)
    # co-rotating X part of the saved vectors, strobed at pericentres
    vecrows = []
    for j, k in enumerate(svk):
        ph = np.exp(-1j * np.angle(chi[ii]))
        dq = (sv[ii, j, 0] + 1j * sv[ii, j, 1]) * ph
        dp = (sv[ii, j, 6] + 1j * sv[ii, j, 7]) * ph
        nrm = np.sqrt(abs(dq) ** 2 + abs(dp) ** 2 / max(k, 0.5) ** 2)
        pfit = np.polyfit(ee, np.log(nrm), 1)
        full = np.linalg.norm(sv[ii, j, :], axis=1)
        pfull = np.polyfit(ee, np.log(full), 1)
        vecrows.append(dict(k=float(k), X_corot_strobe_mu=float(pfit[0]), full12_strobe_mu=float(pfull[0]),
                            rms=float(np.sqrt(np.mean((np.log(nrm) - np.polyval(pfit, ee)) ** 2))),
                            floquet_late=flo.get(round(float(k), 4)), floquet_tau10=flo10.get(round(float(k), 4)),
                            growth_200_800=float(full[-1] / full[0])))
    # (c) effective growth from the unit samples with quadrature eta
    cr = {c['phase_seed_radians']: c for c in rr[tag]['crossings']}
    def eta_at(tau):
        return float(np.interp(tau, t, eta))
    mu_eff = {}
    for name in ('amplitude_unity', 'local_charge_unity'):
        e03 = eta_at(cr[0.03][name]['tau_bracket'][1]); e06 = eta_at(cr[1e-6][name]['tau_bracket'][1])
        e03l = eta_at(cr[0.03][name]['tau_bracket'][0]); e06l = eta_at(cr[1e-6][name]['tau_bracket'][0])
        mu_eff[name] = dict(mu=float(np.log(3e4) / (e06 - e03)),
                            spread=[float(np.log(3e4) / (e06 - e03l)), float(np.log(3e4) / (e06l - e03))],
                            d_eta=float(e06 - e03))
    # positions of the published crossing samples relative to the nearest pericentre
    pos = {}
    for s in (0.03, 1e-6):
        for name in ('amplitude_unity', 'local_charge_unity'):
            tau = cr[s][name]['tau_bracket'][1]
            i = int(np.argmin(abs(t - tau)))
            pn = min(peri, key=lambda p: abs(p['eta_min'] - eta[i]))
            pos[f'{s}_{name}'] = dict(tau=float(t[i]), a=float(a[i]), r=float(r[i]),
                                      r_over_true_rmin=float(r[i] / np.sqrt(pn['r2_min_parab'])),
                                      eta_offset_from_pericentre=float(eta[i] - pn['eta_min']),
                                      k=cr[s][name]['k'], sampled=cr[s][name]['sampled_value'])
    # radial period in tau near the crossings
    pt = np.array([p['tau'] for p in peri]); pe = np.array([p['eta_min'] for p in peri])
    res[tag] = dict(C_end=C_end, tau_star=tau_star, C_late_ptp_rel_tau_ge_80=float(np.ptp(Cs[late]) / C_end),
                    radiation_tail_max_rel_dev_tau_ge_80=tail_dev,
                    eta_quadrature_minus_a_over_C_after80=dict(mean=float(off.mean()), ptp=float(np.ptp(off))),
                    E_rel_ptp_tau_ge_20=float(np.ptp(E[t >= 20]) / E[t >= 20].mean()),
                    L_rel_ptp_tau_ge_20=float(np.ptp(L[t >= 20]) / L[t >= 20].mean()),
                    L_rel_ptp_tau_ge_200=float(np.ptp(L[t >= 200]) / L[t >= 200].mean()),
                    L_change_800_vs_64_80_mean=float(L[-1] / L[(t >= 64) & (t <= 80)].mean() - 1),
                    L_mean_700_800_vs_64_80_mean=float(L[t >= 700].mean() / L[(t >= 64) & (t <= 80)].mean() - 1),
                    L_mean_200_300_vs_64_80_mean=float(L[(t >= 200) & (t <= 300)].mean() / L[(t >= 64) & (t <= 80)].mean() - 1),
                    n_pericentres_200_800=len(sel), radial_period_eta_mean=float(np.mean(np.diff(pe[pt >= 200]))),
                    radial_period_tau_near=dict((str(T0), float(np.interp(T0, pt[1:], np.diff(pt)))) for T0 in (134, 147.5, 772, 782)),
                    strobe_fits=rows, vector_strobe=vecrows, mu_eff_quadrature_eta=mu_eff, crossing_positions=pos)
    print(tag, 'C', C_end, 'tail', tail_dev, 'eta-a/C', off.mean(), np.ptp(off), 'L ptp', res[tag]['L_rel_ptp_tau_ge_20'],
          'E ptp', res[tag]['E_rel_ptp_tau_ge_20'], 'L800', res[tag]['L_change_800_vs_64_80_mean'],
          'Lmean700-800', res[tag]['L_mean_700_800_vs_64_80_mean'])
    print('  mu_eff', mu_eff)
    for q in rows:
        if q['floquet_late'] > 1e-3 or q['k'] <= 0.1:
            print('  k %.2f strobe amp %.4f chg %.4f | Floquet late %.4f tau10 %.4f | rms %.3f' % (
                q['k'], q['amp_strobe_mu'], q['chg_strobe_mu'], q['floquet_late'], q['floquet_tau10'], q['amp_strobe_rms']))
    print('  vectors', vecrows)
    print('  positions', json.dumps(pos))
    print('  period tau', res[tag]['radial_period_tau_near'])
json.dump(res, open(out, 'w'), indent=1)
