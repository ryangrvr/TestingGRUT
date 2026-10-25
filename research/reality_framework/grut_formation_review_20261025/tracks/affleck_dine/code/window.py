#!/usr/bin/env python3
"""Generation-versus-survival window for one Z_N-allowed coefficient c_d (own code).

Reads results/ad_scan.json (ad_condensate.py).  Every input is listed in INPUTS.
Yield formulas (derived in NOTES.md, a=1 at release, n a^3 = q X0^2 w0, w0 = sqrt(2 lX) X0,
h = H_rel/w0, X0 = H_rel/(h sqrt(2 lX))):
  release into radiation:        Y = q h^-3 (H_rel/M)^(3/2) / (2 lX D(g*)),
                                 D = (2 pi^2/45) g* (90/(pi^2 g*))^(3/4)
  release into inflaton matter:  Y = q h^-3 T_R H_rel / (8 lX M^2)       (T_R < T(H_rel))
  eta_d = (c_d/lX) (X0/M)^(d-4).
Bounds on c_d at fixed inputs:
  c_req   : smallest c giving Y = Y_Q at the best scanned phase (and at phi0 = pi/2)
  c_cap   : boundedness/no-runaway of lX|X|^4 + 2|lambda_d| R^d cos at the largest amplitude
            of the history, eta(X_max) <= eta_safe (largest scanned eta with no runaway at any phase)
  c_heavy : angular direction heavier than 3H/2 during inflation (isocurvature-safe DRT branch):
            m_theta^2/H^2 = d^2 eta_inf c_H / 2 > 9/4
  c_surv  : Codex 8986d34 full-response late-survival bound [Q]
"""
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUTS = dict(
    lambda_X=1e-3, g_X=0.1, M_reduced_GeV=2.435e18, Y_Q=4.8294009076988175e-9,
    g_star_AD_epoch=106.75,
    H_inf_max_GeV=4.7e13,  # r<0.036 (BICEP/Keck 2021), P_R=2.1e-9: H < pi M sqrt(r P_R/2)
    c_surv={'6': 765399729417561.2, '5': 0.00012226814094523126},  # Codex CHARGE_LEAKAGE_RESULTS rows (R=40, tol 2e-8)
    beta_iso_max=0.038, P_R=2.1e-9,
    N_CMB_efolds=50.0,  # e-folds between CMB-scale horizon exit and end of inflation
)
lX = INPUTS['lambda_X']; M = INPUTS['M_reduced_GeV']; YQ = INPUTS['Y_Q']; gs = INPUTS['g_star_AD_epoch']
D = (2 * np.pi ** 2 / 45) * gs * (90 / (np.pi ** 2 * gs)) ** 0.75


def T_of_H(H):
    return (90 / (np.pi ** 2 * gs)) ** 0.25 * np.sqrt(H * M)


def tables(scan, d):
    """Per case: eta grid, best |q| over phases, |q| at pi/2, eta_safe, linear coefficient."""
    out = {}
    rows = [r for r in scan['rows'] if r['d'] == d]
    for case in sorted(set(r['case'] for r in rows)):
        rr = [r for r in rows if r['case'] == case]
        etas = sorted(set(r['eta'] for r in rr))
        best, half, safe_ok, eps_best = [], [], [], []
        for e in etas:
            re_ = [r for r in rr if r['eta'] == e]
            ok = [r for r in re_ if not r['ran_away']]
            best.append(max(abs(r['q_inf']) for r in ok) if ok else 0.0)
            eps_best.append(max(r['epsilon'] for r in ok) if ok else 0.0)
            hh = [r for r in re_ if abs(r['phi0'] - np.pi / 2) < 1e-9 and not r['ran_away']]
            half.append(abs(hh[0]['q_inf']) if hh else 0.0)
            safe_ok.append(all(not r['ran_away'] for r in re_))
        etas = np.array(etas); best = np.array(best); half = np.array(half)
        safe = etas[np.array(safe_ok)]
        eta_safe = float(max(e for e, s in zip(etas, safe_ok) if s and all(safe_ok[:list(etas).index(e) + 1])))
        lin = half[0] / etas[0]
        out[case] = dict(eta=etas.tolist(), q_best=best.tolist(), q_half=half.tolist(), eps_best=eps_best,
                         eta_safe=eta_safe, linear_coeff_q_over_eta=float(lin), h=rr[0]['h'], era=rr[0]['era'])
    return out


def eta_for_q(tab, qneed, which='q_best'):
    """Smallest eta <= eta_safe with |q(eta)| >= qneed (log-log interpolation; linear below grid)."""
    e = np.array(tab['eta']); q = np.array(tab[which])
    m = e <= tab['eta_safe'] + 1e-12
    e, q = e[m], q[m]
    if qneed <= q[0]:
        return qneed / (q[0] / e[0])
    # running max so the function is monotone
    qm = np.maximum.accumulate(q)
    if qneed > qm[-1]:
        return None
    i = int(np.searchsorted(qm, qneed))
    if qm[i] == qm[i - 1]:
        return float(e[i])
    t = (np.log(qneed) - np.log(qm[i - 1])) / (np.log(qm[i]) - np.log(qm[i - 1]))
    return float(np.exp(np.log(e[i - 1]) + t * (np.log(e[i]) - np.log(e[i - 1]))))


def c_of_eta(eta, X0, d):
    return eta * lX * (M / X0) ** (d - 4)


def scenario(tab, d, H_rel, era, T_R=None, H_inf=None, cH=None, X_max=None, phase='q_best'):
    h = tab['h']
    w0 = H_rel / h
    X0 = w0 / np.sqrt(2 * lX)
    if era == 'RD':
        F = (H_rel / M) ** 1.5 / (2 * lX * D)
    else:
        F = T_R * H_rel / (8 * lX * M ** 2)
    qneed = YQ * h ** 3 / F
    eta_req = eta_for_q(tab, qneed, phase)
    Xm = X_max if X_max is not None else X0
    c_cap = c_of_eta(tab['eta_safe'], Xm, d)
    c_req = c_of_eta(eta_req, X0, d) if eta_req is not None else None
    res = dict(d=d, era=era, H_rel=H_rel, T_R=T_R, H_inf=H_inf, c_H=cH, h=h, X0_GeV=X0, X_max_GeV=Xm,
               q_needed=qneed, eta_req=eta_req, c_req=c_req, c_cap=c_cap, c_surv=INPUTS['c_surv'][str(d)])
    if cH is not None and H_inf is not None:
        X_inf = np.sqrt(cH / (2 * lX)) * H_inf
        res['c_heavy_angle'] = c_of_eta(9 / (2 * d * d * cH), X_inf, d)  # m_theta > 3H/2 (conservative)
        # isocurvature-based: light-field S ~ d*dtheta (cot phi0 = 1) decays as exp[-(3/2-nu)N]
        S_light = d * np.sqrt(lX / cH) / (2 * np.pi)
        S_max = np.sqrt(INPUTS['beta_iso_max'] / (1 - INPUTS['beta_iso_max']) * INPUTS['P_R'])
        L = np.log(S_light / S_max) / INPUTS['N_CMB_efolds']
        m2 = 9 / 4 - (1.5 - L) ** 2 if L < 1.5 else 9 / 4
        res['m2_over_H2_needed_iso'] = m2
        res['c_heavy_iso'] = c_of_eta(2 * m2 / (d * d * cH), X_inf, d)
        res['radial_m2_over_H2'] = 2 * cH
        res['radial_iso_ok'] = bool(2 * cH >= m2)
    upper = min(c_cap, res['c_surv'])
    lower = c_req if c_req is not None else np.inf
    res['window'] = [lower, upper] if lower <= upper else None
    res['window_decades'] = float(np.log10(upper / lower)) if lower <= upper else None
    if res.get('c_heavy_angle') is not None:
        lo2 = max(lower, res['c_heavy_angle'])
        res['window_isocurvature_safe'] = [lo2, upper] if lo2 <= upper else None
        res['window_isocurvature_safe_decades'] = float(np.log10(upper / lo2)) if lo2 <= upper else None
        lo3 = max(lower, res['c_heavy_iso'])
        res['window_iso_relaxed'] = [lo3, upper] if lo3 <= upper else None
        res['window_iso_relaxed_decades'] = float(np.log10(upper / lo3)) if lo3 <= upper else None
    res['binding_upper'] = 'boundedness' if c_cap < res['c_surv'] else 'late survival'
    return res


def main():
    scan = json.loads((ROOT / 'results' / 'ad_scan.json').read_text())
    tabs = {d: tables(scan, d) for d in (6, 5)}
    out = dict(inputs=INPUTS, D_gstar=D, tables={str(d): t for d, t in tabs.items()}, scenarios=[])
    Hs = [4.7e13, 1e13, 1e12, 1e11, 1e10]
    for d in (6, 5):
        t = tabs[d]
        for cH, key in [(1, 'cH1'), (10, 'cH10'), (100, 'cH100')]:
            for H_inf in Hs:
                # W2: Hubble-induced mass held during inflation, switched off at its end, instant reheating
                r = scenario(t[f'hubble_RD_{key}_b0'], d, H_inf, 'RD', H_inf=H_inf, cH=cH)
                r['scenario'] = 'W2 end-of-inflation release, instant reheating'; out['scenarios'].append(r)
                # W1: same release into inflaton matter domination, T_R supplied
                for T_R in (1e9, 1e12, 1e14):
                    if T_R >= T_of_H(H_inf):
                        continue
                    r = scenario(t[f'hubble_MD_{key}_b0'], d, H_inf, 'MD', T_R=T_R, H_inf=H_inf, cH=cH)
                    r['scenario'] = 'W1 end-of-inflation release into inflaton MD'; out['scenarios'].append(r)
            # W3: Hubble-induced mass persists until reheating (tracking X ~ H), release into RD at H_R
            if key == 'cH1':
                for H_inf in (4.7e13, 1e13):
                    for T_R in (1e15, 1e14, 1e12, 1e9):
                        H_R = np.sqrt(np.pi ** 2 * gs / 90) * T_R ** 2 / M
                        if H_R >= H_inf:
                            continue
                        X_inf = np.sqrt(cH / (2 * lX)) * H_inf
                        r = scenario(t['hubble_RD_cH1_b1.5'], d, H_R, 'RD', T_R=T_R, H_inf=H_inf, cH=cH, X_max=X_inf)
                        r['scenario'] = 'W3 tracking until reheating, release into RD'; out['scenarios'].append(r)
        # W4: frozen light field (no Hubble-induced mass), w0 <= H_inf; isocurvature listed separately
        for w0 in (4.7e13, 1e13, 1e12, 1e11):
            tt = t['frozen_RD']
            H_rel = tt['h'] * w0
            r = scenario(tt, d, H_rel, 'RD')
            X0 = w0 / np.sqrt(2 * lX)
            r['scenario'] = 'W4 frozen light field, onset in RD (H_rel is a dummy >> w0)'
            r['w0'] = w0
            r['radial_isocurvature_dY_over_Y'] = (d - 0.5) * w0 / (2 * np.pi * np.sqrt(2) * X0)
            out['scenarios'].append(r)
    # isocurvature of the light-angle (random phase) branch at the Hubble-induced amplitude
    S_max = np.sqrt(INPUTS['beta_iso_max'] / (1 - INPUTS['beta_iso_max']) * INPUTS['P_R'])
    out['isocurvature'] = dict(S_max=S_max, light_angle={
        f'd{d}_cH{cH}': d * np.sqrt(lX / cH) / (2 * np.pi) for d in (5, 6) for cH in (1, 10, 100)},
        note='S = d*dtheta*cot(phi0), dtheta = H/(2 pi sqrt2 R) = sqrt(lX/cH)/(2 pi); value printed is S/cot(phi0)')
    # minimal H_inf for which a window exists (generation within boundedness and survival; optionally iso-safe)
    hmin = []
    Hgrid = np.geomspace(1e8, INPUTS['H_inf_max_GeV'], 400)
    for d in (6, 5):
        for cH, key in [(1, 'cH1'), (10, 'cH10'), (100, 'cH100')]:
            for era, T_R in [('RD', None), ('MD', 1e9), ('MD', 1e12), ('MD', 1e14)]:
                tab = tabs[d][f'hubble_{era}_{key}_b0']
                first, first_iso = None, None
                for H in Hgrid:
                    if era == 'MD' and T_R >= T_of_H(H):
                        continue
                    r = scenario(tab, d, H, era, T_R=T_R, H_inf=H, cH=cH)
                    if first is None and r['window'] is not None:
                        first = H
                    if first_iso is None and r.get('window_iso_relaxed') is not None:
                        first_iso = H
                hmin.append(dict(d=d, c_H=cH, era=era, T_R=T_R, H_min_window=first, H_min_iso_relaxed=first_iso))
                print('Hmin d', d, 'cH', cH, era, T_R, first, first_iso)
    out['H_min'] = hmin
    (ROOT / 'results' / 'window.json').write_text(json.dumps(out, indent=1, default=float))
    for r in out['scenarios']:
        print(r['scenario'][:3], 'd', r['d'], 'cH', r['c_H'], 'Hrel %.2e' % r['H_rel'], 'TR', r['T_R'],
              'c_req', None if r['c_req'] is None else '%.3g' % r['c_req'], 'c_cap %.3g' % r['c_cap'],
              'c_heavy', '%.3g' % r['c_heavy_angle'] if r.get('c_heavy_angle') else '-',
              'win_dec', r['window_decades'], 'iso-safe dec', r.get('window_isocurvature_safe_decades'),
              'c_heavy_iso', '%.3g' % r['c_heavy_iso'] if r.get('c_heavy_iso') else '-', 'iso-relaxed', r.get('window_iso_relaxed'))


if __name__ == '__main__':
    main()
