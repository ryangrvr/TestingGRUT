#!/usr/bin/env python3
"""e2e_chain_skeptic (Claude): the spectator's inflationary potential including the hidden sector's tree-level response.

e2e_chain's naturalness estimate (points.py) puts s in Starobinsky-Yokoyama equilibrium with V = m^2 s^2/2 only, m^2 the
residual inflationary masses (4.5e-8 .. 3e-5 H^2).  But the holder -(g_h e^{i alpha} s Sigma* Y^2 + h.c.) makes the held
hidden state depend on s: with the heavy hidden fields at their stationary point, s feels
        V_eff(s) = m^2 s^2/2 [+ lambda_s s^4] + V_min(|u(s)|),   u = mu + c G,  c = |g_h| s/|G|,  G = 10i - mu,
V_min(|u|) = the held stationary energy (R = 12, constant-H inflation), from dq_heat skeptic's exact solver fc.stationary
(read-only, via mf_ideal's mfi_common).  mf_hybrid quoted the slope of this term at c = 1 as a "tadpole" (0.6 H^3); here
the whole function is computed.  Outputs (results/sk_natural.json):
  * V_min(c) on a grid, its slope and curvature at M2-X / M2-Q, the induced m^2_eff = V_min''/s_u^2 (s_u = s at c = 1 =
    |G|/|g_h|, g_h fixed by the isocurvature floor or by the bench s_*), drift per e-fold, the runaway time;
  * the large-|u| asymptote V_min/|u|^4 (boundedness: a quartic lambda_s s^4 >= A g_h^4 is needed);
  * with V_s = lambda_s s^4 (Lambda = lambda_s/g_h^4 in units where s is measured through |u| = g_h s ... see code) the
    position of the positive-c local minimum vs Lambda, and whether it is the global minimum (negative-c side).
SY reference: Starobinsky & Yokoyama, arXiv:astro-ph/9407016: P(s) ~ exp(-8 pi^2 V/(3 H^4)).
Labels: [C] numbers on the model; reading [A].  All in H_i units (H = 1 during inflation).
"""
import json, math
import numpy as np
from scipy.optimize import brentq
import sk_common as S

mc = S.E.mc
fc, P, MU = mc.fc, mc.P05, mc.MU
G = 10j - MU
AG = abs(G)


def Vmin_grid(cs):
    out, prev, cprev = [], None, None
    for c in cs:
        u = MU + c * G
        st = fc.stationary(u, P, R=12.0, guess=prev if (cprev is not None and abs(c - cprev) < 0.11) else None)
        if st['grad'] > 1e-6:
            st = fc.stationary(u, P, R=12.0)
        out.append(dict(c=float(c), abs_u=abs(u), V=st['V'], grad=st['grad']))
        prev, cprev = st['z'], c
    return out


def main():
    cs = np.round(np.r_[np.arange(-3.0, 3.0001, 0.025)], 4)
    rows = Vmin_grid(cs)
    c = np.array([r['c'] for r in rows]); V = np.array([r['V'] for r in rows])
    dV = np.gradient(V, c); d2V = np.gradient(dV, c)
    res = dict(method='fc.stationary at R = 12 along u = mu + c(10i - mu); V_min depends on |u| only [C, checked at |u| = 10 for 4 phases]',
               grid=[dict(c=float(a), V=float(b)) for a, b in zip(c, V)])
    # asymptote
    asy = []
    for cc in (5, 10, 20, 40):
        st = fc.stationary(MU + cc * G, P, R=12.0)
        asy.append(dict(c=cc, abs_u=abs(MU + cc * G), V_over_u4=st['V'] / abs(MU + cc * G) ** 4))
    res['asymptote'] = asy
    A = -asy[-1]['V_over_u4']
    res['boundedness'] = dict(A=A, statement='V_min -> -A |u|^4 = -A g_h^4 s^4 at large s: unbounded below along s unless lambda_s s^4 with lambda_s > A g_h^4 is added')
    # points
    pts = {}
    for name, c0, sfl, H in (('M2X', 1.55, 1.1327e16, 1.1229e11), ('M2Q', 1.45, 1.3601e16, 9.8864e10)):
        Vp = float(np.interp(c0, c, dV)); Vpp = float(np.interp(c0, c, d2V)); V0 = float(np.interp(c0, c, V))
        loc = {}
        for lab, s_star_GeV in (('s_at_isocurvature_floor', sfl), ('bench_s_star_7.3e15', 7.305e15)):
            s_star = s_star_GeV / H                     # in H_i
            gh = c0 * AG / s_star                       # |g_h| s_* = c0 |G|
            su = AG / gh                                # s at c = 1
            m2eff = Vpp / su ** 2
            drift = -Vp / su / 3.0                      # ds/dN = -V'(s)/(3H^2), H = 1
            loc[lab] = dict(s_star_over_H=s_star, g_h=gh, s_u_over_H=su, m2_eff_over_H2=m2eff, drift_s_per_efold_over_H=drift,
                            drift_c_per_efold=drift / su, dc_over_c_in_50_efolds=50 * drift / su / c0,
                            runaway_efolds_c_over_cdot=c0 / abs(drift / su),
                            lambda_s_min_boundedness=A * gh ** 4,
                            m2_needed_to_balance_slope=-Vp / (su ** 2 * c0),
                            balance_is_a_maximum_for_any_m2=bool(Vpp < 0),
                            ratio_to_residual_masses={'vs_3e-5': abs(m2eff) / 3e-5, 'vs_2.5e-7': abs(m2eff) / 2.5e-7, 'vs_4.5e-8': abs(m2eff) / 4.5e-8})
        n_loc = math.log(float(np.interp(c0 + 0.05, c, -V)) / float(np.interp(c0 - 0.05, c, -V))) / math.log((c0 + 0.05) / (c0 - 0.05))
        pts[name] = dict(c0=c0, V_min=V0, dV_dc=Vp, d2V_dc2=Vpp, local_power_law_n=n_loc, by_s_star=loc,
                         SY_exponent_8pi2_dV_over_3H4_for_dc_0p01=8 * math.pi ** 2 * abs(Vp) * 0.01 / 3)
    res['points'] = pts
    # quartic-stabilised attractor: V_eff(c) = L |G|^4 c^4 + V_min(c)  (L = lambda_s/g_h^4; s = c |G|/g_h)
    att = []
    for L in np.r_[np.arange(7.5, 40, 0.5), np.arange(40, 400, 10)]:
        Ve = L * AG ** 4 * c ** 4 + V
        dVe = np.gradient(Ve, c)
        mins = [float(c[i]) for i in range(1, len(c) - 1) if Ve[i] < Ve[i - 1] and Ve[i] < Ve[i + 1]]
        pos = [m for m in mins if m > 0.4]; neg = [m for m in mins if m < 0.2]
        vp = float(np.interp(pos[0], c, Ve)) if pos else None
        vn = float(np.interp(neg[0], c, Ve)) if neg else None
        att.append(dict(Lambda=float(L), c_min_pos=pos[0] if pos else None, c_min_neg=neg[0] if neg else None, V_pos=vp, V_neg=vn,
                        global_is_pos=(vp is not None and (vn is None or vp < vn))))
    res['quartic_attractor'] = att
    win = [a for a in att if a['c_min_pos'] and 1.45 <= a['c_min_pos'] <= 1.60]
    res['quartic_attractor_window_M2X'] = dict(Lambda_range=[min(a['Lambda'] for a in win), max(a['Lambda'] for a in win)] if win else None,
                                              global_in_window=[a['global_is_pos'] for a in win])
    S.dump('sk_natural.json', res)
    for nm, p in pts.items():
        print(nm, json.dumps({k: v for k, v in p.items() if k != 'by_s_star'}, default=float))
        for lab, l in p['by_s_star'].items():
            print('   ', lab, json.dumps(l, default=float))
    print('asymptote', asy)
    for a in att[::4]:
        print(a)
    print(res['quartic_attractor_window_M2X'])


if __name__ == '__main__':
    main()
