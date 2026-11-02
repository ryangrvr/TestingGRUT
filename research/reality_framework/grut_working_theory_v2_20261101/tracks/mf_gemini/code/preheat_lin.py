"""mf_gemini: non-perturbative transfer from the quiet, dominant chi (m_chi = 3 H0) with radiatively
allowed couplings [C on a semi-analytic model; E for its validity].

Background: results/waterfall_bg.json (chi(t) about v after the mild waterfall; code units M_P = H0 = 1).

Channel B (bosonic, added field phi, model choice): 1/2 g^2 (chi - v)^2 phi^2, so m_phi vanishes at
chi = v (crossings 1.33 per H_end^-1).  Each crossing is a Landau-Zener scattering with
|beta|^2 = exp(-pi kappa^2), kappa^2 = (k/a)^2/(g |chi_dot|), and between crossings the mode is adiabatic.
Occupations follow the Kofman-Linde-Starobinsky map (hep-ph/9705347):
  n' = e^{-pi k2} + (1 + 2 e^{-pi k2}) n - 2 e^{-pi k2/2} sqrt(1 + e^{-pi k2}) sqrt(n(1+n)) sin(theta)
with theta the accumulated WKB phase; we report (i) the WKB phase from the background, (ii) random phases,
(iii) maximal growth (sin theta = -1), which bounds the transfer from above.  No back-reaction (valid
while rho_phi << rho_chi) and no rescattering; the products are assumed to decay instantly to the SM
("instant preheating"), which is optimistic.

Radiative cap on g: the phi loop V_CW = (g^2 (chi-v)^2)^2/(64 pi^2)[ln - 3/2] is not chi -> -chi even;
after the best polynomial counterterms {1, chi, ..., chi^4}, its residual must stay below 1% of V0 on
chi in [0, 1.5 v] (the same criterion as for the Yukawa in g1g2_checks.py).

Channel F (fermions, the proposal's y_N with the bare mass tuned so that m_N = y (chi - v)): Pauli-limited
Greene-Kofman production, n_k <= 1, energy per crossing f1 = y^{5/2} (Phi/m)^{1/2}/(2 pi^3) at most;
evaluated in g1g2_checks.py.

Run: python3 -B preheat_lin.py -> results/preheat_lin.json
"""
import json
import os
import numpy as np
from scipy.optimize import brentq

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
MP = 2.435e18


def bos_residual(g2, v, lam):
    x = np.linspace(0, 1.5, 4001)
    mm = g2 * (x - 1.0)**2 * v**2             # in M_P^2 (code: hub drops out of the ratio)
    with np.errstate(divide='ignore', invalid='ignore'):
        ff = np.where(mm > 0, mm**2 * np.log(mm / v**2), 0.0)
    Vb = ff / (64 * np.pi**2)
    V0 = lam * v**4 / 4
    A = np.vstack([x**p for p in range(5)]).T * v**4
    c, *_ = np.linalg.lstsq(A / V0, Vb / V0, rcond=None)
    return float(np.max(np.abs(Vb / V0 - (A / V0) @ c)))


def main():
    with open(os.path.join(RES, 'waterfall_bg.json')) as f:
        bg = json.load(f)
    T = bg['table']
    t = np.array(T['t_after_end_H0'])
    N = np.array(T['N_after_end'])
    c = np.array(T['chi'])
    dc = np.array(T['dchi'])
    h = np.array(T['h'])
    hend = bg['H_end_over_H0']
    H0 = 1.3e11
    hub = H0 / MP
    v = np.sqrt(8.0 / 3.0)
    lam_phys = 27.0 / 16.0 * hub**2          # dimensionless lambda
    a = np.exp(N)
    # crossings chi = v
    s = np.sign(c - v)
    idx = np.where(s[1:] != s[:-1])[0]
    tc = []
    for i in idx:
        f = (v - c[i]) / (c[i + 1] - c[i])
        tc.append(dict(t=t[i] + f * (t[i + 1] - t[i]), a=a[i] + f * (a[i + 1] - a[i]),
                       dchi=abs(dc[i] + f * (dc[i + 1] - dc[i]))))
    # radiative cap
    g2cap = float(np.exp(brentq(lambda l: np.log(bos_residual(np.exp(l), v, lam_phys)) - np.log(0.01),
                                np.log(1e-14), np.log(1e-2))))
    out = dict(g2_cap_1pct=g2cap, crossings_first=[dict(t_Hend=x['t'] * hend, a=x['a'], dchi=x['dchi']) for x in tc[:6]])
    rho_chi = 3 * h**2                       # code units (M_P^2 H0^2)
    rng = np.random.default_rng(7)
    res = {}
    for lab, g2 in [('cap', g2cap), ('cap x 10', 10 * g2cap), ('cap x 100', 100 * g2cap)]:
        g = np.sqrt(g2)
        # in H0 units: fields chi_H = chi/hub; g |chi_dot_H| = g dchi/hub  (H0^2)
        kstar0 = np.sqrt(g * tc[0]['dchi'] / hub) * tc[0]['a']     # comoving k at the first crossing scale
        ks = np.linspace(0.01, 4.0, 400) * kstar0
        dk = ks[1] - ks[0]
        modes = {}
        for mode in ['wkb', 'random', 'max']:
            n = np.zeros_like(ks)
            hist = []
            for j, x in enumerate(tc):
                k2 = (ks / x['a'])**2 / (g * x['dchi'] / hub)
                e = np.exp(-np.pi * k2)
                if mode == 'max':
                    sn = -1.0
                elif mode == 'random':
                    sn = np.sin(rng.uniform(0, 2 * np.pi, size=ks.size))
                else:
                    # WKB phase between this crossing and the next: theta = 2 int omega dt (mass-dominated)
                    if j + 1 < len(tc):
                        sel = (t >= x['t']) & (t <= tc[j + 1]['t'])
                        mphi = g * np.abs(c[sel] - v) / hub
                        om = np.sqrt((ks[:, None] / a[sel][None, :])**2 + mphi[None, :]**2)
                        theta = 2 * np.trapezoid(om, t[sel], axis=1)
                    else:
                        theta = np.zeros_like(ks)
                    sn = np.sin(theta)
                n = e + (1 + 2 * e) * n - 2 * np.sqrt(e) * np.sqrt(1 + e) * np.sqrt(n * (1 + n)) * sn
                n = np.clip(n, 0, None)
                # energy at the next turning point (m_phi max): omega ~ g |chi - v|max / hub
                if j + 1 < len(tc):
                    sel = (t >= x['t']) & (t <= tc[j + 1]['t'])
                    i_t = np.argmax(np.abs(c[sel] - v))
                    tt = t[sel][i_t]
                    at = a[sel][i_t]
                    mphi = g * abs(c[sel][i_t] - v) / hub
                    om = np.sqrt((ks / at)**2 + mphi**2)
                    rho_phi = np.sum(ks**2 * om * n) * dk / (2 * np.pi**2 * at**3)     # H0^4
                    rho_tot = 3 * np.interp(tt, t, h)**2 / hub**2                       # H0^4 (M_P^2 H0^2 / hub^2)
                    hist.append(dict(t_Hend=float(tt * hend), frac=float(rho_phi / rho_tot)))
            modes[mode] = hist
        def t_at(hist, f):
            for x in hist:
                if x['frac'] >= f:
                    return x['t_Hend']
            return None
        res[lab] = dict(g2=g2, q_first=float(g2 * (0.147 * v / hub)**2 / (4 * 9.0)),
                        first_crossing_frac=modes['max'][0]['frac'] if modes['max'] else None,
                        t_frac_0p1_Hend={m: t_at(modes[m], 0.1) for m in modes},
                        t_frac_0p5_Hend={m: t_at(modes[m], 0.5) for m in modes},
                        frac_at_2Hend={m: (max([x['frac'] for x in modes[m] if x['t_Hend'] <= 2.0], default=0.0))
                                       for m in modes},
                        history_max=modes['max'][:12], history_wkb=modes['wkb'][:12])
    out['bosonic'] = res
    out['K_if_complete_at_t'] = bg['K_if_sudden_decay_at_t_Hend']
    with open(os.path.join(RES, 'preheat_lin.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: (v if k != 'bosonic' else {kk: {x: vv[x] for x in vv if not x.startswith('history')}
                                                      for kk, vv in v.items()}) for k, v in out.items()}, indent=1))


if __name__ == '__main__':
    main()
