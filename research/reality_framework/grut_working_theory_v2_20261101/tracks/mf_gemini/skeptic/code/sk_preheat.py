"""mf_gemini SKEPTIC: the bosonic crossing-map channel with a CONSISTENT decay of the produced quanta [C on a model].

The track's preheat_lin.py takes the maximal Bose growth of phi (no decay between crossings) AND counts the phi
energy as SM radiation ("instant decay").  Those two are exclusive: a decay between crossings removes the
occupation that the next crossing amplifies.  Here phi decays at a constant width Gamma_phi (a model choice) into
radiation, so
  n_k -> n_k exp(-Gamma_phi dt) between crossings,  d(rho_r a^4) = a^4 Gamma_phi rho_phi dt,
and the crossing map (Kofman-Linde-Starobinsky hep-ph/9705347, maximal-growth phase sin theta = -1 as an upper
bound) acts on what is left.  Background: own integration of chi after eps_H = 1 (no back-reaction), code units
M_P = H0 = 1.  Coupling 1/2 g^2 (chi - v)^2 phi^2 at the track's loop cap g^2 = 4.32e-7 and x10, x100.
Reported: the radiation fraction F_r(t) and the times at which F_r = 0.1, 0.5; for Gamma_phi = 0 the phi fraction
(the track's quantity).
Run: python3 -B sk_preheat.py -> ../results/sk_preheat.json
"""
import json
import os
import numpy as np
from scipy.integrate import solve_ivp

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
MP = 2.435e18
HUB = 1.3e11 / MP
V2 = 8.0 / 3.0
V = np.sqrt(V2)
LAM = 27.0 / 16.0
G2CAP = 4.32224237549203e-07


def background(T=45.0):
    Vp = lambda c: 0.25 * LAM * (c * c - V2)**2
    dVp = lambda c: LAM * (c * c - V2) * c

    def f(t, y):
        H = np.sqrt((0.5 * y[2]**2 + Vp(y[1])) / 3.0)
        return [H, y[2], -3 * H * y[2] - dVp(y[1])]
    s = (-3 + np.sqrt(9 + 4 * LAM * V2)) / 2

    def ev(t, y):
        H = np.sqrt((0.5 * y[2]**2 + Vp(y[1])) / 3.0)
        return 0.5 * y[2]**2 / H**2 - 1
    ev.terminal = True
    ev.direction = 1
    s1 = solve_ivp(f, (0, 400), [0, 1e-6, s * 1e-6], method='LSODA', rtol=1e-11, atol=1e-16, events=ev)
    ye = s1.y_events[0][0]
    He = np.sqrt((0.5 * ye[2]**2 + Vp(ye[1])) / 3.0)
    s2 = solve_ivp(f, (0, T / He), [0.0, ye[1], ye[2]], method='LSODA', rtol=1e-11, atol=1e-16, max_step=0.002,
                   dense_output=True)
    t = np.arange(0, T / He, 0.002)
    N, c, dc = s2.sol(t)
    H = np.sqrt((0.5 * dc**2 + Vp(c)) / 3.0)
    return dict(t=t, a=np.exp(N), c=c, dc=dc, H=H, He=He)


def run(bg, g2, Gphi_Hend, phase='max', seed=5):
    rng = np.random.default_rng(seed)
    g = np.sqrt(g2)
    t, a, c, dc, H, He = bg['t'], bg['a'], bg['c'], bg['dc'], bg['H'], bg['He']
    Gphi = Gphi_Hend * He
    sg = np.sign(c - V)
    idx = np.where(sg[1:] != sg[:-1])[0]
    k0 = np.sqrt(g * abs(dc[idx[0]]) / HUB) * a[idx[0]]
    ks = np.exp(np.linspace(np.log(1e-3), np.log(5.0), 400)) * k0
    w_k = np.gradient(ks)                           # dk weights
    n = np.zeros_like(ks)
    rr_a4 = 0.0                                     # rho_r a^4 (H0^4 units)
    rho_tot = 3 * H**2 / HUB**2
    hist = []
    j_next = 0
    for i in range(len(t) - 1):
        if j_next < len(idx) and i == idx[j_next]:
            k2 = (ks / a[i])**2 / (g * abs(dc[i]) / HUB)
            e = np.exp(-np.pi * k2)
            sn = -1.0 if phase == 'max' else np.sin(rng.uniform(0, 2 * np.pi, ks.size))
            n = e + (1 + 2 * e) * n - 2 * np.sqrt(e * (1 + e)) * np.sqrt(n * (1 + n)) * sn
            n = np.clip(n, 0, None)
            j_next += 1
        dt = t[i + 1] - t[i]
        mphi = g * abs(c[i] - V) / HUB
        om = np.sqrt((ks / a[i])**2 + mphi**2)
        rho_phi = np.sum(ks**2 * om * n * w_k) / (2 * np.pi**2 * a[i]**3)
        if Gphi > 0:
            rr_a4 += a[i]**4 * Gphi * rho_phi * dt
            n = n * np.exp(-Gphi * dt)
        if i % 50 == 0:
            hist.append((float(t[i] * He), float(rho_phi / rho_tot[i]), float(rr_a4 / a[i]**4 / rho_tot[i])))
    hist = np.array(hist)

    def t_at(col, f):
        w = np.where(hist[:, col] >= f)[0]
        return float(hist[w[0], 0]) if w.size else None
    col = 2 if Gphi > 0 else 1
    return dict(g2=g2, Gamma_phi_over_Hend=Gphi_Hend, phase=phase, quantity='radiation' if Gphi > 0 else 'phi (track)',
                t_0p1_Hend=t_at(col, 0.1), t_0p5_Hend=t_at(col, 0.5),
                max_frac_by_2Hend=float(hist[hist[:, 0] <= 2.0, col].max()),
                max_frac_by_10Hend=float(hist[hist[:, 0] <= 10.0, col].max()),
                max_frac_by_40Hend=float(hist[:, col].max()))


def main():
    bg = background()
    out = dict(H_end_over_H0=float(bg['He']), runs=[])
    for mult in [1, 10, 100]:
        for G in [0.0, 0.3, 1.0, 3.0, 10.0]:
            r = run(bg, mult * G2CAP, G)
            r['cap_multiple'] = mult
            out['runs'].append(r)
            print(r, flush=True)
    for G in [0.0, 0.3, 1.0]:
        for mult in [1, 100]:
            r = run(bg, mult * G2CAP, G, phase='random')
            r['cap_multiple'] = mult
            out['runs'].append(r)
            print(r, flush=True)
    with open(os.path.join(RES, 'sk_preheat.json'), 'w') as fh:
        json.dump(out, fh, indent=1)


if __name__ == '__main__':
    main()
