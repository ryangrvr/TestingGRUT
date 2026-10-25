#!/usr/bin/env python3
"""When and at what scale does a (near-circular) AD condensate of X fragment? (own code)

Full two-field linearization of the homogeneous rotating condensate
X = R e^{i w t}, phi = phi0(R) static, in U = lam/4 (phi^2-v^2)^2 + g^2 phi^2 |X|^2 + lX |X|^4.
  R > Xc = v sqrt(lam/(2 g^2)):  phi0 = 0,  w^2 = 2 lX R^2   (phi heavy, condensate stable)
  R < Xc:  phi0^2 = v^2 - 2 g^2 R^2/lam,  w^2 = g^2 phi0^2 + 2 lX R^2
Perturbations delta X = (a + i b) e^{i w t}, delta phi = p, ~ e^{i k x + gamma t}:
  det [[A, -2 w gamma, 2 g^2 phi0 R], [2 w gamma, B, 0], [4 g^2 phi0 R, 0, C]] = 0,
  A = k^2+gamma^2+4 lX R^2, B = k^2+gamma^2, C = k^2+gamma^2+m_p^2,
  m_p^2 = lam(3 phi0^2 - v^2) + 2 g^2 R^2.
(det = A B C + 4 w^2 gamma^2 C - 8 g^4 phi0^2 R^2 B.)  Single-field check: heavy p gives the
band k^2 < 4 a R^2 with a = g^4/lam - lX.
Outputs: gamma_max(k), k_peak, Q_cell = n (2 pi/k_peak)^3 with n = 2 w R^2 (circular, eps=1),
the onset temperature for net yield Y, and a slow-quench onset estimate: R = Xc (1 - delta),
delta = H t, growth exponent N(delta) = int gamma_max/H d delta.  Conditional diagnostic only.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import brentq

OUT = Path(__file__).resolve().parents[1] / 'results'
M = 2.435e18
YQ = 4.8294009076988175e-9
POINTS = {
    'C1-HP-03/v1 (Q=200000)': dict(g=0.1, lam=0.03, lX=0.001, v=0.924089637375362),
    'old benchmark (Q=39944)': dict(g=0.1, lam=0.01, lX=0.001, v=1.0),
}


def gstar(T):  # coarse SM g* (= g*s) table, adequate for O(10%) temperature estimates
    tab = [(1e-3, 10.75), (0.1, 10.75), (0.15, 17.25), (0.3, 61.75), (1.0, 72.25), (5.0, 86.25),
           (100.0, 96.25), (200.0, 106.75), (1e19, 106.75)]
    Ts, gs = zip(*tab)
    return float(np.interp(np.log(T), np.log(Ts), gs))


def background(P, R):
    g, lam, lX, v = P['g'], P['lam'], P['lX'], P['v']
    Xc2 = lam * v * v / (2 * g * g)
    phi0sq = max(v * v - 2 * g * g * R * R / lam, 0.0)
    w2 = g * g * phi0sq + 2 * lX * R * R
    mp2 = lam * (3 * phi0sq - v * v) + 2 * g * g * R * R
    return phi0sq, w2, mp2, np.sqrt(Xc2)


def growth(P, R, k):
    g, lX = P['g'], P['lX']
    phi0sq, w2, mp2, _ = background(P, R)
    c = 8 * g ** 4 * phi0sq * R * R
    # det as cubic in x = gamma^2: A=x+k2+4lX R^2, B=x+k2, C=x+k2+mp2
    k2 = k * k
    a0 = k2 + 4 * lX * R * R
    c0 = k2 + mp2
    pA = np.poly1d([1, a0]); pB = np.poly1d([1, k2]); pC = np.poly1d([1, c0])
    det = pA * pB * pC + np.poly1d([4 * w2, 0]) * pC - c * pB
    roots = det.roots
    pos = [r.real for r in roots if abs(r.imag) < 1e-12 * max(1, abs(r.real)) and r.real > 0]
    return np.sqrt(max(pos)) if pos else 0.0


def peak(P, R):
    phi0sq, w2, mp2, Xc = background(P, R)
    kmax_guess = 2 * np.sqrt(abs(P['g'] ** 4 / P['lam'] - P['lX'])) * R
    ks = np.geomspace(1e-9 * P['v'], 3 * kmax_guess + 1e-12, 1200)
    gs = np.array([growth(P, R, k) for k in ks])
    i = int(np.argmax(gs))
    if gs[i] <= 0:
        return 0.0, 0.0, 0.0
    # refine
    lo, hi = ks[max(i - 1, 0)], ks[min(i + 1, len(ks) - 1)]
    kk = np.geomspace(lo, hi, 200)
    gg = np.array([growth(P, R, k) for k in kk])
    j = int(np.argmax(gg))
    kpos = ks[gs > 0]
    return float(kk[j]), float(gg[j]), float(kpos.max())


def main():
    res = {}
    for name, P in POINTS.items():
        g, lam, lX, v = P['g'], P['lam'], P['lX'], P['v']
        _, _, _, Xc = background(P, 0.0)
        a_eff = g ** 4 / lam - lX
        rows = []
        for frac in [0.9999, 0.999, 0.99, 0.95, 0.9, 0.7, 0.5, 0.3]:
            R = frac * Xc
            phi0sq, w2, mp2, _ = background(P, R)
            kp, gm, kedge = peak(P, R)
            n = 2 * np.sqrt(w2) * R * R
            rows.append(dict(R_over_Xc=frac, R_GeV=R, omega_GeV=np.sqrt(w2), m_phi_GeV=np.sqrt(max(mp2, 0)),
                             k_peak_GeV=kp, gamma_max_GeV=gm, k_band_edge_GeV=kedge,
                             single_field_k_edge_GeV=2 * np.sqrt(a_eff) * R,
                             n_circ_GeV3=n, Q_cell_eps1=n * (2 * np.pi / kp) ** 3 if kp > 0 else None))
        # stability just above Xc (phi=0 phase)
        above = [dict(R_over_Xc=f, gamma_any=max(growth(P, f * Xc, k) for k in np.geomspace(1e-6, 1, 400)))
                 for f in (1.001, 1.2, 2.0, 10.0)]
        # onset temperature: n_circ(Xc) eps = YQ s(T)
        n_c = 2 * np.sqrt(2 * lX) * Xc ** 3
        onset = []
        deltas = np.geomspace(1e-12, 0.1, 140)
        pk = [peak(P, Xc * (1 - dl)) for dl in deltas]
        gam = np.array([x[1] for x in pk])
        for eps in (1.0, 0.5, 0.1):
            f = lambda lT: eps * n_c - YQ * (2 * np.pi ** 2 / 45) * gstar(np.exp(lT)) * np.exp(3 * lT)
            T = np.exp(brentq(f, np.log(1e-3), np.log(1e8)))
            H = np.sqrt(np.pi ** 2 * gstar(T) / 90) * T * T / M
            # slow quench: delta grid, N(delta) = int gamma/H d delta
            N = np.concatenate([[0], np.cumsum(0.5 * (gam[1:] + gam[:-1]) * np.diff(deltas))]) / H
            ent = dict(eps=eps, T_onset_GeV=T, H_GeV=H)
            for Ns in (5.0, 30.0, 60.0):
                idx = np.searchsorted(N, Ns)
                if idx >= len(deltas):
                    ent[f'N{int(Ns)}'] = None; continue
                dl = deltas[idx]; R = Xc * (1 - dl)
                kp, gm, _ = pk[idx]
                _, w2, _, _ = background(P, R)
                n = 2 * np.sqrt(w2) * R * R * eps
                ent[f'N{int(Ns)}'] = dict(delta=float(dl), k_peak_GeV=kp, gamma_max_GeV=gm,
                                         Q_cell=float(n * (2 * np.pi / kp) ** 3))
            onset.append(ent)
        sigma = 2 * np.sqrt(2 * lam) / 3 * v ** 3
        res[name] = dict(params=P, Xc_GeV=Xc, a_eff=a_eff, rows=rows, stability_above_Xc=above,
                         onset=onset, wall_tension_GeV3=sigma, wall_tension_cuberoot_MeV=1e3 * sigma ** (1 / 3),
                         bias_for_annihilation_before_T1MeV_GeV4=sigma * np.sqrt(np.pi ** 2 * 10.75 / 90) * 1e-6 / M,
                         planck_odd_bias_v5_over_M_GeV4=v ** 5 / M)
        print(name, json.dumps(res[name], indent=1, default=float)[:4000])
    OUT.mkdir(exist_ok=True)
    (OUT / 'fragmentation.json').write_text(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
