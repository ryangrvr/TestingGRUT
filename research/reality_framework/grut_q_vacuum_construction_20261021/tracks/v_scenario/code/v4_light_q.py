"""v_scenario / v4: independent recomputation of the light-q relic (scenario NOTES §4.2, verdict row 6).

Equation (my derivation, v1 check C; K0 = M^2, y = (K - K0)/K0, |y| << 1 so the KV quartic terms are dropped
and checked a posteriori):
    y'' + 3 (H/m) y' + y = sigma(tau),   tau = m t,   sigma = S/(3 K0 m^2) = (S/rho) (H/m)^2  (GR background)
with S = rho - 3p = (1 - 3w) rho_pl + rho_d from the Saikawa-Shirai w table (own loader, own download), H from
3 K0 H^2 = rho_pl + rho_d, N = ln(a/a0) from entropy conservation. Integrated in tau with DOP853 from T_i to
m/H = 3000; the Einstein-frame oscillation energy about the quasi-static track
    E = (3/4) K0 m^2 [ (y' - sigma')^2 + (y - sigma)^2 ]         (period-averaged; sigma track good to (H/m)^2)
is then constant x a^3 and is carried to today: Omega_q h^2 = E a^3 / rho_crit,100.
Cases (as in the scenario): quasi-static start at 1 TeV (leading order and O(m^-4)), q = q0 at rest at
T_i = 1, 10, 100, 1000 TeV; trace variants at 100 TeV (S/rho = 0 above 300 GeV, S/rho + 1 sigma above 300 GeV);
mass scan point m_q = 1 eV.
Output: ../results/v4_light_q.json
"""
import json
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from scipy.integrate import solve_ivp

DL = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/qtheory/v_scenario_dl/ss'
M = 2.435e18
K0 = M * M
HBAR = 6.582119569e-25
KB = 8.617333262e-14
H100 = 100.0 / 3.0856775814913673e19 * HBAR
RHO_C100 = 3 * K0 * H100 ** 2
T0 = 2.7255 * KB
OMEGA_C = 0.120


def load(fn):
    rows = [l.split() for l in open(fn) if l.strip() and not l.startswith('#')]
    return np.array([[float(x) for x in r] for r in rows])


G = load(os.path.join(DL, 'x/anc/standardmodel2018.dat'))
W = load(os.path.join(DL, 'eos_w.dat'))
lTg = np.log(G[:, 0])
lTw = np.log(W[:, 0])
s0 = 2 * math.pi ** 2 / 45 * np.interp(math.log(T0), lTg, G[:, 3]) * T0 ** 3
rho_d0 = (0.0224 + OMEGA_C) * RHO_C100

# background grid in ln T (high to low)
LT = np.linspace(math.log(1e7), math.log(1.0), 200001)
T_ = np.exp(LT)
grho = np.interp(LT, lTg, G[:, 1])
gs = np.interp(LT, lTg, G[:, 3])
rho_pl = math.pi ** 2 / 30 * grho * T_ ** 4
rho_d = rho_d0 * (2 * math.pi ** 2 / 45 * gs * T_ ** 3) / s0
NN = -LT + math.log(T0) - np.log(gs / np.interp(math.log(T0), lTg, G[:, 3])) / 3
HH = np.sqrt((rho_pl + rho_d) / (3 * K0))
w_ = np.interp(LT, lTw, W[:, 1])
werr = np.interp(LT, lTw, W[:, 2])


def S_over_rho(variant):
    sr = (1 - 3 * w_) * rho_pl / (rho_pl + rho_d) + rho_d / (rho_pl + rho_d)
    if variant == 'zero_above_300':
        sr = np.where(T_ > 300, rho_d / (rho_pl + rho_d), sr)
    elif variant == 'plus1sigma_above_300':
        sr = np.where(T_ > 300, sr + 3 * werr * rho_pl / (rho_pl + rho_d), sr)
    return sr


def run(m, start, T_i=None, variant='nominal', mH_end=3000.0, rtol=1e-10):
    t0w = time.time()
    SR = S_over_rho(variant)
    # functions of N (N increasing as T decreases): use arrays reversed? NN decreases with LT decreasing -> NN increases
    Ngrid = NN                     # increasing along the array (T decreasing)
    Ht = HH / m                    # H/m
    sig = SR * Ht ** 2             # sigma(N)
    # dsigma/dN for the quasi-static velocity: sigma' (tau) = (H/m) dsigma/dN
    dsig_dN = np.gradient(sig, Ngrid)
    d2sig_dN2 = np.gradient(dsig_dN, Ngrid)
    dHt_dN = np.gradient(Ht, Ngrid)
    Hf = lambda N: np.interp(N, Ngrid, Ht)
    sf = lambda N: np.interp(N, Ngrid, sig)
    if start == 'qs_lo' or start == 'qs_o4':
        N_i = float(np.interp(-math.log(T_i), -LT, Ngrid))
        h = Hf(N_i)
        y = sf(N_i)
        yp = h * np.interp(N_i, Ngrid, dsig_dN)
        if start == 'qs_o4':
            # y2 = -(y1'' + 3 h y1') in tau units, y1 = sigma: y1'' = h (h sigma_NN + h_N sigma_N)
            s1 = np.interp(N_i, Ngrid, dsig_dN)
            s2 = np.interp(N_i, Ngrid, d2sig_dN2)
            hN = np.interp(N_i, Ngrid, dHt_dN)
            y1pp = h * (h * s2 + hN * s1)
            y2 = -(y1pp + 3 * h * yp)
            # derivative of y2 by finite difference along N
            def y2_of(Nv):
                hh = Hf(Nv)
                a1 = np.interp(Nv, Ngrid, dsig_dN)
                a2 = np.interp(Nv, Ngrid, d2sig_dN2)
                b1 = np.interp(Nv, Ngrid, dHt_dN)
                return -(hh * (hh * a2 + b1 * a1) + 3 * hh * hh * a1)
            dN = 1e-3
            y2p = h * (y2_of(N_i + dN) - y2_of(N_i - dN)) / (2 * dN)
            y, yp = y + y2, yp + y2p
    else:   # q at q0 at rest
        N_i = float(np.interp(-math.log(T_i), -LT, Ngrid))
        y, yp = 0.0, 0.0
    # end: m/H = mH_end
    N_end = float(np.interp(-1.0 / mH_end, -Ht, Ngrid))

    def rhs(tau, s):
        yv, ypv, Nv = s
        h = Hf(Nv)
        return [ypv, -3 * h * ypv - yv + sf(Nv), h]
    ev = lambda tau, s: s[2] - N_end
    ev.terminal = True
    sol = solve_ivp(rhs, (0, 1e12), [y, yp, N_i], method='DOP853', rtol=rtol, atol=1e-22, dense_output=True,
                    events=ev, max_step=0.3)
    tau_f = sol.t[-1]
    # E a^3 over late windows (one period each), m/H from 300 to mH_end
    taus = np.linspace(0, tau_f, 200001)
    Ns = sol.sol(taus)[2]
    out = []
    for mH in (300.0, 1000.0, 2000.0, 0.95 * mH_end):
        Nw = float(np.interp(-1.0 / mH, -Ht, Ngrid))
        tc = float(np.interp(Nw, Ns, taus))
        tw = tc + 2 * math.pi * (np.arange(64) / 64 - 0.5)
        ss = sol.sol(tw)
        hs = Hf(ss[2])
        sg = sf(ss[2])
        sgp = hs * np.interp(ss[2], Ngrid, dsig_dN)
        E = 0.75 * K0 * m * m * ((ss[1] - sgp) ** 2 + (ss[0] - sg) ** 2)
        Ea3 = float(np.mean(E * np.exp(3 * ss[2])))
        amp = float(np.sqrt(np.mean((ss[0] - sg) ** 2) * 2))
        out.append(dict(m_over_H=mH, E_a3=Ea3, Omega_h2=Ea3 / RHO_C100, amp_y=amp))
    ymax = float(np.max(np.abs(sol.sol(np.linspace(0, tau_f, 20001))[0])))
    Oh2 = out[-1]['Omega_h2']
    return dict(m_GeV=m, start=start, T_i=T_i, variant=variant, Omega_h2=Oh2, Omega_over_c=Oh2 / OMEGA_C,
                E_a3_rel_spread=float(np.std([o['E_a3'] for o in out]) / np.mean([o['E_a3'] for o in out])),
                windows=out, max_abs_y=ymax, m_over_H_start=float(1 / Hf(N_i)), nfev=int(sol.nfev),
                wall_s=time.time() - t0w)


def job(kw):
    return run(**kw)


if __name__ == '__main__':
    m4 = 4e-12
    jobs = [dict(m=m4, start='q0', T_i=1e6), dict(m=m4, start='q0', T_i=1e5), dict(m=m4, start='q0', T_i=1e4),
            dict(m=m4, start='q0', T_i=1e3), dict(m=m4, start='qs_lo', T_i=1e3), dict(m=m4, start='qs_o4', T_i=1e3),
            dict(m=m4, start='q0', T_i=1e5, variant='zero_above_300'),
            dict(m=m4, start='q0', T_i=1e5, variant='plus1sigma_above_300'),
            dict(m=m4, start='q0', T_i=1e3, rtol=1e-12)]
    # 1 eV: quasi-static at m/H = 3 and q0 at m/H = 0.01 (temperatures from the background)
    m1 = 1e-9
    Tq3 = float(np.exp(np.interp(-1 / 3.0, -HH / m1, LT)))
    Tq001 = float(np.exp(np.interp(-100.0, -HH / m1, LT)))
    jobs += [dict(m=m1, start='qs_lo', T_i=Tq3), dict(m=m1, start='q0', T_i=Tq001)]
    t0 = time.time()
    with ProcessPoolExecutor(4) as ex:
        res = list(ex.map(job, jobs))
    here = os.path.dirname(os.path.abspath(__file__))
    with open(os.path.join(here, '..', 'results', 'v4_light_q.json'), 'w') as f:
        json.dump(dict(script='v4_light_q.py', wall_s=time.time() - t0, runs=res), f, indent=1)
    for r in res:
        print('m=%.3g %-6s T_i=%.3g %-22s Omega h2=%.4g (/Oc %.3g) spread=%.2g |y|max=%.2g m/H_i=%.3g nfev=%d %.0fs'
              % (r['m_GeV'], r['start'], r['T_i'], r['variant'], r['Omega_h2'], r['Omega_over_c'],
                 r['E_a3_rel_spread'], r['max_abs_y'], r['m_over_H_start'], r['nfev'], r['wall_s']))
