"""mf_gemini SKEPTIC: n_s for the hybrid family with the waterfall attached consistently [C/D].

Independent of mf_gemini/code/ns_hybrid.py:
  * N_* by direct horizon matching (k_* = a_* H_*, a matter interval with a^2 H = K after eps_H = 1, then
    radiation with entropy conservation), not the Planck closed form;
  * the log family V = V0 [1 + (A/2) ln((s^2 + sc^2)/sc^2)] (s0 = sc, the shape run in sk_waterfall.py) with
    FINITE sigma_c: N_sigma = (r^2 - 1)/(4 Delta) + ln r/(2 Delta), r = sigma_*/sigma_c, A = 2 Delta sc^2; A_s
    fixes sc = p_* (r^2 + 1)/(2 Delta r);
  * second-order slow roll WITH the xi^2 term (Stewart-Lyth: n_s - 1 = 2 eta - 6 eps + 2[(1/3) eta^2
    + (8C - 1) eps eta - (5/3 + 12 C) eps^2 - (C - 1/3) xi^2], C = -0.7296) -- the track's ns2 multiplies the xi^2
    term by zero;
  * an EXACT check: Mukhanov-Sasaki modes on the numerical sigma background (chi = 0, heavy), n_s and alpha_s
    from d ln P/d ln k at the pivot exit;
  * N_wf(Delta) from this skeptic's own log-family waterfall runs (results/sk_waterfall.json);
  * hilltop: |eta| window by root finding (xi = 0 there).
Run: python3 -B sk_ns.py -> ../results/sk_ns.json
"""
import json
import os
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
MP = 2.435e18
AS = 2.100e-9
C_SL = -0.7296
T0 = 2.7255 * 8.617333e-14        # GeV
MPC_GEV = 1.0 / 1.5637e38          # 1 Mpc^-1 in GeV


def Nstar_direct(H, Hend_over_H=0.7126, K=1.0, k=0.05, g=106.75, gs0=3.909):
    """e-folds between exit of k (a_* H_* = k, H_* = H) and eps_H = 1, by direct matching."""
    Hend = Hend_over_H * H
    a_th = K**2
    H_th = Hend * K**-3
    rho_th = 3 * MP**2 * H_th**2
    T_th = (30 * rho_th / (np.pi**2 * g))**0.25
    a_end_over_a0 = (1.0 / a_th) * (gs0 / g)**(1.0 / 3.0) * T0 / T_th
    return float(np.log(a_end_over_a0 / (k * MPC_GEV / H)))


def ns_sr2(eta, eps, xi2):
    return 1 + 2 * eta - 6 * eps + 2 * (eta**2 / 3 + (8 * C_SL - 1) * eps * eta - (5.0 / 3 + 12 * C_SL) * eps**2
                                         - (C_SL - 1.0 / 3) * xi2)


def ms_ns(fV, dfV, s_star_guess, Nsig, sc):
    """exact n_s, alpha_s: background sigma(N) for V = V0 (1 + f), then MS modes R'' + (3 + e2 - e1) R' + (k/aH)^2 R = 0
    for three k around the mode exiting Nsig e-folds before sigma = sc.  Units M_P = 1, V0 = 3 (H ~ 1)."""
    V = lambda s: 3.0 * (1 + fV(s))
    dV = lambda s: 3.0 * dfV(s)
    # start 8 e-folds before the pivot: integrate forward from a guess with slow-roll IC, locate N(sc)
    s0 = s_star_guess * 1.6 + 1e-12
    def rhs(N, y):
        s, ds = y
        e = 0.5 * ds * ds
        return [ds, -(3 - e) * (ds + dV(s) / V(s))]
    def ev(N, y):
        return y[0] - sc
    ev.terminal = True
    sol = solve_ivp(rhs, (0, 2000), [s0, -dV(s0) / V(s0)], method='LSODA', rtol=1e-12, atol=1e-20, events=ev,
                    dense_output=True)
    Nc = sol.t_events[0][0]
    Npiv = Nc - Nsig
    if Npiv < 8:
        return None
    lnaH = lambda N: N + 0.5 * np.log(V(sol.sol(N)[0]) / (3 - 0.5 * sol.sol(N)[1]**2))
    out = {}
    lnks = lnaH(Npiv) + np.array([-0.5, 0.0, 0.5])
    P = []
    for lnk in lnks:
        k = np.exp(lnk)
        Ni = brentq(lambda N: lnk - lnaH(N) - np.log(200.0), 0.0, Npiv + 2)
        Nf = brentq(lambda N: lnk - lnaH(N) + np.log(1e3), Npiv - 2, Nc)
        s, ds = sol.sol(Ni)
        e1 = 0.5 * ds * ds
        dds = -(3 - e1) * (ds + dV(s) / V(s))
        z = np.exp(Ni) * ds
        x = k / np.exp(lnaH(Ni))
        R0 = 1 / (z * np.sqrt(2 * k))
        dR0 = R0 * (-1j * x - (1 + dds / ds))

        def mrhs(N, y):
            s, ds = sol.sol(N)
            e1 = 0.5 * ds * ds
            dds = -(3 - e1) * (ds + dV(s) / V(s))
            e2 = 2 * dds / ds
            q = (k / np.exp(lnaH(N)))**2
            fr = 3 + e2 - e1
            return [y[2], y[3], -fr * y[2] - q * y[0], -fr * y[3] - q * y[1]]
        m = solve_ivp(mrhs, (Ni, Nf), [R0.real, R0.imag, dR0.real, dR0.imag], method='DOP853', rtol=1e-11,
                      atol=1e-30 * abs(R0))
        P.append(k**3 * (m.y[0, -1]**2 + m.y[1, -1]**2))
    lnP = np.log(P)
    out['n_s'] = float(1 + (lnP[2] - lnP[0]) / 1.0)
    out['alpha_s'] = float((lnP[2] - 2 * lnP[1] + lnP[0]) / 0.25)
    out['A_s_code'] = float(P[1] / (2 * np.pi**2))
    return out


def log_family(Delta, Nsig, pstar, exact=True):
    r = brentq(lambda x: (x * x - 1) / (4 * Delta) + np.log(x) / (2 * Delta) - Nsig, 1.0 + 1e-12, 1e4)
    sc = pstar * (r * r + 1) / (2 * Delta * r)
    A = 2 * Delta * sc**2
    f = lambda s: 0.5 * A * np.log((s * s + sc * sc) / (sc * sc))
    df = lambda s: A * s / (s * s + sc * sc)
    d2f = lambda s: A * (sc * sc - s * s) / (s * s + sc * sc)**2
    d3f = lambda s: 2 * A * s * (s * s - 3 * sc * sc) / (s * s + sc * sc)**3
    s = r * sc
    V = 1 + f(s)
    eps = 0.5 * (df(s) / V)**2
    eta = d2f(s) / V
    xi2 = df(s) * d3f(s) / V**2
    res = dict(n_s_sr2_with_xi=float(ns_sr2(eta, eps, xi2)), n_s_sr2_without_xi=float(ns_sr2(eta, eps, 0.0)),
               n_s_track_formula_1_minus_1_over_Nsigma=float(1 - 1 / Nsig), N_eff=float(Nsig + 1 / (2 * Delta)),
               alpha_s_sr=float(16 * eps * eta - 24 * eps**2 - 2 * xi2), r=float(16 * eps), sigma_c=float(sc),
               sigma_star=float(s), kappa=float(np.sqrt(8 * np.pi**2 * A)), eta=float(eta), xi2=float(xi2))
    if exact:
        ms = ms_ns(f, df, s, Nsig, sc)
        if ms:
            res.update(n_s_MS=ms['n_s'], alpha_s_MS=ms['alpha_s'],
                       A_s_MS_over_target=float(ms['A_s_code'] * (1.3e11 / MP)**2 / AS))
    return res


def main():
    with open(os.path.join(RES, 'sk_waterfall.json')) as fh:
        wf = json.load(fh)
    Nwf_log = {c['Delta']: c['results']['q1']['N_wf_mean'] for c in wf['cases'] if c['kind'] == 'log'}
    out = dict(Nwf_log=Nwf_log)
    for H in [1.3e11, 1.25e11, 1.55e11, 1.81e11]:
        hub = H / MP
        eps = hub**2 / (8 * np.pi**2 * AS)
        pstar = np.sqrt(2 * eps)
        d = dict(eps=eps, r=16 * eps, p_star=pstar)
        for K in [1.0, 2.0, 3.0]:
            d[f'N_star_K{K:g}'] = Nstar_direct(H, K=K)
        f = lambda e, t: ns_sr2(-e, eps, 0.0) - t
        d['hill_abs_eta_window'] = [float(brentq(f, 1e-4, 0.1, args=(0.973,))), float(brentq(f, 1e-4, 0.1, args=(0.967,)))]
        if H == 1.3e11:
            lf = {}
            for K in [1.0, 2.0, 3.0]:
                Ns = d[f'N_star_K{K:g}']
                for D, Nw in sorted(Nwf_log.items()):
                    lf[f'K{K:g}_D{D:g}'] = dict(N_wf=Nw, N_sigma=Ns - Nw, **log_family(D, Ns - Nw, pstar, exact=(K == 2.0)))
            d['log_family'] = lf
            # the track's quoted range: N_wf = 17.4-40.7 with 1 - 1/N_sigma, vs the consistent family at those Delta
            # (Delta = 0.01 is not in this skeptic's waterfall runs; use the track's N_wf = 40.7 there)
            Ns = d['N_star_K2']
            d['log_family_track_Delta0.01'] = dict(N_wf=40.7, **log_family(0.01, Ns - 40.7, pstar, exact=False)) \
                if Ns - 40.7 > (1 - 1) else None
        out[f'H0={H:.3g}'] = d
    l = out['H0=1.3e+11']['log_family']
    vals = [v['n_s_MS'] for k, v in l.items() if k.startswith('K2_') and 'n_s_MS' in v]
    out['log_K2_MS_range'] = dict(n_s=[min(vals), max(vals)],
                                  Planck_pull=[(x - 0.9649) / 0.0042 for x in (min(vals), max(vals))],
                                  ACT_pull=[(x - 0.9743) / 0.0034 for x in (min(vals), max(vals))])
    with open(os.path.join(RES, 'sk_ns.json'), 'w') as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
