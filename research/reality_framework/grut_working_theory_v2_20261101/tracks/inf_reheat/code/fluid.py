#!/usr/bin/env python3
"""inf_reheat (Claude, own code): inflaton-only background for slow reheating, H_i units, M_r = 1.

Stage 1 (tau <= tau_match): the oscillating inflaton released at rest with H = 1 (quadratic, or the RH1 T-model plateau
V = (b^2 w^2/2) tanh^2(I/w)), I'' + (3H + Gamma) I' + V'(I) = 0, rho_r' = -4 H rho_r + Gamma I'^2 (the engine's background
with the hidden sector dropped: eta2 V_h ~ 1e-10, dq_heat SKEPTIC sk_K).
Stage 2: the cycle-averaged fluid in N = ln a with comoving x = rho_I a^3, y = rho_r a^4:
    dx/dN = -(Gamma/H) x,   dy/dN = (Gamma/H) a x,   3 H^2 = x/a^3 + y/a^4      (w_I = 0 for a quadratic minimum)
run until x/a^3 < 1e-14 y/a^4.  Outputs K = a^2 H in radiation (= sqrt(y/3)), D = (y/3)^{3/4} = K^{3/2},
the e-folds to radiation-inflaton equality N_eq, to H = Gamma (N_RH), and T_SM(tau) for instantaneous thermalisation.
Validated against the rprh.c engine runs (results/engine_runs.json, gscan).  Labels: [C].
Usage: fluid.py  -> ../results/fluid.json
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import rh_common as rc


def V_and_dV(I, b, w):
    if w <= 0:
        return 0.5 * b * b * I * I, b * b * I
    t = math.tanh(I / w)
    return 0.5 * b * b * w * w * t * t, b * b * w * t * (1 - t * t)


def stage1(b, G, w=0.0, tau_match=60.0, dense=False):
    if w <= 0:
        I0 = math.sqrt(6.0) / b
    else:
        I0 = w * math.atanh(math.sqrt(6.0) / (b * w))

    def rhs(t, y):
        la, I, P, rr = y
        V, dV = V_and_dV(I, b, w)
        H = math.sqrt(max((0.5 * P * P + V + rr) / 3.0, 0.0))
        return [H, P, -(3 * H + G) * P - dV, -4 * H * rr + G * P * P]

    sol = solve_ivp(rhs, (0, tau_match), [0.0, I0, 0.0, 0.0], method='DOP853', rtol=1e-11, atol=1e-14,
                    max_step=0.02 / max(1.0, b / 3), dense_output=dense)
    return sol


def run(b, G, w=0.0, tau_match=60.0, want_T=False):
    s1 = stage1(b, G, w, tau_match, dense=want_T)
    la, I, P, rr = s1.y[:, -1]
    a = math.exp(la)
    V, _ = V_and_dV(I, b, w)
    rhoI = 0.5 * P * P + V
    x0, y0 = rhoI * a ** 3, rr * a ** 4
    if G == 0:
        return dict(b=b, Gamma=G, w=w, K=None)

    def rhs(N, z):
        x, y = z
        a = math.exp(N)
        H = math.sqrt(max((x / a ** 3 + y / a ** 4) / 3.0, 1e-300))
        return [-(G / H) * x, (G / H) * a * x]

    def ev_done(N, z):
        a = math.exp(N)
        return z[0] / a ** 3 - 1e-14 * z[1] / a ** 4
    ev_done.terminal = True

    def ev_eq(N, z):
        a = math.exp(N)
        return z[0] * a - z[1]          # rho_I = rho_r
    ev_eq.terminal = False

    def ev_HG(N, z):
        a = math.exp(N)
        return math.sqrt((z[0] / a ** 3 + z[1] / a ** 4) / 3.0) - G
    ev_HG.terminal = False

    sol = solve_ivp(rhs, (la, la + 80.0), [x0, y0], method='LSODA', rtol=1e-11, atol=1e-30, events=(ev_done, ev_eq, ev_HG))
    y = sol.y[1, -1]
    K = math.sqrt(y / 3.0)
    out = dict(b=b, Gamma=G, w=w, K=K, D=K ** 1.5, N_eq=float(sol.t_events[1][0]) if len(sol.t_events[1]) else None,
               N_RH=float(sol.t_events[2][0]) if len(sol.t_events[2]) else None, N_end=float(sol.t[-1]),
               a_match=a, radfrac_match=rr / (rhoI + rr))
    if want_T:
        out['s1'] = s1
    return out


def T_SM_early(b, G, taus, w=0.0, H_GeV=rc.H_CAL_FREEX):
    """SM temperature (GeV) for instantaneous thermalisation of direct SM decay products, at physical times taus."""
    s1 = stage1(b, G, w, max(taus) + 1.0, dense=True)
    out = []
    for t in taus:
        la, I, P, rr = s1.sol(t)
        rho = rr * rc.MP ** 2 * H_GeV ** 2
        out.append(rc.T_of_rho(rho) if rho > 0 else 0.0)
    return out


def main():
    res = dict(method=__doc__)
    eng = json.load(open(rc.RES / 'engine_runs.json'))
    val = []
    for row in eng['gscan']:
        f = run(3.0, row['Gamma'])
        val.append(dict(Gamma=row['Gamma'], K_engine=row['K'], K_fluid=f['K'], rel=f['K'] / row['K'] - 1, D_engine=row['D'], D_fluid=f['D']))
        print('validate G=%.3g  K engine %.6f fluid %.6f  rel %.2e' % (row['Gamma'], row['K'], f['K'], f['K'] / row['K'] - 1))
    res['validation_b3'] = val
    # tau_match sensitivity
    res['tau_match_check'] = [dict(tau_match=tm, K=run(3.0, 1e-3, tau_match=tm)['K']) for tm in (30.0, 60.0, 120.0)]
    print('tau_match', res['tau_match_check'])
    grid = []
    for b in (1.5, 2.0, 2.5, 3.0, 4.0, 5.0):
        for G in (0.5, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12):
            f = run(b, G, tau_match=60.0 if G < 0.05 else 300.0)
            f.pop('s1', None)
            grid.append(f)
    res['grid'] = grid
    # power-law fit K = c (Gamma)^(-1/3) for small Gamma, each b
    fits = {}
    for b in (1.5, 2.0, 2.5, 3.0, 4.0, 5.0):
        rows = [g for g in grid if g['b'] == b and g['Gamma'] <= 1e-6]
        c = [g['K'] * g['Gamma'] ** (1 / 3) for g in rows]
        fits[str(b)] = dict(c_K=float(np.mean(c)), spread=float(np.ptp(c) / np.mean(c)),
                            N_RH_minus_two_thirds_ln=float(np.mean([g['N_RH'] - (2 / 3) * math.log(1 / g['Gamma']) for g in rows])))
        print('b=%g  K Gamma^(1/3) = %.5f (spread %.1e)  N_RH - (2/3) ln(1/G) = %.4f' % (b, fits[str(b)]['c_K'], fits[str(b)]['spread'],
                                                                              fits[str(b)]['N_RH_minus_two_thirds_ln']))
    res['fits'] = fits
    # T-model plateau variants at b = 3 (RH1), slow decay
    res['tmodel'] = [{k: v for k, v in run(3.0, 1e-6, w=w).items() if k != 's1'} for w in (1.0, 1.5, 2.5)]
    json.dump(res, open(rc.RES / 'fluid.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
