#!/usr/bin/env python3
"""inf_chain SKEPTIC (Claude, own code): re-derivation of the chain arithmetic from the track's stored release/vector
counts, with my own deposits, my own yield, the analytic large-W limit of the S/Q closed form, the leak estimate with
benchmark-specific logarithms, the fixed-point structure of the one-step iteration, and two relic diagnostics.

Inputs read (data only): tracks/inf_chain/results/{release,vectors,yy_born,chain,calib,pmap_kin,leak_estimate}_*.json.
Nothing of the track's code is imported.
 (1) S/Q closed form 0.190255 (W - 0.97494)^0.75013 - 0.00428: its W -> infinity limit is the ideal massless gas of
     X (complex, 2 dof) + phi (1 dof) at fixed charge, with E_c the circular-orbit energy of V = lam_X |X|^4:
       S/Q -> (pi^2/15)(30/pi^2)^{3/4} lam^{1/4} W^{3/4}/sqrt2 ,   N_c/E^{3/4} = 2 sqrt2 /(3^{3/4} lam^{1/4})   [D]
 (2) ledger: R_Y = (N_YL + N_YH) m_L <a_L>/E_X, R_casc = (N_YH + 2 N_sigma) m_A <a_A>/E_X, R_YY from the Born counts,
     R_vec from the vector file; <a_L> mine: radiation (pi K/2 Gamma_L)^{1/2} for Gamma_I = 0.5, own matter+radiation fluid
     integral for TB; Y/Y_Q from N_X and K mine; H_i,req = H_M1 (0.9374/(M Y/Y_Q))^{2/3} with the track's p.
 (3) leak: c(M_P) = ln(M_P/M_Y)/(4 pi^2) at each benchmark's own M_Y = 30 H_i (the track used M1r's 13.31 throughout).
 (4) iteration: secant (log-log) fixed point of H_req(H_pot) from the two points per family.
 (5) relic: entropy per final particle S/Q/(1 + 2p); asymmetric share of Omega_DM = 1/M; n_Xbar/s_hidden = p/(S/Q).
Output ../results/sk_ledger.json.  Labels: [D] (1), [C] (2)-(5) on the track's stored engine/vector counts.
"""
import json, math, pathlib
import numpy as np
from scipy.integrate import solve_ivp

HERE = pathlib.Path(__file__).resolve().parent
RES = HERE.parent.parent / 'results'
OUT = {}
MP = 2.435e18; H_M1 = 1.342842735e11; YQ = 5.07941423e-9; GST = 106.75
M_L, G_L = 22.24859546128699, 1.005903273746766e-07
MU, MSIG = 7.159456, 40.0
LX = 1e-3


def SQ(W):
    return 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428


def asymptote():
    coef = (math.pi ** 2 / 15) * (30 / math.pi ** 2) ** 0.75 * LX ** 0.25 / math.sqrt(2)
    chart = 2 * math.sqrt(2) / (3 ** 0.75 * LX ** 0.25)
    Ws = [1e3, 1e4, 1e5, 1e6]
    return dict(coef_ideal_gas=coef, coef_fit=0.190255, rel=coef / 0.190255 - 1, C_HART_analytic=chart, C_HART_track=6.977568,
                ratio_closed_over_ideal={str(W): SQ(W) / (coef * W ** 0.75) for W in Ws})


def a_table_fluid(x0, gam, tmax):
    def rhs(t, y):
        a, rm, rr = y
        H = math.sqrt((rm + rr) / 3)
        return [a * H, -3 * H * rm - gam * rm, -4 * H * rr + gam * rm]
    s = solve_ivp(rhs, (0, tmax), [1.0, x0, 0.0], method='LSODA', rtol=1e-10, atol=1e-30, dense_output=True)
    return s


def mean_a_fluid(x0, gam, GL):
    tmax = 60 / GL
    s = a_table_fluid(x0, gam, tmax)
    t = np.concatenate([np.linspace(0, 10 / GL, 200001)])
    a = s.sol(t)[0]
    w = GL * np.exp(-GL * t)
    inside = np.trapezoid(w * a, t)
    # tail beyond 10/GL: radiation a ~ a_end (t/t_end)^1/2 if decayed
    return float(inside)


def yield_YQ(NX, K):
    rr = 3 * K * K * MP ** 2 * H_M1 ** 2
    Ta = (30 * rr / (math.pi ** 2 * GST)) ** 0.25
    s = 2 * math.pi ** 2 / 45 * GST * Ta ** 3
    return NX * H_M1 ** 3 / s / YQ


def ledger():
    tags = {'HT6t': ('HT6t', None), 'T1': ('T1', None), 'HT6t_Send': ('HT6t_Send', None), 'HT6t_cal1': ('HT6t_cal1', None),
            'T1_cal1': ('T1_cal1', None), 'TB1': ('TB1', 7.512e-09), 'TB2': ('TB2', 2.113e-09)}
    rows = {}
    for key, (rk, g) in tags.items():
        vt = rk + ('' if g is None else '_g%.3e' % g)
        rel = json.load(open(RES / ('release_%s.json' % rk)))['result']
        vec = json.load(open(RES / ('vectors_%s.json' % vt)))
        chn = json.load(open(RES / ('chain_%s.json' % vt)))
        v6 = [x for x in vec['by_e'] if abs(x['e'] - 1e-6) < 1e-12][0]
        K = vec['background']['K']
        EX, NX, eps = rel['E_X'], abs(rel['N_X']), rel['eps']
        if g is None:
            aL = math.sqrt(math.pi * K / (2 * G_L))
        else:
            aL = mean_a_fluid(rel['x0'], g, G_L)
        DL = M_L * aL
        # deposit per vector: m_A <a_A>, <a_A> = (pi K/2 Gamma_A)^1/2, Gamma_A = e^2 m_A/48 pi, m_A = e F_v
        e = 1e-6; mA = e * vec['background']['Fv']; GA = e * e * mA / (48 * math.pi)
        DA = mA * math.sqrt(math.pi * K / (2 * GA))
        RY = (rel['N_YL'] + rel['N_YH']) * DL / EX
        Rc = (rel['N_YH'] + 2 * rel['N_sigma']) * DA / EX
        yy = chn.get('yy')
        Ryy = 0.0
        if yy:
            Ryy = ((yy['N_YL'] + yy['N_YH'] + yy['tail']) * DL + (yy['N_YH'] + 0.5 * yy['tail']) * DA) / EX
        Rv = v6['N_A'] * DA / EX
        R = Rv + RY + Rc + Ryy
        W = eps ** (-4 / 3) * (1 + R)
        sq = SQ(W)
        RYinf = (rel['N_YL'] + rel['N_YH'] + ((yy['N_YL'] + yy['N_YH'] + yy['tail']) if yy else 0.0)) * DL / EX
        sqinf = SQ(eps ** (-4 / 3) * (1 + RYinf))
        yyq = yield_YQ(NX, K)
        p = chn['relic_e1e6']['p']; M = 1 + 2 * p
        Hreq = H_M1 * (0.9374 / (M * yyq)) ** (2 / 3)
        r = dict(K=K, aL_mine=aL, aL_track=chn['a_mean_L'], R_vec=Rv, R_Y=RY, R_casc=Rc, R_YY=Ryy, R_total=R, W=W, S_per_Q=sq,
                 S_per_Q_track=chn['e1e6']['S_per_Q'], S_per_Q_einf=sqinf, S_per_Q_einf_track=chn['S_per_Q_e_inf'],
                 Y_over_YQ_mine=yyq, Y_over_YQ_track=chn['Y_over_YQ_at_HM1'], p_track=p, M=M, H_req_mine=Hreq,
                 H_req_track=chn['relic_e1e6']['H_i_req_GeV'], entropy_per_final_particle=sq / M, asym_share_of_DM=1 / M,
                 nXbar_over_s_hidden=p / sq)
        r['rel'] = dict(S_per_Q=sq / r['S_per_Q_track'] - 1, S_per_Q_einf=sqinf / r['S_per_Q_einf_track'] - 1,
                        YYQ=yyq / r['Y_over_YQ_track'] - 1, H_req=Hreq / r['H_req_track'] - 1, aL=aL / r['aL_track'] - 1)
        # Y_L-only threshold: R_Y needed for S/Q(e->inf) = 3.4716 (bulk)
        Wb = ((3.4716 + 0.00428) / 0.190255) ** (1 / 0.75013) + 0.97494
        r['R_Y_threshold_bulk'] = Wb * eps ** (4 / 3) - 1
        rows[key] = r
        print(key, {k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items() if k != 'rel'})
        print('    rel', r['rel'])
    return rows


def leak(lrows):
    d = json.load(open(RES / 'leak_estimate.json'))['rows']
    Hs = {'HT6t': 1.29e11, 'HT6t_cal1': 2.805e10, 'T1': 1.29e11, 'T1_cal1': 1.043e10, 'TB1': 7.466e12, 'TB2': 1.0605e13}
    out = {}
    for key, r in d.items():
        L = math.log(MP / (30 * Hs[key]))
        c = L / (4 * math.pi ** 2)
        tr = r['c_MP']
        scale = (c / tr['c']) ** 2
        Rl = tr['R_leak_e1e6'] * scale
        lr = lrows.get(key)
        sq = None
        if lr:
            eps_W = lr['W'] / (1 + lr['R_total'])
            sq = SQ(eps_W * (1 + lr['R_total'] + Rl))
        # sensitivity of the near-resonant HT6t to m_sigma (+-2%)
        sens = None
        if key.startswith('HT6t'):
            mI = r['m_I']
            sens = {str(f): ((MSIG ** 2 - mI ** 2) / ((f * MSIG) ** 2 - mI ** 2)) ** 2 for f in (0.98, 1.02)}
        out[key] = dict(ln_MP_over_MY=L, c_MP_own=c, c_MP_track=tr['c'], R_leak_MP_track=tr['R_leak_e1e6'], R_leak_MP_own=Rl,
                        factor=scale, S_per_Q_with_leak_MP_own=sq, S_per_Q_with_leak_MP_track=tr['S_per_Q_e1e6'],
                        msigma_pm2pct_factor=sens)
        print('leak', key, out[key])
    return out


def iteration():
    fam = {'HT6t': [(1.29e11, 2.6985e10), (2.805e10, 1.0434e10)], 'T1': [(1.29e11, 9.8635e9), (1.043e10, 1.1945e10)]}
    out = {}
    for k, ((x0, y0), (x1, y1)) in fam.items():
        b = (math.log(y1) - math.log(y0)) / (math.log(x1) - math.log(x0))
        a = math.log(y0) - b * math.log(x0)
        Hfix = math.exp(a / (1 - b))
        mis0, mis1 = math.log(y0 / x0), math.log(y1 / x1)
        out[k] = dict(slope_dlnHreq_dlnHpot=b, H_fixed_point_secant=Hfix, log_mismatch=[mis0, mis1], contraction=abs(mis1 / mis0),
                      r_at_fixed=2 * Hfix ** 2 / (math.pi ** 2 * 2.1e-9 * MP ** 2))
        print('iteration', k, out[k])
    return out


def main():
    OUT['asymptote'] = asymptote(); print(OUT['asymptote'])
    OUT['ledger'] = ledger()
    OUT['leak'] = leak(OUT['ledger'])
    OUT['iteration'] = iteration()
    pk = json.load(open(RES / 'pmap_kin.json'))
    OUT['pmap_entropy_per_particle'] = {k: v['S_per_Q_closed_form'] / (1 + 2 * v['p']) for k, v in pk['rows'].items() if 'p' in v}
    print('entropy per particle in the p map:', OUT['pmap_entropy_per_particle'])
    json.dump(OUT, open(HERE.parent / 'results' / 'sk_ledger.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
