#!/usr/bin/env python3
"""mf_ideal item 3 (Claude, own code): heat ledger, S/Q(e), branch, e_crit, relic and H_i calibration for each idealised
history, on its own release, own vectors and own K.

Ledger (dq_heat's, as inf_chain chain.py; every quantity in units of the history's own E_X, a = 1 at the slice)
[C arithmetic on C inputs; the ledger itself is dq_heat's, E]:
  R_vec(e)  : vectors_mfi.py, log-log interpolated in e (computed 1e-7..1e-4), extrapolated with the last slope;
  R_Y       : (N_YL + N_YH) m_L <a_L>/E_X, <a_L> = (pi K/2 Gamma_L)^{1/2} (Gamma_L = 1.006e-7 H_i decays deep in radiation);
  R_casc(e) : (N_YH + 2 N_sigma) m_A <a_A>/E_X;
  R_YY      : 0 (the holder sits on s with m_s <= 4 H_i < 2 m_L = 44.5 H_i: s -> Y_L Y_L is closed; s is not the background).
S/Q: W = (E/E_c)(1 + R_total), S/Q = 0.190255 (W - 0.97494)^0.75013 - 0.00428 [Q dq_entropy; checked to W = 1.3e5 by
     inf_chain, and here by direct Hartree solves (hartree_check)].
Branch: free X if S/Q > 3.4716 (bulk); Q = 2e5 objects need S/Q < 2.9284 [Q fx_chain skeptic].
Relic: free X -> p(S/Q) from wt_critic's gas map (S/Q <= 6.74) and inf_chain's case-L runs (17-1308), with inf_chain's
     self-consistent-yield correction; M = 1 + 2p.  Q-matter -> M = 0.79-0.82 m_X [Q dq_critic / qm_micro], <Q> ~ 1e20.
H_i: Y/Y_Q(H_M1) = |N_X| H_M1^3/(a^3 s) with s from rho_r a^4 = 3 K^2 (g* = 106.75); H_i,req = H_M1 (0.9374/(M Y/Y_Q))^{2/3};
     r = 2 H_i^2/(pi^2 A_s M_P^2) (constant-H inflation: H_* = H_i).
Isocurvature: Omega_DM ∝ M |N_X| K^{-3/2} at fixed H_i; d ln Omega/d ln s_* from the c0 = 1 +- 0.01 releases (fixed g_s,
     preparation re-solved) with their own vectors.  Planck 2018 uncorrelated CDM isocurvature beta_iso < 0.038 (k = 0.05/Mpc)
     [R, arXiv:1807.06211] -> s_* > |D| H_i/(2 pi sqrt(beta A_s/(1 - beta))).
Usage: python3 -B chain_mfi.py SET|TAG ...   (SET = P, D, I, or tags)   -> ../results/chain/chain_<tag>.json
"""
import json, math, sys
import numpy as np
import mfi_common as mc

SQB, SQ2 = 3.4716191256166837, 2.9283778774213727
M_L, G_L = 22.24859546128699, 1.005903273746766e-07
AS = 2.100e-9
BETA_ISO = 0.038
OUT = mc.RES / 'chain'


def SQ(W):
    return 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428


def W_of_SQ(s):
    return ((s + 0.00428) / 0.190255) ** (1 / 0.75013) + 0.97494


# ------------------------------------------------------------------ pair fraction (read-only data of wt_critic, inf_chain)
def p_rows():
    rows = []
    cc = json.load(open(mc.TRK / 'wt_critic' / 'results' / 'chain_consistency.json'))['p_scan_SQ']
    rows += [(r['S_per_Q'], r['p']) for r in cc]
    pk = json.load(open(mc.TRK / 'inf_chain' / 'results' / 'pmap_kin.json'))
    rows += [(r['S_per_Q_closed_form'], r['p']) for r in pk['rows'].values() if 'p' in r]
    rows.sort()
    sc = sorted((float(k), v['p_self_consistent'] / v['p_fixed_yield']) for k, v in pk['self_consistent'].items())
    return rows, sc


ROWS, SCR = p_rows()


def p_of_SQ(s):
    x = np.log([r[0] for r in ROWS]); y = np.log([r[1] for r in ROWS]); lx = math.log(s)
    if lx <= x[0]:
        return math.exp(y[0] + (y[1] - y[0]) / (x[1] - x[0]) * (lx - x[0]))
    if lx >= x[-1]:
        return math.exp(y[-1] + (y[-1] - y[-2]) / (x[-1] - x[-2]) * (lx - x[-1]))
    p0 = math.exp(float(np.interp(lx, x, y)))
    # self-consistent-yield correction (inf_chain pmap_kin), -> 1 below S/Q = 6.74
    xs = [math.log(6.742)] + [math.log(a) for a, _ in SCR]; ys = [1.0] + [b for _, b in SCR]
    f = 1.0 if lx <= xs[0] else float(np.interp(lx, xs, ys))
    return p0 * f


def M_qmatter(s):
    S = [1.5, 2.5, 3.0, 3.2]
    lo = [0.7883, 0.7974, 0.8023, 0.8048]; hi = [0.7944, 0.8076, 0.8152, 0.8182]
    return float(np.interp(s, S, lo)), float(np.interp(s, S, hi))


def relic(s):
    if s > SQB:
        p = p_of_SQ(s)
        return dict(branch='free X', p=p, M=1 + 2 * p, asym_share=1 / (1 + 2 * p))
    lo, hi = M_qmatter(s)
    return dict(branch='Q-matter' + ('' if s < SQ2 else ' (bulk only; no Q = 2e5 objects)'), M=0.5 * (lo + hi), M_range=[lo, hi],
                Q_mean='~1e20 (0.4-1.5e20) [Q dq_critic/qm_micro]', vapour='0.12-0.175 of the charge [Q]')


def yield_YQ(NX, K):
    rho_a4 = 3 * K ** 2 * mc.MP ** 2 * mc.H_M1 ** 2
    T_a = (30 * rho_a4 / (math.pi ** 2 * mc.GSTAR)) ** 0.25
    s_a3 = 2 * math.pi ** 2 / 45 * mc.GSTAR * T_a ** 3
    return abs(NX) * mc.H_M1 ** 3 / s_a3 / mc.YQ_TARGET


def H_of(M, yyq):
    return mc.H_M1 * (0.9374 / (M * yyq)) ** (2.0 / 3.0)


def r_of(H):
    return 2 * H * H / (math.pi ** 2 * AS * mc.MP ** 2)


# ------------------------------------------------------------------ ledger
def ledger(tag):
    rel = json.load(open(mc.RES / 'release' / ('release_%s.json' % tag)))
    vec = json.load(open(mc.RES / 'vectors' / ('vectors_%s.json' % tag)))
    r = rel['result']
    K = vec['background']['K']
    NX, EX, eps = abs(r['N_X']), r['E_X'], r['eps']
    E_Ec = eps ** (-4 / 3)
    aL = math.sqrt(math.pi * K / (2 * G_L))
    DL = M_L * aL
    es = np.array([x['e'] for x in vec['by_e']]); Rv = np.array([x['R_vec'] for x in vec['by_e']])
    dA = np.array([x['deposit_per_vector'] for x in vec['by_e']])
    le, lR, lD = np.log(es), np.log(Rv), np.log(dA)

    def interp(lx, X, Y):
        if lx <= X[0]:
            return Y[0] + (Y[1] - Y[0]) / (X[1] - X[0]) * (lx - X[0])
        if lx >= X[-1]:
            return Y[-1] + (Y[-1] - Y[-2]) / (X[-1] - X[-2]) * (lx - X[-1])
        return float(np.interp(lx, X, Y))

    def comps(e):
        if not np.isfinite(e):
            RY = (r['N_YL'] + r['N_YH']) * DL / EX
            return dict(R_vec=0.0, R_Y=RY, R_casc=0.0, R_total=RY)
        Rvec = math.exp(interp(math.log(e), le, lR)); DA = math.exp(interp(math.log(e), le, lD))
        RY = (r['N_YL'] + r['N_YH']) * DL / EX
        Rc = (r['N_YH'] + 2 * r['N_sigma']) * DA / EX
        return dict(R_vec=Rvec, R_Y=RY, R_casc=Rc, R_total=Rvec + RY + Rc)
    return dict(rel=rel, vec=vec, K=K, NX=NX, EX=EX, eps=eps, E_Ec=E_Ec, aL=aL, DL=DL, comps=comps,
                N_A_1e6=[x['N_A'] for x in vec['by_e'] if abs(x['e'] - 1e-6) < 1e-12][0])


def analyse(tag):
    L = ledger(tag)
    egrid = np.geomspace(1e-7, 1e-2, 181)
    rows = []
    for e in egrid:
        c = L['comps'](e); W = L['E_Ec'] * (1 + c['R_total'])
        rows.append(dict(e=float(e), W=W, S_per_Q=SQ(W), **c))

    def ecrit(thr):
        for r1, r2 in zip(rows, rows[1:]):
            if (r1['S_per_Q'] - thr) * (r2['S_per_Q'] - thr) < 0:
                return math.exp(math.log(r1['e']) + (math.log(r2['e']) - math.log(r1['e'])) * (r1['S_per_Q'] - thr) / (r1['S_per_Q'] - r2['S_per_Q']))
        return None
    cinf = L['comps'](float('inf')); SQinf = SQ(L['E_Ec'] * (1 + cinf['R_total']))
    yyq = yield_YQ(L['NX'], L['K'])
    out = dict(tag=tag, spec=L['rel']['spec'], K=L['K'], N_X=-L['NX'], E_X=L['EX'], eps=L['eps'], E_over_Ec=L['E_Ec'],
               N_YL=L['rel']['result']['N_YL'], N_YH=L['rel']['result']['N_YH'], N_sigma=L['rel']['result']['N_sigma'],
               checks=dict(gauss=L['rel']['result']['gauss'], hidden_work=L['rel']['result']['hidden_work'],
                           H_consistency=L['rel']['result']['H_consistency'], N_X_err=L['rel']['result']['N_X_err'],
                           charge_balance=L['rel']['result']['charge_balance'], min_mL2=L['rel']['result']['min_mL2']),
               N_A_1e6=L['N_A_1e6'], deposit_per_YL_over_EX=L['DL'] / L['EX'], Y_over_YQ_at_HM1=yyq,
               S_per_Q_e_inf=SQinf, R_Y=cinf['R_Y'], e_crit_bulk=ecrit(SQB), e_crit_Q2e5=ecrit(SQ2), rows=rows[::6])
    for lab, e in (('e1e-6', 1e-6), ('e3e-6', 3e-6)):
        c = L['comps'](e); W = L['E_Ec'] * (1 + c['R_total']); s = SQ(W)
        rr = relic(s)
        rr['H_i_req_GeV'] = H_of(rr['M'], yyq); rr['r'] = r_of(rr['H_i_req_GeV'])
        out[lab] = dict(**c, W=W, S_per_Q=s, margin=s / SQB, **rr)
    return out


def hartree_check(NX, EX, Ws):
    """direct relaxed-Hartree solves (dq_entropy's copied solver, read-only import) against the closed form."""
    sys.path.insert(0, str(mc.TRK / 'dq_entropy' / 'code'))
    import hartree_dq as hq
    hq.set_release(NX, EX)
    out = []; g = None; Wp = None
    E_Ec = EX / (abs(NX) / mc.C_HART) ** (4 / 3)
    for W in sorted(Ws):
        try:
            if Wp is not None:
                for Wm in np.geomspace(Wp, W, int(math.ceil(math.log(W / Wp) / math.log(1.3))) + 1)[1:-1]:
                    gg = list(g); gg[0] += 0.25 * math.log(Wm / Wp)
                    st, g = hq.solve(Wm / E_Ec - 1, True, 0.0, gg); Wp = Wm
                g = list(g); g[0] += 0.25 * math.log(W / Wp)
            st, g = hq.solve(W / E_Ec - 1, True, 0.0, g); Wp = W
            sq = st['S'] / abs(NX)
            out.append(dict(W=W, S_per_Q_hartree=sq, S_per_Q_closed=SQ(W), rel=sq / SQ(W) - 1, cond_frac=st['cond_frac'],
                            Xbar_per_Q=st['n_minus'] / abs(NX), resid=max(abs(x) for x in st['resid'])))
        except Exception as ex:
            out.append(dict(W=W, error=repr(ex)))
    return out


def isocurvature(K, ton, ms=3.0):
    """D = d ln Omega_DM / d ln s_* at e = 1e-6 and 3e-6 (fixed H_i, fixed g_s), from c0 = 0.99, 1.01."""
    from grid import tag as gtag
    res = {}
    for lab in ('e1e-6', 'e3e-6'):
        vals = {}
        for c0 in (0.99, 1.01):
            o = analyse(gtag(K, ms, ton, c0=c0))
            M = o[lab]['M']
            vals[c0] = dict(lnOmega=math.log(M * abs(o['N_X']) * o['K'] ** -1.5), S_per_Q=o[lab]['S_per_Q'], branch=o[lab]['branch'],
                            N_X=o['N_X'], M=M)
        D = (vals[1.01]['lnOmega'] - vals[0.99]['lnOmega']) / (math.log(1.01) - math.log(0.99))
        DN = (math.log(abs(vals[1.01]['N_X'])) - math.log(abs(vals[0.99]['N_X']))) / (math.log(1.01) - math.log(0.99))
        amp = math.sqrt(BETA_ISO * AS / (1 - BETA_ISO))
        res[lab] = dict(D=D, D_charge_only=DN, pm=vals,
                        s_star_min_GeV_at_Hi=dict(note='s_* > |D| H_i/(2 pi sqrt(beta A_s/(1-beta))), beta = 0.038',
                                                  per_H_i=abs(D) / (2 * math.pi * amp)))
    return res


def main(args):
    from grid import grid
    tags = []
    for a in args:
        tags += [j['tag'] for j in grid(a)] if a in ('P', 'D', 'I', 'F', 'X') else [a]
    OUT.mkdir(parents=True, exist_ok=True)
    for t in tags:
        o = analyse(t)
        json.dump(o, open(OUT / ('chain_%s.json' % t), 'w'), indent=1, default=float)
        a, b = o['e1e-6'], o['e3e-6']
        print('%-18s K=%.3g eps=%.4f E/Ec=%6.2f N_YL=%7.4g | 1e-6: Rv=%6.3g RY=%6.3g Rc=%6.2g W=%7.1f S/Q=%6.3f %-9s | 3e-6: S/Q=%6.3f %-9s | inf %6.3f | ecrit %s / %s | M %.3g H_i %.3g' % (
            t, o['K'], o['eps'], o['E_over_Ec'], o['N_YL'], a['R_vec'], a['R_Y'], a['R_casc'], a['W'], a['S_per_Q'], a['branch'][:9],
            b['S_per_Q'], b['branch'][:9], o['S_per_Q_e_inf'],
            ('%.3g' % o['e_crit_bulk']) if o['e_crit_bulk'] else '-', ('%.3g' % o['e_crit_Q2e5']) if o['e_crit_Q2e5'] else '-',
            a['M'], a['H_i_req_GeV']), flush=True)


if __name__ == '__main__':
    main(sys.argv[1:])
