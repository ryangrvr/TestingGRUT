#!/usr/bin/env python3
"""inf_chain item 3 (Claude, own code): heat ledger, S/Q(e), branch, e_crit, relic, H_i calibration and r per benchmark.

Ledger (dq_heat's, on the benchmark's own release and own K; every quantity in units of the benchmark's own E_X,
a = 1 at the slice) [C arithmetic on C inputs; the ledger itself is dq_heat's, E]:
  R_vec(e)    : vectors.py (own background, exit factor along the benchmark path), log-log interpolated in e;
  R_Y         : (N_YL + N_YH) m_L <a_L>/E_X, each Y_H giving one Y_L parent (Y_H -> Y_L + A_L), <a> = int Gamma e^{-Gamma t} a dt
                on the benchmark's own a(t) table (radiation (pi K/2 Gamma)^{1/2} when the decay is late);
  R_casc(e)   : (N_YH + 2 N_sigma) m_A <a_A>/E_X (cascade and sigma -> A_L A_L vectors decay with the primordial ones);
  R_YY(e)     : optional, the holder-induced phi_I -> YY Born counts of inf_potential (yy_harmonics.json, Gamma = 0.5,
                g_I at S1; converted e^3 from its eps_H = 1 normalisation) with the same deposits [Q inf_potential, E].
S/Q: W = (E/E_c)(1 + R_total), E/E_c = eps^{-4/3}; S/Q = 0.190255 (W - 0.97494)^0.75013 - 0.00428 [Q dq_entropy, fitted on
     W = 13-140], checked here against the copied relaxed-Hartree solver up to W = 3e4 (hartree_check()).
Branch: free X if S/Q > 3.4716 (bulk), objects of Q = 2e5 only if S/Q < 2.9284 [Q fx_chain skeptic thresholds].
Relic: free X -> pair fraction p(S/Q) from (i) wt_critic's gas-only map (S/Q <= 6.74) and (ii) our runs of the reviewed
     case-L moment network at large heat (pmap_kin.py, S/Q = 17-1308), log-log interpolated; mass per charge M = 1 + 2p.
     Q-matter -> M from qm_micro SKEPTIC via dq_critic (0.79-0.82), <Q> ~ 1e20, f_v 0.12-0.175 [Q].
H_i: Y/Y_Q at H_M1 = |N_X| H_M1^3/(a^3 s) (s from the radiation, g* = 106.75) [fx_holder convention]; Y ∝ H_i^{3/2} at fixed
     dimensionless history: H_i,req = H_M1 (0.9374/(M Y/Y_Q))^{2/3}.  r = 2 H_*^2/(pi^2 A_s M_P^2), H_* from the potential.
Usage: python3 -B chain.py KEY [gam]   (TB*: gam = Gamma_I/H_i of the fluid continuation)   -> ../results/chain_<tag>.json
"""
import json, math, os, sys
import numpy as np
import ch_common as ch

SQB, SQ2 = 3.4716191256166837, 2.9283778774213727
M_L, G_L = 22.24859546128699, 1.005903273746766e-07
M_H, G_H = 36.124783736376884, 0.0350420463551497
AS = 2.100e-9


def SQ(W):
    return 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428


def W_of_SQ(s):
    return ((s + 0.00428) / 0.190255) ** (1 / 0.75013) + 0.97494


def table(path):
    with open(path) as fh:
        n = int(fh.readline())
        d = np.loadtxt(fh, max_rows=n)
    return d[:, 0], d[:, 1]


def mean_a(gam, t, la):
    """<a> = int_0^inf gam e^{-gam t} a(t) dt on the table, radiation (a ∝ t^{1/2}) beyond its end."""
    a = np.exp(la)
    tt = np.concatenate([[0.0], t]); aa = np.concatenate([[1.0], a])
    w = np.exp(-gam * tt)
    # integral of a dW with W = 1 - e^{-gam t}
    W = 1 - w
    inside = float(np.sum(0.5 * (aa[1:] + aa[:-1]) * np.diff(W)))
    tail_t, tail_a = tt[-1], aa[-1]
    # beyond: a = tail_a (t/tail_t)^{1/2}: int_T^inf gam e^{-gam t} a dt
    from scipy.special import gammaincc, gamma
    x = gam * tail_t
    tail = tail_a / math.sqrt(tail_t) * gam ** -0.5 * gamma(1.5) * gammaincc(1.5, x)
    return inside + tail


def p_map():
    rows = []
    cc = json.load(open(ch.TRK / 'wt_critic' / 'results' / 'chain_consistency.json'))['p_scan_SQ']
    rows += [(r['S_per_Q'], r['p'], 'wt_critic gas-only map') for r in cc]
    pk = json.load(open(ch.RES / 'pmap_kin.json'))['rows']
    rows += [(r['S_per_Q_closed_form'], r['p'], 'case-L network, a_ref start (pmap_kin)') for r in pk.values() if 'p' in r]
    rows.sort()
    return rows


def p_of_SQ(s, rows=None):
    rows = rows or p_map()
    x = np.log([r[0] for r in rows]); y = np.log([r[1] for r in rows])
    lx = math.log(s)
    if lx <= x[0]:
        sl = (y[1] - y[0]) / (x[1] - x[0]); return math.exp(y[0] + sl * (lx - x[0]))
    if lx >= x[-1]:
        sl = (y[-1] - y[-2]) / (x[-1] - x[-2]); return math.exp(y[-1] + sl * (lx - x[-1]))
    return math.exp(float(np.interp(lx, x, y)))


def sc_ratio(s):
    """p_self-consistent / p_fixed-yield from pmap_kin.py's self-consistent runs (the network's visible clock run with the
    calibrated charge yield Y_Q 0.9374/M); -> 1 where p << 1 (anchored at the map's S/Q = 6.74)."""
    d = json.load(open(ch.RES / 'pmap_kin.json')).get('self_consistent', {})
    pts = sorted((float(k), v['p_self_consistent'] / v['p_fixed_yield']) for k, v in d.items())
    x = [math.log(6.742)] + [math.log(a) for a, _ in pts]; y = [1.0] + [b for _, b in pts]
    lx = math.log(s)
    if lx <= x[0]:
        return 1.0
    if lx >= x[-1]:
        return y[-1] + (y[-1] - y[-2]) / (x[-1] - x[-2]) * (lx - x[-1])
    return float(np.interp(lx, x, y))


def M_qmatter(s):
    S = [1.5, 2.5, 3.0, 3.2]
    lo = [0.7883, 0.7974, 0.8023, 0.8048]; hi = [0.7944, 0.8076, 0.8152, 0.8182]
    return float(np.interp(s, S, lo)), float(np.interp(s, S, hi))


def H_of(M, yyq):
    return ch.H_M1 * (0.9374 / (M * yyq)) ** (2.0 / 3.0)


def r_of(H):
    return 2 * H * H / (math.pi ** 2 * AS * ch.MP ** 2)


def hartree_check(NX, EX, Ws):
    """direct relaxed-Hartree solve (copied solver dq_entropy/code/hartree_dq.py, read-only import) vs the closed form."""
    sys.path.insert(0, str(ch.TRK / 'dq_entropy' / 'code'))
    import hartree_dq as hq
    hq.set_release(NX, EX)
    out = []; g = None; Wp = None
    E_Ec = EX / (abs(NX) / ch.C_HART) ** (4 / 3)
    for W in sorted(Ws):
        R = W / E_Ec - 1
        try:
            # continuation in W (steps <= x1.3, guess T scaled as W^{1/4}, the other unknowns held)
            if Wp is not None:
                for Wm in np.geomspace(Wp, W, int(math.ceil(math.log(W / Wp) / math.log(1.3))) + 1)[1:-1]:
                    gg = list(g); gg[0] += 0.25 * math.log(Wm / Wp)
                    st, g = hq.solve(Wm / E_Ec - 1, True, 0.0, gg); Wp = Wm
                g = list(g); g[0] += 0.25 * math.log(W / Wp)
            st, g = hq.solve(R, True, 0.0, g); Wp = W
            sq = st['S'] / abs(NX)
            out.append(dict(W=W, R=R, S_per_Q_hartree=sq, S_per_Q_closed=SQ(W), rel=sq / SQ(W) - 1, Xbar_per_Q=st['n_minus'] / abs(NX),
                            cond_frac=st['cond_frac'], resid=max(abs(x) for x in st['resid']), success=st['success']))
        except Exception as e:
            out.append(dict(W=W, error=repr(e)))
    return out


def ledger(key, gam=None, with_yy=True):
    rel = json.load(open(ch.RES / ('release_%s.json' % key)))
    res = rel['result']
    vtag = key + ('' if gam is None else '_g%.3e' % gam)
    vec = json.load(open(ch.RES / ('vectors_%s.json' % vtag)))
    bgv = vec['background']
    K = bgv['K']
    NX, EX, eps = abs(res['N_X']), res['E_X'], res['eps']
    E_Ec = eps ** (-4 / 3)
    t, la = table(ch.SP / ('table_%s.txt' % vtag))
    aL = mean_a(G_L, t, la)
    DL = M_L * aL
    NYL, NYH, NS = res['N_YL'], res['N_YH'], res['N_sigma']
    yy = None
    fy = ch.RES / ('yy_born_%s.json' % key)
    if with_yy and fy.exists():
        r = json.load(open(fy))
        yy = dict(N_YL=r['N_YL'], N_YH=r['N_YH'], tail=r['tail_NY_estimate'] or 0.0, source='own Born count yy_born.py (slice normalisation)')
    elif with_yy and key in ('HT6t', 'T1', 'HT6t_Send'):
        base = 'HT6t' if key.startswith('HT6t') else key
        conv = 'Send' if key == 'HT6t_Send' else 'S1'
        r = [x for x in json.load(open(ch.TRK / 'inf_potential' / 'results' / 'yy_harmonics.json'))['runs']
             if x['benchmark'] == base and x['gamma_over_Hi'] == 0.5][0]['res'][conv]
        n3 = math.e ** 3                                     # eps_H = 1 -> S1 normalisation of comoving numbers
        tail = r['tail_NY_estimate'] or 0.0
        yy = dict(N_YL=r['N_YL'] * n3, N_YH=r['N_YH'] * n3, tail=tail * n3, source='inf_potential yy_harmonics.json (%s, Gamma = 0.5), x e^3' % conv)
    es = np.array([x['e'] for x in vec['by_e']]); Rv = np.array([x['R_vec'] for x in vec['by_e']])
    dA = np.array([x['deposit_per_vector'] for x in vec['by_e']])

    def comps(e, yyon=True):
        Rvec = math.exp(np.interp(math.log(e), np.log(es), np.log(Rv)))
        DA = math.exp(np.interp(math.log(e), np.log(es), np.log(dA)))
        RY = (NYL + NYH) * DL / EX
        Rc = (NYH + 2 * NS) * DA / EX
        Ryy = 0.0
        if yy and yyon:
            Ryy = ((yy['N_YL'] + yy['N_YH'] + yy['tail']) * DL + (yy['N_YH'] + 0.5 * yy['tail']) * DA) / EX
        return dict(R_vec=Rvec, R_Y=RY, R_casc=Rc, R_YY=Ryy, R_total=Rvec + RY + Rc + Ryy)
    return dict(rel=rel, vec=vec, K=K, NX=NX, EX=EX, eps=eps, E_Ec=E_Ec, aL=aL, DL=DL, yy=yy, comps=comps, es=es, vtag=vtag)


def yield_YQ(L, gam=None):
    """Y/Y_Q at H_M1 from N_X and the late radiation constant Rr = rho_r a^4 (= 3 K^2 in radiation)."""
    res = L['rel']['result']
    if gam is None and 'Y_over_YQ_at_HM1_frozen' in res:
        return res['Y_over_YQ_at_HM1_frozen']
    Rr = 3 * L['K'] ** 2
    rho_rad_a4 = Rr * ch.MP ** 2 * ch.H_M1 ** 2              # GeV^4 x a^4
    T_a = (30 * rho_rad_a4 / (math.pi ** 2 * ch.GSTAR)) ** 0.25
    s_a3 = 2 * math.pi ** 2 / 45 * ch.GSTAR * T_a ** 3
    return L['NX'] * ch.H_M1 ** 3 / s_a3 / ch.YQ_TARGET


def analyse(key, gam=None, with_yy=True, hcheck=False):
    L = ledger(key, gam, with_yy)
    egrid = np.geomspace(1e-7, L['es'][-1], 161)
    rows = []
    for e in egrid:
        c = L['comps'](e)
        W = L['E_Ec'] * (1 + c['R_total'])
        c0 = L['comps'](e, yyon=False); W0 = L['E_Ec'] * (1 + c0['R_total'])
        rows.append(dict(e=float(e), W=W, S_per_Q=SQ(W), S_per_Q_noYY=SQ(W0), **c))

    def ecrit(thr, k='S_per_Q'):
        for r1, r2 in zip(rows, rows[1:]):
            if (r1[k] - thr) * (r2[k] - thr) < 0:
                x = math.log(r1['e']) + (math.log(r2['e']) - math.log(r1['e'])) * (r1[k] - thr) / (r1[k] - r2[k])
                return math.exp(x)
        return None
    c_inf = L['comps'](1e12) if False else None
    # e -> infinity: vectors and A_L deposits vanish
    RYinf = (L['rel']['result']['N_YL'] + L['rel']['result']['N_YH']) * L['DL'] / L['EX']
    if L['yy'] and with_yy:
        RYinf += (L['yy']['N_YL'] + L['yy']['N_YH'] + L['yy']['tail']) * L['DL'] / L['EX']
    SQinf = SQ(L['E_Ec'] * (1 + RYinf))
    r6 = [r for r in rows if abs(r['e'] - 1e-6) / 1e-6 < 0.03]
    c6 = L['comps'](1e-6); W6 = L['E_Ec'] * (1 + c6['R_total']); s6 = SQ(W6)
    c60 = L['comps'](1e-6, yyon=False); s60 = SQ(L['E_Ec'] * (1 + c60['R_total']))
    out = dict(key=key, gam=gam, vtag=L['vtag'], K=L['K'], N_X=L['NX'], E_X=L['EX'], eps=L['eps'], E_over_Ec=L['E_Ec'],
               a_mean_L=L['aL'], deposit_per_YL=L['DL'], deposit_per_YL_over_EX=L['DL'] / L['EX'], yy=L['yy'],
               e1e6=dict(**c6, W=W6, S_per_Q=s6, S_per_Q_noYY=s60), S_per_Q_e_inf=SQinf,
               e_crit_bulk=ecrit(SQB), e_crit_Q2e5=ecrit(SQ2), e_crit_bulk_noYY=ecrit(SQB, 'S_per_Q_noYY'),
               rows=rows[::4])
    # relic at e = 1e-6
    rowsP = p_map()
    def relic(s):
        if s > SQB:
            p0 = p_of_SQ(s, rowsP); p = p0 * sc_ratio(s)
            return dict(branch='free X', p=p, p_fixed_yield=p0, M=1 + 2 * p, M_fixed_yield=1 + 2 * p0)
        lo, hi = M_qmatter(s)
        return dict(branch='Q-matter', M=[lo, hi], Q_mean='~1e20 (0.4-1.5e20) [Q dq_critic / qm_micro]', f_v='0.12-0.175 (charge) [Q]')
    yyq = yield_YQ(L, gam)
    out['Y_over_YQ_at_HM1'] = yyq
    for lab, s in (('relic_e1e6', s6), ('relic_e1e6_noYY', s60)):
        rr = relic(s)
        M = rr['M'] if isinstance(rr['M'], float) else 0.5 * sum(rr['M'])
        rr['H_i_req_GeV'] = H_of(M, yyq)
        if 'M_fixed_yield' in rr:
            rr['H_i_req_fixed_yield_GeV'] = H_of(rr['M_fixed_yield'], yyq)
            rr['H_i_req_pmap_x2_GeV'] = [H_of(1 + 2 * rr['p'] * 2, yyq), H_of(1 + 2 * rr['p'] / 2, yyq)]
        out[lab] = rr
    if hcheck:
        out['hartree_check'] = hartree_check(L['NX'], L['EX'], [L['E_Ec'] * 1.01, 50.0, 200.0, 1000.0, 5000.0, min(3e4, max(W6, 6000.0))])
    return out


def main(key, gam=None):
    o = analyse(key, gam, hcheck=False)      # closed form checked separately (sq_check.py)
    tag = o['vtag']
    json.dump(o, open(ch.RES / ('chain_%s.json' % tag), 'w'), indent=1, default=float)
    c = o['e1e6']
    print('%s: K=%.4g eps=%.4f E/Ec=%.3f | e=1e-6: R_vec=%.4g R_Y=%.4g R_casc=%.4g R_YY=%.4g R=%.4g W=%.4g S/Q=%.4g (noYY %.4g) | e->inf S/Q=%.4g | e_crit bulk=%s Q2e5=%s (noYY %s)' % (
        tag, o['K'], o['eps'], o['E_over_Ec'], c['R_vec'], c['R_Y'], c['R_casc'], c['R_YY'], c['R_total'], c['W'], c['S_per_Q'],
        c['S_per_Q_noYY'], o['S_per_Q_e_inf'], o['e_crit_bulk'], o['e_crit_Q2e5'], o['e_crit_bulk_noYY']))
    print('   relic(1e-6): %s | noYY: %s | Y/YQ(HM1)=%.5g' % (o['relic_e1e6'], o['relic_e1e6_noYY'], o['Y_over_YQ_at_HM1']))
    if 'hartree_check' in o:
        for h in o['hartree_check']:
            print('   hartree W=%.4g: %s' % (h['W'], ('S/Q %.5g closed %.5g rel %.2e Xbar/Q %.3f resid %.1e' % (h['S_per_Q_hartree'], h['S_per_Q_closed'], h['rel'], h['Xbar_per_Q'], h['resid'])) if 'rel' in h else h.get('error')))


if __name__ == '__main__':
    k = sys.argv[1]
    g = float(sys.argv[2]) if len(sys.argv) > 2 else None
    main(k, g)
