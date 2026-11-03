#!/usr/bin/env python3
"""e2e_chain item 2 (Claude): the full chain on a spliced history at holder strength c0 and phase delta, with the
corrected relic step and the hidden-heat items.

Per (history, c0, delta):
  1. release: mf_ideal release_mfi.release on the MODE 5 table of the history (bench_run.py builds the table at c0 = 1;
     other c0 reuse it, MFI_SM2TAB = 1, as mf_critic hyb_c0.py / the mf_chain skeptic's sk_c0chain.py).  The holder at
     the slice is u0 = mu + c0 e^{i delta} (10i - mu): delta rotates the holder coupling's phase (mfi_common.U_PREP is
     patched to mu + e^{i delta}(10i - mu) for the call and restored) [model choice].
  2. vectors (vectors_mfi) and the dq_heat ledger + closed-form S/Q (chain_mfi.analyse), read-only.
  3. corrected relic [e2e_relic, read-only]: free X -> M = 1 + 2p with p from pmap.p_selfconsistent (relaxed start,
     abundance-calibrated yield); Q-matter -> relic.relic(S/Q) (abundance-calibrated <Q>, f_v, M_central).  H_i,req =
     chain_mfi.H_of(M, Y/Y_Q / D_s) (Y ~ H_i^{3/2} at fixed dimensionless history), D_s the entropy ratio of the s decay.
  4. hidden heat (mf_chain post.py formulae, re-typed; [E]): the s leak with A = 2 mu Re(G_eff) c_L, G_eff = c0 e^{i delta}
     (10i - mu); Gamma_leak = (A/s_*)^2 m_s^3/[32 pi (m_sig^2 - m_s^2)^2], n_s a^3 = 3 mP^2 K^2 (f_s/a)/m_s, BR = Gamma_leak/
     Gamma_s, R_leak = 2 n_s a^3 BR dep_A / E_X; gravitational production at the abrupt end from mf_chain's grav run on the
     Pauli history (grav_itHFa_P0_it2, same waterfall) [Q], scaled to this run's E_X.  c_L is stated with every number:
     c_L = 0 (portals vanish at M_Y, no threshold piece: the declared condition) and 1/16 pi^2 (a threshold-size variant).
Labels: [C] given the model choices; the ledger and closed form are dq_heat / dq_entropy's [Q].
"""
import json, math, cmath
import numpy as np
import e2c_common as E
import pmap as PM                                         # e2e_relic, read-only
import relic as RL                                        # e2e_relic, read-only

cm, rm, vm, mc = E.cm, E.rm, E.vm, E.mc
MU = mc.MU
MSIG = 40.0
GRAV = json.load(open(E.TRK / 'mf_chain' / 'results' / 'grav' / 'grav_itHFa_P0_it2.json'))['fields']
CLS = {'cL0_vanish_at_MY': 0.0, 'cL_1/16pi^2': 1 / (16 * math.pi ** 2)}


def col(h, k):
    return np.array([np.nan if v is None else float(v) for v in h[k]])


def fs_over_a(d):
    h = d['history']
    t, fs, a, fr = col(h, 't'), col(h, 'f_s'), col(h, 'a'), col(h, 'f_rad')
    m = (t > 100) & (fr > 0.999) & np.isfinite(fs) & (fs > 0)
    if m.sum() < 3:
        m = (t > 100) & np.isfinite(fs) & (fs > 0)
    return float(np.median(fs[m] / a[m]))


def release_c0(base_tag, tag, c0, delta_deg=0.0):
    rel0 = json.load(open(E.ENG / 'release' / ('release_%s.json' % base_tag)))
    sp = rel0['spec']
    orig_env, orig_U = rm.spec_env, mc.U_PREP

    def spec_env(*a, **kw):
        env, gam = orig_env(*a, **kw)
        env['MFI_SM2TAB'] = '1'
        return env, gam
    rm.spec_env = spec_env
    mc.U_PREP = MU + cmath.exp(1j * math.radians(delta_deg)) * (10j - MU)
    try:
        rm.release(tag, sp['K_target'], sp['m_s'], 0.0, gs=sp['Gamma_s'], sstar=sp['s_star_GeV'], c0=c0, bg=sp['bg'],
                   Tl=rel0['result']['T_long'], bcap=sp['step_cap_b'])
    finally:
        rm.spec_env, mc.U_PREP = orig_env, orig_U
    vm.main(tag)
    o = cm.analyse(tag)
    E.dump(E.ENG / 'chain' / ('chain_%s.json' % tag), o)
    return o


def relic_corr(sq, H_i):
    if sq > cm.SQB:
        p, f = PM.p_selfconsistent(sq)
        return dict(branch='free X', p=p, M=1 + 2 * p, asym_share=1 / (1 + 2 * p))
    r = RL.relic(sq, H_i=H_i)
    return dict(branch=r['branch'] + ('' if sq < cm.SQ2 else ' (bulk only; no Q = 2e5 objects)'), M=r['M_central'],
                M_range=r.get('M_range'), Q_mean=r.get('Q_mean'), f_v=r.get('f_v'))


def full(tag, bench_path, c0=1.0, delta_deg=0.0, E_X_c1=None):
    """corrected relic, hidden heat, H_req for one engine run (chain_<tag>.json must exist)."""
    ch = json.load(open(E.ENG / 'chain' / ('chain_%s.json' % tag)))
    vec = json.load(open(E.ENG / 'vectors' / ('vectors_%s.json' % tag)))
    d = json.load(open(bench_path))
    Hi = float(d['params']['H_i_GeV'])
    mP = E.MP_GEV / Hi
    K, EX, EEc = ch['K'], ch['E_X'], ch['E_over_Ec']
    ms = 3.0
    sstar = float(d['params']['s_star_MP']) * mP
    fsa = fs_over_a(d)
    ns_a3 = 3 * mP ** 2 * K ** 2 * fsa / ms
    Gs = float(d['spectator']['Gamma_s_assumed'])
    a_dec = math.sqrt(2 * K / Gs)
    f_dec = fsa * a_dec
    Dent = (1 + 0.8862 * f_dec) ** 0.75
    yyq = ch['Y_over_YQ_at_HM1']
    Geff = c0 * cmath.exp(1j * math.radians(delta_deg)) * (10j - MU)
    u0 = MU + Geff
    ge = np.array([x['e'] for x in vec['by_e']]); gd = np.array([x['deposit_per_vector'] for x in vec['by_e']])
    gscale = (E_X_c1 / EX) if E_X_c1 else 1.0
    gX = GRAV['X_m2_20.2']['E_total_over_EX'] * gscale
    dNYL, dNYH, dNs = GRAV['Y_L']['N_total'], GRAV['Y_H']['N_total'], GRAV['sigma_Sigma']['N_total']
    out = dict(tag=tag, c0=c0, delta_deg=delta_deg, mu_eff=[u0.real, u0.imag], abs_mu_eff=abs(u0), arg_mu_eff_deg=math.degrees(cmath.phase(u0)),
               H_i_GeV=Hi, K=K, eps=ch['eps'], E_over_Ec=EEc, N_X=ch['N_X'], E_X=EX, N_YL=ch['N_YL'], N_YH=ch['N_YH'],
               N_sigma=ch['N_sigma'], R_Y=ch['R_Y'], S_per_Q_e_inf=ch['S_per_Q_e_inf'], checks=ch['checks'],
               spectator=dict(f_s_over_a=fsa, n_s_a3=ns_a3, Gamma_s=Gs, a_decay=a_dec, f_s_at_decay=f_dec, entropy_ratio_D=Dent,
                              s_star_over_Hi=sstar),
               grav=dict(R_X=gX, dN_YL=dNYL, dN_YH=dNYH, dN_sigmaSigma=dNs, source='mf_chain grav_itHFa_P0_it2 [Q], scaled to E_X'),
               Y_over_YQ_at_HM1=yyq, by_e={})
    egrid = np.geomspace(1e-7, 1e-3, 121)
    L = cm.ledger(tag)
    for nm, cl in CLS.items():
        A = 2 * MU * Geff.real * cl
        Gl = (A / sstar) ** 2 * ms ** 3 / (32 * math.pi * (MSIG ** 2 - ms ** 2) ** 2)
        BR = Gl / Gs
        rows = []
        for e in egrid:
            c = L['comps'](e)
            DA = math.exp(float(np.interp(math.log(e), np.log(ge), np.log(gd))))
            RY_per = c['R_Y'] / max(ch['N_YL'] + ch['N_YH'], 1e-300)
            Rc_per = c['R_casc'] / max(ch['N_YH'] + 2 * ch['N_sigma'], 1e-300)
            dRg = gX + dNYL * RY_per + dNYH * (RY_per + Rc_per) + 2 * dNs * Rc_per
            Rl = 2 * ns_a3 * BR * DA / EX
            Rt = c['R_total'] + dRg + Rl
            rows.append(dict(c, e=float(e), S_per_Q=cm.SQ(EEc * (1 + Rt)), R_total=Rt, R_leak=Rl, dR_grav=dRg))

        def ecrit(thr):
            for r1, r2 in zip(rows, rows[1:]):
                if (r1['S_per_Q'] - thr) * (r2['S_per_Q'] - thr) < 0:
                    return math.exp(math.log(r1['e']) + (math.log(r2['e']) - math.log(r1['e'])) * (r1['S_per_Q'] - thr) / (r1['S_per_Q'] - r2['S_per_Q']))
            return None
        blk = dict(c_L=cl, A_Hi2=A, Gamma_leak=Gl, BR=BR, e_crit_bulk=ecrit(cm.SQB), e_crit_Q2e5=ecrit(cm.SQ2),
                   SQ_of_e=[dict(e=r['e'], S_per_Q=r['S_per_Q']) for r in rows[::10]])
        for lab, e in (('e1e-6', 1e-6), ('e3e-6', 3e-6)):
            r = min(rows, key=lambda x: abs(math.log(x['e'] / e)))
            c = L['comps'](e)
            DA = math.exp(float(np.interp(math.log(e), np.log(ge), np.log(gd))))
            RY_per = c['R_Y'] / max(ch['N_YL'] + ch['N_YH'], 1e-300)
            Rc_per = c['R_casc'] / max(ch['N_YH'] + 2 * ch['N_sigma'], 1e-300)
            dRg = gX + dNYL * RY_per + dNYH * (RY_per + Rc_per) + 2 * dNs * Rc_per
            Rl = 2 * ns_a3 * BR * DA / EX
            Rt = c['R_total'] + dRg + Rl
            sq = cm.SQ(EEc * (1 + Rt))
            rr = relic_corr(sq, Hi)
            Hreq = cm.H_of(rr['M'], yyq / Dent)
            blk[lab] = dict(R_vec=c['R_vec'], R_Y=c['R_Y'], R_casc=c['R_casc'], dR_grav=dRg, R_leak=Rl, R_total=Rt, S_per_Q=sq,
                            S_per_Q_engine=ch[lab]['S_per_Q'], M_inherited_map=ch[lab]['M'], H_req_inherited_map=ch[lab]['H_i_req_GeV'],
                            H_i_req_GeV=Hreq, H_req_over_H_i=Hreq / Hi, r_at_req=cm.r_of(Hreq), **rr)
        out['by_e'][nm] = blk
    return out


def line(o):
    b = o['by_e']['cL0_vanish_at_MY']
    a, c = b['e1e-6'], b['e3e-6']
    return ('%s c0 %.3f d %+.0f |u| %.2f arg %.1f: K %.3f eps %.4f E/Ec %.3f N_YL %.4g N_YH %.4g | 1e-6: Rv %.3g RY %.3g Rc %.3g '
            'S/Q %.4f %s M %.4f H_req %.4e (x%.3f) | 3e-6: S/Q %.4f %s | inf %.3f | e_crit %s' % (
                o['tag'], o['c0'], o['delta_deg'], o['abs_mu_eff'], o['arg_mu_eff_deg'], o['K'], o['eps'], o['E_over_Ec'], o['N_YL'],
                o['N_YH'], a['R_vec'], a['R_Y'], a['R_casc'], a['S_per_Q'], a['branch'][:8], a['M'], a['H_i_req_GeV'],
                a['H_req_over_H_i'], c['S_per_Q'], c['branch'][:8], o['S_per_Q_e_inf'],
                ('%.3g' % b['e_crit_bulk']) if b['e_crit_bulk'] else '-'))
