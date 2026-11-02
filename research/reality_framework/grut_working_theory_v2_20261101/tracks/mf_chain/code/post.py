#!/usr/bin/env python3
"""mf_chain: the final ledger on each benchmark history, with the items bench_run.py does not do.

For each case (exported bench at the design H_i, and the self-consistently calibrated histories):
  1. base chain [C]: mf_ideal's release, vectors, dq_heat ledger, closed-form S/Q (results/engine/chain).
  2. hidden heat from s decays [E]: the induced holder portal A (s/s_*)|Sigma|^2, A = 2 mu Re(G) c_L with G = 10i - mu
     (the |u|^2 - mu^2 structure of dq_entropy's induced portals), mixes s with the Sigma radial mode (40 H_i), and
     sigma (d pi)^2/h sends s -> A_L A_L:  Gamma_leak = (A/s_*)^2 m_s^3 / [32 pi (m_sig^2 - m_s^2)^2]  (derived; the
     inf_chain/inf_reheat-skeptic form with the canonical normalisation made explicit).  c_L brackets: 0 (portals
     vanish at M_Y, threshold pieces dropped), 1/16 pi^2, 1/4 pi^2 (threshold size), ln(M_P/M_Y)/4 pi^2 (vanish at M_P).
     s quanta: n_s a^3 = 3 mP^2 K^2 (f_s/a)/m_s from the history's own late f_s/a (radiation era); BR = Gamma_leak/Gamma_s;
     each A_L deposits the vector deposit m_A <a_A> of the history's own vector run.  R_leak = 2 n_s BR dep/E_X.
     Gamma_s = 1.3e-6 H_i x (H_i/1.3e11 GeV) (mf_hybrid's y_s s Nbar_s N_s at the light-s limit; y_s^2 ~ H_i/M_P under the
     rescaling rule) [model choice].
  3. s decay into the SM (N_s -> L H, a model choice; otherwise N_s is dark radiation) [D/E]: entropy ratio
     D = (1 + 0.8862 f_s(Gamma_s t = 1))^{3/4}, f_s(t) = (f_s/a) sqrt(2 K t) (subdominant matter decaying in radiation);
     the yield falls by D, so H_i,req rises by D^{2/3}.
  4. sigma and chi decays: no hidden channel (no sigma/chi-hidden portal, sequestered by choice; the N singlets carry
     no hidden charge) [A]; chi -> s s (HFa, Gamma = 9.6e-14 H_i) adds n_s ~ 1e3 against the condensate's ~1e10 [E].
  5. gravitational production at the abrupt end [C/E, grav.py]: finite-k X (heat E/E_X), Y_L, Y_H, sigma quanta.
  6. relic: free X -> inherited p(S/Q) map, M = 1 + 2p; Q-matter -> qm_micro with its skeptic's corrections (<Q> x0.54-0.70,
     M +0.4-1.0%) -- not reached on any history here (S/Q >= 10.8 even at e -> infinity).
  7. isocurvature [C]: D_iso from iso.py; beta_iso = (H_i/2 pi s_*)^2 D^2/A_s at the history's s_* (3e-3 M_P).
  8. mf_ideal's map at the same K and m_s = 3 (t_on = 0.1 and 1; linear in K, extrapolated beyond K = 3 and flagged).
Output: results/final_ledger.json and results/summary.json.
"""
import json, math
import numpy as np
import mfc_common as M
import chain_mfi as cm

MU = 7.159456
G_EFF = complex(-MU, 10.0)
MSIG = 40.0
AS = 2.1e-9
HSP = M.SP / 'hist'

CASES = [
    # key, label, engine tag, history json, grav tag, iso tag, iterate json
    ('HFa_exported', 'HFa as exported (unblocked axial reheating), design H_i', 'HFa_raw', M.BENCH_DIR / 'bench_HFa.json', 'HFa_raw', None, None),
    ('HFb_exported', 'HFb as exported (unblocked), design H_i', 'HFb_raw', M.BENCH_DIR / 'bench_HFb.json', 'HFa_raw', None, None),
    ('HFa_unblocked_sc', 'HFa unblocked, self-consistent H_i', 'itHFa_track_it1', HSP / 'bench_itHFa_track_it1.json', 'HFa_raw',
     'itHFa_track_it1', 'iterate_itHFa_track.json'),
    ('HFb_unblocked_sc', 'HFb unblocked, self-consistent H_i', 'itHFb_track_it0', HSP / 'bench_itHFb_track_it0.json', 'HFa_raw',
     'itHFb_track_it0', 'iterate_itHFb_track.json'),
    ('HFa_pauli_sc', 'HFa Pauli-blocked (Gamma_th = 0, tails free), self-consistent H_i', 'itHFa_P0_it2', HSP / 'bench_itHFa_P0_it2.json',
     'itHFa_P0_it2', 'itHFa_P0_it2', 'iterate_itHFa_P0.json'),
    ('HFb_pauli_sc', 'HFb Pauli-blocked (Gamma_th = 0, tails free), self-consistent H_i', 'itHFb_P0_it0', HSP / 'bench_itHFb_P0_it0.json',
     'itHFa_P0_it2', 'itHFb_P0_it0', 'iterate_itHFb_P0.json'),
    ('HFa_pauli_Gth10_sc', 'HFa Pauli-blocked, Gamma_th = 10 H_i, self-consistent H_i', 'itHFa_P10_it3', HSP / 'bench_itHFa_P10_it3.json',
     'itHFa_P0_it2', None, 'iterate_itHFa_P10.json'),
    ('HFa_pauli_band_sc', 'HFa Pauli-blocked, tails blocked while in band, self-consistent H_i', 'itHFa_PBb_it1',
     HSP / 'bench_itHFa_PBb_it1.json', 'itHFa_P0_it2', None, 'iterate_itHFa_PBb.json'),
]


def col(h, k):
    return np.array([np.nan if v is None else float(v) for v in h[k]])


def fs_over_a(d):
    h = d['history']
    t, fs, a, fr = col(h, 't'), col(h, 'f_s'), col(h, 'a'), col(h, 'f_rad')
    m = (t > 100) & (fr > 0.999) & np.isfinite(fs) & (fs > 0)
    if m.sum() < 3:
        m = (t > 100) & np.isfinite(fs) & (fs > 0)
    return float(np.median(fs[m] / a[m]))


def mfi_map(K):
    sq = json.load(open(M.TRK / 'mf_ideal' / 'results' / 'sq_map.json'))['groups']['primary']
    out = {}
    for ton in (0.1, 1.0):
        rows = sorted((r['K'], r['S_per_Q']['1e-06']) for r in sq if abs(r['m_s'] - 3.0) < 1e-9 and abs(r['t_on'] - ton) < 1e-9)
        Ks = np.array([x[0] for x in rows]); S = np.array([x[1] for x in rows])
        if K <= Ks[-1]:
            v = float(np.interp(K, Ks, S)); ext = False
        else:
            v = float(S[-1] + (S[-1] - S[-2]) / (Ks[-1] - Ks[-2]) * (K - Ks[-1])); ext = True
        out['t_on_%g' % ton] = dict(S_per_Q_1e6=v, extrapolated=ext)
    return out


def case(key, label, tag, hpath, gtag, itag, ipath):
    ch = json.load(open(M.ENG / 'chain' / ('chain_%s.json' % tag)))
    vec = json.load(open(M.ENG / 'vectors' / ('vectors_%s.json' % tag)))
    rel = json.load(open(M.ENG / 'release' / ('release_%s.json' % tag)))
    d = json.load(open(hpath))
    Hi = float(d['params'].get('H_i_GeV', d.get('H_i_GeV')))
    mP = M.MP_GEV / Hi
    K, EX, EEc = ch['K'], ch['E_X'], ch['E_over_Ec']
    ms = 3.0
    sstar = float(d['params']['s_star_MP']) * mP
    fsa = fs_over_a(d)
    ns_a3 = 3 * mP ** 2 * K ** 2 * fsa / ms
    Gs = 1.3e-6 * Hi / 1.3e11
    a_dec = math.sqrt(2 * K / Gs)
    f_dec = fsa * a_dec
    Dent = (1 + 0.8862 * f_dec) ** 0.75
    dep = {x['e']: x['deposit_per_vector'] for x in vec['by_e']}
    LMP = math.log(M.MP_GEV / (30 * Hi))
    CL = {'vanish_at_MY_no_threshold': 0.0, 'threshold_1/16pi^2': 1 / (16 * math.pi ** 2), 'threshold_1/4pi^2': 1 / (4 * math.pi ** 2),
          'vanish_at_MP': LMP / (4 * math.pi ** 2)}
    # gravitational production
    g = json.load(open(M.RES / 'grav' / ('grav_%s.json' % gtag)))['fields']
    gX = g['X_m2_20.2']['E_total_over_EX']              # E/E_X of the grav run's own release (HFb uses HFa's: same R(t))
    dNYL, dNYH, dNs = g['Y_L']['N_total'], g['Y_H']['N_total'], g['sigma_Sigma']['N_total']
    yyq = ch['Y_over_YQ_at_HM1']
    out = dict(key=key, label=label, engine_tag=tag, history=str(hpath), H_i_GeV=Hi, K=K, eps=ch['eps'], E_over_Ec=EEc,
               N_X=ch['N_X'], E_X=EX, N_YL=ch['N_YL'], N_YH=ch['N_YH'], N_sigma=ch['N_sigma'], N_A_1e6=ch['N_A_1e6'],
               S_per_Q_e_inf=ch['S_per_Q_e_inf'], e_crit_bulk=ch['e_crit_bulk'], e_crit_Q2e5=ch['e_crit_Q2e5'],
               reheating=d.get('reheating', {}), spectator=dict(f_s_over_a=fsa, n_s_a3=ns_a3, Gamma_s=Gs, a_decay=a_dec,
                                                                f_s_at_decay=f_dec, entropy_ratio_D=Dent, s_star_over_Hi=sstar),
               grav=dict(R_X=gX, dN_YL=dNYL, dN_YH=dNYH, dN_sigmaSigma=dNs), by_e={})
    for lab, e in (('e1e-6', 1e-6), ('e3e-6', 3e-6)):
        c = ch[lab]
        DA = dep[e]
        RY_per = c['R_Y'] / max(ch['N_YL'] + ch['N_YH'], 1e-300)
        Rc_per = c['R_casc'] / max(ch['N_YH'] + 2 * ch['N_sigma'], 1e-300)
        dR_grav = gX + dNYL * RY_per + dNYH * (RY_per + Rc_per) + 2 * dNs * Rc_per
        base = dict(R_vec=c['R_vec'], R_Y=c['R_Y'], R_casc=c['R_casc'], R_total=c['R_total'], S_per_Q=c['S_per_Q'], branch=c['branch'],
                    M=c['M'], p=c.get('p'), H_i_req_GeV=c['H_i_req_GeV'], r_at_req=c['r'], deposit_A_over_EX=DA / EX, dR_grav=dR_grav)
        leak = {}
        for nm, cl in CL.items():
            A = 2 * MU * G_EFF.real * cl
            Gl = (A / sstar) ** 2 * ms ** 3 / (32 * math.pi * (MSIG ** 2 - ms ** 2) ** 2)
            BR = Gl / Gs
            Rl = 2 * ns_a3 * BR * DA / EX
            Rt = c['R_total'] + dR_grav + Rl
            sq = cm.SQ(EEc * (1 + Rt))
            rr = cm.relic(sq)
            Hreq = cm.H_of(rr['M'], yyq / Dent)
            leak[nm] = dict(c_L=cl, A_Hi2=A, Gamma_leak=Gl, BR=BR, R_leak=Rl, R_total=Rt, S_per_Q=sq, branch=rr['branch'], M=rr['M'],
                            H_i_req_GeV=Hreq, H_req_over_H_i=Hreq / Hi)
        base['with_hidden_heat'] = leak
        out['by_e'][lab] = base
    # isocurvature
    if itag:
        iso = json.load(open(M.ENG / 'chain' / ('iso_%s.json' % itag)))
        amp2 = (Hi / (2 * math.pi * float(d['params']['s_star_MP']) * M.MP_GEV)) ** 2 / AS
        out['isocurvature'] = {lab: dict(D=v['D'], D_charge_only=v['D_charge_only'], beta_iso=amp2 * v['D'] ** 2 / (1 + amp2 * v['D'] ** 2),
                                         beta_iso_charge_only=amp2 * v['D_charge_only'] ** 2 / (1 + amp2 * v['D_charge_only'] ** 2),
                                         s_star_min_GeV=v['s_star_min_per_H_i'] * Hi, s_star_GeV=float(d['params']['s_star_MP']) * M.MP_GEV)
                               for lab, v in iso.items()}
    if ipath:
        it = json.load(open(M.RES / ipath))
        out['iteration'] = [dict(step=i, H_i_GeV=s['H_i_GeV'], K=s['K'], S_per_Q_1e6=s['SQ_1e6'], M_1e6=s['M_1e6'], H_req_GeV=s['H_req_GeV'],
                                 ratio=s['ratio']) for i, s in enumerate(it['steps'])]
        out['iteration_converged'] = it.get('converged')
    out['mf_ideal_map_same_K_ms3'] = mfi_map(K)
    return out


def main():
    rows = [case(*c) for c in CASES]
    M.dump(M.RES / 'final_ledger.json', dict(method=__doc__, cases=rows))
    for r in rows:
        a = r['by_e']['e1e-6']; b = r['by_e']['e3e-6']
        L = a['with_hidden_heat']
        print('%-20s H_i %.3e K %.3f | 1e-6: S/Q %.2f (%s, M %.2f, H_req/H_i %.3f) leak[1/16pi2, 1/4pi2, MP] R %.2g/%.2g/%.3g -> S/Q %.1f/%.1f/%.0f '
              '| 3e-6: S/Q %.2f %s | inf %.2f | D_s %.4f f_s,dec %.4f | grav dR %.2g | map(t_on .1/1) %.1f/%.1f' % (
                  r['key'], r['H_i_GeV'], r['K'], a['S_per_Q'], a['branch'], a['M'], a['H_i_req_GeV'] / r['H_i_GeV'],
                  L['threshold_1/16pi^2']['R_leak'], L['threshold_1/4pi^2']['R_leak'], L['vanish_at_MP']['R_leak'],
                  L['threshold_1/16pi^2']['S_per_Q'], L['threshold_1/4pi^2']['S_per_Q'], L['vanish_at_MP']['S_per_Q'],
                  b['S_per_Q'], b['branch'], r['S_per_Q_e_inf'], r['spectator']['entropy_ratio_D'], r['spectator']['f_s_at_decay'],
                  a['dR_grav'], r['mf_ideal_map_same_K_ms3']['t_on_0.1']['S_per_Q_1e6'], r['mf_ideal_map_same_K_ms3']['t_on_1']['S_per_Q_1e6']))
        if 'isocurvature' in r:
            for lab, v in r['isocurvature'].items():
                print('     iso %s: D %.3f (charge %.3f) beta %.2e (charge-only %.2e), s_* min %.2e GeV vs %.2e' % (
                    lab, v['D'], v['D_charge_only'], v['beta_iso'], v['beta_iso_charge_only'], v['s_star_min_GeV'], v['s_star_GeV']))
    return rows


if __name__ == '__main__':
    main()
