#!/usr/bin/env python3
"""mf_chain SKEPTIC: independent ledger, S/Q, relic, H_i calibration and isocurvature from the stored release and vector
outputs (mf_chain's and the skeptic's), with
  (a) my own yield/entropy/abundance arithmetic (no chain_mfi function is used for it):
        rho_r a^4 = (pi^2/30) g T^4 a^4 = 3 K^2 M_P^2 H_i^2,  s a^3 = (2 pi^2/45) g (T a)^3,  Y = |N_X| H_i^3/(s a^3),
        rho_DM/s = M m_X Y = 4.4e-10 GeV  (BRIEF_DECISIVE target; Planck 2018 Omega_c h^2 = 0.120 gives 4.37e-10) -> H_req;
  (b) the pair fraction p(S/Q) from the case-L network computed INSIDE the inherited map's gap (sk_pmap.json) instead of
      the log-linear bridge between two different recipes;
  (c) D = d ln Omega/d ln s_* recomputed from the c0 = 0.99 / 1.01 releases with (b), and beta_iso at the H_i of the
      self-consistent point that actually belongs to each e (the track evaluated the e = 3e-6 derivative at the 1e-6 H_i).
[C] arithmetic on [C] inputs; ledger (dq_heat), closed form (dq_entropy) and the case-L recipe are inherited [Q].
Usage: python3 -B sk_ledger.py -> ../results/sk_ledger.json
"""
import json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TRK = HERE.parents[2]
MFC = TRK / 'mf_chain' / 'results' / 'engine'
SKE = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_chain_skeptic/engine')
OUT = HERE.parent / 'results' / 'sk_ledger.json'
MP = 2.435e18
MX = 0.0924089637
DTARGET = 4.4e-10
G = 106.75
AS = 2.100e-9
BETA = 0.038
SQB = 3.4716191256166837
M_L, G_L = 22.24859546128699, 1.005903273746766e-07


def SQ(W):
    return 0.190255 * (W - 0.97494) ** 0.75013 - 0.00428


# ------------------------------------------------------------------ pair fraction maps
def load_maps():
    cc = json.load(open(TRK / 'wt_critic' / 'results' / 'chain_consistency.json'))['p_scan_SQ']
    wt = sorted((r['S_per_Q'], r['p']) for r in cc)
    pk = json.load(open(TRK / 'inf_chain' / 'results' / 'pmap_kin.json'))
    ic = sorted((r['S_per_Q_closed_form'], r['p']) for r in pk['rows'].values() if 'p' in r)
    sc = sorted((float(k), v['p_self_consistent'] / v['p_fixed_yield']) for k, v in pk['self_consistent'].items())
    sk = json.load(open(HERE.parent / 'results' / 'sk_pmap.json'))
    mine = sorted((r['S_per_Q'], r['p']) for r in sk['rows'].values())
    return wt, ic, sc, mine


WT, IC, SCR, MINE = load_maps()


def _loglog(rows, s):
    x = np.log([r[0] for r in rows]); y = np.log([r[1] for r in rows]); lx = math.log(s)
    if lx <= x[0]:
        return math.exp(y[0] + (y[1] - y[0]) / (x[1] - x[0]) * (lx - x[0]))
    if lx >= x[-1]:
        return math.exp(y[-1] + (y[-1] - y[-2]) / (x[-1] - x[-2]) * (lx - x[-1]))
    return math.exp(float(np.interp(lx, x, y)))


def sc_factor(s):
    xs = [math.log(6.742)] + [math.log(a) for a, _ in SCR]; ys = [1.0] + [b for _, b in SCR]
    return 1.0 if math.log(s) <= xs[0] else float(np.interp(math.log(s), xs, ys))


def p_track(s):
    """the inherited map as chain_mfi uses it: wt_critic (<= 6.74) + inf_chain case-L (>= 17.1), log-linear bridge."""
    return _loglog(sorted(WT + IC), s) * sc_factor(s)


def p_skeptic(s):
    """case-L network everywhere it was computed (mine 6.74-21.8, inf_chain above); wt_critic below 6.74."""
    rows = sorted([r for r in WT if r[0] < 6.7] + MINE + [r for r in IC if r[0] > 22.0])
    return _loglog(rows, s) * sc_factor(s)


def M_relic(s, pf):
    """mass per net charge (m_X): free X 1 + 2p; Q-matter: qm_micro / dq_critic table [Q] (mid of its lo/hi band)."""
    if s > SQB:
        return 1 + 2 * pf(s)
    S = [1.5, 2.5, 3.0, 3.2]; lo = [0.7883, 0.7974, 0.8023, 0.8048]; hi = [0.7944, 0.8076, 0.8152, 0.8182]
    return 0.5 * (float(np.interp(s, S, lo)) + float(np.interp(s, S, hi)))


# ------------------------------------------------------------------ independent ledger and calibration
def ledger(rel, vec):
    r = rel['result']
    K = vec['background']['K']
    NX, EX, eps = abs(r['N_X']), r['E_X'], r['eps']
    E_Ec = eps ** (-4 / 3)
    aL = math.sqrt(math.pi * K / (2 * G_L))
    RY = (r['N_YL'] + r['N_YH']) * M_L * aL / EX
    out = dict(K=K, N_X=NX, E_X=EX, eps=eps, E_Ec=E_Ec, R_Y=RY, by_e={})
    es = [x['e'] for x in vec['by_e']]
    le = np.log(es); lR = np.log([x['R_vec'] for x in vec['by_e']]); lD = np.log([x['deposit_per_vector'] for x in vec['by_e']])
    for e in (1e-6, 3e-6):
        Rv = math.exp(float(np.interp(math.log(e), le, lR))); DA = math.exp(float(np.interp(math.log(e), le, lD)))
        Rc = (r['N_YH'] + 2 * r['N_sigma']) * DA / EX
        W = E_Ec * (1 + Rv + RY + Rc)
        out['by_e']['%g' % e] = dict(R_vec=Rv, R_casc=Rc, W=W, S_per_Q=SQ(W))
    out['S_per_Q_inf'] = SQ(E_Ec * (1 + RY))
    return out


def H_req(NX, K, M):
    """my own abundance arithmetic (GeV units): Y(H) = |N_X| H^3/(s a^3) ~ H^{3/2}; M m_X Y = D_target."""
    c = (2 * math.pi ** 2 / 45) * G * (90 * K * K * MP * MP / (math.pi ** 2 * G)) ** 0.75   # s a^3 / H^{3/2}
    return (DTARGET * c / (M * MX * NX)) ** (2.0 / 3.0)


def beta_iso(H, sstar, D):
    x = (H / (2 * math.pi * sstar)) ** 2 * D * D / AS
    return x / (1 + x)


def load(tag, eng, need_chain=True):
    rel = json.load(open(eng / 'release' / ('release_%s.json' % tag)))
    vec = json.load(open(eng / 'vectors' / ('vectors_%s.json' % tag)))
    ch = json.load(open(eng / 'chain' / ('chain_%s.json' % tag))) if need_chain else None
    return rel, vec, ch


CASES = [  # (label, base tag, engine dir, has iso variants, H_i of that history [GeV])
    ('HFa exported (design 1.30e11)', 'HFa_raw', MFC, False),
    ('HFa unblocked sc 1e-6', 'itHFa_track_it1', MFC, True),
    ('HFa unblocked sc 3e-6', 'e3U', SKE, True),
    ('HFa Pauli P0 sc 1e-6', 'itHFa_P0_it2', MFC, True),
    ('HFb Pauli P0 sc 1e-6', 'itHFb_P0_it0', MFC, True),
    ('HFa Pauli P0 sc 3e-6', 'e3P', SKE, True),
    ('HFa unblocked at 1.12e11 (skeptic-map iteration step)', 'skUH112', SKE, False),
    ('HFa Pauli P0 at 8.40e10, c0 = 1.40 (skeptic)', 'skP0H084_C1.400', SKE, False),
    ('HFa Pauli P0 at 8.40e10, c0 = 1.36 (skeptic)', 'skP0H084_C1.360', SKE, False),
]


def main():
    res = dict(method=__doc__, p_map=dict(track={str(s): p_track(s) for s in (6.742, 8, 10, 12, 14, 15, 16, 17.13, 21.8)},
                                          skeptic={str(s): p_skeptic(s) for s in (6.742, 8, 10, 12, 14, 15, 16, 17.13, 21.8)},
                                          slope_track_14_16=math.log(p_track(16) / p_track(14)) / math.log(16 / 14),
                                          slope_skeptic_14_16=math.log(p_skeptic(16) / p_skeptic(14)) / math.log(16 / 14),
                                          wt_vs_caseL_at_6p742=[p_track(6.742), MINE[0][1]]), cases={})
    for lab, tag, eng, iso in CASES:
        try:
            rel, vec, ch = load(tag, eng)
        except FileNotFoundError as e:
            res['cases'][lab] = dict(missing=str(e)); continue
        L = ledger(rel, vec)
        row = dict(tag=tag, K=L['K'], N_X=L['N_X'], E_Ec=L['E_Ec'], R_Y=L['R_Y'], S_per_Q_inf_mine=L['S_per_Q_inf'],
                   S_per_Q_inf_track=ch['S_per_Q_e_inf'], by_e={})
        for e in ('1e-06', '3e-06'):
            lab_e = 'e1e-6' if e == '1e-06' else 'e3e-6'
            b = L['by_e'][e]
            s = b['S_per_Q']
            pt, ps = p_track(s), p_skeptic(s)
            Mt, Ms = M_relic(s, p_track), M_relic(s, p_skeptic)
            row['by_e'][e] = dict(S_per_Q_mine=s, S_per_Q_track=ch[lab_e]['S_per_Q'], branch='free X' if s > SQB else 'Q-matter',
                                  p_track_map=pt, p_skeptic_map=ps, M_track=ch[lab_e]['M'], M_track_recomputed=Mt, M_skeptic=Ms,
                                  asym_share_skeptic=1 / Ms,
                                  H_req_track_GeV=ch[lab_e]['H_i_req_GeV'], H_req_mine_trackmap_GeV=H_req(L['N_X'], L['K'], Mt),
                                  H_req_mine_skepticmap_GeV=H_req(L['N_X'], L['K'], Ms))
        if iso:
            v = {}
            for c0 in ('0.99', '1.01'):
                try:
                    rr, vv, cc = load('%s_c%s' % (tag, c0), eng, need_chain=False)
                except FileNotFoundError:
                    v = None; break
                LL = ledger(rr, vv)
                v[c0] = (LL, cc)
            if v:
                isd = {}
                for e in ('1e-06', '3e-06'):
                    lab_e = 'e1e-6' if e == '1e-06' else 'e3e-6'
                    lo = {}
                    for mapname, pf in (('track', p_track), ('skeptic', p_skeptic)):
                        lnO = {c0: math.log(M_relic(LL['by_e'][e]['S_per_Q'], pf) * LL['N_X'] * LL['K'] ** -1.5) for c0, (LL, cc) in v.items()}
                        lo[mapname] = (lnO['1.01'] - lnO['0.99']) / (math.log(1.01) - math.log(0.99))
                    DN = (math.log(v['1.01'][0]['N_X']) - math.log(v['0.99'][0]['N_X'])) / (math.log(1.01) - math.log(0.99))
                    isd[e] = dict(D_trackmap=lo['track'], D_skepticmap=lo['skeptic'], D_charge_only=DN)
                row['iso'] = isd
        res['cases'][lab] = row
    json.dump(res, open(OUT, 'w'), indent=1, default=float)
    return res


if __name__ == '__main__':
    r = main()
    print(json.dumps(r['p_map'], indent=1))
    for lab, row in r['cases'].items():
        if 'missing' in row:
            print(lab, 'MISSING', row['missing']); continue
        print('== %s (%s): K %.3f E/Ec %.2f R_Y %.3f S/Q(inf) %.3f/%.3f' % (lab, row['tag'], row['K'], row['E_Ec'], row['R_Y'],
              row['S_per_Q_inf_mine'], row['S_per_Q_inf_track']))
        for e, b in row['by_e'].items():
            print('   e=%s S/Q %.4f (track %.4f) %s | M track %.4f (recomp %.4f) skeptic %.4f | H_req track %.4e mine(trackmap) %.4e mine(skmap) %.4e' % (
                e, b['S_per_Q_mine'], b['S_per_Q_track'], b['branch'], b['M_track'], b['M_track_recomputed'], b['M_skeptic'],
                b['H_req_track_GeV'], b['H_req_mine_trackmap_GeV'], b['H_req_mine_skepticmap_GeV']))
        if 'iso' in row:
            for e, x in row['iso'].items():
                print('   iso e=%s: D trackmap %.3f skepticmap %.3f charge-only %.3f' % (e, x['D_trackmap'], x['D_skepticmap'], x['D_charge_only']))
