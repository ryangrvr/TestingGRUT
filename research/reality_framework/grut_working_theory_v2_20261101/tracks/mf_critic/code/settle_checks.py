#!/usr/bin/env python3
"""mf_critic (Claude): quick checks K1-K9 that settle the contradictions between the multi-field tracks.

K1  G1/G2 algebra of the proposal as written (v, lam, amplitude, blocking, CW ratio, V0^{1/4})            [D]
K2  the perturbative bound for a field holding V0 (mf_hybrid S2)                                          [D]
K3  the log-family tilt with the waterfall attached (mf_gemini NS vs its skeptic)                        [D]
K4  beta_iso = x/(1+x), x = (H/2 pi s_*)^2 D^2/A_s for every quoted D (mf_hybrid, mf_chain, its skeptic)   [D on Q]
K5  r, V0^{1/4} and the axial EFT ratio at every candidate H_i                                           [D on Q]
K6  gravitational-only M_t,crit at the strong-holder H_i, interpolated in ln H_i from mf_higgs H7          [E on Q]
K7  T_RH and related temperatures from every history file (read-only)                                     [C]
K8  e_crit(c0) on M1r's own history and the c0 that puts the threshold at e = 1e-6 (from m1r_c0.json)     [C]
K9  self-consistent strong-holder points on both ends of the Pauli bracket (skeptic grid + hyb_c0.json)   [C/E]
Usage: python3 -B settle_checks.py  -> results/settle_checks.json
Every choice (criteria, histories, rescaling rule, holder amplitude) is a model choice.
"""
import json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TRACK = HERE.parent
TRK = HERE.parents[1]
RES = TRACK / 'results'
SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2')
MP = 2.435e18
AS = 2.10e-9
G = 106.75
out = {}

# ---------------------------------------------------------------- K1
H0 = 1.3e11
v = math.sqrt(8 / 3)                                     # M_P units, from m = 3H and lam v^4/4 = 3 H^2 (M_P = 1)
lam = 27 / 16 * (H0 / MP) ** 2
m = 3 * H0 / MP
out['K1'] = dict(v_MP=v, lam=lam, check_m=math.sqrt(2 * lam) * v / (H0 / MP),
                 Phi_MP={str(k): math.sqrt(6) / k for k in (2, 3, 4)},
                 blocking_y2=2 * v / (m / 2), CW_over_tree_y2=2 ** 4 / (16 * math.pi ** 2 * lam),
                 V0_quarter_GeV=(3 * MP ** 2 * H0 ** 2) ** 0.25,
                 y_for_Gamma_half_H_Majorana=math.sqrt(0.5 * 16 * math.pi / 3.0))

# ---------------------------------------------------------------- K2
out['K2'] = dict(bound='m >= (384 pi)^(1/5) (M_P/H)^(2/5) H  [amplitude sqrt6 H M_P/m, y A < m/2, Gamma = y^2 m/16 pi >= H]',
                 m_min_over_H={'%.3g' % h: (384 * math.pi) ** 0.2 * (MP / h) ** 0.4 for h in (0.5e11, 0.84e11, 1.3e11, 2.4e11)})

# ---------------------------------------------------------------- K3
# pure log V = V0 (1 + A ln sigma): eta = -A/sigma^2, N = (sigma^2 - sigma_c^2)/2A, Delta = A/sigma_c^2
# -> n_s - 1 = -1/(N_sigma + 1/(2 Delta)) (leading order); track used -1/N_sigma; skeptic's regulated log in sk_ns.json
nsk = json.load(open(TRK / 'mf_gemini' / 'skeptic' / 'results' / 'sk_ns.json'))['H0=1.3e+11']['log_family']
k3 = {}
for D in (1.0, 0.3, 0.1, 0.03):
    r = nsk['K1_D%g' % D]
    Ns = r['N_sigma']
    k3['Delta=%g' % D] = dict(N_wf=r['N_wf'], N_sigma=Ns, pure_log_consistent=1 - 1 / (Ns + 1 / (2 * D)),
                             track_formula=1 - 1 / Ns, skeptic_regulated_log=r['n_s_sr2_with_xi'])
out['K3'] = dict(rows=k3, reading='both consistent attachments are bluer for slower crossings; the regulated-log shape of the '
                 'skeptic gives larger n_s at small Delta than the pure log (shape-dependent), the track formula is wrong in sign')

# ---------------------------------------------------------------- K4
def beta(H, sstar, D):
    x = (H / (2 * math.pi * sstar)) ** 2 * D * D / AS
    return x / (1 + x)


def smin(H, D, b=0.038):
    x = b / (1 - b)
    return abs(D) * H / (2 * math.pi * math.sqrt(x * AS))


k4 = {}
for lab, H, D in (('HFa unblocked 1e-6 (mf_chain map)', 1.28e11, 1.36), ('HFa unblocked 1e-6 (gap-filled map)', 1.11e11, 2.49),
                  ('HFa unblocked 3e-6 own fixed point (track map)', 1.752e11, 2.40), ('HFa unblocked 3e-6 (gap-filled)', 1.752e11, 2.68),
                  ('HFa Pauli 1e-6', 2.40e11, -0.07), ('strong holder Q-matter 1e-6 (c0 1.40)', 8.40e10, 9.245),
                  ('charge-only p = 4.5 (asymmetric relic)', 1.28e11, 4.54), ('mf_hybrid p = 5.1 (mf_ideal)', 1.30e11, 5.1)):
    k4[lab] = dict(H_i=H, D=D, beta_at_7p3e15=beta(H, 7.305e15, D), s_star_min_GeV=smin(H, D))
out['K4'] = dict(bound='beta_iso < 0.038 (uncorrelated CDI, n = 1; Planck 2015 XX arXiv:1502.02114, Planck 2018 X arXiv:1807.06211)',
                 rows=k4)

# ---------------------------------------------------------------- K5
def r_of(H):
    return 2 * H * H / (math.pi ** 2 * AS * MP * MP)


k5 = {}
for lab, H in (('design', 1.30e11), ('unblocked sc (track map)', 1.28e11), ('unblocked sc (gap-filled)', 1.11e11),
               ('Pauli P0 sc', 2.40e11), ('strong holder P0', 8.40e10), ('strong holder P0 c0 1.36', 8.85e10),
               ('strong holder Gamma_th = 10', 5.2e10), ('M1r c0 = 1', 1.292e11), ('M1r c0 = 1.1', 9.56e10)):
    k5[lab] = dict(H_i=H, r=r_of(H), V0_quarter_GeV=(3 * MP ** 2 * H ** 2) ** 0.25,
                   axial_EFT_ratio=[0.42 * 1.30e11 / H, 0.49 * 1.30e11 / H])
out['K5'] = dict(note='axial EFT ratio mu5,max/4 pi f under mf_chain\'s rescaling rule (f fixed in H_i units) [Q mf_chain skeptic]',
                 rows=k5)

# ---------------------------------------------------------------- K6
tab = [(0.6e11, 173.12, 173.16), (1.0e11, 172.96, 172.99), (1.3e11, 172.88, 172.92), (2.1e11, 172.74, 172.77),
       (3.0e11, 172.64, 172.68), (5.0e11, 172.54, 172.56)]          # mf_higgs H7 [Q], gravitational only, xi_EW = 0.1
lh = np.log([x[0] for x in tab])
k6 = {}
for H in (5.2e10, 8.4e10, 8.85e10, 1.11e11, 1.28e11, 2.40e11):
    lo = float(np.interp(math.log(H), lh, [x[1] for x in tab], left=None))
    hi = float(np.interp(math.log(H), lh, [x[2] for x in tab]))
    if H < tab[0][0]:                                    # linear extrapolation in ln H below the table
        s1 = (tab[1][1] - tab[0][1]) / (lh[1] - lh[0]); s2 = (tab[1][2] - tab[0][2]) / (lh[1] - lh[0])
        lo = tab[0][1] + s1 * (math.log(H) - lh[0]); hi = tab[0][2] + s2 * (math.log(H) - lh[0])
    k6['%.3g' % H] = [lo, hi]
out['K6'] = dict(note='M_t,crit (GeV), gravitational reheating only, interpolated in ln H_i across histories (K-dependence ignored) [E]; '
                 'with the N-gas and SM hard-gas medium masses every history is safe to >= 174.6 GeV [Q mf_higgs + skeptic, E]',
                 Mt_crit=k6)

# ---------------------------------------------------------------- K7
def temps(path):
    d = json.load(open(path))
    h = d['history']
    Hi = float(d['params'].get('H_i_GeV', d.get('H_i_GeV')))
    t = np.array(h['t'], float); H = np.array(h['H'], float)
    fr = np.array([np.nan if x is None else x for x in h['f_rad']], float)
    fsm = np.array([np.nan if x is None else x for x in h['f_SM']], float)
    a = np.array(h['a'], float) if 'a' in h else np.exp(np.array(h['N'], float))
    K = float(d['reheating']['K'])
    T = lambda i, f: (30 * f * 3 * MP ** 2 * (H[i] * Hi) ** 2 / (math.pi ** 2 * G)) ** 0.25
    post = t >= 0
    res = dict(H_i_GeV=Hi, K=K, T_K_GeV=(90 * K * K * MP ** 2 * Hi ** 2 / (math.pi ** 2 * G)) ** 0.25)
    for lab, arr, thr in (('rad50', fr, 0.5), ('rad99', fr, 0.99), ('SM50', fsm, 0.5)):
        idx = np.where(post & (arr >= thr))[0]
        if len(idx):
            # last entry into the region (radiation fraction can dip after a first rise)
            ok = post & (arr >= thr)
            j = idx[0]
            k = np.where(post & ~ok & (np.arange(len(t)) > j))[0]
            if len(k):
                later = idx[idx > k[-1]]
                j = later[0] if len(later) else j
            res[lab] = dict(t=float(t[j]), a=float(a[j]), T_GeV=float(T(j, arr[j])))
        else:
            res[lab] = None
    ipk = np.nanargmax(np.where(post, fr * H ** 2, np.nan))
    res['T_max_rad_GeV'] = float(T(ipk, fr[ipk]))
    return res


hists = {'HFa exported (unblocked, c0-independent background)': TRK / 'mf_hybrid' / 'results' / 'bench_HFa.json',
         'HFb exported': TRK / 'mf_hybrid' / 'results' / 'bench_HFb.json',
         'HFa unblocked sc 1.284e11 (mf_chain)': TRK / 'mf_chain' / 'results' / 'histories' / 'bench_itHFa_track_it1.json',
         'HFa Pauli P0 sc 2.40e11 (mf_chain)': TRK / 'mf_chain' / 'results' / 'histories' / 'bench_itHFa_P0_it2.json',
         'HFb Pauli P0 2.40e11 (mf_chain)': TRK / 'mf_chain' / 'results' / 'histories' / 'bench_itHFb_P0_it0.json',
         'HFa Pauli P0 8.40e10 (rebuilt here)': SCR / 'mf_critic' / 'hist' / 'bench_crP0_H084.json',
         'HFa Pauli Gamma_th=10 8.40e10 (here)': SCR / 'mf_critic' / 'hist' / 'bench_crP10_H084.json',
         'HFa Pauli Gamma_th=10 5.0e10 (here)': SCR / 'mf_critic' / 'hist' / 'bench_crP10_H050.json'}
out['K7'] = {k: temps(p) for k, p in hists.items() if p.exists()}
Km = 1.98581; Hm = 1.292e11
out['K7']['M1r (3, 0.5) supplied'] = dict(H_i_GeV=Hm, K=Km, T_K_GeV=(90 * Km * Km * MP ** 2 * Hm ** 2 / (math.pi ** 2 * G)) ** 0.25,
                                          note='oscillator decaying at Gamma_I = 0.5 H_i; radiation dominates at t ~ 2-3 H_i^-1')

# ---------------------------------------------------------------- K8
m1 = json.load(open(RES / 'm1r_c0.json'))
rows = []
for k in sorted(m1, key=float):
    s = m1[k]
    rows.append(dict(c0=s['c0'], eps=s['eps'], E_over_Ec=s['E_over_Ec'], N_YL=s['N_YL'], N_YH=s['N_YH'],
                     SQ_1e6=s['e1e-6']['S_per_Q'], SQ_3e6=s['e3e-6']['S_per_Q'], SQ_inf=s['S_per_Q_e_inf'],
                     branch_1e6=s['e1e-6']['branch'], branch_3e6=s['e3e-6']['branch'], e_crit_bulk=s['e_crit_bulk'],
                     e_crit_Q2e5=s['e_crit_Q2e5'], M_1e6=s['e1e-6']['M'], H_req_1e6=s['e1e-6']['H_i_req_GeV'],
                     r=s['e1e-6']['H_i_req_GeV'] and r_of(s['e1e-6']['H_i_req_GeV'])))
c = np.array([r_['c0'] for r_ in rows]); sq = np.array([r_['SQ_1e6'] for r_ in rows])
ec = [(r_['c0'], r_['e_crit_bulk']) for r_ in rows if r_['e_crit_bulk']]
i = int(np.where((sq[:-1] > 3.4716) & (sq[1:] <= 3.4716))[0][0])
c_thr = float(np.exp(np.interp(math.log(3.4716), [math.log(sq[i + 1]), math.log(sq[i])], [math.log(c[i + 1]), math.log(c[i])])))
ce = np.array([x[0] for x in ec]); ee = np.array([x[1] for x in ec])
j = int(np.argmin(abs(ce - 1.0)))
slope = float(np.polyfit(np.log(ce[(ce >= 0.9) & (ce <= 1.1)]), np.log(ee[(ce >= 0.9) & (ce <= 1.1)]), 1)[0])
out['K8'] = dict(rows=rows, c0_threshold_at_1e6=c_thr, dln_ecrit_dln_c0_near_1=slope,
                 reading='on M1r\'s own supplied history the holder amplitude moves e_crit as ~c0^%.1f; a holder %.0f%% stronger than '
                         'the inherited 10i calibration puts the bulk threshold at e = 1e-6' % (slope, 100 * (c_thr - 1)))

# ---------------------------------------------------------------- K9
sk = json.load(open(TRK / 'mf_chain' / 'skeptic' / 'results' / 'sk_summary.json'))['full_chain_runs']
hy = json.load(open(RES / 'hyb_c0.json'))
pts = []
for tag, H in (('skP0H084', 8.40e10), ('skP0H110', 1.10e11), ('skP0', 2.40e11)):
    for k, s in sk.items():
        if k == tag or k.startswith(tag + '_C'):
            pts.append(dict(src='mf_chain skeptic', bracket='Gamma_th = 0', H_i=H, K=s['K'], c0=s['c0'],
                            SQ_1e6=s['e1e-6']['S_per_Q'], SQ_3e6=s['e3e-6']['S_per_Q'], e_crit=s['e_crit_bulk'],
                            H_req=s['e1e-6']['H_i_req_GeV'], ratio=s['e1e-6']['H_i_req_GeV'] / H, branch=s['e1e-6']['branch'],
                            M=s['e1e-6']['M']))
for tag, blk in hy.items():
    H = blk['meta']['H_i_GeV']
    for k, s in blk['runs'].items():
        pts.append(dict(src='mf_critic', bracket='Gamma_th = %g' % blk['meta']['reheating']['Gamma_th'], H_i=H, K=s['K'], c0=s['c0'],
                        SQ_1e6=s['e1e-6']['S_per_Q'], SQ_3e6=s['e3e-6']['S_per_Q'], e_crit=s['e_crit_bulk'],
                        H_req=s['e1e-6']['H_i_req_GeV'], ratio=s['e1e-6']['H_i_req_GeV'] / H, branch=s['e1e-6']['branch'],
                        M=s['e1e-6']['M']))
pts.sort(key=lambda p: (p['bracket'], p['H_i'], p['c0']))
# regression: my rebuilt P0 8.4e10 history at c0 = 1.40 against the skeptic's
reg = [p for p in pts if p['src'] == 'mf_critic' and p['bracket'] == 'Gamma_th = 0' and abs(p['c0'] - 1.4) < 1e-9]
skr = sk['skP0H084_C1.400']['e1e-6']
out['K9'] = dict(points=pts,
                 regression_c0_1p40_P0=dict(mine=reg[0]['SQ_1e6'] if reg else None, skeptic=skr['S_per_Q'],
                                            H_req_mine=reg[0]['H_req'] if reg else None, H_req_skeptic=skr['H_i_req_GeV']),
                 reading=('Gamma_th = 0 (K = 6.5 at 8.4e10): self-consistent within 5-9% for c0 = 1.36-1.40 at H_i = 8.4-9.2e10 '
                          '(c0 1.36: free X 4.16, e_crit 2.2e-6; 1.40: Q-matter 2.88).  Gamma_th = 10 (K = 3.2 at 8.4e10, 2.5 at 5e10): '
                          'at 8.4e10 the self-consistent c0 ~ 1.24 gives free X with no e_crit; the M1r-like point moves to c0 ~ 1.25 at '
                          'H_i ~ 5.2e10 (S/Q 3.64, e_crit 1.3e-6) and Q-matter at 1e-6 needs c0 ~ 1.27-1.30 at H_i ~ 4.3-5e10, where the '
                          'axial EFT ratio is 1.1-1.5 under the rescaling rule.'))

json.dump(out, open(RES / 'settle_checks.json', 'w'), indent=1, default=float)
print(json.dumps({k: (v if k in ('K1', 'K2', 'K5', 'K6') else '...') for k, v in out.items()}, indent=1, default=float)[:4000])
print('K3', json.dumps(out['K3']['rows'], default=float))
print('K4', json.dumps({k: (round(v['beta_at_7p3e15'], 4), '%.3g' % v['s_star_min_GeV']) for k, v in out['K4']['rows'].items()}))
print('K7', json.dumps({k: (v.get('K'), '%.3g' % v['T_K_GeV'], v.get('rad50') and '%.3g@t%.3g' % (v['rad50']['T_GeV'], v['rad50']['t']),
                            v.get('rad99') and '%.3g@t%.3g' % (v['rad99']['T_GeV'], v['rad99']['t']),
                            v.get('SM50') and '%.3g@t%.3g' % (v['SM50']['T_GeV'], v['SM50']['t'])) for k, v in out['K7'].items()}, indent=0))
print('K8 c0_thr %.4f slope %.2f' % (out['K8']['c0_threshold_at_1e6'], out['K8']['dln_ecrit_dln_c0_near_1']))
print('K9 regression', out['K9']['regression_c0_1p40_P0'])
for p in pts:
    print('  %-16s %-13s H %.3g K %.2f c0 %.3f S/Q %.3f / %.3f e_crit %s ratio %.3f %s' % (p['src'], p['bracket'], p['H_i'], p['K'], p['c0'],
          p['SQ_1e6'], p['SQ_3e6'], p['e_crit'] and '%.2e' % p['e_crit'], p['ratio'], p['branch']))
