"""Verifier v_vacua, part 3: the Standard Model relation of track `vacua` (NOTES.md section 5), recomputed.

Independent of the track's code: RG running and beta functions are those of SMDR v1.3 (Martin-Robertson,
arXiv:1907.02500; full 3-loop SM betas from Bednyakov-Pikelner-Velizhanin / Chetyrkin-Zoller, plus higher-order
pieces at loopOrder 4, 5), called through our own driver smdr_driver.c.  Two matchings are used:
  (B) the NNLO fit formulas of Buttazzo et al. 1307.3536 sec. 3 (transcribed here from the TeX, as the track did);
  (S) SMDR's own full fit from pole masses (2-loop + leading 3-loop EW/QCD, 4-loop QCD in Mt <-> yt), cached in
      ../results/v3a_smdr_fit_grid.json by v3a_smdr_fit_grid.py.
The track's own beta functions (tracks/vacua/code/sm_rge.py) are imported ONLY in part A, to audit their
transcription against SMDR; no number of parts B-F uses the track's code.

Run (after sh v3_build_smdr.sh and python3 v3a_smdr_fit_grid.py):  python3 v3b_sm_check.py
-> ../results/v3b_sm_check.json
"""
import importlib.util
import json
import os
import sys
import time

import numpy as np
from scipy.optimize import brentq, fsolve
from scipy.stats import norm

import smdr_common as S

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results', 'v3b_sm_check.json')
GRID = os.path.join(HERE, '..', 'results', 'v3a_smdr_fit_grid.json')
TRACK_RGE = os.path.join(HERE, '..', '..', 'vacua', 'code', 'sm_rge.py')
MPL, MPL_RED = S.MPL, S.MPL_RED
K16 = 16 * np.pi ** 2
out = {}
t_start = time.time()

# ================================================================ A. audit of the track's beta functions vs SMDR
sys.dont_write_bytecode = True        # never write into the audited track's directory
spec = importlib.util.spec_from_file_location('track_sm_rge', TRACK_RGE)
TR = importlib.util.module_from_spec(spec)
spec.loader.exec_module(TR)


def smdr_in_track_units(L, g3, g, gp, yt, k):
    """SMDR betas (Q d/dQ) converted to the track's d/d ln mu^2 of (G1, G2, G3, T, lam), and gamma = d ln v/d ln Q."""
    b = S.betas(L, g3, g, gp, yt, 0.0, 0.0, k)
    return np.array([5 / 3 * gp * b['beta_gp'], g * b['beta_g'], g3 * b['beta_g3'], yt * b['beta_yt'],
                     b['beta_k'] / 2]), b['beta_v'] / 246.0


rng = np.random.default_rng(1)
pts = [(1.1666, 0.64779, 0.35830, 0.93690, 0.12604), (0.4873, 0.5055, 0.4767, 0.3825, -0.0143),
       (0.70, 0.58, 0.39, 0.60, 0.012)]
pts += [tuple(rng.uniform([0.4, 0.45, 0.3, 0.3, -0.03], [1.2, 0.7, 0.5, 1.0, 0.15])) for _ in range(3)]
audit = []
names = ['G1', 'G2', 'G3', 'T', 'lam']
for (g3, g, gp, yt, k) in pts:
    G1 = 5 / 3 * gp ** 2
    y = np.array([G1, g ** 2, g3 ** 2, yt ** 2, 0.0, 0.0, k])
    row = {'point(g3,g,gp,yt,lam)': [g3, g, gp, yt, k]}
    for L in (1, 2, 3):
        tb = TR.betas(y, loops=L)[[0, 1, 2, 3, 6]]
        if L == 3:   # remove the 4-loop pure-QCD term the track adds at loops=3, compare it separately
            tb = tb.copy()
            tb[2] -= g3 ** 10 / K16 ** 4 * (-2472.28)
        sb, gam_s = smdr_in_track_units(L, g3, g, gp, yt, k)
        # loop-L pieces only (differences of cumulative orders), to make the comparison sharp
        if L > 1:
            tb_prev = TR.betas(y, loops=L - 1)[[0, 1, 2, 3, 6]]
            sb_prev, gam_prev = smdr_in_track_units(L - 1, g3, g, gp, yt, k)
            dt, ds = tb - tb_prev, sb - sb_prev
            gt = TR.gamma_h(y, loops=L) - TR.gamma_h(y, loops=L - 1)
            gs = gam_s - gam_prev
        else:
            dt, ds = tb, sb
            gt, gs = TR.gamma_h(y, loops=1), gam_s
        row[f'{L}-loop piece rel. diff (track-SMDR)/SMDR'] = {n: float((a - b_) / b_) if b_ != 0 else None
                                                               for n, a, b_ in zip(names, dt, ds)}
        row[f'{L}-loop gamma: track, SMDR(d ln v/d ln Q)'] = [float(gt), float(gs)]
    # 4-loop pure QCD g3 term: track's coefficient vs SMDR loop-4 minus loop-3 in the pure-QCD limit
    s3 = S.betas(3, g3, 1e-9, 1e-9, 1e-9, 0, 0, 0)['beta_g3']
    s4 = S.betas(4, g3, 1e-9, 1e-9, 1e-9, 0, 0, 0)['beta_g3']
    row['4-loop QCD dG3/dlnmu2: track, SMDR'] = [g3 ** 10 / K16 ** 4 * (-2472.28), float(g3 * (s4 - s3))]
    audit.append(row)
out['A_beta_audit'] = audit


# ================================================================ B. Buttazzo matching (own transcription) + SMDR running
def matching_B(Mt, Mh, als, MW=80.384, dlam=0.0, dyt=0.0):
    """1307.3536 sec. 3 (TeX artv2.tex lines 1183-1225): MSbar couplings at mu = Mt."""
    a = (als - 0.1184) / 0.0007
    w = (MW - 80.384) / 0.014
    return {'k': 0.12604 + 0.00206 * (Mh - 125.15) - 0.00004 * (Mt - 173.34) + dlam,
            'yt': 0.93690 + 0.00556 * (Mt - 173.34) - 0.00042 * a + dyt,
            'g3': 1.1666 + 0.00314 * a - 0.00046 * (Mt - 173.34),
            'g': 0.64779 + 0.00004 * (Mt - 173.34) + 0.00011 * w,
            'gp': 0.35830 + 0.00011 * (Mt - 173.34) - 0.00020 * w,
            'yb': 0.0156, 'ytau': 0.0100}


def lam_min_B(Mt, Mh, als, L=3, Qcut=MPL, MW=80.3692, **kw):
    return S.min_lambda(L, Mt, matching_B(Mt, Mh, als, MW, **kw), Qcut=Qcut)


def Mtcrit_B(Mh, als, L=3, Qcut=MPL, MW=80.3692, **kw):
    f = lambda Mt: lam_min_B(Mt, Mh, als, L, Qcut, MW, **kw)[0]
    Mt = brentq(f, 166.0, 178.0, xtol=1e-6)
    return Mt, lam_min_B(Mt, Mh, als, L, Qcut, MW, **kw)[1]


B = {}
# validation at Buttazzo's reference point (their MW = 80.384)
c_ref = matching_B(173.34, 125.15, 0.1184)
r = S.run(3, 173.34, c_ref, N=3, Qmin=1e8, Qmax=MPL)
B['ref_point_at_MPl_3loop'] = {'g1': float(np.sqrt(5 / 3) * r['gp'][-1]), 'g2': float(r['g'][-1]),
                               'g3': float(r['g3'][-1]), 'yt': float(r['yt'][-1]), 'lam': float(r['k'][-1]),
                               'track': {'g1': 0.61537, 'g2': 0.50546, 'g3': 0.48731, 'yt': 0.38248, 'lam': -0.01437},
                               'Buttazzo': {'g1': 0.6154, 'g2': 0.5055, 'g3': 0.4873, 'yt': 0.3825, 'lam': -0.0143}}
lamP = lambda Mt, Mh, als, L=3: float(S.run(L, Mt, matching_B(Mt, Mh, als), N=2, Qmin=1e8, Qmax=MPL)['k'][-1])
l0 = lamP(173.34, 125.15, 0.1184)
B['dlam(MPl)'] = {'per GeV Mt': lamP(174.34, 125.15, 0.1184) - l0, 'per 0.0007 as': lamP(173.34, 125.15, 0.1191) - l0,
                  'per GeV Mh': lamP(173.34, 126.15, 0.1184) - l0,
                  'track': [-0.00652, 0.00185, 0.00286], 'Buttazzo': [-0.0066, 0.0018, 0.0029]}
B['lam(MPl) by loop order 1/2/3/4/5'] = [lamP(173.34, 125.15, 0.1184, L) for L in (1, 2, 3, 4, 5)]
B['track lam(MPl) 1/2/3'] = [-0.03371, -0.01508, -0.01437]
rr = S.run(3, 173.34, c_ref, N=2001, Qmin=1e8, Qmax=1e12)
B['lambda=0 crossing at ref (GeV)'] = float(np.exp(np.interp(0.0, -rr['k'], np.log(rr['Q']))))
B['Mtcrit_ref(125.15,0.1184) msbar'] = Mtcrit_B(125.15, 0.1184, MW=80.384)
B['Mtcrit_Bednyakov_inputs(125.7,0.1185) msbar'] = Mtcrit_B(125.7, 0.1185, MW=80.385)
B['Buttazzo linear formula at Mh=125.15'] = 173.34 - (129.6 - 125.15) / 2.0

# PDG 2025 point
Mh0, sMh, a0, sa, MW0 = 125.20, 0.11, 0.1180, 0.0009, 80.3692
Mtc, mus = Mtcrit_B(Mh0, a0)
dMh = (Mtcrit_B(Mh0 + 0.5, a0)[0] - Mtcrit_B(Mh0 - 0.5, a0)[0]) / 1.0
das = (Mtcrit_B(Mh0, a0 + 0.0009)[0] - Mtcrit_B(Mh0, a0 - 0.0009)[0]) / 2.0
B['PDG2025'] = {'Mt_crit': Mtc, 'mu_star': mus, 'dMtcrit/dMh': dMh, 'dMtcrit per 0.0009 as': das,
                'track': {'Mt_crit': 170.9366, 'mu_star': 4.742e17, 'dMh': 0.4947, 'das': 0.3535},
                'Mt_crit_cut_at_reduced_MPl': Mtcrit_B(Mh0, a0, Qcut=MPL_RED)[0]}
B['shifts'] = {'dlam +0.0003': Mtcrit_B(Mh0, a0, dlam=0.0003)[0] - Mtc,
               'dyt +0.0005': Mtcrit_B(Mh0, a0, dyt=0.0005)[0] - Mtc,
               'RGE 2-loop minus 3-loop': Mtcrit_B(Mh0, a0, L=2)[0] - Mtc,
               'RGE loopOrder4 minus 3': Mtcrit_B(Mh0, a0, L=4)[0] - Mtc,
               'RGE loopOrder5 minus 3': Mtcrit_B(Mh0, a0, L=5)[0] - Mtc,
               'track': {'dlam': 0.0720, 'dyt': -0.0867, '2-3 loop': -0.1092}}
B['Mh_crit_at_Mt=172.56'] = brentq(lambda Mh: lam_min_B(172.56, Mh, a0)[0], 120, 140, xtol=1e-6)
B['Mh_crit_at_Mt=172.4'] = brentq(lambda Mh: lam_min_B(172.4, Mh, a0)[0], 120, 140, xtol=1e-6)


# FN double condition lambda(mu*) = beta_lambda(mu*) = 0 at fixed mu*
def fn_double(mu_star, als=a0):
    def eqs(xv):
        Mt, Mh = xv
        rr_ = S.run(3, Mt, matching_B(Mt, Mh, als, MW0), N=2, Qmin=1e3, Qmax=mu_star)
        return [rr_['k'][-1] * 100, rr_['betak'][-1] * 1e3]
    xv = fsolve(eqs, [172.0, 128.0], xtol=1e-12)
    return [float(xv[0]), float(xv[1])]


B['FN_double'] = {'1.2209e19': fn_double(MPL), '2.435e18': fn_double(MPL_RED), '1e18': fn_double(1e18),
                  'track': {'1.2209e19': [174.20, 132.03], '2.435e18': [172.61, 128.64], '1e18': [171.71, 126.77]}}

# potential at PDG central values (Mt = 172.56): lambda, 1-loop lambda_eff (own), V = e^{4Gamma} lambda_eff h^4/4
c_pdg = matching_B(172.56, Mh0, a0, MW0)
rr = S.run(3, 172.56, c_pdg, N=4001, Qmin=172.56, Qmax=1.3e19)
lnv0 = np.log(rr['v'][0])


def lam_eff_1loop(i):
    g, gp, yt, k = rr['g'][i], rr['gp'][i], rr['yt'][i], rr['k'][i]
    Gam = np.log(rr['v'][i]) - lnv0           # Gamma = int gamma dln mu = ln v(mu)/v(Mt)  (Landau gauge)
    e2G = np.exp(2 * Gam)
    terms = [(6, g ** 2 / 4, 5 / 6), (3, (g ** 2 + gp ** 2) / 4, 5 / 6), (-12, yt ** 2 / 2, 3 / 2),
             (1, 3 * k, 3 / 2), (3, k, 3 / 2)]
    s = sum(n * kap ** 2 * (np.log(abs(kap) * e2G) - cc) for n, kap, cc in terms)
    return k + s / K16, Gam


pot = {}
for h in (1e10, 1e12, 1e15, 1e17, MPL_RED, MPL):
    i = int(np.argmin(np.abs(np.log(rr['Q'] / h))))
    le, Gam = lam_eff_1loop(i)
    pot[f'{h:.4g}'] = {'lambda': float(rr['k'][i]), 'lambda_eff_1loop': float(le), 'Gamma': float(Gam),
                       'V_eff_GeV4': float(np.exp(4 * Gam) * le * h ** 4 / 4)}
pot['lambda zero crossing GeV'] = float(np.exp(np.interp(0.0, -rr['k'][:2000], np.log(rr['Q'][:2000]))))
pot['track'] = {'1e17': [-0.01047, -0.00974, -2.03e65], '2.435e18': [-0.01043, -0.00999, -7.48e70],
                '1.221e19': [-0.01022, -0.00991, -4.75e73], 'crossing': 4.21e10}
pot['|V(2.435e18)|/rho_Lambda'] = abs(pot[f'{MPL_RED:.4g}']['V_eff_GeV4']) / 2.52e-47
B['potential_PDG_central'] = pot
out['B_buttazzo_matching_smdr_running'] = B

# ================================================================ C. SMDR full matching (fit grid) + SMDR running
grid = json.load(open(GRID))['fits']
keys = ['g3', 'g', 'gp', 'yt', 'yb', 'ytau', 'k']
base = [f for f in grid if abs(f['inputs'][1] - 125.20) < 1e-9 and abs(f['inputs'][2] - 0.1180) < 1e-9]
Mts = np.array([f['inputs'][0] for f in base])
coef = {kk: np.polyfit(Mts, [f[kk] for f in base], 2) for kk in keys}
get = lambda Mt, Mh, als: [f for f in grid if np.allclose(f['inputs'], [Mt, Mh, als])][0]
dMh_ = {kk: (get(171.2, 125.70, 0.1180)[kk] - get(171.2, 124.70, 0.1180)[kk]) / 1.0 for kk in keys}
das_ = {kk: (get(171.2, 125.20, 0.1189)[kk] - get(171.2, 125.20, 0.1171)[kk]) / 0.0018 for kk in keys}


def matching_S(Mt, Mh, als):
    return {kk: float(np.polyval(coef[kk], Mt) + dMh_[kk] * (Mh - 125.20) + das_[kk] * (als - 0.1180)) for kk in keys}


interp_check = {}
for f in grid:
    p = f['inputs']
    m = matching_S(*p)
    interp_check[str(p)] = {kk: m[kk] - f[kk] for kk in ('g3', 'yt', 'k')}
C = {'interpolation_residuals(model - fit)': interp_check}
C['couplings_at_Mt: SMDR minus Buttazzo-fit'] = {}
for p in ((173.34, 125.15, 0.1184), (172.56, 125.20, 0.1180), (171.2, 125.20, 0.1180)):
    fS, fB = get(*p), matching_B(*p, MW=80.3692)
    C['couplings_at_Mt: SMDR minus Buttazzo-fit'][str(p)] = {kk: fS[kk] - fB[kk] for kk in ('g3', 'g', 'gp', 'yt', 'k')}


def Mtcrit_S(Mh, als, L=3, Qcut=MPL):
    f = lambda Mt: S.min_lambda(L, Mt, matching_S(Mt, Mh, als), Qcut=Qcut)[0]
    Mt = brentq(f, 168.0, 175.5, xtol=1e-6)
    return Mt, S.min_lambda(L, Mt, matching_S(Mt, Mh, als), Qcut=Qcut)[1]


C['Mtcrit_PDG2025'] = {f'loopOrder{L}': Mtcrit_S(Mh0, a0, L) for L in (2, 3, 4, 5)}
C['dMtcrit/dMh'] = (Mtcrit_S(Mh0 + 0.5, a0)[0] - Mtcrit_S(Mh0 - 0.5, a0)[0])
C['dMtcrit per 0.0009 as'] = (Mtcrit_S(Mh0, a0 + 0.0009)[0] - Mtcrit_S(Mh0, a0 - 0.0009)[0]) / 2
C['Mtcrit_ref(125.15,0.1184)'] = Mtcrit_S(125.15, 0.1184)
C['Mtcrit_Bednyakov(125.7,0.1185) loop3/loop5'] = [Mtcrit_S(125.7, 0.1185, 3), Mtcrit_S(125.7, 0.1185, 5)]
C['Mh_crit_at_Mt=172.56'] = brentq(lambda Mh: S.min_lambda(3, 172.56, matching_S(172.56, Mh, a0))[0], 120, 140,
                                   xtol=1e-6)
C['literature'] = {'Hiller et al. 2401.08811 (PDG 2024, Alam et al./SMDR-type matching, lambda>0 below MPl)': 171.04,
                   'Hiller eff': 171.10, 'Bednyakov et al. 1507.08833 at (125.7, 0.1185)': 171.44}
out['C_smdr_matching_smdr_running'] = C

# ================================================================ D. tension and held-out projections
def tension(Mtc, dMh, das, s_match, s_crit=0.032):
    rows = {}
    s_par = np.hypot(dMh * sMh, das * sa / 0.0009)
    for lab, (Mt, sMt) in (('direct 172.56+-0.31', (172.56, 0.31)), ('xsec 172.4+-0.7', (172.4, 0.7))):
        for blab, extra in (('exp+param+matching', 0.0), ('+criterion', s_crit), ('+criterion+0.3 MC', np.hypot(s_crit, 0.3))):
            st = np.sqrt(sMt ** 2 + s_par ** 2 + s_match ** 2 + extra ** 2)
            z = (Mt - Mtc) / st
            rows[f'{lab} | {blab}'] = {'Delta': Mt - Mtc, 'sigma': st, 'z': z, 'p_one_sided': norm.sf(z)}
    return rows


Bp = B['PDG2025']
s_match_B = float(np.hypot(B['shifts']['dlam +0.0003'], B['shifts']['dyt +0.0005']))
D = {'B (Buttazzo matching)': tension(Bp['Mt_crit'], Bp['dMtcrit/dMh'], Bp['dMtcrit per 0.0009 as'], s_match_B)}
McS = C['Mtcrit_PDG2025']['loopOrder3'][0]
D['S (SMDR matching), same error budget'] = tension(McS, C['dMtcrit/dMh'], C['dMtcrit per 0.0009 as'], s_match_B)
# matching systematic = |S - B| added in quadrature
s_sys = abs(McS - Bp['Mt_crit'])
D['S with |S-B| matching difference added'] = tension(McS, C['dMtcrit/dMh'], C['dMtcrit per 0.0009 as'],
                                                       float(np.hypot(s_match_B, s_sys)))
D['sigma_match_B'] = s_match_B
D['matching difference S - B (GeV)'] = McS - Bp['Mt_crit']


def heldout(Mtc, dMh, das, sm, sMt=0.04, sas=0.00012, sMh_=0.004, s_crit=0.032):
    st = np.sqrt(sMt ** 2 + (dMh * sMh_) ** 2 + (das * sas / 0.0009) ** 2 + sm ** 2 + s_crit ** 2)
    return {'sigma': st, 'z(direct)': (172.56 - Mtc) / st, 'z(xsec)': (172.4 - Mtc) / st}


D['heldout FCC-ee (track inputs), B'] = heldout(Bp['Mt_crit'], Bp['dMtcrit/dMh'], Bp['dMtcrit per 0.0009 as'], s_match_B)
D['heldout FCC-ee, S'] = heldout(McS, C['dMtcrit/dMh'], C['dMtcrit per 0.0009 as'], s_match_B)
D['heldout FCC-ee, matching halved, B'] = heldout(Bp['Mt_crit'], Bp['dMtcrit/dMh'], Bp['dMtcrit per 0.0009 as'],
                                                  s_match_B / 2)
D['heldout HL-LHC (0.2, 0.00045, 0.02), B'] = heldout(Bp['Mt_crit'], Bp['dMtcrit/dMh'], Bp['dMtcrit per 0.0009 as'],
                                                     s_match_B, 0.2, 0.00045, 0.02)
D['heldout HL-LHC, S'] = heldout(McS, C['dMtcrit/dMh'], C['dMtcrit per 0.0009 as'], s_match_B, 0.2, 0.00045, 0.02)
D['track'] = {'z': [3.34, 2.84, 1.84], 'FCC sigma': 0.1325, 'FCC z': 12.25, 'HL z': 5.57}
out['D_tension_heldout'] = D


# ================================================================ E. CDL quench: where does eps reach the threshold?
def cdl(dMt, Mtc=Bp['Mt_crit']):
    Mt = Mtc + dMt
    rr_ = S.run(3, Mt, matching_B(Mt, Mh0, a0, MW0), N=6001, Qmin=1e3, Qmax=MPL)
    h = rr_['Q']
    V = (rr_['v'] / rr_['v'][0]) ** 4 * rr_['k'] * h ** 4 / 4
    i = int(np.argmin(V))
    eps = -V[i]
    if eps <= 0:
        return {'dMt': dMt, 'eps': 0.0}
    mask = (V > 0) & (h < h[i])
    sig_pos = np.trapezoid(np.sqrt(2 * V[mask]), h[mask])
    # alternative: tension of the degenerate potential, int sqrt(2|V|) over the whole path 0..h_min
    sig_all = np.trapezoid(np.sqrt(2 * np.abs(V[:i + 1])), h[:i + 1])
    M = MPL_RED
    return {'dMt': dMt, 'h_min': float(h[i]), 'eps': float(eps), 'sigma_pos': float(sig_pos),
            'sigma_all': float(sig_all), 'eps/thr(pos)': float(eps / (0.75 * sig_pos ** 2 / M ** 2)),
            'eps/thr(all)': float(eps / (0.75 * sig_all ** 2 / M ** 2))}


out['E_CDL'] = [cdl(d) for d in (1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.05)]
out['E_CDL_track'] = {'0.001': 5934, '0.01': 7.6e8, '0.05': 2.9e17}

# ================================================================ F. decomposition of the matching difference; 1-loop 'eff'
def Mtcrit_mixed(swap):
    """Buttazzo matching with the couplings in `swap` replaced by SMDR's (same Mt, Mh, as)."""
    def cpl(Mt):
        c = matching_B(Mt, Mh0, a0, MW0)
        cs = matching_S(Mt, Mh0, a0)
        for kk in swap:
            c[kk] = cs[kk]
        return c
    f = lambda Mt: S.min_lambda(3, Mt, cpl(Mt))[0]
    return brentq(f, 166.0, 178.0, xtol=1e-6)


F = {'B (all Buttazzo)': Bp['Mt_crit']}
for sw in (['yt'], ['g3'], ['k'], ['g', 'gp'], ['yt', 'g3', 'k', 'g', 'gp', 'yb', 'ytau']):
    F['B with SMDR ' + '+'.join(sw)] = Mtcrit_mixed(sw)


def lam_eff_min(Mt, cfun, L=3):
    rr_ = S.run(L, Mt, cfun(Mt), N=401, Qmin=1e8, Qmax=MPL)
    lnv0_ = np.log(rr_['v'][0])
    vals = []
    for i in range(len(rr_['Q'])):
        g, gp, yt, k = rr_['g'][i], rr_['gp'][i], rr_['yt'][i], rr_['k'][i]
        e2G = np.exp(2 * (np.log(rr_['v'][i]) - lnv0_))
        terms = [(6, g ** 2 / 4, 5 / 6), (3, (g ** 2 + gp ** 2) / 4, 5 / 6), (-12, yt ** 2 / 2, 3 / 2),
                 (1, 3 * k, 3 / 2), (3, k, 3 / 2)]
        vals.append(k + sum(n * kap ** 2 * (np.log(abs(kap) * e2G) - cc) for n, kap, cc in terms) / K16)
    return float(np.min(vals))


cB = lambda Mt: matching_B(Mt, Mh0, a0, MW0)
cS = lambda Mt: matching_S(Mt, Mh0, a0)
F['eff_1loop_B'] = brentq(lambda Mt: lam_eff_min(Mt, cB), 168, 175, xtol=1e-5)
F['eff_1loop_S'] = brentq(lambda Mt: lam_eff_min(Mt, cS), 168, 175, xtol=1e-5)
F['eff_1loop_minus_msbar_B'] = F['eff_1loop_B'] - Bp['Mt_crit']
F['track eff(2-loop) - msbar'] = 0.0641
F['as needed for stability at Mt=172.4 (S)'] = 0.1180 + (172.4 - McS) / C['dMtcrit per 0.0009 as'] * 0.0009
F['Hiller: as needed at Mt=172.4'] = 0.1215
F['Mh_crit at Mt=172.4 (S)'] = brentq(lambda Mh: S.min_lambda(3, 172.4, matching_S(172.4, Mh, a0))[0], 120, 140,
                                      xtol=1e-6)
F['Hiller: Mh needed at Mt=172.4'] = 127.97
out['F_decomposition_eff'] = F

out['runtime_s'] = time.time() - t_start
with open(OUT, 'w') as f:
    json.dump(out, f, indent=1, default=float)
print(json.dumps(out, indent=1, default=float))
