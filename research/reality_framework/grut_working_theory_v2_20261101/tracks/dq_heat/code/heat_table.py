#!/usr/bin/env python3
"""dq_heat (Claude, own code), item 3: the hidden heat budget R(e) of the actual M1r release and its timing.

R = (deposited conformal energy)/E_X, E_X the X release energy of the same run (comoving H_i^4, a = 1 at the slice):
  R_vec     : inflationary (piece A, GMR-type) + release-era (piece B: UV modes and separate-universe angle) longitudinal
              vectors, N_A from vector_modes.json / sep_universe.json, deposit m_A <a_A> per vector (A -> X Xbar),
              <a_A> = (pi K / 2 Gamma_A)^{1/2}, Gamma_A = e^2 m_A beta^3/(48 pi), m_A = 6 e h = 240 e H_i (beta = 1 to 1e-30);
  R_Y,hom   : homogeneous Y_L + Y_H (Y_H cascades to Y_L first): (N_YL + N_YH) m_L <a_L>, Gamma_L = 3 kap^2 m_L/(256 pi^3);
  R_Y,k     : finite-k Y_L from the release (ymodes YS/XYS): direct k < 15 (cascade-free) and the subtraction estimate;
  R_cascade : Y_H -> Y_L + A_L vectors, N_YH m_A <a_A>  (scales e^{-1/2} at fixed N_YH); sigma -> A_L A_L adds 2 N_sigma.
Timing against the phi onset: a_on from the relaxed Hartree plasma (copied solver, release-state invariant
W = |eps|^{-4/3}(1 + R)), converted to this frame and H_i: a_on = a_on,spot (N_X/N_spot)^{1/3} (H_i/H_spot).
S/Q from the same solver (wt_critic convention), against 3.4716 (bulk) and 2.9284 (Q = 2e5).
Both renormalisation conditions of the induced portals: vanish at M_Y (benchmark) and at M_P.
"""
import json, math, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hartree_copy as HI

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
WT = HERE.parents[1]
KAP, MY, HH, MU = 0.0034586847211042523, 30.0, 40.0, 7.159456159513793
E_TABLE = [1e-6, 1.5e-6, 2e-6, 3e-6, 5e-6, 1e-5]
MRED, GSTAR = 2.435e18, 106.75
H_CAL = {'freeX_1.29e11': 1.29e11, 'Qmatter_1.51e11': 1.51e11}
SIG = json.load(open(WT / 'wt_critic' / 'results' / 'chain_consistency.json'))['sigma_min']
W_THR = json.load(open(WT / 'wt_critic' / 'results' / 'chain_consistency.json'))['W_thresholds']
EPS_SPOT = 0.14571744282991028          # kin_source spot point |eps| (frame of the copied Hartree solver)


def mediators(K):
    mL = math.sqrt(MY ** 2 - math.sqrt(2) * MU * HH); mH = math.sqrt(MY ** 2 + math.sqrt(2) * MU * HH)
    GL = 3 * KAP ** 2 * mL / (256 * math.pi ** 3)
    GH = MU ** 2 * (mH ** 2 - mL ** 2) / (8 * math.pi * mH ** 3)        # Y_H -> Y_L + A_L (exact Goldstone vertex sqrt2 mu)
    msig = math.sqrt(2 * 0.5 * HH ** 2); Gsig = 0.5 ** 2 * HH ** 2 / (8 * math.pi * msig)   # sigma -> G G (lambda_S h sigma G^2)
    am = lambda G: math.sqrt(math.pi * K / (2 * G))
    aend = lambda G: math.sqrt(2 * K * math.log(1e6) / G)                 # 1 - 1e-6 decayed (Rayleigh law in radiation)
    return dict(m_L=mL, m_H=mH, Gamma_L=GL, Gamma_H=GH, m_sigma=msig, Gamma_sigma=Gsig, a_mean_L=am(GL), a_end_L=aend(GL),
                a_mean_H=am(GH), a_GammaH_eq_H=math.sqrt(K / GH), a_sigma=math.sqrt(K / Gsig), p_star_cascade=(mH ** 2 - mL ** 2) / (2 * mH))


def vec_decay(e, K):
    m = 240.0 * e; G = e * e * m / (48 * math.pi)
    return dict(m_A=m, Gamma_A=G, a_mean=math.sqrt(math.pi * K / (2 * G)), a_end=math.sqrt(2 * K * math.log(1e6) / G))


def T_vis(a, K, H_GeV):
    """visible temperature at comoving a in radiation, H = K H_i / a^2, g* = 106.75 (all decays have T_vis > TeV)."""
    Hp = K * H_GeV / a ** 2
    return math.sqrt(Hp * MRED / math.sqrt(math.pi ** 2 * GSTAR / 90))


def gmr_NA(e, Rr, H_GeV=1.29e11):
    """[R] Graham-Mardon-Rajendran: Omega_A/Omega_DM = sqrt(m/6e-6 eV)(H_I/1e14 GeV)^2 -> n_A/s, times this frame's a^3 s/H_i^3
    (late radiation rho_r a^4 = Rr in M_r^2 H_i^2 units); H_i-independent.  D = 4.4e-10 GeV per unit entropy."""
    m = 240.0 * e * H_GeV
    Om = math.sqrt(m / 6e-15) * (H_GeV / 1e14) ** 2
    nA_s = 4.4e-10 * Om / m
    T = (30 * Rr * MRED ** 2 * H_GeV ** 2 / (math.pi ** 2 * GSTAR)) ** 0.25      # a T at a = 1 units
    a3s = 2 * math.pi ** 2 / 45 * GSTAR * T ** 3 / H_GeV ** 3
    return nA_s * a3s


def finite_k():
    out = {}
    for tag in ('YS_r11', 'XYS_r11'):
        d = json.load(open(RES / ('ymodes_%s.json' % tag)))
        K = np.array(d['kgrid']); lk = np.log(K)
        T = d['times'][-1]
        iL = [i for i, c in enumerate(T['eig']) if c['Y'] > 0.5 and 400 < c['m2'] < 600][0]
        iH = [i for i, c in enumerate(T['eig']) if c['Y'] > 0.5 and 1200 < c['m2'] < 1400][0]
        iS = [i for i, c in enumerate(T['eig']) if c['S'] > 0.5 and 1500 < c['m2'] < 1700][0]
        res = dict(tau=T['tau'], a=T['a'])
        for nm, i in (('YL', iL), ('YH', iH), ('sigma', iS)):
            n = np.array(T['n_k'][i]); integ = K ** 3 * n / (2 * np.pi ** 2)
            sel = K <= 15.0
            res['N_%s_all_k' % nm] = float(np.trapezoid(integ, lk))
            res['N_%s_k_le_15' % nm] = float(np.trapezoid(integ[sel], lk[sel]))
        # plateau history of N_YL (all k) for the cascade subtraction
        res['N_YL_all_k_history'] = [[t['tau'], t['N'][iL]] for t in d['times'] if t['tau'] >= 10]
        out[tag] = res
    return out


def SQ_of_W(W, cache=[None]):
    W0s = EPS_SPOT ** (-4.0 / 3.0)
    Reff = W / W0s - 1.0
    st, g = HI.solve(Reff, True, 0.0, cache[0])
    if (not st['success']) or max(abs(x) for x in st['resid']) > 1e-9:
        st, g = HI.solve(Reff, True, 0.0, None)
    cache[0] = g
    a = HI.onset_a(st); mc = HI.MX * a / HI.H
    return dict(W=W, S_per_Q=st['S'] / HI.NX, a_on_spot_frame=a, T_on_over_mX=st['T'] / mc, nQ_on_mX3=HI.NX / mc ** 3,
                Xbar_per_Q=st['n_minus'] / HI.NX, cond_frac=st['cond_frac'], resid=max(abs(x) for x in st['resid']))


def budget(bgfile, vmfile, su, fk, label):
    bg = json.load(open(RES / bgfile)); V = bg['validation']
    K = bg['late']['K_end']; EX = V['EX']; NX = abs(V['NX_400_1000']); eps = abs(V['eps'])
    W0 = eps ** (-4.0 / 3.0)
    md = mediators(K)
    dL = md['m_L'] * md['a_mean_L'] / EX
    vm = {c['e']: c for c in json.load(open(RES / vmfile))['cases'] if c['variant'] == 'full'}
    # dense resolution of the release-era band k >= 1 (e-independent, vector_band.py); benchmark background only
    band_file = RES / ('vector_band.json' if vmfile == 'vector_modes.json' else 'vector_band_planck.json')
    NA_band = None
    if band_file.exists():
        NA_band = json.load(open(band_file))['sets']['full_1e-06']['N_A_k_ge_1']
    suB = {r['e']: r['N_A_pieceB_SU'] for r in su['per_e']}
    NYLk_direct = fk['YS_r11']['N_YL_k_le_15']
    # subtraction estimate: plateau (all k <= 80, YS, tau >= 14) minus the linear-theory cascade born before a = 80/p*
    hist = np.array(fk['YS_r11']['N_YL_all_k_history']); plateau = float(np.mean(hist[hist[:, 0] >= 20, 1]))
    bgm = json.load(open(RES / 'background.json'))
    tab = bgm['F_table']
    a_cut = 80.0 / md['p_star_cascade']
    t_cut = float(np.interp(a_cut, [r['a'] for r in tab], [r['tau'] for r in tab]))
    casc_lin = md['Gamma_H'] * bgm['validation']['NYH'] * (t_cut - 0.7)
    NYLk_sub = plateau - casc_lin
    NYHk = max(fk['YS_r11']['N_YH_k_le_15'], 0.0); Nsigk = max(fk['YS_r11']['N_sigma_k_le_15'], 0.0)
    rows = []
    for e in sorted(vm):
        c = vm[e]; vd = vec_decay(e, K); dA = vd['m_A'] * vd['a_mean'] / EX
        NA_uv = NA_band if NA_band is not None else c['N_A_k_gt_1']
        NA_tot = c['N_A_k_le_1'] + NA_uv
        R_vecA = NA_tot * dA; R_vecB = suB.get(e, 0.0) * dA
        R_vec = R_vecA + R_vecB
        R_Yh = (V['NYL'] + V['NYH']) * dL
        NYLk_xys = fk['XYS_r11']['N_YL_k_le_15']
        R_Yk_lo = (NYLk_direct + NYHk) * dL                                   # YS, k < 15 (cascade-free), direct
        R_Yk = (max(NYLk_sub, NYLk_xys, NYLk_direct) + NYHk) * dL              # table entry: the larger estimate
        R_Yk_hi = (fk['YS_r11']['N_YL_all_k'] + NYHk) * dL                   # hard upper: every k <= 80 incl. cascade products
        R_casc = V['NYH'] * dA + 2 * V['Nsig'] * dA + NYHk * dA + 2 * Nsigk * dA
        R_tot = R_vec + R_Yh + R_Yk + R_casc
        row = dict(e=e, R_vec=R_vec, R_vec_modes=R_vecA, R_vec_k_le_1=c['N_A_k_le_1'] * dA, R_vec_k_ge_1_release_era=NA_uv * dA,
                   R_vec_pieceB_SU=R_vecB, N_A_k_le_1=c['N_A_k_le_1'], N_A_k_ge_1=NA_uv, N_A_k_ge_1_source='dense band' if NA_band is not None else 'coarse grid',
                   R_Y_hom=R_Yh, R_Y_k=R_Yk, R_Y_k_low=R_Yk_lo, R_Y_k_hard_upper=R_Yk_hi, R_cascade=R_casc, R_total=R_tot,
                   R_total_upper=R_tot - R_Yk + R_Yk_hi, N_A=NA_tot, N_A_pieceB_SU=suB.get(e, 0.0), N_A_per_net_charge=NA_tot / NX,
                   N_A_GMR_formula=gmr_NA(e, bg['validation']['Rr_end']), ratio_to_GMR_formula=NA_tot / gmr_NA(e, bg['validation']['Rr_end']),
                   m_A_over_H=vd['m_A'], Gamma_A_over_H=vd['Gamma_A'], a_mean_vec=vd['a_mean'], a_end_vec=vd['a_end'],
                   deposit_per_vector_EX=dA, W=W0 * (1 + R_tot), W_vec=W0 * R_vec)
        rows.append(row)
    return dict(label=label, K=K, EX=EX, NX=NX, eps=eps, W0=W0, NYL_hom=V['NYL'], NYH_hom=V['NYH'], Nsig_hom=V['Nsig'],
                mediators=md, deposit_per_YL_EX=dL, N_YL_k_direct_k_le_15=NYLk_direct, N_YL_k_subtraction=NYLk_sub,
                N_YL_k80_plateau=plateau, cascade_linear_before_cut=casc_lin, a_cut_cascade=a_cut, tau_cut_cascade=t_cut,
                N_YH_k=NYHk, N_sigma_k=Nsigk, rows=rows)


def main():
    t0 = time.time()
    su = json.load(open(RES / 'sep_universe.json'))
    fk = finite_k()
    out = dict(method=__doc__, finite_k=fk, budgets={})
    for bgf, vmf, lab in (('background.json', 'vector_modes.json', 'portals_vanish_at_MY_benchmark'),
                          ('background_planck.json', 'vector_modes_planck.json', 'portals_vanish_at_MP')):
        b = budget(bgf, vmf, su, fk, lab)
        for r in b['rows']:
            s = SQ_of_W(r['W'])
            r['S_per_Q'] = s['S_per_Q']; r['margin_bulk'] = s['S_per_Q'] / SIG['bulk_Q_inf']; r['margin_Q2e5'] = s['S_per_Q'] / SIG['Q_2e5']
            r['Xbar_per_Q_preonset'] = s['Xbar_per_Q']; r['T_on_over_mX'] = s['T_on_over_mX']; r['nQ_on_mX3'] = s['nQ_on_mX3']
            for hk, Hv in H_CAL.items():
                a_on = s['a_on_spot_frame'] * (b['NX'] / HI.NX) ** (1.0 / 3.0) * (Hv / HI.H)
                md = b['mediators']
                r['timing_' + hk] = dict(a_on=a_on, a_on_over_a_vec_mean=a_on / r['a_mean_vec'], a_on_over_a_vec_end=a_on / r['a_end_vec'],
                                         a_on_over_a_YL_mean=a_on / md['a_mean_L'], a_on_over_a_YL_end=a_on / md['a_end_L'],
                                         a_on_over_a_cascade_emission=a_on / md['a_mean_H'], a_on_over_a_sigma=a_on / md['a_sigma'],
                                         T_vis_vec_decay_GeV=T_vis(r['a_mean_vec'], b['K'], Hv), T_vis_YL_decay_GeV=T_vis(md['a_mean_L'], b['K'], Hv),
                                         T_vis_cascade_GeV=T_vis(md['a_mean_H'], b['K'], Hv), T_vis_onset_GeV=T_vis(a_on, b['K'], Hv))
            print('%-32s e=%.2e R_vec=%.4f (k<=1 %.4f, k>=1 %.4f, SU %.1e) R_Yhom=%.4f R_Yk=%.1e R_casc=%.4f R_tot=%.4f W=%.2f S/Q=%.3f (x%.3f) a_on(1.29)=%.3e a_on/<a_A>=%.2e' % (
                b['label'], r['e'], r['R_vec'], r['R_vec_k_le_1'], r['R_vec_k_ge_1_release_era'], r['R_vec_pieceB_SU'], r['R_Y_hom'], r['R_Y_k'],
                r['R_cascade'], r['R_total'], r['W'], r['S_per_Q'], r['margin_bulk'], r['timing_freeX_1.29e11']['a_on'],
                r['timing_freeX_1.29e11']['a_on_over_a_vec_mean']), flush=True)
        # critical e for survival (bulk, Q = 2e5): W(e) = W_thr; R_vec ~ e^{-1} interpolated in log-log, cascade e^{-1/2}
        rr = b['rows']; le = np.log([r['e'] for r in rr]); lW = np.log([r['W'] for r in rr])
        crit = {}
        for nm in ('bulk_Q_inf', 'Q_2e5', 'Q_bind_71263'):
            Wt = W_THR[nm]
            f = lambda x: np.interp(x, le, lW) - math.log(Wt)
            if f(le[0]) > 0 > f(le[-1]):
                x = brentq(f, le[0], le[-1]); crit[nm] = dict(e_crit=math.exp(x), W_threshold=Wt)
            else:
                crit[nm] = dict(e_crit=None, W_threshold=Wt, W_range=[float(np.exp(lW[0])), float(np.exp(lW[-1]))])
        b['e_crit'] = crit
        # check: direct S/Q at the bulk e_crit
        if crit['bulk_Q_inf']['e_crit']:
            ec = crit['bulk_Q_inf']['e_crit']
            Wc = math.exp(np.interp(math.log(ec), le, lW))
            crit['bulk_Q_inf']['S_per_Q_check'] = SQ_of_W(Wc)['S_per_Q']
        print(b['label'], 'e_crit', json.dumps(crit, default=float), flush=True)
        out['budgets'][b['label']] = b
    # e = 1e-7 reference row (benchmark) and the old spot comparison
    vm = {c['e']: c for c in json.load(open(RES / 'vector_modes.json'))['cases'] if c['variant'] == 'full'}
    if 1e-7 in vm:
        b = out['budgets']['portals_vanish_at_MY_benchmark']
        out['reference_e1e-7'] = dict(R_vec=vm[1e-7]['R_vec'], N_A=vm[1e-7]['N_A'])
    out['seconds'] = time.time() - t0
    json.dump(out, open(RES / 'heat_table.json', 'w'), indent=1, default=float)
    print('seconds', out['seconds'])


if __name__ == '__main__':
    main()
