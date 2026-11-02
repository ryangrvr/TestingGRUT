#!/usr/bin/env python3
"""Track mf_higgs: controls, the non-equilibrium inflationary tail, R < 0 windows and lines, the G6 portal window,
the plasma-modified mean-field upper edge, induced Higgs masses, and the thermal-only plasma reach (own code).

  python3 -B mh_extra.py controls     -> results/controls.json
  python3 -B mh_extra.py tail         -> results/inflation_tail.json
  python3 -B mh_extra.py windows      -> results/windows.json
  python3 -B mh_extra.py g6 <hist>    -> results/g6_portal_<hist>.json
  python3 -B mh_extra.py upperpl      -> results/upper_plasma.json
  python3 -B mh_extra.py induced      -> results/induced.json
  python3 -B mh_extra.py thermal      -> results/thermal_reach.json
Every portal, renormalisation condition and criterion here is a model choice.
"""
import sys
sys.dont_write_bytecode = True
import json
import math
import time

import numpy as np

import mh_common as C
import mh_modes as M
import mh_run as RN
import mh_inflation as INF

LN10 = math.log(10.0)


def cmd_controls():
    name, hist = 'HFa', 'EXP'
    bg = RN.bgc(name, hist)
    xi = RN.xi_at(bg.Hi_GeV)
    base = RN.run_cached(name, hist, xi)
    Nob = base['res']['N_obs']
    Mts = (172.86, 172.88, 173.0)
    ref = {f'{m}': RN.D_of(base, m, Nob)['D'] for m in Mts}
    out = dict(reference=ref, variants={})
    jobs = [('R lines off after the join (mean R)', dict(osc=False)), ('t0 = -5', dict(t0=-5.0)),
            ('n = 400 modes', dict(n=400)), ('cK = 0.025', dict(cK=0.025)), ('Kmax = 15', dict(Kmax=15.0)),
            ('Hartree at M_t = 172.92', dict(sm_mt=172.92))]
    for lab, kw in jobs:
        R = RN.run_cached(name, hist, xi, **kw)
        out['variants'][lab] = {f'{m}': RN.D_of(R, m, Nob)['D'] - ref[f'{m}'] for m in Mts}
        out['variants'][lab]['wronskian_err_end'] = R['res']['wronskian_err_end']
        print(lab, out['variants'][lab], flush=True)
    R = M.run(bg, xi, RN.smt(172.57), hartree=False)
    R['res']['Hi_GeV'] = bg.Hi_GeV
    out['variants']['Hartree off'] = {f'{m}': RN.D_of(R, m, Nob)['D'] - ref[f'{m}'] for m in Mts}
    print('Hartree off', out['variants']['Hartree off'], flush=True)
    # PAULI0 R lines off
    bgp = RN.bgc(name, 'PAULI0')
    Rp = RN.run_cached(name, 'PAULI0', xi)
    Rq = RN.run_cached(name, 'PAULI0', xi, osc=False)
    out['variants']['PAULI0: R lines off'] = {f'{m}': RN.D_of(Rq, m, Rp['res']['N_obs'])['D'] - RN.D_of(Rp, m, Rp['res']['N_obs'])['D']
                                             for m in Mts}
    print('PAULI0 R lines off', out['variants']['PAULI0: R lines off'], flush=True)
    C.jdump(out, C.RES / 'controls.json')


def cmd_tail():
    """non-equilibrium Gaussian count against the instantaneous tree de Sitter barrier, from the EXP mode run."""
    out = {}
    for name in ('HFa', 'HFb'):
        bg = RN.bgc(name, 'EXP')
        xi = RN.xi_at(bg.Hi_GeV)
        R = RN.run_cached(name, 'EXP', xi)
        rows, spec, ks, wq = R['rows'], R['spec'], R['ks'], R['wq']
        te = R['res']['tau_endinf']
        d = C.load_bench(name)
        W = INF.waterfall_table(name)
        N0 = R['res']['N_obs'] + W['N_end']
        Ne = W['N_end']
        o = {}
        for Mt in (172.57, 172.86, 173.45, 174.05, 174.6):
            sm = INF.SMfull(Mt)
            best = -np.inf
            rec = []
            for j in range(len(rows)):
                t, a, H, Rc = rows[j, 0], rows[j, 1], rows[j, 2], rows[j, 3]
                if t > te:
                    break
                if t < -2.5 or Rc <= 0:
                    continue
                Nb = Ne - math.log(a)
                m = ks < a * H
                s2 = float(np.sum(wq[m] * spec[m, j]))
                Bv, hb = INF.IHI.barrier_B(sm, H * bg.Hi_GeV, Rc * bg.Hi_GeV ** 2, 0.1, False)
                if hb is None or not np.isfinite(Bv):
                    continue
                phib = math.exp(float(sm.cpl(hb)['Gam'])) * hb / bg.Hi_GeV
                x = phib ** 2 / (2 * s2)
                lp = (3 * (N0 - Nb) + math.log1p(x) - x) / LN10
                rec.append((t, Nb, Rc / H ** 2, s2, phib, x, lp))
                best = max(best, lp)
            o[f'{Mt}'] = dict(best_log10_patches=best, rows=rec[::6] + rec[-6:])
            print(name, Mt, 'best log10 patches %.1f' % best, flush=True)
        out[name] = o
        C.jdump(out, C.RES / 'inflation_tail.json')


def cmd_windows():
    """R < 0 windows after eps_H = 1 (start, duration, R_min/H^2, tachyonic exponent at xi(H_i)) and the Mathieu q of
    the 2 m_chi and 2 m_sigma R lines for every history."""
    out = {}
    for hist in ('EXP', 'PAULI0', 'PAULIb'):
        bg = RN.bgc('HFa', hist)
        xi = RN.xi_at(bg.Hi_GeV)
        te = bg.tau_endinf
        t = np.arange(te - 0.02, 6.0, 2e-5)
        aHR = np.array([bg.a_H_R(x) for x in t])
        H, R = aHR[:, 1], aHR[:, 2]
        neg = R < 0
        wins = []
        i = 0
        while i < len(t) and len(wins) < 400:
            if neg[i]:
                j = i
                while j < len(t) and neg[j]:
                    j += 1
                seg = slice(i, j)
                k = i + int(np.argmin(R[seg]))
                wins.append(dict(t_start=float(t[i]), dur_H=float((t[min(j, len(t) - 1)] - t[i]) * H[k]),
                                 R_min_over_H2=float(R[k] / H[k] ** 2), R_min=float(R[k]),
                                 exponent=float(np.trapezoid(np.sqrt(np.maximum(-xi * R[seg], 0)), t[seg]))))
                i = j
            else:
                i += 1
        expo = np.array([w['exponent'] for w in wins])
        frac = float(np.mean(neg))
        # lines: q = 2 xi R_osc / omega^2 (h'' + [k^2/a^2 + xi Rbar + xi R_osc cos(omega t)] h = 0, Mathieu q = 2 Delta/omega^2)
        lines = {}
        for x in (0.3, 0.6, 1.0, 3.0):
            Hx = float(np.interp(x, bg.tt, bg.HH))
            fc = float(np.interp(x, bg.tt, bg.fchi)); fs = float(np.interp(x, bg.tt, bg.fsig))
            lines[f'{x}'] = dict(q_chi=2 * xi * 9 * Hx ** 2 * fc / C.OMEGA_CHI ** 2,
                                 q_sig=2 * xi * 9 * Hx ** 2 * fs / C.OMEGA_SIG ** 2,
                                 omega_chi_over_H=C.OMEGA_CHI / Hx, R_osc_chi_over_H2=9 * fc, R_osc_sig_over_H2=9 * fs)
        out[hist] = dict(n_windows=len(wins), first=wins[:5], R_min_over_H2_min=float(min(w['R_min_over_H2'] for w in wins)),
                         max_exponent=float(expo.max()), sum_exponent_to_t6=float(expo.sum()),
                         fraction_time_R_negative=frac, lines=lines, xi_H=xi)
        print(hist, len(wins), 'first', wins[:2], 'max expo %.4f sum %.3f frac %.3f' % (expo.max(), expo.sum(), frac),
              lines['0.3'], flush=True)
    C.jdump(out, C.RES / 'windows.json')


# ------------------------------------------------------------------------------------------- G6 portal
def chi_factory(bg, name):
    """1 - chi(t)^2/v^2 on the history: exact bench chi where homogeneous, rebuilt oscillation after the join:
    chi = v (1 + delta), 1 - chi^2/v^2 = -2 (A/v) cos(m_chi t) + (A/v)^2 (mean 2<delta^2>/v^2 from the virial of the
    quartic valley), A = sqrt(2 rho_chi)/m_chi."""
    d = C.load_bench(name)
    h = d['h']
    v = d['params']['v_MP']
    mchi = d['params']['m_chi_over_Hi']
    mP = C.mP_REF
    tb, chib = h['t'], h['chi_MP']
    tj = bg.t_join

    def f(t):
        if t <= tj or not np.isfinite(np.interp(t, tb, chib)) and t <= 0.6:
            c = float(np.interp(t, tb, chib))
            return 1.0 - (c / v) ** 2
        Hx = float(np.interp(t, bg.tt, bg.HH))
        fc = float(np.interp(t, bg.tt, bg.fchi))
        rho = fc * 3 * Hx ** 2 * mP ** 2
        Av = math.sqrt(2 * rho) / mchi / (v * mP)
        return -2 * Av * math.cos(mchi * (t - tj)) + Av * Av
    return f


def cmd_g6(hist='EXP'):
    """the G6 portal lam_chiH (chi^2 - v^2)|H|^2 with lam_chiH < 0: m_H^2(t) = Lam H_i^2 (1 - chi^2/v^2),
    Lam = |lam_chiH| v^2 / H_i^2.  Inflation: HM with m2x = Lam H_i^2; after: modes with the oscillating mass,
    D (grav only, period-mean mass in the criterion), the produced variance, runaway."""
    name = 'HFa'
    path = C.RES / f'g6_portal_{hist}.json'
    try:
        out = json.load(open(path))
    except Exception:
        out = {}
    bg = RN.bgc(name, hist)
    d = C.load_bench(name)
    vH = d['params']['v_MP'] * C.mP_REF
    out['conversion'] = dict(lam_chiH_per_Lam=1.0 / vH ** 2, v_over_Hi=vH)
    f1 = chi_factory(bg, name)
    xi = RN.xi_at(bg.Hi_GeV)
    W = INF.waterfall_table(name)
    for Lam in (0.3, 1.0, 1.5, 3.0, 30.0, 300.0, 3000.0, 1e4, 3e4):
        k = f'{Lam:g}'
        if k in out:
            continue
        t0w = time.time()
        m2x = lambda t, L=Lam: L * f1(t)
        Kmax = 12.0 if Lam < 100 else 60.0 if Lam < 1000 else 300.0
        R = RN.run_cached(name, hist, xi, t0=-1.0 if Lam >= 3 else -3.0, Kmax=Kmax, n=300, m2x=m2x, tag=f'_g6L{Lam:g}')
        rows = R['rows']
        # period mean of the extra mass at the saves: window = exactly two periods of the rebuilt m_chi line
        wwin = 2 * 2 * math.pi / d['params']['m_chi_over_Hi']
        tt = np.linspace(rows[0, 0], rows[-1, 0], 400001)
        mm = np.array([m2x(x) for x in tt])
        cm = np.concatenate([[0.0], np.cumsum(0.5 * (mm[1:] + mm[:-1]) * np.diff(tt))])
        mbar = np.empty(len(rows))
        for i, x in enumerate(rows[:, 0]):
            lo, hi = max(x - wwin / 2, tt[0]), min(x + wwin / 2, tt[-1])
            mbar[i] = (np.interp(hi, tt, cm) - np.interp(lo, tt, cm)) / (hi - lo) if x > bg.tau_endinf else m2x(x)
        Nob = R['res']['N_obs']
        ie = int(np.searchsorted(rows[:, 0], bg.tau_endinf))
        row = dict(Lam=Lam, lam_chiH=-Lam / vH ** 2, runaway=R['res']['runaway'], t_stop=R['res']['t_end'],
                   max_sigma2=R['res']['max_sigma2'], sigma2_at_endinf=float(rows[ie, 6]),
                   max_nk=float(np.max(R['res']['nk_end'])), seconds=time.time() - t0w,
                   mean_m2x_over_H2_at={f'{x}': float(np.interp(x, rows[:, 0], mbar) / np.interp(x, rows[:, 0], rows[:, 2]) ** 2)
                                        for x in (0.3, 1.0, 3.0, 10.0, 30.0)})
        for Mt in (172.57, 172.86, 173.45, 174.05, 174.6):
            s = RN.smt(Mt)
            dd = M.danger(R, s, xi, Nob, m2xbar=mbar)
            row[f'{Mt}'] = dict(D_grav=dd['D'], note=dd['note'], tau_c=dd.get('tau_c'), hbar_over_hmax=dd.get('hbar_over_hmax'),
                                x=dd.get('x'), D_plasma=M.danger(R, s, xi, Nob, m2xbar=mbar, mp2='bg')['D'])
            sm = INF.SMfull(Mt)
            N0 = Nob + W['N_end']
            row[f'{Mt}']['inflation_log10N_tree'] = INF.log10N_const(sm, bg.Hi_GeV, 0.1, N0, False, m2x=Lam * bg.Hi_GeV ** 2)[0]
        out[k] = row
        print(hist, 'Lam', Lam, json.dumps({kk: (vv if not isinstance(vv, dict) else {a: (round(b, 2) if isinstance(b, float) else b)
                                                                                    for a, b in vv.items()})
                                           for kk, vv in row.items()}, default=float), flush=True)
        C.jdump(out, path)


# ------------------------------------------------------------------------------------------- upper edge with plasma
def plasma_barrier(smt, H_GeV, m2_Hi2):
    """canonical barrier position (H_i units) of V = lam_eff phi^4/4 + m^2 phi^2/2 (flat 2-loop lam_eff), or inf."""
    eG = math.exp(smt.G(H_GeV))
    hs = np.exp(np.linspace(math.log(smt.flat['h_max'] * 0.5), math.log(1e17), 1500))
    lam = smt.lam(hs)
    phi = eG * hs / H_GeV
    V = lam * phi ** 4 / 4 / eG ** 4 + 0.5 * m2_Hi2 * phi ** 2     # lam_eff includes e^{4 Gamma}: V(h0) = lam h0^4/4
    k = int(np.argmax(V))
    return float(phi[k]) if k < len(hs) - 1 else np.inf


def cmd_upperpl():
    out = {}
    for hist in ('EXP', 'PAULI0'):
        bg = RN.bgc('HFa', hist)
        o = {}
        for xi in (30.0, 100.0, 300.0, 1000.0):
            n = 320 if xi <= 100 else 400
            R = RN.run_cached('HFa', hist, xi, t0=-1.0, t_end=120.0, n=n)
            rows = R['rows']
            te = R['res']['tau_endinf']
            nk = np.maximum(np.array(R['res']['nk_end']), 0.0)
            res = {}
            for Mt in (172.57, 172.86, 173.45):
                s = RN.smt(Mt)
                hmc = math.exp(s.G(bg.Hi_GeV)) * s.flat['h_max'] / bg.Hi_GeV
                worst = np.inf
                tw = None
                for i in np.where(rows[:, 0] > te + 0.05)[0][::3]:
                    a, H, Rb, mp2 = rows[i, 1], rows[i, 2], rows[i, 4], rows[i, 8]
                    w = np.sqrt(R['ks'] ** 2 / a ** 2 + max(xi * Rb, 0.0) + H ** 2)
                    sp = float(np.sum(R['wq'] * nk / (a ** 3 * w)))
                    hb = plasma_barrier(s, bg.Hi_GeV, mp2 + max(xi * Rb, 0.0))
                    r = hb / (2 * math.sqrt(sp)) if sp > 0 else np.inf
                    if r < worst:
                        worst, tw = r, float(rows[i, 0])
                res[f'{Mt}'] = dict(min_hbarrier_over_2sigma=worst, at_t=tw, h_max_canon=hmc)
            o[f'{xi:g}'] = res
            print(hist, xi, res, flush=True)
        out[hist] = o
    C.jdump(out, C.RES / 'upper_plasma.json')


# ------------------------------------------------------------------------------------------- induced masses
def cmd_induced():
    """[E] sizes of loop-induced Higgs masses during inflation on the benchmark.
    (i) axial N_j: |m_N| = M_j constant -> the one-loop Higgs mass is field independent (absorbed in the EW tuning).
    (ii) two-loop phase mismatch between N1 (sigma phase) and N2-4 (chi phase): ~ Re[(y^dag y)_{1j}^2 e^{i dtheta}] M1 Mj
         /(16 pi^2)^2 (one power of log omitted), frozen during inflation, rapidly rotating after the waterfall.
    (iii) the spectator's N_s IF it is given the seesaw Yukawa the mf_hybrid skeptic says it needs (N8): field-dependent
         M_s(s) = M_soft + y_s s, so delta m_h^2(s) = -(y_nus^2/8 pi^2) M_s^2 [ln(M_s^2/mu^2) - 1] (inf_higgs skeptic N4
         form), with y_nus^2 = mt M_soft / v^2."""
    d = C.load_bench('HFa')
    p = d['params']
    H = d['H_i_GeV']
    Ms = [p['axial_sigma'][0]['M_over_Hi']] + [c['M_over_Hi'] for c in p['axial_chi']]
    ynu2 = [C.MNU * (M * H) / C.VEW ** 2 for M in Ms]
    oneloop_const = [y2 * M ** 2 / (8 * math.pi ** 2) for y2, M in zip(ynu2, Ms)]
    mism = [math.sqrt(ynu2[0] * ynu2[j]) ** 2 * Ms[0] * Ms[j] / (16 * math.pi ** 2) ** 2 for j in (1, 2, 3)]
    # spectator N_s
    ys = 4.6e-3
    sstar = p['s_star_MP'] * C.mP_REF
    MsI = ys * sstar
    rows = {}
    for Msoft in (0.3, 1.0, 1.5):
        y2 = C.MNU * Msoft * H / C.VEW ** 2
        cM2 = y2 * MsI ** 2 / (8 * math.pi ** 2)
        rows[f'M_soft={Msoft}'] = dict(y_nus2=y2, coeff_Hi2=cM2,
                                       mu_eq_Msoft=-cM2 * (math.log(MsI ** 2 / Msoft ** 2) - 1),
                                       mu_eq_MsI=cM2,
                                       mu_eq_MP=-cM2 * (math.log(MsI ** 2 / C.mP_REF ** 2) - 1),
                                       mu_min_for_tolerance_Hi=MsI * math.exp(-0.5 * (0.66 / cM2 + 1)))
    out = dict(axial_N=dict(M_over_Hi=Ms, y_nu2=ynu2, oneloop_constant_Hi2=oneloop_const,
                            note='field independent (|m_N| = M for every phi): absorbed in the electroweak tuning; no '
                                 'inflationary or post-inflation Higgs mass from it [D]'),
               phase_mismatch_twoloop_Hi2=mism,
               spectator_Ns=dict(M_s_inflation_over_Hi=MsI, rows=rows,
                                 tolerance_note='-0.66 H_i^2: the negative Higgs mass^2 during inflation that the '
                                                'constant-H HM count tolerates at xi_EW = 0.1, M_t = 172.86 (tree): '
                                                '12 [xi(H_i) - xi_min(H_i)] H^2 = 12 (0.0877 - 0.0326)'),
               label='E')
    C.jdump(out, C.RES / 'induced.json')
    print(json.dumps(out, indent=1, default=float))


def cmd_thermal():
    """reach if only a thermal bath counted: sharp onset at the pessimistic (LPM-type) thermalisation time against the
    gravitational-only deadline curve D(tau_th) = 0."""
    P = json.load(open(C.RES / 'plasma.json'))
    Dd = json.load(open(C.RES / 'post_D.json'))
    out = {}
    for hist, plab in (('EXP', 'UNBL'), ('UNBL', 'UNBL'), ('PAULI0', 'PAULI0'), ('PAULI10', 'PAULI10'), ('PAULIb', 'PAULIb')):
        o = Dd[f'HFa|{hist}']
        tth = P['HFa'][plab]['t_thermal_pessimistic']['t']
        Mts = sorted(float(m) for m in o['deadlines'])
        dl = [o['deadlines'][f'{m}']['deadline'] for m in Mts]
        reach = None
        for i in range(len(Mts) - 1):
            a, b = dl[i], dl[i + 1]
            a = 1e9 if a is None else a
            if a >= tth > b:
                reach = Mts[i] + (a - tth) * (Mts[i + 1] - Mts[i]) / (a - b) if a < 1e8 else Mts[i]
        out[hist] = dict(t_thermal=tth, deadlines=dict(zip(map(str, Mts), dl)), Mt_reach_thermal_only=reach)
        print(hist, out[hist], flush=True)
    C.jdump(out, C.RES / 'thermal_reach.json')


if __name__ == '__main__':
    c = sys.argv[1]
    if c == 'g6':
        cmd_g6(sys.argv[2] if len(sys.argv) > 2 else 'EXP')
    else:
        {'controls': cmd_controls, 'tail': cmd_tail, 'windows': cmd_windows, 'upperpl': cmd_upperpl,
         'induced': cmd_induced, 'thermal': cmd_thermal}[c]()
