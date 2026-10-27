#!/usr/bin/env python3
"""st_clocks (Claude, own code): independent decay-clock and kinematics checks of Codex b73de64
(grut_source_timing_20261026), issue #13.

Independent: widths re-derived from the couplings (own formulas, own three-body phase-space quadrature), background
a(tau) from OUR reheating engine (background800.py, rp_reheating/rpreheat.c), own spline quadrature, own tail
integral (mpmath incomplete gamma AND direct mpmath quadrature), own N_X.  Codex numbers are read from its
SOURCE_CLOCKS.json / RESPONSE_RESULTS.json as DATA ONLY (comparison targets, crossing tau samples); no Codex code is
imported or run.

Usage: python3 -B clocks.py   (after background800.py)
"""
import json, math
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.integrate import dblquad

HERE = Path(__file__).resolve().parent
TR = HERE.parents[1]
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')
CODEX = SP / 'b73/research/reality_framework/grut_source_timing_20261026'
RUNS = SP / 'review/st_clocks/rp/runs'
mp.mp.dps = 40

# ---------------------------------------------------------------- couplings (spot point, H_i = 1 units) [Q own tracks]
SRCI = json.loads((TR / 'kin_source/results/source_interface.json').read_text())
FP = SRCI['meta']['fixed_point']
H_GEV = FP['H_GeV']                                   # 133349619398.33566 GeV
KAP, MY, HH, MU, LS = 0.0034586847211042523, 30.0, 40.0, 7.159456159513794, 1.0   # rp_reheating/rp.py P, hw_gates SPOT
MX_GEV = FP['mX_GeV']
COLS = ['tau', 'a', 'H', 'R', 'NX', 'NY', 'NS', 'Gauss', 'EX', 'eps', 'NYH', 'NYL', 'Nsig', 'absX', 'absY', 'S_over_vev',
        'rho_r', 'rho_I', 'radfrac', 'Hcheck', 'Wcheck', 'w_I', 'm2curv', 'KV', 'Haux', 'rho']
CI = {k: i for i, k in enumerate(COLS)}


def phase_space_3body_massless(m, mx=0.0):
    """Lorentz-invariant 3-body phase space Phi_3 = 1/(128 pi^3 m^2) * Dalitz area (own quadrature)."""
    def lim(s12):
        e2 = s12 / (2 * math.sqrt(s12)); e3 = (m * m - s12) / (2 * math.sqrt(s12))
        p2 = math.sqrt(max(e2 * e2 - mx * mx, 0)); p3 = math.sqrt(max(e3 * e3 - mx * mx, 0))
        return (e2 + e3) ** 2 - (p2 + p3) ** 2, (e2 + e3) ** 2 - (p2 - p3) ** 2
    area, _ = dblquad(lambda s23, s12: 1.0, (2 * mx) ** 2, (m - mx) ** 2, lambda s: lim(s)[0], lambda s: lim(s)[1],
                      epsabs=0, epsrel=1e-12)
    return area / (128 * math.pi ** 3 * m * m)


def widths():
    mL = math.sqrt(MY ** 2 - math.sqrt(2) * MU * HH); mH = math.sqrt(MY ** 2 + math.sqrt(2) * MU * HH)
    # Y = (Y1 + i Y2)/sqrt2; 2Re(kap Y* X^3) gives each real eigenmode (kap/sqrt2)(X^3 + X*^3) (up to a phase):
    # |M|^2 = (3! kap/sqrt2)^2 = 18 kap^2 per channel (XXX and XbarXbarXbar), identical-particle factor 1/3!.
    phi3 = phase_space_3body_massless(mL)
    G_ch = 18 * KAP ** 2 * phi3 / (2 * mL) / 6.0
    G_L_own = 2 * G_ch
    G_L_formula = 3 * KAP ** 2 * mL / (256 * math.pi ** 3)
    # vector: m_A = q_Sigma e h with q_Sigma = 6 (Sigma = h/sqrt2), Gamma = e^2 m_A beta^3/(48 pi) for charge-1 complex X
    vec = {}
    for e in (1e-7, 1e-6):
        mA = 6 * e * HH
        mx = MX_GEV / H_GEV
        beta = math.sqrt(1 - 4 * mx * mx / (mA * mA))
        # own: |M|^2 summed over pols / 3 = e^2 (4 p*^2)/3 ... Gamma = p*/(8 pi mA^2) * (4/3) e^2 p*^2
        pstar = mA * beta / 2
        G_own = pstar / (8 * math.pi * mA ** 2) * (4.0 / 3.0) * e * e * pstar ** 2
        vec['e%g' % e] = dict(mA_over_Hi=mA, mA_GeV=mA * H_GEV, beta=beta, gamma_own=G_own,
                              gamma_formula=e * e * mA / (48 * math.pi) * beta ** 3, Gamma_GeV=G_own * H_GEV,
                              e_cubed_law=e ** 3 * HH / (8 * math.pi))
    # prescribed cascade and Sigma radial decays (Goldstone equivalence, E_A >> m_A)
    G_casc = MU ** 2 * (mH ** 2 - mL ** 2) / (8 * math.pi * mH ** 3)        # coupling sqrt2 mu G Y1 Y2
    msig = math.sqrt(2 * LS) * HH
    G_sig_AA = LS ** 2 * HH ** 2 / (8 * math.pi * msig)                      # coupling lambda h sigma G^2
    p_sYY = math.sqrt(msig ** 2 / 4 - mL ** 2)
    G_sig_YLYL = MU ** 2 * p_sYY / (8 * math.pi * msig ** 2)                 # coupling (mu/sqrt2) sigma Y_L^2
    return dict(mL=mL, mH=mH, phi3_numeric_over_formula=phi3 / (mL ** 2 / (256 * math.pi ** 3)),
                gamma_L_own=G_L_own, gamma_L_formula=G_L_formula, Gamma_L_GeV=G_L_own * H_GEV, vectors=vec,
                gamma_cascade=G_casc, Gamma_cascade_GeV=G_casc * H_GEV, E_A_cascade=(mH ** 2 - mL ** 2) / (2 * mH),
                m_sigma=msig, gamma_sigma_AA=G_sig_AA, gamma_sigma_YLYL=G_sig_YLYL)


# ---------------------------------------------------------------- background clock
class Background:
    def __init__(self, path):
        S = np.load(path)
        self.S = S
        self.tau = S[:, CI['tau']]; self.a = S[:, CI['a']]; self.H = S[:, CI['H']]
        self.C_series = self.a ** 2 * self.H
        self.tm = self.tau[-1]
        self.C = float(self.C_series[-1])
        self.ts = float(self.tm - self.a[-1] ** 2 / (2 * self.C))           # matched offset tau_*
        sel = self.tau >= 64
        A = np.vstack([self.tau[sel], np.ones(sel.sum())]).T
        slope, icpt = np.linalg.lstsq(A, self.a[sel] ** 2, rcond=None)[0]
        self.C_fit, self.ts_fit = slope / 2, -icpt / slope
        resid = self.a[sel] ** 2 - 2 * self.C * (self.tau[sel] - self.ts)
        self.match_resid = float(np.max(np.abs(resid / self.a[sel] ** 2)))
        self.spl_a = CubicSpline(self.tau, self.a)

    def a_of(self, t):
        return float(self.spl_a(t)) if t <= self.tm else math.sqrt(2 * self.C * (t - self.ts))

    def early(self, gam, t):
        """int_0^t gam a(s) e^{-gam s} ds on our own 0.05 grid: cubic spline of the integrand, exact integral."""
        g = gam * self.a * np.exp(-gam * self.tau)
        return float(CubicSpline(self.tau, g).integrate(0.0, min(t, self.tm)))

    def early_simpson(self, gam, t):
        n = int(round(t / 0.05)); assert abs(n * 0.05 - t) < 1e-9 and n % 2 == 0
        g = gam * self.a[:n + 1] * np.exp(-gam * self.tau[:n + 1])
        return float(0.05 / 3 * (g[0] + g[-1] + 4 * g[1:-1:2].sum() + 2 * g[2:-1:2].sum()))

    def tail_closed(self, gam, t):
        u = mp.mpf(gam) * (mp.mpf(t) - mp.mpf(self.ts))
        return float(mp.sqrt(2 * mp.mpf(self.C) / gam) * mp.e ** (-mp.mpf(gam) * self.ts) * mp.gammainc(1.5, u, mp.inf))

    def tail_quad(self, gam, t):
        C, ts, g = mp.mpf(self.C), mp.mpf(self.ts), mp.mpf(gam)
        f = lambda s: g * mp.sqrt(2 * C * (s - ts)) * mp.e ** (-g * s)
        t = mp.mpf(t)
        return float(mp.quad(f, [t, t + 0.1 / g, t + 1 / g, t + 5 / g, t + 30 / g, mp.inf]))

    def total(self, gam):
        return self.early(gam, self.tm) + self.tail_closed(gam, self.tm)

    def energy_fraction(self, gam, t, tot=None):
        tot = self.total(gam) if tot is None else tot
        if t <= self.tm:
            return self.early(gam, t) / tot
        return (tot - self.tail_closed(gam, t)) / tot


def main():
    W = widths()
    gam = {'YL': W['gamma_L_own'], 'A_e1e-7': W['vectors']['e1e-07']['gamma_own'],
           'A_e1e-6': W['vectors']['e1e-06']['gamma_own']}
    # source-dictionary values (our kin_source) and Codex values (data)
    src = {'YL': SRCI['mediator_debris']['YL_three_body']['parents']['Gamma_L_GeV'],
           'A_e1e-7': SRCI['vector_sources']['e1e-7_band']['parent']['Gamma_A_GeV'],
           'A_e1e-6': SRCI['vector_sources']['e1e-6_band']['parent']['Gamma_A_GeV']}
    cx = json.loads((CODEX / 'SOURCE_CLOCKS.json').read_text())
    cxh = {h['tag']: h for h in cx['histories']}
    cr = {r['tag']: r for r in json.loads((CODEX / 'RESPONSE_RESULTS.json').read_text())['rows']}
    out = dict(widths=W, width_checks={}, medians={}, histories={})
    for k in gam:
        out['width_checks'][k] = dict(Gamma_GeV_own=gam[k] * H_GEV, Gamma_GeV_source_dictionary=src[k],
                                      rel_own_vs_dictionary=gam[k] * H_GEV / src[k] - 1,
                                      Gamma_GeV_codex=cxh['quiet3']['sources'][k]['Gamma_GeV'],
                                      rel_own_vs_codex=gam[k] * H_GEV / cxh['quiet3']['sources'][k]['Gamma_GeV'] - 1)
        out['medians'][k] = dict(tau_half_own=math.log(2) / gam[k],
                                 tau_half_codex=cxh['quiet3']['sources'][k]['number_quantiles'][2]['tau'])
    out['vector_median_ratio_e1e-7_over_e1e-6'] = gam['A_e1e-6'] / gam['A_e1e-7']
    for tag, f in [('quiet3', 'q3_nominal_800.npy'), ('quiet5', 'q5_nominal_800.npy')]:
        bg = Background(RUNS / f)
        bgt = Background(RUNS / f.replace('nominal', 'tight'))
        h = dict(C_own=bg.C, C_codex=cxh[tag]['late_C'], C_rel=bg.C / cxh[tag]['late_C'] - 1,
                 C_relative_range_64_800=float((bg.C_series[bg.tau >= 64].max() - bg.C_series[bg.tau >= 64].min()) / bg.C),
                 C_series_samples={str(t): float(bg.C_series[np.argmin(abs(bg.tau - t))]) for t in (20, 40, 64, 80, 134, 300, 500, 800)},
                 C_fit_64_800=float(bg.C_fit), tau_star_fit=float(bg.ts_fit),
                 tau_star_matched_own=bg.ts, tau_star_codex=cxh[tag]['matched_radiation_tau_offset'],
                 a2_linear_max_rel_residual_64_800=bg.match_resid,
                 C_tight=bgt.C, a800_own=float(bg.a[-1]), a800_codex=cr[tag]['final_a'], a800_tight=float(bgt.a[-1]))
        # where does C settle (inflaton gone)?
        rel = np.abs(bg.C_series / bg.C - 1)
        bad = np.flatnonzero(rel > 1e-6)
        h['tau_after_which_C_within_1e-6'] = float(bg.tau[bad[-1] + 1]) if len(bad) else 0.0
        bad = np.flatnonzero(rel > 1e-3)
        h['tau_after_which_C_within_1e-3'] = float(bg.tau[bad[-1] + 1]) if len(bad) else 0.0
        h['sources'] = {}
        for k, g in gam.items():
            tot = bg.total(g); tot_t = bgt.total(g)
            tail = bg.tail_closed(g, bg.tm); tailq = bg.tail_quad(g, bg.tm)
            s = dict(gamma=g, mean_a_total=tot, mean_a_codex=cxh[tag]['sources'][k]['mean_a_per_decay'],
                     mean_a_rel=tot / cxh[tag]['sources'][k]['mean_a_per_decay'] - 1,
                     late_limit_sqrt_piC_over_2gamma=math.sqrt(math.pi * bg.C / (2 * g)),
                     late_limit_rel_to_total=math.sqrt(math.pi * bg.C / (2 * g)) / tot - 1,
                     tail_closed=tail, tail_quad=tailq, tail_closed_vs_quad=tail / tailq - 1,
                     tight_total_rel=tot_t / tot - 1)
            nq = []
            for fr, cq in zip([.01, .1, .5, .9, .99], cxh[tag]['sources'][k]['number_quantiles']):
                t = -math.log1p(-fr) / g
                nq.append(dict(fraction=fr, tau=t, a_own=bg.a_of(t), a_codex=cq['a'], a_rel=bg.a_of(t) / cq['a'] - 1))
            s['number_quantiles'] = nq
            eq = []
            for fr, cq in zip([.1, .5, .9], cxh[tag]['sources'][k]['energy_quantiles']):
                # F_E(t) = 1 - tail(t)/tot for t > 800; solve in u = g (t - ts)
                u = brentq(lambda u: 1 - bg.tail_closed(g, bg.ts + u / g) / tot - fr, 1e-6, 40, xtol=1e-14)
                t = bg.ts + u / g
                eq.append(dict(fraction=fr, tau=t, tau_codex=cq['tau'], tau_rel=t / cq['tau'] - 1, a_own=bg.a_of(t), a_codex=cq['a']))
            s['energy_quantiles'] = eq
            xs = []
            for c in cxh[tag]['sources'][k]['at_response_crossings']:
                t = c['tau']
                nf = -math.expm1(-g * t)
                ef = bg.energy_fraction(g, t, tot)
                row = dict(seed=c['phase_seed'], kind=c['kind'], tau=t, a_own=bg.a_of(t), a_codex=c['a'],
                           number_fraction_own=nf, number_fraction_codex=c['decayed_number_fraction_rest_clock'],
                           energy_fraction_own=ef,
                           energy_fraction_codex=c['cumulative_conformal_energy_fraction_prescribed_source'],
                           energy_fraction_rel=ef / c['cumulative_conformal_energy_fraction_prescribed_source'] - 1,
                           energy_fraction_tight=bgt.energy_fraction(g, t, tot_t))
                if abs(t / 0.05 - round(t / 0.05)) < 1e-9 and round(t / 0.05) % 2 == 0:
                    row['energy_fraction_simpson'] = bg.early_simpson(g, t) / tot
                xs.append(row)
            s['at_crossings'] = xs
            h['sources'][k] = s
        # kinematic wavelength cube from OUR N_X
        NX800 = float(bg.S[-1, CI['NX']])
        rows = {r['tag']: r for r in json.loads((TR / 'rp_reheating/results/extra.json').read_text())['rows']}
        NX80 = rows['yr_b3_g0.5' if tag == 'quiet3' else 'yr_b5_g1']['late']['NX']     # tau-weighted mean 64-80
        h['NX_own_mean_64_80'] = NX80; h['NX_codex_norm_64_80'] = cr[tag]['late_NX_normalization_through80']
        h['NX_own_800'] = NX800; h['NX_codex_800'] = cr[tag]['final_homogeneous_NX']
        h['wavelength_cube'] = {}
        for kk in (2.13, 5.8):
            cube = (2 * math.pi / kk) ** 3
            h['wavelength_cube'][str(kk)] = dict(Q_full_wavelength_own_64_80=abs(NX80) * cube, Q_full_wavelength_own_800=abs(NX800) * cube,
                                                 Q_codex=[w for w in cxh[tag]['wavelength_diagnostics'] if w['k_over_Hi'] == kk][0]['Q_in_full_wavelength_cube'],
                                                 Q_half_wavelength_cube=abs(NX80) * cube / 8,
                                                 Q_sphere_diameter_lambda=abs(NX80) * cube * math.pi / 6)
        # other prescribed sources on this background
        gc = W['gamma_cascade']
        t_half_c = math.log(2) / gc
        # Gamma_casc = H on our background
        iH = np.flatnonzero(bg.H <= gc)
        h['cascade_YH_to_YL_AL'] = dict(gamma=gc, tau_half=t_half_c, a_at_tau_half=bg.a_of(t_half_c),
                                         tau_99=math.log(100) / gc, a_at_tau_99=bg.a_of(math.log(100) / gc),
                                         tau_where_Gamma_eq_H=float(bg.tau[iH[0]]), a_where_Gamma_eq_H=float(bg.a[iH[0]]),
                                         fraction_done_by_134=-math.expm1(-gc * 134.0), fraction_done_by_147_5=-math.expm1(-gc * 147.5),
                                         NYH_late_own=float(bg.S[-1, CI['NYH']]), NYH_dictionary=FP['N_YH'],
                                         NYL_late_own=float(np.mean(bg.S[bg.tau >= 640, CI['NYL']])), NYL_dictionary=FP['N_YL'],
                                         EX_late_own=float(np.mean(bg.S[bg.tau >= 640, CI['EX']])), EX_dictionary=FP['E_X'])
        # [E] condensate-assisted two-body Y_L -> XX / XbarXbar (not in the dictionary): Gamma_2 = 9 kap^2 |X|^2/(8 pi m_L)
        c2 = 9 * KAP ** 2 / (8 * math.pi * W['mL'])
        X2 = bg.S[:, CI['absX']] ** 2
        cumX2 = CubicSpline(bg.tau, X2).antiderivative()
        h['condensate_assisted_YL_E'] = dict(coefficient=c2, gamma2_at_tau80=c2 * float(np.mean(X2[(bg.tau > 70) & (bg.tau < 90)])),
                                              gamma2_at_tau134=c2 * float(np.mean(X2[(bg.tau > 124) & (bg.tau < 144)])),
                                              integrated_fraction_by_134=c2 * float(cumX2(134.0) - cumX2(1.0)),
                                              integrated_fraction_by_800=c2 * float(cumX2(800.0) - cumX2(1.0)),
                                              three_body_fraction_by_134=gam['YL'] * 134.0,
                                              scope='[E] tree-level, condensate leg |X| from our homogeneous run, massless X, no Bose/thermal factors; from tau=1')
        # [C] what the homogeneous classical Y counts of THIS history would give (Codex does not transfer populations;
        # quantum Y production is not in a homogeneous run, so this is a diagnostic, not a source normalisation)
        cc = h['cascade_YH_to_YL_AL']
        NL_q = cc['NYL_late_own'] + cc['NYH_late_own']
        mA = h['sources']['YL']['mean_a_total']
        h['YL_energy_on_this_history_C'] = dict(
            N_L_homogeneous=NL_q, R_Y_homogeneous=NL_q * W['mL'] * mA / cc['EX_late_own'],
            R_Y_dictionary_check=(FP['N_YH'] + FP['N_YL']) * W['mL'] * SRCI['mediator_debris']['YL_three_body']['parents']['a_mean_decay'] / FP['E_X'],
            injected_over_EX_by_unit_sample=NL_q * W['mL'] * mA * h['sources']['YL']['at_crossings'][1]['energy_fraction_own'] / cc['EX_late_own'],
            injected_over_EX_by_unit_sample_dictionary_population=(FP['N_YH'] + FP['N_YL']) * W['mL'] * mA *
            h['sources']['YL']['at_crossings'][1]['energy_fraction_own'] / FP['E_X'])
        out['histories'][tag] = h
    (HERE.parent / 'results' / 'clocks.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=1)[:200])


if __name__ == '__main__':
    main()
