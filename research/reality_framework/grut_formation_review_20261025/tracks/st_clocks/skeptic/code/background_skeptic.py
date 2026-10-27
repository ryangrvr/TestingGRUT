#!/usr/bin/env python3
"""st_clocks_skeptic (Claude, own code): independent radiation clock and source-clock fractions for Codex b73de64.

DIFFERENT METHOD from st_clocks (which extended our C engine rpreheat.c, DP5, full 19-variable system, then splined a(tau)
on a 0.05 grid):
* Only the gravitationally relevant sector: quadratic inflaton phi (M_r = H_i = 1 units), friction Gamma, radiation.
  The hidden X/Y/Sigma sector enters the Friedmann equation with eta^2 = (H_i/M_r)^2 = 3.0e-15 and shifts f = 1 - 2 xi
  eta^2 |X|^2 and rho_I(0) only at ~1e-11 relative [A]; that is the expected size of any C difference to the full models.
* scipy DOP853 with rtol 1e-13 (control 1e-11); state (ln a, phi, phidot, E_r = rho_r a^4) plus quadrature states
  J_g(tau) = int_0^tau g a(s) e^{-g s} ds for each width g and A(tau) = int_0^tau a ds (no interpolation at all).
* After tau_1 (inflaton energy < 1e-30 of radiation) the radiation era is exact: C = sqrt(E_r/3), a^2 = a_1^2 + 2C(tau-tau_1);
  every later number (a at medians, J beyond tau_1, total) is closed form: int_t^inf g sqrt(2C(s-ts)) e^{-gs} ds =
  sqrt(2C/g) e^{-g ts} Gamma(3/2, g (t - ts)), evaluated with mpmath at 50 digits and checked by mpmath quadrature.
* Also re-evaluates Codex's energy fraction from Codex's OWN saved a(tau) (npz data) with (i) PCHIP + Gauss (Codex's
  interpolant), (ii) a not-a-knot cubic spline on the same grid, to test st_clocks' claim that Codex's 2e-7 offset is
  PCHIP interpolation error.
Codex JSON/npz files are read as DATA only; no Codex or st_clocks code is imported.
Output: ../results/background_skeptic.json
"""
import json, math
from pathlib import Path
import numpy as np
import mpmath as mp
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator, CubicSpline
from scipy.optimize import brentq

mp.mp.dps = 50
HERE = Path(__file__).resolve().parent
OUT = HERE.parent / 'results' / 'background_skeptic.json'
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')
CODEX = SP / 'b73/research/reality_framework/grut_source_timing_20261026'
STC = HERE.parents[1] / 'results' / 'clocks.json'
W = json.loads((HERE.parent / 'results' / 'widths_skeptic.json').read_text())
GAM = {'YL': W['scalar']['theta_minus_pi_over_2']['gamma_L'], 'A_e1e-7': W['vector']['e1e-07']['gamma'],
       'A_e1e-6': W['vector']['e1e-06']['gamma']}
KEYS = list(GAM)
CASES = {'quiet3': (3.0, 0.5), 'quiet5': (5.0, 1.0)}


def rhs_factory(b, Gam):
    gs = np.array([GAM[k] for k in KEYS])

    def rhs(t, y):
        lna, phi, v, Er = y[0], y[1], y[2], y[3]
        a = math.exp(lna)
        rhoI = 0.5 * v * v + 0.5 * b * b * phi * phi
        H = math.sqrt((rhoI + Er / a ** 4) / 3.0)
        d = np.empty_like(y)
        d[0] = H; d[1] = v; d[2] = -(3 * H + Gam) * v - b * b * phi; d[3] = Gam * v * v * a ** 4
        d[4] = a                                  # A(tau) = int a ds
        d[5:] = gs * a * np.exp(-gs * t)          # J_g
        return d
    return rhs


def run(b, Gam, rtol, tau1, t_eval):
    y0 = np.zeros(5 + len(KEYS)); y0[1] = math.sqrt(6.0) / b          # rho_I(0) = 3 => H(0) = 1, phi at rest
    sol = solve_ivp(rhs_factory(b, Gam), (0.0, tau1), y0, method='DOP853', rtol=rtol, atol=1e-30, t_eval=t_eval,
                    dense_output=True)
    assert sol.success
    return sol


class Clock:
    def __init__(self, b, Gam, rtol=1e-13, tau1=160.0):
        self.grid = np.round(np.arange(0, tau1 + 1e-9, 0.05), 10)
        self.sol = run(b, Gam, rtol, tau1, self.grid)
        y1 = self.sol.y[:, -1]
        self.tau1 = tau1
        self.a1 = math.exp(y1[0])
        rhoI1 = 0.5 * y1[2] ** 2 + 0.5 * b * b * y1[1] ** 2
        self.inflaton_fraction_tau1 = rhoI1 * self.a1 ** 4 / y1[3]
        self.C = math.sqrt(y1[3] / 3.0)
        self.ts = tau1 - self.a1 ** 2 / (2 * self.C)
        self.J1 = {k: y1[5 + i] for i, k in enumerate(KEYS)}
        a = np.exp(self.sol.y[0]); lna = self.sol.y[0]
        rhoI = 0.5 * self.sol.y[2] ** 2 + 0.5 * b * b * self.sol.y[1] ** 2
        H = np.sqrt((rhoI + self.sol.y[3] / a ** 4) / 3)
        self.Cser = a * a * H
        self.a_grid = a

    def a_of(self, t):
        if t <= self.tau1:
            return float(math.exp(self.sol.sol(t)[0]))
        return math.sqrt(2 * self.C * (t - self.ts))

    def tail(self, k, t):           # int_t^inf g a e^{-g s} ds in the exact radiation era (t >= tau1)
        g = mp.mpf(GAM[k]); C = mp.mpf(self.C); ts = mp.mpf(self.ts)
        return mp.sqrt(2 * C / g) * mp.e ** (-g * ts) * mp.gammainc(mp.mpf(3) / 2, g * (mp.mpf(t) - ts), mp.inf)

    def tail_quad(self, k, t):
        g = mp.mpf(GAM[k]); C = mp.mpf(self.C); ts = mp.mpf(self.ts); t = mp.mpf(t)
        f = lambda u: mp.sqrt(2 * C * (t - ts + u / g)) * mp.e ** (-g * t - u)      # s = t + u/g
        return mp.quad(f, [0, 1, 5, 20, 60, mp.inf])

    def total(self, k):
        return mp.mpf(self.J1[k]) + self.tail(k, self.tau1)

    def J(self, k, t):
        if t <= self.tau1:
            return mp.mpf(float(self.sol.sol(t)[5 + KEYS.index(k)]))
        return mp.mpf(self.J1[k]) + self.tail(k, self.tau1) - self.tail(k, t)

    def F_E(self, k, t):
        return float(self.J(k, t) / self.total(k))


def codex_interp_check(tag, k, t, total_ref):
    """Codex's own a(tau) data, two interpolants, high-order Gauss on each knot interval (own quadrature)."""
    z = np.load(CODEX / (tag + '.npz'))
    times = z['times']; a = np.exp(z['bg'][:, 15]); g = GAM[k]
    xg, wg = np.polynomial.legendre.leggauss(12)
    out = {}
    for name, f in [('pchip', PchipInterpolator(times, a)), ('cubic_spline', CubicSpline(times, a))]:
        edges = np.r_[times[times < t], t]
        half = np.diff(edges) / 2; mid = (edges[1:] + edges[:-1]) / 2
        tt = mid[:, None] + half[:, None] * xg[None, :]
        val = float(np.sum(half[:, None] * wg[None, :] * g * f(tt) * np.exp(-g * tt)))
        out[name] = val
    out['max_codex_step_below_t'] = float(np.max(np.diff(times[times <= t])))
    return out


def main():
    cx = json.loads((CODEX / 'SOURCE_CLOCKS.json').read_text())
    cxh = {h['tag']: h for h in cx['histories']}
    cr = {r['tag']: r for r in json.loads((CODEX / 'RESPONSE_RESULTS.json').read_text())['rows']}
    st = json.loads(STC.read_text())['histories']
    res = dict(gammas=GAM, histories={})
    for tag, (b, Gam) in CASES.items():
        ck = Clock(b, Gam, 1e-13)
        ck2 = Clock(b, Gam, 1e-11)
        h = dict(b=b, Gamma_I=Gam, tau1=ck.tau1, inflaton_fraction_at_tau1=ck.inflaton_fraction_tau1,
                 C_own=ck.C, C_own_rtol1e_11=ck2.C, C_codex=cxh[tag]['late_C'], C_st_clocks=st[tag]['C_own'],
                 C_rel_vs_codex=ck.C / cxh[tag]['late_C'] - 1, C_rel_vs_st_clocks=ck.C / st[tag]['C_own'] - 1,
                 C_codex_gauge_run=cr['gauge' + tag[-1]]['late_C'],
                 tau_star_own=ck.ts, tau_star_codex=cxh[tag]['matched_radiation_tau_offset'],
                 tau_star_st_clocks=st[tag]['tau_star_matched_own'],
                 a800_own=ck.a_of(800.0), a800_codex=cr[tag]['final_a'], a800_st_clocks=st[tag]['a800_own'])
        h['a800_rel_vs_codex'] = h['a800_own'] / h['a800_codex'] - 1
        h['a800_rel_vs_st_clocks'] = h['a800_own'] / h['a800_st_clocks'] - 1
        rel = np.abs(ck.Cser / ck.C - 1)
        for thr in (1e-3, 1e-6):
            bad = np.flatnonzero(rel > thr)
            h['tau_after_which_C_within_%g' % thr] = float(ck.grid[bad[-1] + 1]) if len(bad) else 0.0
        h['C_at_20'] = float(ck.Cser[np.argmin(abs(ck.grid - 20))])
        h['C_drift_64_160'] = float(np.ptp(ck.Cser[ck.grid >= 64]) / ck.C)
        # gamma_cascade = H crossing
        gc = W['scalar']['theta_minus_pi_over_2']['gamma_cascade']
        Hgrid = ck.Cser / ck.a_grid ** 2
        i = int(np.flatnonzero(Hgrid <= gc)[0])
        h['cascade'] = dict(tau_Gamma_eq_H=float(ck.grid[i]), a_Gamma_eq_H=float(ck.a_grid[i]),
                            done_by_134=-math.expm1(-gc * 134.0), done_by_147_5=-math.expm1(-gc * 147.5),
                            tau_half=math.log(2) / gc)
        h['sources'] = {}
        for k in KEYS:
            g = GAM[k]
            tot = ck.total(k)
            s = dict(mean_a_own=float(tot), mean_a_codex=cxh[tag]['sources'][k]['mean_a_per_decay'],
                     mean_a_rel=float(tot) / cxh[tag]['sources'][k]['mean_a_per_decay'] - 1,
                     late_limit_rel=math.sqrt(math.pi * ck.C / (2 * g)) / float(tot) - 1,
                     tail_closed_vs_quad=float(ck.tail(k, ck.tau1) / ck.tail_quad(k, ck.tau1) - 1),
                     tail_closed_vs_quad_800=float(ck.tail(k, 800) / ck.tail_quad(k, 800) - 1),
                     total_rtol_change=float(ck2.total(k) / tot - 1))
            s['number_quantiles'] = []
            for fr, cq in zip([.01, .1, .5, .9, .99], cxh[tag]['sources'][k]['number_quantiles']):
                t = -math.log1p(-fr) / g
                s['number_quantiles'].append(dict(fraction=fr, tau=t, a_own=ck.a_of(t), a_codex=cq['a'],
                                                  a_rel=ck.a_of(t) / cq['a'] - 1))
            s['energy_quantiles'] = []
            for fr, cq in zip([.1, .5, .9], cxh[tag]['sources'][k]['energy_quantiles']):
                # all energy quantiles lie in the exact radiation era: solve tail(t)/total = 1 - fr in u = g (t - ts)
                fn = lambda u: float(ck.tail(k, ck.ts + u / g) / tot) - (1 - fr)
                u = brentq(fn, max(1e-12, g * (ck.tau1 - ck.ts)), 60, xtol=1e-15, rtol=1e-15)
                t = ck.ts + u / g
                s['energy_quantiles'].append(dict(fraction=fr, tau=t, tau_codex=cq['tau'], tau_rel=t / cq['tau'] - 1,
                                                  a_own=ck.a_of(t), a_codex=cq['a']))
            s['at_crossings'] = []
            for c in cxh[tag]['sources'][k]['at_response_crossings']:
                t = c['tau']
                fe = ck.F_E(k, t)
                row = dict(seed=c['phase_seed'], kind=c['kind'], tau=t, a_own=ck.a_of(t), a_codex=c['a'],
                           number_fraction_own=-math.expm1(-g * t),
                           number_fraction_codex=c['decayed_number_fraction_rest_clock'],
                           energy_fraction_own=fe,
                           energy_fraction_codex=c['cumulative_conformal_energy_fraction_prescribed_source'],
                           energy_rel_vs_codex=fe / c['cumulative_conformal_energy_fraction_prescribed_source'] - 1)
                if k == 'YL':
                    ci = codex_interp_check(tag, k, t, float(tot))
                    row['codex_data_pchip_over_total_own'] = ci['pchip'] / float(tot)
                    row['codex_data_spline_over_total_own'] = ci['cubic_spline'] / float(tot)
                    row['codex_pchip_rel_vs_own'] = ci['pchip'] / float(ck.J(k, t)) - 1 if t <= ck.tau1 else None
                    row['codex_spline_rel_vs_own'] = ci['cubic_spline'] / float(ck.J(k, t)) - 1 if t <= ck.tau1 else None
                    row['codex_max_step'] = ci['max_codex_step_below_t']
                s['at_crossings'].append(row)
            h['sources'][k] = s
        res['histories'][tag] = h
        print(tag, json.dumps({kk: vv for kk, vv in h.items() if kk != 'sources'}), flush=True)
        for k in KEYS:
            s = h['sources'][k]
            print('  ', k, 'mean_a', s['mean_a_own'], s['mean_a_rel'], 'median a', s['number_quantiles'][2]['a_own'],
                  s['number_quantiles'][2]['a_rel'])
            for r in s['at_crossings'][:3]:
                print('     ', {kk: vv for kk, vv in r.items() if kk not in ('seed',)})
    OUT.write_text(json.dumps(res, indent=2) + '\n')


if __name__ == '__main__':
    main()
