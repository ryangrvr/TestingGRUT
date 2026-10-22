"""Bubble nucleation and the 'final vacuum' (track `vacua`, NOTES.md section 3).

Past-light-cone 4-volume V4(t) = (4 pi/3) int_0^eta a^4 (eta - eta')^3 deta'  (flat FRW, conformal time).
Probability that no bubble has nucleated in the past light cone of a point: exp(-Gamma V4).
Comoving false-vacuum fraction in de Sitter (bubble walls at c): f(t) = exp(-(4 pi/3) Gamma t / H^3);
physical false-vacuum volume ~ a^3 f = exp[(3 - (4 pi/3) Gamma/H^4) H t]  ->  grows iff Gamma/H^4 < 9/(4 pi).
Run: python3 nucleation.py  -> ../results/nucleation.json
"""
import json, os
import numpy as np
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results', 'nucleation.json')

# units and inputs
GEV_INV_S = 6.582119569e-25            # 1 GeV^-1 in s
GYR = 3.15576e16 / GEV_INV_S           # 1 Gyr in GeV^-1
GPC = 3.0857e25 * 5.067731e15          # 1 Gpc in GeV^-1
H0 = 67.4 / 3.0857e19 * GEV_INV_S      # GeV (Planck 2018 H0 = 67.4 km/s/Mpc) [R, as in the brief's inputs]
OM, OR = 0.315, 9.0e-5                 # [R]
OL = 1 - OM - OR
M_RED = 2.435e18
K = M_RED ** 2
RHO_L = 2.52e-47                       # brief


def Hn(N):
    a = np.exp(N)
    return np.sqrt(OR / a ** 4 + OM / a ** 3 + OL)


def V4_at_N(N_end, N_start=np.log(1e-9), n=400001):
    """LambdaCDM, units H0 = 1. V4 = (4 pi/3) int dN a^3 r(N)^3 / H, r(N) = int_N^{N_end} dN'/(a H)."""
    from scipy.integrate import cumulative_trapezoid
    N = np.linspace(N_start, N_end, n)
    a = np.exp(N)
    H = Hn(N)
    g = 1.0 / (a * H)
    # r(N) = int_N^{N_end} g dN', accumulated from the end (no cancellation at late times)
    r = -cumulative_trapezoid(g[::-1], N[::-1], initial=0.0)[::-1]
    integrand = a ** 3 * r ** 3 / H
    return 4 * np.pi / 3 * np.trapezoid(integrand, N)


def t_of_N(N_end, N_start=np.log(1e-9), n=400001):
    N = np.linspace(N_start, N_end, n)
    return np.trapezoid(1.0 / Hn(N), N)


def main():
    out = {'inputs': {'H0_GeV': H0, 'Omega_m': OM, 'Omega_r': OR, 'Omega_L': OL, 'rho_Lambda_GeV4': RHO_L,
                      'M_reduced_GeV': M_RED}}
    t0 = t_of_N(0.0)
    out['t0_H0units'] = t0
    out['t0_Gyr'] = t0 / H0 / GYR
    out['V4_past_light_cone_today_H0^-4'] = V4_at_N(0.0)
    out['Buttazzo_quoted'] = '0.15/H0^4 (1307.3536 eq. tunrate)'
    HL = np.sqrt(OL)
    NA, NB = 30.0, 40.0
    slope = (V4_at_N(NB) - V4_at_N(NA)) / (t_of_N(NB) - t_of_N(NA))
    out['dV4/dt_late_H0units'] = slope
    out['expected_(4pi/3)/H_L^3'] = 4 * np.pi / 3 / HL ** 3
    # analytic forms
    out['analytic'] = {
        'Minkowski_from_t=0': 'V4 = pi t^4 / 3',
        'matter_dominated': 'V4 = (16 pi/1485) H^-4 = (81 pi/1485) t^4 = 0.1714 t^4',
        'de_Sitter_late': 'V4 -> (4 pi/(3 H^4)) (H t - 11/6) (pure dS from a = 1; 1307.3536)',
    }
    # numerical check of the MD coefficient
    md = 4 * np.pi / 3 * quad(lambda e: e ** 8 * (1 - e) ** 3, 0, 1)[0]
    # with a = eta^2 (units where eta0 = 1 and t = eta^3/3): V4 = (4pi/3) int eta^8 (1-eta)^3 = (4pi/3)(1/1980)
    t_md = 1.0 / 3.0
    out['MD_check_V4/t^4'] = md / t_md ** 4
    out['MD_expected_81pi/1485'] = 81 * np.pi / 1485

    # SM decay rates (quoted) and consequences
    rates = {'Chigusa-Moroi-Shoji_best': -577, 'CMS_fast_end': -295, 'CMS_slow_end': -1465}
    G_unit = 1.0 / (GYR * GPC ** 3)        # 1 Gyr^-1 Gpc^-3 in GeV^4
    rr = {}
    V4_today = out['V4_past_light_cone_today_H0^-4'] / H0 ** 4
    HLg = H0 * HL
    for k, lg in rates.items():
        lgG = lg + np.log10(G_unit)                   # log10 Gamma[GeV^4]
        lgP = lgG + np.log10(V4_today)
        lg_eps = lgG - 4 * np.log10(HLg)              # Gamma / H_L^4
        # expected number of bubbles in the past cone N(t) = Gamma (4 pi/3) t / H^3 = 1  ->  t = 3 H^3/(4 pi Gamma)
        lg_tdec_dS = np.log10(3.0 / (4 * np.pi)) + 3 * np.log10(HLg) - lgG     # GeV^-1
        lg_tdec_MD = 0.25 * (np.log10(1485 / (81 * np.pi)) - lgG)           # GeV^-1
        lg_tdec_Mink = 0.25 * (np.log10(3 / np.pi) - lgG)
        rr[k] = {'log10_Gamma_GeV4': lgG, 'log10_P_decay_in_past_light_cone': lgP,
                 'log10_Gamma/H_L^4': lg_eps, 'eternal_false_vacuum_volume (Gamma/H^4 < 9/4pi)': bool(lg_eps < np.log10(9 / (4 * np.pi))),
                 'log10_t_worldline_decay_dS_Gyr': lg_tdec_dS - np.log10(GYR),
                 'log10_t_decay_if_asymptotically_MD_Gyr': lg_tdec_MD - np.log10(GYR),
                 'log10_t_decay_Minkowski_Gyr': lg_tdec_Mink - np.log10(GYR)}
    out['SM_rates'] = rr
    out['rate_sources'] = ('Chigusa, Moroi, Shoji, PRL 119 (2017) 211801, arXiv:1707.09301 (abstract): '
                           '~1e-577 Gyr^-1 Gpc^-3, range 1e-295 .. 1e-1465 from Mt; Andreassen, Frost, Schwartz, '
                           'PRD 97 (2018) 056006, arXiv:1707.08124 (abstract): lifetime 1e139 yr, >1e58 yr at 95%')
    # AdS crunch time inside a bubble with Lambda_in < 0 (pure negative Lambda, open slicing)
    crunch = {}
    for lab, L in (('|V(2.4e18 GeV)| at PDG central', 7.48e70), ('|V(1e17 GeV)|', 2.03e65),
                   ('rho_Lambda (degenerate within rho_L)', RHO_L)):
        ell = np.sqrt(3 * K / L)
        crunch[lab] = {'Lambda_in_GeV4': -L, 't_crunch_GeV^-1': np.pi * ell, 't_crunch_s': np.pi * ell * GEV_INV_S}
    out['AdS_crunch'] = crunch
    out['notes'] = ['every worldline: P(no bubble) = exp(-Gamma V4) -> 0 as V4 -> infinity, in Minkowski (V4 ~ t^4), '
                    'matter domination (t^4) and de Sitter (linear in t); only Gamma = 0 avoids it',
                    'de Sitter volume: false vacuum persists (eternally) iff Gamma/H^4 < 9/(4 pi) = 0.716']
    with open(OUT, 'w') as f:
        json.dump(out, f, indent=1, default=float)
    return out


if __name__ == '__main__':
    print(json.dumps(main(), indent=1, default=float))
