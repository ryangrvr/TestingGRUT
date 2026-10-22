"""Verifier v_vacua, part 2: nucleation, light-cone volumes, rates, crunch times (track `vacua` NOTES.md section 2).

Independent code. Run: python3 v2_nucleation.py -> ../results/v2_nucleation.json
"""
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results', 'v2_nucleation.json')
res = {}

# ---------------------------------------------------------------- constants (CODATA/IAU; [R])
HBAR_GeVs = 6.582119569e-25
HBARC_GeVm = 1.973269804e-16
GYR_s = 3.15576e16            # Julian Gyr
GPC_m = 3.0856775814913673e25
MPC_m = GPC_m / 1e3
H0_kms = 67.4
H0_GeV = H0_kms * 1e3 / MPC_m * HBAR_GeVs        # H0 in GeV
M = 2.435e18
K = M ** 2
rhoL = 2.52e-47

# ---------------------------------------------------------------- analytic V4
tt, tp, eta, etap, Hs = sp.symbols('t tp eta etap H', positive=True)
V4_mink = sp.simplify(sp.Rational(4, 3) * sp.pi * sp.integrate((tt - tp) ** 3, (tp, 0, tt)))
# matter domination: a = c eta^2, t = c eta^3/3
c = sp.symbols('c', positive=True)
V4_md = sp.Rational(4, 3) * sp.pi * sp.integrate((c * etap ** 2) ** 4 * (eta - etap) ** 3, (etap, 0, eta))
t_md = c * eta ** 3 / 3
ratio_md = sp.simplify(V4_md / t_md ** 4)
# H = 2/(3t): V4 H^4
ratio_md_H = sp.simplify(ratio_md * (sp.Rational(3, 2)) ** 4)   # V4 = ratio * t^4 = ratio*(3/(2H))^4 ... H^4 V4
# de Sitter from a = 1 at t = 0: a = 1/(1 - H eta'), eta' in [0, eta), t = -ln(1 - H eta)/H
ad = 1 / (1 - Hs * etap)
V4_ds = sp.Rational(4, 3) * sp.pi * sp.integrate(ad ** 4 * (eta - etap) ** 3, (etap, 0, eta))
u = sp.symbols('u', positive=True)  # u = H t, eta = (1 - e^{-u})/H
V4_ds_u = sp.simplify(V4_ds.subs(eta, (1 - sp.exp(-u)) / Hs) * Hs ** 4 * 3 / (4 * sp.pi))
asym = sp.limit(V4_ds_u - u, u, sp.oo)
res['analytic'] = {'Minkowski V4': str(V4_mink), 'MD V4/t^4': str(ratio_md), 'MD V4/t^4 numeric': float(ratio_md),
                   'MD V4 H^4': str(sp.nsimplify(ratio_md * sp.Rational(16, 81))),
                   'track_MD': '(81 pi/1485) t^4 = (16 pi/1485) H^-4 = 0.1714 t^4',
                   'dS: (3H^4/4pi) V4 - H t as t->inf': str(asym), 'track_dS': 'H t - 11/6 (Buttazzo)'}

# ---------------------------------------------------------------- LambdaCDM light-cone volume, my own integration
def lcdm(Om=0.315, Or=9e-5):
    OL = 1 - Om - Or
    E = lambda a: np.sqrt(Or / a ** 4 + Om / a ** 3 + OL)
    # conformal time eta(a) = int_0^a da'/(a'^2 E)  (H0 = 1 units)
    eta_of = lambda a: quad(lambda s: 1 / (s ** 2 * E(s)), 0, a, epsabs=1e-13, epsrel=1e-12, limit=200)[0]
    eta0 = eta_of(1.0)
    t0 = quad(lambda s: 1 / (s * E(s)), 0, 1, epsabs=1e-13, epsrel=1e-12, limit=200)[0]
    # V4 = (4pi/3) int a^4 (eta0 - eta)^3 d eta = (4pi/3) int_0^1 a^2 (eta0 - eta(a))^3 / E(a) da
    f = lambda a: a ** 2 * (eta0 - eta_of(a)) ** 3 / E(a)
    V4 = 4 * np.pi / 3 * quad(f, 0, 1, epsabs=1e-12, epsrel=1e-10, limit=200)[0]
    return {'Om': Om, 'Or': Or, 'eta0_H0': eta0, 't0_H0': t0, 't0_Gyr': t0 / H0_GeV * HBAR_GeVs / GYR_s,
            'V4_H0^4': V4}


res['LCDM'] = [lcdm(0.315, 9e-5), lcdm(0.315, 0.0), lcdm(0.30, 0.0), lcdm(0.27, 0.0), lcdm(1.0, 0.0)]
res['LCDM_note'] = ('Buttazzo et al. 1307.3536 eq. (tunrate) quote ~0.15/H0^4 (internal note in their TeX: 0.148); '
                    'their Omega_m is not stated at that line. Om = 1 gives 16 pi/1485 = 0.0338.')


# de Sitter late slope dV4/dt -> (4pi/3)/H_L^3 for LCDM; check numerically by V4(t) at two late times
def lcdm_late_slope(Om=0.315, Or=9e-5):
    OL = 1 - Om - Or
    HL = np.sqrt(OL)
    E = lambda a: np.sqrt(Or / a ** 4 + Om / a ** 3 + OL)
    eta_of = lambda a: quad(lambda s: 1 / (s ** 2 * E(s)), 0, a, epsabs=1e-14, epsrel=1e-12, limit=400)[0]
    t_of = lambda a: quad(lambda s: 1 / (s * E(s)), 0, a, epsabs=1e-14, epsrel=1e-12, limit=400)[0]

    def V4(aobs):
        e0 = eta_of(aobs)
        return 4 * np.pi / 3 * quad(lambda a: a ** 2 * (e0 - eta_of(a)) ** 3 / E(a), 0, aobs, epsrel=1e-11,
                                    limit=400)[0]
    a1, a2 = 1e3, 2e3
    slope = (V4(a2) - V4(a1)) / (t_of(a2) - t_of(a1))
    return {'slope_numeric': slope, 'expected_(4pi/3)/H_L^3': 4 * np.pi / 3 / HL ** 3}


res['dS_late_slope'] = lcdm_late_slope()

# ---------------------------------------------------------------- volume criterion (de Sitter)
G4 = sp.symbols('Gamma_over_H4', positive=True)
growth = sp.diff(3 * u - sp.Rational(4, 3) * sp.pi * G4 * u, u)
res['volume_criterion'] = {'d ln(a^3 f)/d(Ht)': str(growth), 'threshold': str(sp.solve(growth, G4)),
                           'numeric': float(sp.Rational(9, 4) / sp.pi)}

# ---------------------------------------------------------------- SM rates (CMS 1707.09301 abstract numbers)
unit = (HBAR_GeVs / GYR_s) * (HBARC_GeVm / GPC_m) ** 3      # 1 Gyr^-1 Gpc^-3 in GeV^4
Om, Or = 0.315, 9e-5
OL = 1 - Om - Or
HL = H0_GeV * np.sqrt(OL)
V4_today = res['LCDM'][0]['V4_H0^4'] / H0_GeV ** 4
rows = {}
for lab, lg in (('1e-577', -577.0), ('1e-295', -295.0), ('1e-1465', -1465.0)):
    lgG = lg + np.log10(unit)
    t_ds = np.log10(3 / (4 * np.pi)) + 3 * np.log10(HL) - lgG          # GeV^-1
    t_mink = (np.log10(3 / np.pi) - lgG) / 4
    t_md = (np.log10(1 / float(ratio_md)) - lgG) / 4
    conv = np.log10(HBAR_GeVs / GYR_s)
    rows[lab] = {'log10_Gamma_GeV4': lgG, 'log10_P_past_cone': lgG + np.log10(V4_today),
                 'log10_Gamma/H_L^4': lgG - 4 * np.log10(HL),
                 'log10_t_dS_worldline_Gyr': t_ds + conv, 'log10_t_Minkowski_Gyr': t_mink + conv,
                 'log10_t_MD_Gyr': t_md + conv}
res['rates'] = {'unit_GeV4_per_Gyr^-1Gpc^-3': unit, 'H0_GeV': H0_GeV, 'H_L_GeV': HL, 'rows': rows}

# ---------------------------------------------------------------- AdS crunch
ts, l = sp.symbols('t l', positive=True)
aa = l * sp.sin(ts / l)
# open FRW (k = -1) with Lambda_in < 0: (a'/a)^2 = -1/l^2 + 1/a^2, l^2 = 3K/|Lambda|
fr = sp.simplify((sp.diff(aa, ts) / aa) ** 2 - (-1 / l ** 2 + 1 / aa ** 2))
cr = {}
for lab, LamIn in (('7.5e70', 7.5e70), ('2e65', 2.0e65), ('rho_Lambda', rhoL)):
    lval = np.sqrt(3 * K / LamIn)
    cr[lab] = {'t_c_s': float(np.pi * lval * HBAR_GeVs), 't_c_Gyr': float(np.pi * lval * HBAR_GeVs / GYR_s)}
res['AdS_crunch'] = {'friedmann_residual': str(fr), 'times': cr}

# ---------------------------------------------------------------- CDL threshold normalization check
# thin wall: decay of a zero-energy false vacuum to negative energy forbidden if eps <= (3/4) kappa sigma^2, kappa=1/M^2
res['CDL'] = 'eps_crit = 3 sigma^2/(4 M^2) with reduced M ([R] Coleman-De Luccia 1980, kappa = 8 pi G)'

with open(OUT, 'w') as f:
    json.dump(res, f, indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
