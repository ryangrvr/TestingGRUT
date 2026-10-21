"""Source-faithfulness checks for the Klinkhamer-Volovik q-theory extraction (track: literature).

Everything here re-derives or re-computes statements quoted in NOTES.md; nothing is fitted.
Run:  python3 -I kv_checks.py   (writes ../results/kv_checks.json)

Sections
  A  sympy: KV 0806.2805 model functions (eqs 4.7, 4.8, 4.11, 5.3) and the brief-convention
     translation of the field equations (FRW 00 component, vacuum-energy balance).
  B  sympy: adiabatic electroweak-kick amplitude rho_V0 = (9/8) kappa^2 H^4 K'^2/eps''
     and its reduction to KV 0905.1919 eq (3.4a) and 1001.1939 eq (3.18b).
  C  numerics: integrate KV 0806.2805 eqs (5.2a,b) [same ODEs as 4.15 for sigma=1] with the
     boundary conditions of KV Figs 1/2 and 5, and check the asymptotic law (5.15):
     tau*<h> -> 2/3, tau^2*<r_V> -> (2/9) d_M^2.
  D  numbers: KV (5.16) estimate, the full-coefficient late-time value (2/3) d_M^2 M^2/t^2,
     KV 0905.1919 eq (1.4) numbers, AHHKM eq (3), and the E_ew needed for rho_Lambda,obs.
"""
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

out = {}

# ---------------------------------------------------------------- A
f, y, u, chi, F0, F = sp.symbols('f y u chi F0 F', real=True)
eps = (1 / (2 * chi)) * (-F**2 / F0**2 + F**4 / (3 * F0**4))          # KV 0806.2805 (4.7)
mu0 = sp.simplify(sp.diff(eps, F).subs(F, F0))                           # should be -1/(3 chi F0)  (4.8)
rhoV_at_F0 = sp.simplify(eps.subs(F, F0) - mu0 * F0)                     # should be 0
chi_at_F0 = sp.simplify(1 / (F**2 * sp.diff(eps, F, 2)).subs(F, F0))     # should be chi
rV_f = sp.Rational(1, 2) * (-f**2 + f**4 / 3) - u * f                    # (4.11)
rV_y = sp.expand(rV_f.subs({u: sp.Rational(-1, 3), f: 1 + y}))           # (5.3): y^2/2 + 2y^3/3 + y^4/6
out['A_model'] = {
    'mu0': str(mu0),
    'rhoV(F0)': str(rhoV_at_F0),
    'chi(F0)': str(chi_at_F0),
    'rV(y) at u0=-1/3': str(rV_y),
    'matches_KV_5.3': bool(sp.simplify(rV_y - (y**2 / 2 + sp.Rational(2, 3) * y**3 + y**4 / 6)) == 0),
}

# brief conventions, flat FRW, signature (-+++), MTW curvature: R = 6(Hdot + 2H^2)
t = sp.symbols('t', real=True)
a = sp.Function('a')(t)
K = sp.Function('K')(t)
H = sp.diff(a, t) / a
# 00 component of  K G_mn = T_mn - rho_V g_mn - (g_mn Box - nabla_m nabla_n) K
# G_00 = 3H^2 ; g_00 = -1 ; Box K = -(Kddot + 3 H Kdot) ; nabla_0 nabla_0 K = Kddot
Kd, Kdd = sp.diff(K, t), sp.diff(K, t, 2)
lhs00 = K * 3 * H**2
box = -(Kdd + 3 * H * Kd)
rhs00_extra = -((-1) * box - Kdd)                                          # -(g_00 Box K - nabla0nabla0 K)
out['A_FRW00'] = {
    'extra_term_from_K': str(sp.simplify(rhs00_extra)),                    # expect -3 H Kdot
    'statement': '3 K H^2 = rho_m + rho_V - 3 H Kdot  <=>  KV (4.2a) H^2 = (8pi/3) G rho_tot - H G dG^-1/dt with K = 1/(8 pi G)',
}

# ---------------------------------------------------------------- B
q, q0, Kp, epp, kap, Hs = sp.symbols('q q0 Kp epp kappa H', positive=True)
Kb = sp.symbols('K0', positive=True)
# q-equation near q0 at mu = mu0:  eps''(q0) dq - (1/2) K' R = 0 ;  R = kappa * 3 H^2 (radiation era, K ~ const)
R = kap * 3 * Hs**2
dq = Kp * R / (2 * epp)
rho_V0 = sp.simplify(epp * dq**2 / 2)
out['B_kick'] = {
    'rho_V0_general': str(rho_V0),                                         # 9/8 kappa^2 H^4 Kp^2/epp
    'note': 'brief K (=1/(8 pi G)); trace: -K R = T - 4 rho_V => R = (rho_m - 3P_m)/K = kappa_M * 3 H^2 for rho_m = 3 K H^2',
}
# KV 0806.2805 units: K' = sigma/3 with sigma=1 ; eps'' = 1/(chi F0^2) ; H = h/sqrt(chi F0) ; rho = r/chi
h_ = sp.symbols('h', positive=True)
rV_kv = sp.simplify(rho_V0.subs({Kp: sp.Rational(1, 3), epp: 1 / (chi * F0**2), Hs: h_ / sp.sqrt(chi * F0)}) * chi)
out['B_kick']['KV_0905.1919_eq3.4a_check'] = str(rV_kv)                    # expect kappa^2 h^4 / 8
# 1001.1939: K_KV = q/2 (K_brief = q, K'=1), rho_V = (q-q0)^2/2 (eps''=1), units t_ew -> (9/8) kappa^2 h^4
out['B_kick']['Klinkhamer_1001.1939_eq3.18b_check'] = str(sp.simplify(rho_V0.subs({Kp: 1, epp: 1})))

# ---------------------------------------------------------------- C
def rv(yv):
    return 0.5 * yv**2 + (2.0 / 3.0) * yv**3 + yv**4 / 6.0


def drv(yv):
    return yv + 2.0 * yv**2 + (2.0 / 3.0) * yv**3


def rhs(tau, s, w):
    yv, yd, hv = s
    hd = drv(yv) - 2.0 * hv**2                                             # (5.2b)
    rM = yd * hv + (1.0 + yv) * hv**2 - rv(yv)                             # from (4.15b), sigma=1
    ydd = yd * hv - 2.0 * (1.0 + yv) * hd - 3.0 * (1.0 + w) * rM           # (5.1)
    return [yd, ydd, hd]


def run(w, h1, f1, rM1, tau_end):
    y1 = f1 - 1.0
    fdot = (rv(y1) + rM1 - f1 * h1**2) / h1                                # (4.15b) at tau=1
    sol = solve_ivp(rhs, (1.0, tau_end), [y1, fdot, h1], args=(w,), method='DOP853',
                    rtol=1e-11, atol=1e-13, dense_output=True)
    T = np.linspace(0.6 * tau_end, tau_end, 400001)
    Y = sol.sol(T)
    hv = Y[2]
    r_V = rv(Y[0])
    rM = Y[1] * hv + (1 + Y[0]) * hv**2 - r_V
    return {
        'tau_window': [float(T[0]), float(T[-1])],
        'mean_tau_h': float(np.mean(T * hv)),
        'mean_tau2_rV': float(np.mean(T**2 * r_V)),
        'mean_tau2_rM_end': float(np.mean((T**2 * rM)[-50000:])),
        'min_h_in_window': float(hv.min()),
        'success': bool(sol.success),
    }


resC = {}
resC['fig1_2_w=1/3'] = run(1.0 / 3.0, 0.5, 0.5, 1.0 / 20.0, 3000.0)
resC['fig1_2_w=1/3']['KV_prediction'] = {'mean_tau_h': 2 / 3, 'mean_tau2_rV': 2 / 9}
resC['fig3_4_w=0'] = run(0.0, 0.5, 0.5, 1.0 / 200.0, 3000.0)
resC['fig5_w=0_high'] = run(0.0, 0.5, 0.5, 0.5, 3000.0)
for key in ('fig3_4_w=0', 'fig5_w=0_high'):
    rinf = resC[key]['mean_tau2_rM_end']
    dM = float(np.sqrt(max(0.0, 1 - 2.25 * rinf)))
    resC[key]['rM_inf_measured'] = rinf
    resC[key]['d_M_from_5.14b'] = dM
    resC[key]['KV_prediction_mean_tau2_rV'] = 2 / 9 * dM**2
out['C_numerics_0806.2805'] = resC

# ---------------------------------------------------------------- D
GeV_to_eV = 1e9
hbar_GeV_s = 6.582119569e-25
yr_s = 3.15576e7
M_red = 2.435e18                     # GeV, brief
E_P_nonred = 1.220890e19             # GeV
rho_obs = 2.52e-47                   # GeV^4, brief (Planck 2018 inputs)
t0 = 13.787e9 * yr_s                 # s, Planck 2018 age (exposed datum)


def per_t2(Esq_GeV2, t_s):           # E^2 / t^2 in GeV^4
    return Esq_GeV2 * (hbar_GeV_s / t_s) ** 2


def mev4(rho):                       # (x meV)^4 -> x
    return (rho ** 0.25) * GeV_to_eV * 1e3


D = {}
# KV (5.16): <rho_V> ~ E_Planck^2 / t^2 with nonreduced E_P, t = 1e10 yr
D['KV_5.16_EP^2/t^2_t=1e10yr_meV'] = mev4(per_t2(E_P_nonred**2, 1e10 * yr_s))
D['KV_5.16_quoted_meV'] = 4.0
# full coefficient from (5.15c): <rho_V> = (2/3) d_M^2 M^2 / t^2   (sigma=1, |F0| = 3 M^2)
for tt, lab in ((1e10 * yr_s, 't=1e10yr'), (t0, 't=13.787Gyr')):
    val = (2.0 / 3.0) * per_t2(M_red**2, tt)
    D[f'(2/3)M^2/t^2_{lab}_GeV4'] = val
    D[f'(2/3)M^2/t^2_{lab}_meV'] = mev4(val)
    D[f'(2/3)M^2/t^2_{lab}_over_rho_obs'] = val / rho_obs
# 0905.1919 (1.4): (E_ew^2/E_P)^4 with E_ew = 1e3 GeV, reduced E_P = 2.44e18 (and the rounded 1e18)
for EP, lab in ((2.44e18, 'EP=2.44e18'), (1e18, 'EP=1e18')):
    val = (1e3**2 / EP) ** 4
    D[f'0905.1919_eq1.4_{lab}_GeV4'] = val
    D[f'0905.1919_eq1.4_{lab}_meV'] = mev4(val)
    D[f'0905.1919_eq1.4_{lab}_rho_obs_over_pred'] = rho_obs / val
# E_ew such that (E_ew^2/M)^4 = rho_obs exactly (unit coefficient), and with KV's r_V^num
D['E_ew_unit_coeff_TeV'] = float(np.sqrt(M_red * rho_obs ** 0.25)) / 1e3
for rnum in (0.051, 2.4e-3, 2.0e-5, 1.0):
    D[f'E_ew_for_rVnum={rnum}_TeV(1001.1939 eq5.2, E_P=2.44e18, Lambda=rho_obs)'] = float(
        (rho_obs / rnum) ** 0.125 * np.sqrt(2.44e18)) / 1e3
# radiation-era H at temperature T: H = sqrt(pi^2 g*/90) T^2/M ; T giving H^4 = rho_obs (kappa=1, unit coefficient)
gstar = 106.75
cH = np.sqrt(np.pi**2 * gstar / 90.0)
D['H_over_T2/M_coefficient_gstar=106.75'] = float(cH)
D['T_with_H^4=rho_obs_TeV'] = float(np.sqrt(rho_obs ** 0.25 * M_red / cH)) / 1e3
D['H^4_at_T=159.5GeV_over_rho_obs'] = float((cH * 159.5**2 / M_red) ** 4 / rho_obs)
out['D_numbers'] = D

os.makedirs(os.path.join(os.path.dirname(__file__), '..', 'results'), exist_ok=True)
path = os.path.join(os.path.dirname(__file__), '..', 'results', 'kv_checks.json')
with open(path, 'w') as fh:
    json.dump(out, fh, indent=2)
print(json.dumps(out, indent=2))
