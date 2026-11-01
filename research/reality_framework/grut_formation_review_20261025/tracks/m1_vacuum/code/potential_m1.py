"""Track m1_vacuum, item 1 and item 4: algebra of the M1 one-vacuum potential and the
homogeneous (phi, |X|) sector at fixed charge density.  Independent code (Claude, 1 Nov 2026).

Conventions (read from Codex d4747cb DERIVATION.md / deformation.py, re-derived here):
  V_phi = lambda/4 (phi^2 - v^2)^2 - eps v (phi - v)^3,  lambda = 0.03, eps = 0.002
  canonical complex X with V ⊃ g^2 phi^2 |X|^2 + lambda_X |X|^4   (g = 0.1, lambda_X = 1e-3)
  X = v F e^{-i omega t}/sqrt2, phi = v S, rho = g v r, Omega = omega/m_X, m_X = g v.
  Dimensionless U = (b/4)(S^2-1)^2 - c(S-1)^3 + S^2 F^2/2 + a F^4/4, a = 0.1, b = 3, c = 0.2
  (energy density in units g^2 v^4, i.e. 100 m_X^4 because v = 10 m_X).
  Charge density n = omega v^2 F^2 -> n/m_X^3 = 100 Omega F^2;  Qbar = g^2 Q = Omega * int F^2 d^3rho.
"""
import json, math, sys
import numpy as np
from scipy.optimize import brentq, minimize_scalar, minimize

OUT = sys.argv[1] if len(sys.argv) > 1 else 'potential.json'
LAM, EPS, G, LX = 0.03, 0.002, 0.1, 1e-3
a, b, c = LX / G**2, LAM / G**2, EPS / G**2
V_OVER_M = 1.0 / G          # v / m_X = 10

res = {}

# ---------------------------------------------------------------- 1. one-vacuum algebra
# V'/v^3 in u = phi/v: lambda(u^2-1)u - 3 eps (u-1)^2 = (u-1)[lam u^2 + (lam-3eps) u + 3 eps]
disc = (LAM - 3 * EPS)**2 - 12 * LAM * EPS
# threshold eps where disc = 0: lam^2 - 18 lam eps + 9 eps^2 = 0
eps_lo = LAM * (1 - 2 * math.sqrt(2) / 3)
eps_hi = LAM * (1 + 2 * math.sqrt(2) / 3)
# symbolic check of the factorisation on random u
rng = np.random.default_rng(1)
u = rng.uniform(-5, 5, 2000)
lhs = LAM * (u * u - 1) * u - 3 * EPS * (u - 1)**2
rhs = (u - 1) * (LAM * u * u + (LAM - 3 * EPS) * u + 3 * EPS)
fact_err = float(np.max(np.abs(lhs - rhs)) / np.max(np.abs(lhs)))
# V/v^4 = (u-1)^2 [lam (u+1)^2/4 - eps (u-1)]; bracket discriminant eps(eps - 2 lam)
Vn = LAM / 4 * (u * u - 1)**2 - EPS * (u - 1)**3
Vf = (u - 1)**2 * (LAM * (u + 1)**2 / 4 - EPS * (u - 1))
pos_fact_err = float(np.max(np.abs(Vn - Vf)) / np.max(np.abs(Vn)))
brk_disc = EPS * (EPS - 2 * LAM)
# global minimum of V on a wide grid
ug = np.linspace(-50, 50, 2000001)
Vg = LAM / 4 * (ug * ug - 1)**2 - EPS * (ug - 1)**3
imin = int(np.argmin(Vg))
# V, V', V'' at v ; V'''(v)
V_at_v = LAM / 4 * 0 - EPS * 0
Vp_at_v = 0.0
Vpp_at_v = 2 * LAM            # /v^2
Vppp_at_v = 6 * (LAM - EPS)   # /v
# inflection points V''=0: 3 lam u^2 - 6 eps u - lam + 6 eps = 0
infl = np.roots([3 * LAM, -6 * EPS, -LAM + 6 * EPS])
res['one_vacuum'] = dict(
    discriminant=disc, discriminant_codex=-0.000144,
    eps_threshold_low=eps_lo, eps_threshold_high=eps_hi, eps_threshold_codex=0.001715728752538096,
    derivative_factorisation_max_rel_err=fact_err,
    positivity_factorisation_max_rel_err=pos_fact_err, bracket_discriminant=brk_disc,
    V_min_over_v4_on_grid=float(Vg[imin]), u_at_min=float(ug[imin]),
    V_v=V_at_v, Vp_v=Vp_at_v, Vpp_v_over_v2=Vpp_at_v, Vppp_v_over_v=Vppp_at_v,
    cubic_reduction_fraction=EPS / LAM,
    V_minus_v_over_v4=8 * EPS, V_zero_over_v4=LAM / 4 + EPS, V_zero_symmetric_over_v4=LAM / 4,
    V_minus_v_mX4=8 * EPS * V_OVER_M**4, V_zero_mX4=(LAM / 4 + EPS) * V_OVER_M**4,
    V_zero_symmetric_mX4=LAM / 4 * V_OVER_M**4,
    slope_at_minus_v_over_v3=-12 * EPS,
    inflection_u=sorted([float(x) for x in infl]),
    monotone_below_v=bool(np.all(np.diff(Vg[ug < 1]) < 0)),
    monotone_above_v=bool(np.all(np.diff(Vg[ug > 1]) > 0)),
)

# ---------------------------------------------------------------- 2. hidden tree lower bound (algebra)
kap, lamY, lamS = 0.0034586847211042523, 1.0, 1.0
MY, hh, mu = 30.0, 40.0, 7.159456159513794   # H_i units; ratios are rescaled together by Codex
lamY_eff = lamY - 27 * kap**4 / (16 * LX**3)
# check the x-minimisation numerically for y = 1
xx = np.linspace(0, 20, 400001)
num_min = float(np.min(LX * xx**4 - 2 * kap * xx**3))
s_c = MY**2 / (2 * mu)
r2 = hh**2 / (2 * s_c**2)
margin = (1 - r2)**2 * s_c**2 - mu**2 / lamY_eff
res['hidden_tree_bound'] = dict(lamY_eff=lamY_eff, lamY_eff_codex=0.758516376953125,
                                x_min_numeric_y1=num_min, x_min_formula_y1=-27 * kap**4 / (16 * LX**3),
                                s_c=s_c, r2=r2, large_Sigma_margin=margin, margin_codex=2445.0408755072076)

# ---------------------------------------------------------------- 3. homogeneous (S, F) sector at fixed density
def Uphi(S):
    return b / 4 * (S * S - 1)**2 - c * (S - 1)**3

def Uphi_p(S):
    return b * (S * S - 1) * S - 3 * c * (S - 1)**2

def Uphi_pp(S):
    return b * (3 * S * S - 1) - 6 * c * (S - 1)

def U(S, F):
    return Uphi(S) + 0.5 * S * S * F * F + a / 4 * F**4

def S_of_F(F):
    """Unique root of Uphi'(S) + F^2 S = 0 in (0,1)."""
    if F == 0:
        return 1.0
    return brentq(lambda s: Uphi_p(s) + F * F * s, 1e-14, 1.0, xtol=1e-15, rtol=1e-15)

# (3a) count roots of h(S) = Uphi'(S) + F^2 S over S in [-3, 3] for many F (sign changes on a fine grid)
Sg = np.linspace(-3, 3, 600001)
root_counts = {}
for F in [0, 0.1, 0.5, 1.0, 1.5, 2.0, 2.4065, 2.63, 3.0, 5.0, 10.0]:
    hv = Uphi_p(Sg) + F * F * Sg
    nroots = int(np.sum(np.sign(hv[1:]) != np.sign(hv[:-1])))
    root_counts[str(F)] = nroots
res['phi_equation_root_count_S_in_[-3,3]'] = root_counts
# analytic: dF^2/dS = -2bS + 3c(1 - 1/S^2) < 0 on (0,1): check numerically
Ss = np.linspace(0.01, 0.999, 2000)
res['dF2_dS_max_on_(0,1)'] = float(np.max(-2 * b * Ss + 3 * c * (1 - 1 / Ss**2)))

# (3b) energy per volume at fixed charge density nbar = Omega F^2 (units g^2 v^4 and (g v)^3):
#     e(F) = nbar^2/(2F^2) + U(S(F), F);  stationary in F: nbar^2/F^4 = S^2 + a F^2.
def e_of_F(F, nbar):
    return nbar**2 / (2 * F * F) + U(S_of_F(F), F)

def count_stationary(nbar, Fmax=12.0, n=24000):
    Fg = np.linspace(0.02, Fmax, n)
    gg = np.array([nbar**2 / F**4 - S_of_F(F)**2 - a * F * F for F in Fg])
    idx = np.where(np.sign(gg[1:]) != np.sign(gg[:-1]))[0]
    roots = [brentq(lambda F: nbar**2 / F**4 - S_of_F(F)**2 - a * F * F, Fg[i], Fg[i + 1], xtol=1e-14) for i in idx]
    return roots

# monotonicity of G(F) = F^4 (S^2 + a F^2): if monotone increasing, exactly one stationary F for every nbar
Fg = np.linspace(0.005, 15, 30000)
GF = np.array([F**4 * (S_of_F(F)**2 + a * F * F) for F in Fg])
res['G_monotone_increasing'] = bool(np.all(np.diff(GF) > 0))
res['G_min_slope'] = float(np.min(np.diff(GF) / np.diff(Fg)))

# bulk point: Omega^2 = min_{S,F} 2U/F^2 ; at fixed S: F^2 = 2 sqrt(Uphi/a), Omega^2 = S^2 + 2 sqrt(a Uphi)
def Om2_of_S(S):
    return S * S + 2 * math.sqrt(max(a * Uphi(S), 0.0))

opt = minimize_scalar(Om2_of_S, bounds=(0, 1), method='bounded', options={'xatol': 1e-14})
S_bulk = float(opt.x)
# my own quartic from squaring 2 S sqrt(Uphi) = -sqrt(a) Uphi' and dividing by (S-1)^2:
#   4 S^2 [ (b/4)(S+1)^2 - c (S-1) ] - a [ b S (S+1) - 3 c (S-1) ]^2 = 0
P1 = np.polymul([4, 0, 0], np.polyadd(b / 4 * np.array([1, 2, 1]), [0, -c, c]))
inner = np.polyadd(b * np.array([1, 1, 0]), [0, -3 * c, 3 * c])
P2 = a * np.polymul(inner, inner)
quart = np.polysub(P1, P2)
qroots = np.roots(quart)
S_bulk_quartic = float([z.real for z in qroots if abs(z.imag) < 1e-12 and 0 < z.real < 1][0])
F2_bulk = 2 * math.sqrt(Uphi(S_bulk_quartic) / a)
F_bulk = math.sqrt(F2_bulk)
Om_bulk = math.sqrt(Om2_of_S(S_bulk_quartic))
nbar_bulk = Om_bulk * F2_bulk
# verify it is a stationary point of the three equations and a minimum (Hessian of e(S,F) at fixed nbar)
eqs = [Om_bulk**2 - S_bulk_quartic**2 - a * F2_bulk, Uphi_p(S_bulk_quartic) + F2_bulk * S_bulk_quartic,
       Uphi(S_bulk_quartic) - a * F2_bulk**2 / 4]
def hess_fixed_n(S, F, nbar):
    eSS = Uphi_pp(S) + F * F
    eSF = 2 * S * F
    eFF = 3 * nbar**2 / F**4 + S * S + 3 * a * F * F
    return np.array([[eSS, eSF], [eSF, eFF]])
Hb = hess_fixed_n(S_bulk_quartic, F_bulk, nbar_bulk)
# old symmetric bulk (c=0): Omega0^2 = 2 sqrt(a b/4) at S=0
Om_bulk_sym = math.sqrt(2 * math.sqrt(a * b / 4))
# also a global 2D minimisation of 2U/F^2 over S in [-3,3], F>0 to look for any other minimum (negative S)
best = []
for S0 in np.linspace(-2.5, 1.5, 41):
    for F0 in (0.3, 1.0, 2.5, 5.0):
        r = minimize(lambda z: 2 * U(z[0], z[1]) / z[1]**2 if z[1] > 1e-6 else 1e9, [S0, F0], method='Nelder-Mead',
                     options={'xatol': 1e-12, 'fatol': 1e-14, 'maxiter': 20000})
        best.append((float(r.fun), float(r.x[0]), float(r.x[1])))
best.sort()
distinct = []
for f, s, F in best:
    if all(abs(s - d[1]) > 1e-3 or abs(F - d[2]) > 1e-3 for d in distinct):
        distinct.append((f, s, F))
res['bulk'] = dict(S=S_bulk_quartic, S_from_bounded_minimiser=S_bulk, S_codex=0.1511721078446304,
                   F=F_bulk, F_codex=2.406476756777611, Omega=Om_bulk, Omega_codex=0.7758647074595397,
                   charge_density_over_mX3=nbar_bulk / G**2, density_codex=449.3133678845478,
                   my_quartic_descending=[float(x) for x in quart], codex_quartic=[2.1, 3.76, 2.864, -0.288, -0.036],
                   quartic_roots=[[float(z.real), float(z.imag)] for z in qroots],
                   three_equation_residuals=[float(x) for x in eqs],
                   hessian_fixed_n_eigs=[float(x) for x in np.linalg.eigvalsh(Hb)],
                   Omega_bulk_symmetric_eps0=Om_bulk_sym,
                   U_phi_at_bulk_S_mX4=Uphi(S_bulk_quartic) * 100,
                   distinct_local_minima_of_2U_over_F2=distinct[:6])

# (3c) homogeneous rotating states vs density: Omega(n), e/n, pressure, dOmega/dn sign, counts of stationary points
rows = []
for F in [0.1, 0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.1, 2.2, 2.3, F_bulk, 2.5, 2.63, 3.0, 4.0]:
    S = S_of_F(F)
    Om = math.sqrt(S * S + a * F * F)
    nbar = Om * F * F
    e = 0.5 * Om**2 * F * F + U(S, F)
    P = nbar * Om - e
    # dOmega/dn numerically along the family parametrised by F
    dF = 1e-5
    def OmF(Fx):
        Sx = S_of_F(Fx); return math.sqrt(Sx * Sx + a * Fx * Fx)
    dOm_dF = (OmF(F + dF) - OmF(F - dF)) / (2 * dF)
    dn_dF = (OmF(F + dF) * (F + dF)**2 - OmF(F - dF) * (F - dF)**2) / (2 * dF)
    D = Uphi_pp(S) + F * F
    stat = count_stationary(nbar)
    rows.append(dict(F=F, S=S, Omega=Om, n_over_mX3=nbar / G**2, E_over_Q=e / nbar, pressure=P,
                     dOmega_dn=dOm_dF / dn_dF, longwave_unstable_codex_criterion=bool(a * D < 2 * S * S),
                     D=D, n_stationary_F_at_this_density=len(stat), stationary_F_values=stat,
                     hess_eigs=[float(x) for x in np.linalg.eigvalsh(hess_fixed_n(S, F, nbar))]))
res['homogeneous_states'] = rows
# spinodal: dOmega/dn = 0
Fs = np.linspace(1.5, 2.6, 2201)
Oms = np.array([math.sqrt(S_of_F(F)**2 + a * F * F) for F in Fs])
ns = Oms * Fs**2
dOm = np.gradient(Oms, ns)
isp = int(np.argmin(np.abs(dOm[(Fs > 1.8) & (Fs < 2.5)]))) + int(np.sum(Fs <= 1.8))
res['spinodal'] = dict(F=float(Fs[isp]), Omega=float(Oms[isp]), n_over_mX3=float(ns[isp] / G**2),
                       Omega_min_over_family=float(Oms.min()))

# ---------------------------------------------------------------- 4. walls / negative-phi energetics (m_X units, v = 10)
v = V_OVER_M
Vphys = lambda phi: LAM / 4 * (phi**2 - v**2)**2 - EPS * v * (phi - v)**3
Vp = lambda phi: LAM * (phi**2 - v**2) * phi - 3 * EPS * v * (phi - v)**2
# old symmetric kink tension sigma0 = (2 sqrt2/3) sqrt(lambda) v^3
sigma0 = 2 * math.sqrt(2) / 3 * math.sqrt(LAM) * v**3
# numerical check: sigma = int sqrt(2 V_sym) dphi from -v to v
pg = np.linspace(-v, v, 200001)
sigma0_num = float(np.trapezoid(np.sqrt(2 * LAM / 4 * (pg**2 - v**2)**2), pg))
dV = Vphys(-v)
# roll-down time from rest at phi=-v to phi=+v (1D ODE, no gradients): t = int dphi / sqrt(2 (V(-v) - V(phi)))
pg2 = np.linspace(-v, v, 400001)[1:-1]
integrand = 1 / np.sqrt(np.maximum(2 * (dV - Vphys(pg2)), 1e-300))
t_roll = float(np.trapezoid(integrand, pg2))
# does any static 1D solution other than phi=v exist with phi->v at +inf? energy integral: phi'^2/2 = V(phi) >= 0
#   turning point needs V(phi_t)=0 -> only phi_t = v. So none. (algebraic; recorded as statement)
# amplitude of oscillation about +v that reaches phi=-v: energy density 8 eps v^4
m_phi = math.sqrt(2 * LAM) * v
res['walls'] = dict(V_minus_v_mX4=dV, V_zero_mX4=Vphys(0.0), sigma0_symmetric_kink_mX3=sigma0,
                    sigma0_numeric=sigma0_num, wall_thickness_1_over_mphi=1 / m_phi, m_phi_over_mX=m_phi,
                    critical_patch_radius_sigma_over_dV=sigma0 / dV,
                    V_minus_v_times_thickness_mX3=dV / m_phi,
                    roll_time_minus_v_to_plus_v_over_mX_inverse=t_roll,
                    roll_time_times_mphi=t_roll * m_phi,
                    force_at_minus_v=-Vp(-v),
                    static_wall_solutions="none: with phi->v at infinity, (dphi/dx)^2/2 = V(phi) and a turning point needs V=0, which occurs only at phi=v",
                    energy_density_needed_to_reach_minus_v_mX4=dV)

json.dump(res, open(OUT, 'w'), indent=1, default=float)
print(json.dumps({k: res[k] for k in ('one_vacuum', 'hidden_tree_bound', 'bulk', 'spinodal', 'walls')}, indent=1, default=float))
print('G monotone:', res['G_monotone_increasing'], res['G_min_slope'])
print('root counts', root_counts, 'dF2/dS max', res['dF2_dS_max_on_(0,1)'])
for r in rows:
    print('F=%.4f S=%.5f Om=%.6f n=%.2f E/Q=%.6f P=%.4f dOm/dn=%.3e crit=%s nstat=%d heig=%s' % (
        r['F'], r['S'], r['Omega'], r['n_over_mX3'], r['E_over_Q'], r['pressure'], r['dOmega_dn'],
        r['longwave_unstable_codex_criterion'], r['n_stationary_F_at_this_density'], r['hess_eigs']))
