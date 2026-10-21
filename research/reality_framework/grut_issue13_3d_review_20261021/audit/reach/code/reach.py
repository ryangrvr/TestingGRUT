"""Empirical reach of unchanged Candidate 1 (gravity-only coupling) in the classes that pass the wall tests.

Track `reach`, audit of Codex a69bd86 vs Claude c7b8480 (issue #13). Read-only inputs:
  - results/qball_family.json (this directory; computed by qball_family.py from the archived soliton)
  - Codex 479f111 production endpoint component lists (extracted read-only into the scratchpad; summary numbers are
    copied below with their provenance; path given by env CODEX479_RAW if a recompute is wanted)
  - Claude c7b8480 module numbers quoted with [Q].
Labels: [C] computed here, [D] closed form, [E] estimate, [Q] quoted, [R] recalled literature value (not fetched in this
session; order-of-magnitude use only), [P] preparation input, [RETRO] retrospective datum (known before the analysis).
"""
import json, os, csv
import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')

# ---------------------------------------------------------------- constants
MP = 2.435323e18            # reduced Planck mass [GeV]
GEV_CM = 1.97327e-14        # 1 GeV^-1 in cm
GEV_G = 1.78266e-24         # 1 GeV in g
MSUN_GEV = 1.11542e57       # solar mass in GeV
MPC_CM = 3.08568e24
PC_CM = 3.08568e18
T0 = 2.7255 * 8.617333e-14  # GeV  [RETRO] FIRAS T_CMB
S0 = 2891.2                 # cm^-3  [Q relic_yield]
RHOC_H2 = 1.05375e-5        # GeV cm^-3 [Q]
OMEGA_C_H2 = 0.120          # [RETRO] Planck 2018
OMEGA_B_H2 = 0.02233        # [Q] BBN (Mossa et al.)
YB = 8.674e-11              # [Q] relic_yield / charge_phase
OMEGA_R_H2 = 4.18e-5        # photons + 3 massless nu [R]
H0_KMSMPC = 67.4            # [RETRO] only for the background integral of free streaming
h = H0_KMSMPC / 100
OMEGA_M = 0.315; OMEGA_L = 1 - OMEGA_M

# benchmark action
v_b = 1000.0; gX = 0.1; lam_phi = 0.01; lam_X = 0.001
a_ = gX**2 / lam_phi; b_ = lam_X / lam_phi
mX = gX * v_b; mh = np.sqrt(2 * lam_phi) * v_b
mstar_over_mX = (b_) ** 0.25 / np.sqrt(a_)          # 0.5623
sigma_wall = (2 * np.sqrt(2 * lam_phi) / 3) * v_b**3  # 9.428e7 GeV^3 [Q]

out = {}

# ---------------------------------------------------------------- g*(T)
_gtab = np.array([  # T [GeV], g*, g*s   (approximate SM step table) [R]
    [1e-6, 3.36, 3.91], [1e-4, 3.36, 3.91], [2e-4, 3.40, 3.95], [5e-4, 7.4, 7.8], [1e-3, 10.75, 10.75],
    [0.08, 10.75, 10.75], [0.12, 14.0, 14.0], [0.15, 17.25, 17.25], [0.2, 50.0, 50.0], [0.4, 61.75, 61.75],
    [1.0, 61.75, 61.75], [2.0, 72.25, 72.25], [5.0, 75.75, 75.75], [20, 86.25, 86.25], [90, 86.25, 86.25],
    [150, 96.25, 96.25], [200, 106.75, 106.75], [1e20, 106.75, 106.75]])


def gstar(T):
    return np.interp(np.log(T), np.log(_gtab[:, 0]), _gtab[:, 1])


def gstars(T):
    return np.interp(np.log(T), np.log(_gtab[:, 0]), _gtab[:, 2])


def Hrad(T):
    return np.sqrt(np.pi**2 * gstar(T) / 90) * T**2 / MP


def entropy(T):
    return 2 * np.pi**2 / 45 * gstars(T) * T**3


def a_of_T(T):
    return (T0 / T) * (3.91 / gstars(T)) ** (1 / 3)


H0_GEV = 2.13312e-42 * h  # GeV


def H_of_a(aa):
    # radiation with g* steps, matter, Lambda; T from a by fixed-point (2 iterations)
    T = T0 / aa
    for _ in range(3):
        T = T0 / aa * (3.91 / gstars(T)) ** (1 / 3)
    om_r = OMEGA_R_H2 / h**2 * (gstar(T) / 3.36) * (3.91 / gstars(T)) ** (4 / 3)
    return H0_GEV * np.sqrt(om_r / aa**4 + OMEGA_M / aa**3 + OMEGA_L)


# ================================================================ 1. calibration: hidden charge = all DM
rho_over_s_DM = OMEGA_C_H2 * RHOC_H2 / S0     # GeV  (0.437 eV)
coef = YB * S0 / RHOC_H2                       # Omega h^2 per (GeV of energy per unit charge * r)
egv_r = OMEGA_C_H2 / coef                      # (E/Q)*r in GeV
fam = json.load(open(os.path.join(RES, 'qball_family.json')))['family']
Qf = np.array([o['Q_phys'] for o in fam]); ef = np.array([o['E_over_mXQ'] for o in fam])
R90f = np.array([o['R90_energy_ell'] for o in fam]); Rhf = np.array([o['Rhalf_charge_ell'] for o in fam])
order = np.argsort(Qf); Qf, ef, R90f, Rhf = Qf[order], ef[order], R90f[order], Rhf[order]
# binding threshold (E/Q = m_X) on the computed branch
Qbind = float(np.interp(1.0, ef[::-1], Qf[::-1]))
# large-Q asymptotic fit e(Q) = m* + c1 Q^-1/3 + c2 Q^-2/3 (fit on Q >= 8e4)
sel = Qf >= 8e4
A = np.vstack([Qf[sel] ** (-1 / 3), Qf[sel] ** (-2 / 3)]).T
c12, *_ = np.linalg.lstsq(A, ef[sel] - mstar_over_mX, rcond=None)


def e_of_Q(Q):
    Q = np.asarray(Q, float)
    inside = np.interp(np.log(Q), np.log(Qf), ef)
    tw = mstar_over_mX + c12[0] * Q ** (-1 / 3) + c12[1] * Q ** (-2 / 3)
    return np.where(Q <= Qf[-1], inside, tw)


def R_of_Q(Q, which='R90'):
    Q = np.asarray(Q, float)
    arr = R90f if which == 'R90' else Rhf
    nstar = 1.778279  # Qbar per ell^3
    Rtw = (3 * (Q * lam_phi) / (4 * np.pi * nstar)) ** (1 / 3)
    # beyond the computed family: thin-wall radius plus the computed offset at the last point
    off = arr[-1] - (3 * (Qf[-1] * lam_phi) / (4 * np.pi * nstar)) ** (1 / 3)
    return np.where(Q <= Qf[-1], np.interp(np.log(Q), np.log(Qf), arr), Rtw + off)


cal = {}
for label, e in [('Q-matter floor e=m*/m_X (relic_yield r_max convention)', mstar_over_mX),
                 ('archived soliton Q=39,944, e=0.91109', 0.91108652),
                 ('Codex production endpoint, charge-weighted Q ~ 4-6e5, e(Q)', float(e_of_Q(5e5))),
                 ('free X quanta e=1', 1.0)]:
    cal[label] = dict(e=e, r_allDM_at_v1TeV=egv_r / (e * gX * v_b))
out['1_calibration_allDM'] = dict(
    label='[C] with [RETRO] omega_c = 0.120 (Planck 2018) and [Q] Y_B = 8.674e-11 (BBN)',
    omega_per_GeV_r=coef, product_e_gX_v_r_GeV=egv_r,
    statement='All-DM calibration fixes only the product (E/Q) * r = e * gX * v * r = %.4f GeV; v and r separately free' % egv_r,
    by_E_over_Q=cal,
    relic_yield_rmax_check=dict(quoted=0.0897, recomputed=egv_r / (mstar_over_mX * gX * v_b),
                                omega_asym_per_r=coef * mstar_over_mX * gX * v_b),
    codex_maser_product_check=dict(bound_GeV=63.93, mX_r_allDM_GeV=[mX * egv_r / (1.0 * mX), mX * egv_r / (mstar_over_mX * mX)],
                                   H0_floor_kms_Mpc=100 * np.sqrt(OMEGA_B_H2 + OMEGA_C_H2),
                                   verdict='satisfied by construction; the floor is an upper bound on omega_h and is saturated only far above omega_c'),
    T_n_threshold_TeV=[198.2 * (egv_r / (e * gX * v_b)) ** (-1 / 3) for e in (1.0, mstar_over_mX)],
)

# ================================================================ 2. Q-ball population in units of v
codex_endpoint = {  # [Q] Codex 479f111 raw/production_*.csv.components.csv (extracted read only); Qbar = Q_phys/100
    'source': 'git show 479f111:.../grut_3d_empirical_20261020/raw_diagnostics.tar.gz -> raw/production_{coarse,medium,fine}.csv.components.csv',
    'box_Qbar_total': 0.2200852108 * 64**3,
    'fine_10pct': dict(components=39, charge_fraction=0.680, top_Qbar=[12460.1, 5775.0, 4589.2, 3454.5, 3043.8], charge_weighted_mean_Qbar=6185.3, n_above_binding=15),
    'fine_5pct': dict(components=111, charge_fraction=0.744, top_Qbar=[15199.1, 14336.6, 3744.7, 3253.9, 1516.6], charge_weighted_mean_Qbar=10721.7),
    'medium_10pct': dict(components=36, charge_weighted_mean_Qbar=3831.6), 'coarse_10pct': dict(components=29, charge_weighted_mean_Qbar=3646.6),
    'fine_2pct': 'one connected network holds 85 % of the charge: not a particle population',
    'scope': 'Codex: threshold-dependent concentrations at tau=240 of an unexpanding replay, not a converged soliton count or mass function',
}
cr = os.environ.get('CODEX479_RAW')
if cr and os.path.isdir(cr):
    rows = list(csv.DictReader(open(os.path.join(cr, 'production_fine.csv.components.csv'))))
    Q10 = np.array([float(r['Q']) for r in rows if r['threshold'] == '0.1'])
    codex_endpoint['recomputed_fine_10pct'] = dict(n=len(Q10), cw_mean=float((Q10**2).sum() / Q10.sum()))


def qball_props(Q, v=v_b):
    e = float(e_of_Q(Q)); M = e * gX * v * Q
    R90 = float(R_of_Q(Q, 'R90')) / (np.sqrt(lam_phi) * v)  # GeV^-1
    Rh = float(R_of_Q(Q, 'half')) / (np.sqrt(lam_phi) * v)
    sig_hi = np.pi * (2 * R90) ** 2 * GEV_CM**2   # contact at R90+R90
    sig_lo = np.pi * (2 * Rh) ** 2 * GEV_CM**2    # contact at half-charge radii
    return dict(Q=Q, e=e, M_GeV=M, M_g=M * GEV_G, M_Msun=M / MSUN_GEV, R90_cm=R90 * GEV_CM,
                sigma_over_m_cm2_g=[sig_lo / (M * GEV_G), sig_hi / (M * GEV_G)])


pop = {}
for lab, Q in [('binding threshold', Qbind), ('archived', 39944.0), ('per fastest wavelength at Codex density (estimate)', None),
               ('Codex endpoint charge-weighted (fine 10%)', 6185.3 * 100), ('Codex endpoint largest (fine 5%)', 15199.1 * 100)]:
    if Q is None:
        continue
    pop[lab] = qball_props(Q)

# ================================================================ 3. fragmentation: dispersion relation (charge_phase cubic)


def growth(x, k):
    S0sq = 1 - a_ * x; Om2 = a_ - (a_**2 - b_) * x; ms2 = 3 * S0sq - 1 + a_ * x
    # (A + 2bx) A (A + ms2) + 4 Om2 (A - k^2)(A + ms2) - 4 a^2 x S0^2 A = 0
    c3 = 1.0
    c2 = 2 * b_ * x + ms2 + 4 * Om2
    c1 = 2 * b_ * x * ms2 + 4 * Om2 * (ms2 - k**2) - 4 * a_**2 * x * S0sq
    c0 = -4 * Om2 * k**2 * ms2
    roots = np.roots([c3, c2, c1, c0])
    lam2 = [r.real - k**2 for r in roots if abs(r.imag) < 1e-12 * max(1, abs(r.real))]
    lam2 = [l for l in lam2 if l > 0]
    return np.sqrt(max(lam2)) if lam2 else 0.0


def gmax(x):
    kk = np.sqrt(x) * np.linspace(0.02, 3.0, 600)
    g = np.array([growth(x, k) for k in kk])
    i = int(np.argmax(g))
    res = minimize_scalar(lambda k: -growth(x, k), bounds=(kk[max(i - 1, 0)], kk[min(i + 1, len(kk) - 1)]), method='bounded')
    return -res.fun, res.x


disp = []
for x in [1e-8, 1e-6, 1e-4, 1e-3, 1e-2, 0.1, 0.25, 0.5]:
    g, k = gmax(x)
    disp.append(dict(x=x, gamma_max=g, k_max=k, gamma_over_x=g / x, k_over_sqrtx=k / np.sqrt(x)))
small = disp[0]
cg, ck = small['gamma_over_x'], small['k_over_sqrtx']
out['3_fragmentation_dispersion'] = dict(
    label='[C] cubic of charge_phase section 2, branch 1; units sqrt(lam_phi) v = 100 GeV (v/TeV) for rates, 1/ell for k',
    check_charge_phase=dict(quoted='gamma_max 0.0045/0.045/0.112/0.226 at x=0.01/0.1/0.25/0.5; k=0.404 at x=0.25',
                            got=[(d['x'], round(d['gamma_max'], 4), round(d['k_max'], 4)) for d in disp if d['x'] >= 0.01]),
    rows=disp, small_x_law=dict(gamma_max='%.4f x' % cg, k_max='%.4f sqrt(x)' % ck))


def nbar_of_T(T, r, v):
    # n_Q = r Y_B s(T); nbar = n / (sqrt(lam_phi) v^3); x ~ nbar / Omega ~ nbar (small x)
    return r * YB * entropy(T) / (np.sqrt(lam_phi) * v**3)


def Qfrag_of_x(x):
    g, k = gmax(x) if x > 1e-8 else (cg * x, ck * np.sqrt(x))
    return x * (2 * np.pi / k) ** 3 / lam_phi   # physical charge per fastest-growing cube (Omega~1) [E]


def T_floor(r, v, Nefold=30.0):
    # lowest supply temperature at which the homogeneous condensate fragments within a Hubble time: gamma_max = N H
    f = lambda lT: np.log(cg * nbar_of_T(np.exp(lT), r, v) * np.sqrt(lam_phi) * v / (Nefold * Hrad(np.exp(lT))))
    lo, hi = np.log(1e-9), np.log(1e8)   # 1 eV .. 1e8 GeV
    if f(lo) > 0:
        return 1e-9   # fragments for every supply temperature above 1 eV (reported as the 1 eV edge)
    return float(np.exp(brentq(f, lo, hi)))


frag_rows = []
for v in [1000.0, 100.0, 10.0, 1.0, 0.3, 0.1]:
    for e_lab, e in [('e=1', 1.0), ('e=m*', mstar_over_mX)]:
        r = egv_r / (e * gX * v)
        Tn = 198.2e3 * r ** (-1 / 3) * (v / 1000.0)  # GeV, n_c threshold [Q] (g*s = 106.75 convention)
        Tf1 = T_floor(r, v, 1.0); Tf30 = T_floor(r, v, 30.0)
        x_f = nbar_of_T(Tf30, r, v)
        frag_rows.append(dict(v_GeV=v, e_assumed=e_lab, r_allDM=r, T_n_GeV=Tn, T_floor_gamma_eq_H_GeV=Tf1, floor_at_1eV_edge=bool(Tf30 <= 1e-9),
                              T_floor_gamma_eq_30H_GeV=Tf30, x_at_floor=x_f, Q_frag_at_floor=Qfrag_of_x(x_f),
                              Q_frag_near_Tn=Qfrag_of_x(0.9 * 0.316)))
out['3_fragmentation_window'] = dict(
    label='[E] homogeneous sub-critical condensate supplied at T_Q (form [P]); all-DM r; fragments into Q-balls only if gamma_max >> H',
    rows=frag_rows,
    note='Below T_floor the supplied condensate stays homogeneous (behaves as CDM); gamma/H falls as 1/a so it never fragments later by this instability.')

# ================================================================ 4. free streaming


def lambda_fs(T_i, beta_i, M=None):
    """comoving free-streaming length [Mpc] for a particle with velocity beta_i at T_i (NR mass M irrelevant if beta_i<1;
    beta_i>=1 treated as p/M = beta_i (momentum over mass))."""
    a_i = a_of_T(T_i)
    p_over_m_i = beta_i / np.sqrt(max(1 - beta_i**2, 1e-30)) if beta_i < 1 else beta_i
    def integrand(lna):
        aa = np.exp(lna)
        pm = p_over_m_i * a_i / aa
        vel = pm / np.sqrt(1 + pm**2)
        return vel / (aa * H_of_a(aa))  # d(comoving)/dln a in GeV^-1
    val, _ = quad(integrand, np.log(a_i), 0.0, limit=400)
    return val * GEV_CM / MPC_CM


fs = []
for T in [4.4e5, 1e3, 10.0, 0.3, 1e-2, 1e-3, 1e-4]:
    for beta in [0.01, 0.5, 0.9]:
        fs.append(dict(T_Q_GeV=T, beta_i=beta, lambda_fs_Mpc=lambda_fs(T, beta)))
lam_floor_v1 = [dict(T=row['T_floor_gamma_eq_30H_GeV'], v=row['v_GeV'], e=row['e_assumed'],
                     lam_fs_Mpc_beta0p5=lambda_fs(row['T_floor_gamma_eq_30H_GeV'], 0.5),
                     lam_fs_Mpc_beta0p9=lambda_fs(row['T_floor_gamma_eq_30H_GeV'], 0.9)) for row in frag_rows]
# debris of a biased network decaying at T_safe with momentum p/m = 1/f
deb = [dict(f=f, T_ann_GeV=Ta, p_over_m=1 / f, lambda_fs_Mpc=lambda_fs(Ta, 1 / f)) for f, Ta in [(1.0, 0.0306), (0.1, 3.3e-3), (0.01, 3.1e-4)]]
rho_m0 = OMEGA_M * 2.775e11 * h**2  # Msun / Mpc^3
# scale at which the latest-supply corner (T_Q at the gamma = 30 H floor) reaches lambda_fs = 0.03 Mpc
def lam_floor_of_v(v, beta, e=0.75):
    r = egv_r / (e * gX * v)
    return lambda_fs(T_floor(r, v, 30.0), beta)
v_lya = {}
for beta in (0.01, 0.5, 0.9):
    g = lambda lv: np.log(lam_floor_of_v(np.exp(lv), beta) / 0.03)
    try:
        v_lya['beta=%g' % beta] = float(np.exp(brentq(g, np.log(2.0), np.log(1000.0))))
    except ValueError:
        v_lya['beta=%g' % beta] = None
out['4_free_streaming'] = dict(
    label='[C] numerical integral of v dt/a with g*(T) steps; formation velocities not computed [E]: beta_i = 0.01 / 0.5 / 0.9 bracket (energy bound below)',
    grid=fs, at_fragmentation_floor=lam_floor_v1, v_GeV_where_floor_corner_reaches_0p03Mpc=v_lya, wall_debris_biased_class=deb,
    Mcut_Msun={'lambda=%g Mpc' % L: (4 * np.pi / 3) * rho_m0 * (0.5 * L) ** 3 for L in (1e-12, 1e-9, 6e-6, 0.03, 0.1)},
    beta_energy_bound=dict(e_hom_codex_density=0.9442, e_Q_1e5=float(e_of_Q(1e5)), e_Q_6e5=float(e_of_Q(6.185e5)),
                           beta_max=[float(np.sqrt(2 * (0.9442 - e_of_Q(Q)) / e_of_Q(Q))) for Q in (1e5, 6.185e5)],
                           note='[E] upper bound if all released binding energy went into Q-ball bulk motion; most goes to X/phi waves; not computed'),
    lyman_alpha_reference='[R] thermal-WDM bounds m_WDM >~ 3-5 keV correspond to lambda_fs ~ 0.03-0.1 Mpc; halo-scale probes reach ~1e7-1e8 Msun')

# ================================================================ 5. self-interaction
# elementary threshold amplitudes [D]
M_XX = abs(-4 * lam_X + 4 * gX**4 / lam_phi)
sig_XX = M_XX**2 / (32 * np.pi * 4 * mX**2) * GEV_CM**2
M_hh = abs(-6 * lam_phi + 30 * lam_phi)
sig_hh = M_hh**2 / (32 * np.pi * 4 * mh**2) * GEV_CM**2
elem = dict(XX_to_XX=dict(M_threshold=M_XX, sigma_cm2=sig_XX, sigma_over_m=sig_XX / (mX * GEV_G)),
            hh_to_hh=dict(M_threshold=M_hh, sigma_cm2=sig_hh, sigma_over_m=sig_hh / (mh * GEV_G)))
# v needed for SIDM-relevant geometric Q-ball cross-sections (sigma/m ~ v^-3 at fixed Q)
sidm = {}
for lab, Q in [('archived Q=39,944', 39944.0), ('Codex endpoint cw mean Q=6.2e5', 6.185e5), ('binding threshold', Qbind)]:
    som = qball_props(Q)['sigma_over_m_cm2_g']
    sidm[lab] = dict(sigma_over_m_v1TeV=som,
                     v_for_1cm2g_GeV=[v_b * (s / 1.0) ** (1 / 3) for s in som],
                     v_for_0p1cm2g_GeV=[v_b * (s / 0.1) ** (1 / 3) for s in som])
v_sidm = 0.15
r_sidm = egv_r / (0.9 * gX * v_sidm)
out['5_self_interaction'] = dict(
    label='[D/C/E] geometric Q-ball contact cross-section pi(R1+R2)^2 (band: half-charge to R90 radii); elementary threshold Born amplitudes',
    qball_population_v1TeV=pop, elementary=elem, sidm_scale_change=sidm,
    sidm_point_consistency=dict(v_GeV=v_sidm, r_allDM=r_sidm, T_n_GeV=198.2e3 * r_sidm ** (-1 / 3) * (v_sidm / 1000)),
    sidm_reference='[R] cluster/Bullet bounds sigma/m <~ 0.1-1 cm^2/g; dwarf-scale preferences up to O(1-10) cm^2/g; forecasts ~0.01-0.1 cm^2/g')

# ================================================================ 6. compact objects / lensing
M_arch = qball_props(39944.0)['M_GeV']
causal = []
for T in [4.4e5, 1e3, 1.0]:
    r = egv_r / (0.75 * gX * v_b)
    Qmax = r * YB * entropy(T) * (4 * np.pi / 3) / Hrad(T) ** 3
    causal.append(dict(T_Q_GeV=T, Q_horizon=Qmax, M_horizon_Msun=Qmax * mstar_over_mX * mX / MSUN_GEV))
out['6_lensing'] = dict(
    label='[C/E]',
    qball_mass_range_v1TeV_Msun=[qball_props(Qbind)['M_Msun'], qball_props(1e14)['M_Msun']],
    archived_mass_Msun=M_arch / MSUN_GEV,
    causal_upper_bound_all_DM_in_one_horizon=causal,
    lensing_reference='[R] microlensing (HSC, EROS/OGLE) probes ~1e-11 to 10 Msun; femtolensing claims ~1e-16-1e-13 Msun (disputed)',
    gap_decades=float(np.log10(1e-16 / qball_props(1e14)['M_Msun'])))

# ================================================================ 7. N_eff and BBN
neff = {}
for T in [4.4e5, 1e4, 1e3]:
    frac = (4 / 3) * egv_r * YB / T * gstars(T) / gstar(T)  # rho_h/rho_rad = (4/3)(E/Q) r Y_B (g*s/g*) / T  [D]
    # if a fraction eps of the radiation energy at T became dark radiation: DeltaN_eff = eps (g*(T)/2)(3.91/g*s(T))^(4/3) / ((7/8)(4/11)^(4/3))
    dneff = frac * (gstar(T) / 2) * (3.91 / gstars(T)) ** (4 / 3) / ((7 / 8) * (4 / 11) ** (4 / 3))
    neff[T] = dict(rho_h_over_rho_rad=frac, DeltaNeff_if_all_hidden_energy_became_radiation_at_T=dneff)
bbn = dict(rho_DM_over_rho_rad_at_1MeV=rho_over_s_DM * entropy(1e-3) / (np.pi**2 / 30 * gstar(1e-3) * 1e-12),
           note='identical to any CDM: the hidden matter is pressureless at BBN with gravity-only coupling; no injection, no portal')
sym = dict(xi_kin_band=[4e-5, 7e-5], Omega_sym_h2=[0.272 * (xi / 1e-3) ** 3 for xi in (4e-5, 7e-5)],
           DeltaNeff_hidden_bath_at_BBN_upper='~ (4/7)(g_h) xi^4 ... < 1e-15 (all hidden species massive >= 100 GeV (v/TeV); none relativistic at BBN)')
out['7_Neff_BBN'] = dict(label='[D/E]', cold_class=neff, BBN=bbn, biased_class_symmetric_relic=sym,
                         reference='[R/RETRO] Planck N_eff = 2.99 +- 0.17; CMB-S4 forecast sigma ~ 0.03')

# ================================================================ 8. gravitational waves
eps_gw = 0.7; A = 0.8


def gw_today_from(T, Om_ann):
    return OMEGA_R_H2 * (gstar(T) / 3.36) * (3.91 / gstars(T)) ** (4 / 3) * Om_ann


def f_peak_today(T):
    return a_of_T(T) * Hrad(T) / 6.58212e-25  # Hz (a/a0 * H, HKS convention)


def hks_today(sigma, T):
    # radiation-era scaling network: Omega_gw(t) = 8 pi eps G^2 A^2 sigma^2 / (3 H^2), G = 1/(8 pi MP^2)
    G = 1 / (8 * np.pi * MP**2)
    Om = 8 * np.pi * eps_gw * G**2 * A**2 * sigma**2 / (3 * Hrad(T) ** 2)
    return gw_today_from(T, Om), Om


gw = {}
# (a) relic-safe bound, any v: rho_w(T_ann)/s <= rho_DM/s / f  ->  rho_w/rho_r <= (rho_DM/s)/f * (4/3) gs/g / T
for f in (1.0, 0.1, 0.01):
    rows = []
    for T in [3e-4, 1e-3, 1e-2, 0.1, 1.0, 100.0, 1e4]:
        x = rho_over_s_DM / f * (4 / 3) * gstars(T) / gstar(T) / T
        Om_ann = 3 * eps_gw / (32 * np.pi) * x**2
        rows.append(dict(T_ann_GeV=T, rho_w_over_rho_r_max=x, Omega_gw_h2_today_max=gw_today_from(T, Om_ann), f_peak_Hz=f_peak_today(T)))
    gw['relic_safe_bound_f=%g' % f] = rows
# (b) benchmark tension at T_safe band and at formation-scale decays
bench = []
for lab, T in [('T_safe f=1 (30.6 MeV) [Q]', 0.0306), ('T_safe f=0.1 (3.3 MeV) [Q]', 3.3e-3), ('T_safe f=0.01 (0.31 MeV) [Q]', 3.1e-4),
               ('charge-forced exit ~ T_n(r=0.06) ~ 500 TeV', 5e5), ('free-streaming gas exit > 2.6e7 GeV', 2.6e7)]:
    g0, Oa = hks_today(sigma_wall, T)
    bench.append(dict(case=lab, T_ann_GeV=T, Omega_gw_h2_today=g0, f_peak_Hz=f_peak_today(T),
                      note='decays at formation are sub-horizon: Omega lower by (H R)^2-ish, frequency higher by 1/(H R)' if T > 1e3 else ''))
gw['benchmark_v1TeV_HKS'] = bench
gw['scaling'] = '[D] Omega_gw(T_ann) ~ (3 eps/(32 pi)) (rho_w/rho_r)^2; relic safety caps rho_w/rho_r independent of v, so Omega_gw h^2 <~ 1.6e-19 f^-2 (MeV/T_ann)^2 (g* factors aside)'
gw['detector_reference'] = '[R] NANOGrav-15 signal ~1e-9-1e-8 at 2e-9-3e-8 Hz; SKA-PTA ~1e-13-1e-12; LISA ~1e-13-1e-12 at mHz; ET/CE ~1e-13-1e-12 at 10-100 Hz; BBO/DECIGO ~1e-17-1e-16 at 0.1-1 Hz'
gw['fragmentation_GW_cold_class'] = dict(rho_h_over_rho_r_at_T_Q=neff[4.4e5]['rho_h_over_rho_rad'],
                                         Omega_gw_upper='< (rho_h/rho_r)^2 ~ 1e-25 at emission, times (H/m)^2 suppression: unobservable [E]')
out['8_gravitational_waves'] = gw

# ================================================================ 9. PBH / late minority domains and H_inf calibration
As = 2.1e-9  # [RETRO]
pbh = dict(
    quoted_minority_relic='[Q] late minority domains re-entering after T_dom = 7.97 keV dominate their patch: black-hole/baby-universe candidates of >= 4 pi M_P^4/sigma = 4.2e9 Msun (v/TeV)^-3; none per Hubble volume if H_inf <~ 206-258 GeV; mass density <= Omega_DM if <~ 301-533 GeV; cold class allows H_inf up to 314-560 GeV (f=1)',
    M_min_Msun_formula=4 * np.pi * MP**4 / sigma_wall / MSUN_GEV,
    tensor_to_scalar_at_H_inf={str(H): 2 * H**2 / (np.pi**2 * MP**2 * As) for H in (150.0, 250.0, 560.0)},
    note='Abundance changes ~1 decade per 1-3 % in H_inf [Q]; H_inf cannot be calibrated independently (r_t ~ 1e-25); BH fraction not computed; signal generic to any exact-Z2 wall sector with the same sigma')
out['9_PBH_minority_domains'] = pbh

# ================================================================ 10. portal thermalization tension (for the extension discussion)
Tc = 1.56 * v_b
kappa_th = np.sqrt(Hrad(Tc) / (2e-3 * Tc))
out['10_portal_estimate'] = dict(
    label='[E] Higgs portal kappa |H|^2 phi^2 (a versioned change of a matching choice, NOT adopted): SM -> hidden rate ~ 2e-3 kappa^2 T',
    kappa_thermalizes_at_Tc=kappa_th,
    statement='A portal large enough for direct/collider detection (kappa ~ 1e-3-1e-2 [R]) thermalizes the hidden sector above T_c if T_RH > T_c ~ 1.56 v, i.e. produces the rejected H-star-type history unless reheating stays below T_c')

json.dump(out, open(os.path.join(RES, 'reach.json'), 'w'), indent=1, default=float)
print(json.dumps({k: (v if k in ('1_calibration_allDM',) else '...') for k, v in out.items()}, indent=1, default=float))
