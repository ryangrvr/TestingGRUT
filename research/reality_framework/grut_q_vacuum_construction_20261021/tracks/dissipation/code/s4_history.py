"""Track `dissipation`, script 4: late-time behaviour of the KV dissipative attractor and compatibility with a
standard thermal history; residual vacuum energy; late-time observables.  Physical units (GeV), brief conventions.

Law used (KV 1601.00601 eq (5.7), constant G, translated; s1 checks B, E):
    d rho_V/dt = -(gamma/36) |H| R^2,  R = (rho_m - 3 p_m + 4 rho_V)/M^2,  the energy goes into massless quanta.
    => per e-fold N = ln a:  d rho_V/dN = -(gamma/36) R^2   (H > 0).
Inputs (provenance in ../NOTES.md section 7): M, rho_Lambda from BRIEF.md; gamma = 1/(32 pi^2) from 1601.00601 (5.5);
H0, Omega_m, t0, T_CMB, Omega_gamma h^2, g* values recalled (Planck 2018 / PDG) [R]; Delta_EW, Delta_QCD from
../../equations/NOTES.md section 7; FIRAS and N_eff limits from ../../ledger/results/ledger.json.
Output: ../results/history.json.  Run: python3 -I s4_history.py  (seconds)
"""
import json, os
import numpy as np
from scipy.integrate import solve_ivp

out = {'labels': '[C] computed here; [D] derived (formula in NOTES); inputs [Q]/[R] as listed under inputs'}
# ------------------------------------------------------------------ inputs
M = 2.435e18                    # GeV, reduced Planck mass (BRIEF)
rhoL = 2.52e-47                 # GeV^4 (BRIEF)
gam = 1.0 / (32 * np.pi**2)     # 1601.00601 eq (5.5)
hbar = 6.582119569e-25          # GeV s [R]
Gyr = 3.15576e16                # s
Mpc_km = 3.0856775814913673e19  # km per Mpc
H0 = 67.4 / Mpc_km * hbar     # GeV  (H0 [1/s] * hbar), Planck 2018 [R]
Om, OL = 0.315, 0.685           # Planck 2018 [R]
t0 = 13.787 * Gyr / hbar        # GeV^-1  (Planck 2018 age [R])
T_cmb = 2.7255 * 8.617333262e-14  # GeV [R]
Og_h2 = 2.473e-5                # Omega_gamma h^2 [R]
h_ = 0.674
D_EW = 1.19e8                   # GeV^4 (equations track section 7: m_H^2 v^2/8)
D_QCD = (4.4e-4, 3.0e-3)        # GeV^4 (equations track section 7)
Estar = np.sqrt(3) * M
inputs = {'M_GeV': M, 'rho_Lambda_GeV4': rhoL, 'gamma': gam, 'H0_GeV': H0, 'Omega_m': Om, 't0_GeVinv': t0,
          'T_CMB_GeV': T_cmb, 'Omega_gamma_h2': Og_h2, 'Delta_EW_GeV4': D_EW, 'Delta_QCD_GeV4': D_QCD,
          'E_star=sqrt(3)M_GeV': Estar}
out['inputs'] = inputs
rho_crit = 3 * M**2 * H0**2


def grho(T, g):
    return np.pi**2 / 30 * g * T**4


def Hrad(T, g):
    return np.sqrt(grho(T, g) / (3 * M**2))


# ------------------------------------------------------------------ (a) exact radiation relaxation law
N_needed = 9 * M**4 / (4 * gam * rhoL)
gs0 = 3.91  # entropy dof today [R]
N_std = np.log((M * (106.75 / gs0)**(1 / 3)) / T_cmb)   # e-folds from T = M (g*s=106.75) to T_CMB today
floor_std = 9 * M**4 / (4 * gam * N_std)
out['a_radiation_law_[D]'] = {
    'formula': 'rho_V(N) = rho_Vi / (1 + (4 gamma/9)(rho_Vi/M^4) N)  (exact for radiation + vacuum, constant G)',
    'floor_formula': 'rho_V(N) >= (1/2) min(rho_Vi, 9 M^4/(4 gamma N))',
    'efolds_needed_to_reach_rho_Lambda': N_needed,
    'efolds_T=M_to_T_CMB_standard_history': N_std,
    'floor_after_standard_history_GeV4': floor_std, 'floor_in_M^4': floor_std / M**4,
    'floor_in_rho_Lambda': floor_std / rhoL,
    'relative_change_if_prepared_at_rho_r(1MeV)': (4 * gam / 9) * grho(1e-3, 10.75) / M**4 * N_std,
}

# ------------------------------------------------------------------ (b) numerical history from T = 1e18 GeV
def gstar(T):   # crude step model [R]
    return 106.75 if T > 100 else (61.75 if T > 0.2 else (17.25 if T > 0.1 else (10.75 if T > 5e-4 else 3.36)))


def gstar_of_rr(rr):
    T = (max(rr, 1e-300) / (np.pi**2 / 30 * 106.75))**0.25
    for _ in range(3):
        T = (max(rr, 1e-300) / (np.pi**2 / 30 * gstar(T)))**0.25
    return gstar(T)


def T_of_rr(rr):
    return (max(rr, 1e-300) / (np.pi**2 / 30 * gstar_of_rr(rr)))**0.25


def Vm(T):
    """SM vacuum energy relative to its final T=0 value (matter sector): EW and QCD plateaux, smoothed."""
    def step(x):   # C^1 step, exactly 0 for x < -1 and 1 for x > 1 (no exponential tails)
        s_ = min(max(0.5 * (x + 1.0), 0.0), 1.0)
        return s_ * s_ * (3.0 - 2.0 * s_)
    return D_EW * step((T - 160.0) / 30.0) + D_QCD[0] * step((T - 0.155) / 0.03)


def history(rhoVi, kappa_hot, Ti=1e18, T_nr=100.0):
    """Integrate in N = ln a: y = [rho_V (q-sector), rho_r (SM radiation + produced quanta), rho_d (dust), t].
    Dust (CDM + baryons, Omega_m) appears at T_nr with comoving-entropy normalisation; kappa_hot = trace fraction of
    the plasma for T > 1 MeV.  The latent heat of V_m goes into radiation.  Stops at turnaround (rho_tot <= 0)."""
    gi = 106.75
    rr0 = grho(Ti, gi)
    rd_nr = Om * rho_crit * (106.75 / gs0) * (T_nr / T_cmb)**3

    def rhs(N, y):
        rV, rr, rd, t = y
        T = T_of_rr(rr)
        kap = kappa_hot if T > 1e-3 else 0.0
        vm = Vm(T)
        rho_tot = rV + rr + rd + vm
        H = np.sqrt(max(rho_tot, 1e-300) / (3 * M**2))
        R = (kap * rr + rd + 4 * (rV + vm)) / M**2
        drain = gam / 36 * R * R
        dVm = (Vm(T * (1 + 1e-6)) - Vm(T)) / (T * 1e-6) * (-T)
        return [-drain, -4 * rr + drain - dVm, -3 * rd, 1.0 / H]

    def e_turn(N, y):   # H^2 has fallen to 1e-4 of the radiation term: turnaround imminent
        return y[0] + y[1] + y[2] + Vm(T_of_rr(y[1])) - 1e-4 * y[1]
    e_turn.terminal = True; e_turn.direction = -1

    def e_nr(N, y):
        return y[1] - grho(T_nr, 106.75)
    e_nr.terminal = True; e_nr.direction = -1

    Hi = np.sqrt((rhoVi + rr0 + Vm(Ti)) / (3 * M**2))
    y0 = [rhoVi, rr0, 0.0, 1.0 / (2.0 * Hi)]
    Nmax = np.log(Ti / T_cmb) + 3
    kw = dict(method='LSODA', rtol=1e-9, atol=[1e-300, 1e-300, 1e-300, 1e-60], dense_output=True)
    segs = []
    s1 = solve_ivp(rhs, (0, Nmax), y0, events=[e_turn, e_nr], **kw)
    segs.append(s1)
    reason = 'reached N_max'
    if s1.status == 1 and len(s1.t_events[1]):
        y1 = s1.y[:, -1].copy(); y1[2] = rd_nr
        s2 = solve_ivp(rhs, (s1.t[-1], Nmax), y1, events=[e_turn], **kw)
        segs.append(s2)
        if s2.status == 1:
            reason = 'turnaround (crunch follows)'
    elif s1.status == 1:
        reason = 'turnaround (crunch follows)'
    return segs, reason


def sample(segs, Nq):
    for s in segs:
        if s.t[0] <= Nq <= s.t[-1]:
            return s.sol(Nq)
    return None


hist = {}
cases = (('formed? rhoVi=M^4, kappa=0', M**4, 0.0),
         ('formed? rhoVi=M^4, kappa=1e-2', M**4, 1e-2),
         ('rhoVi=1e-60 M^4, kappa=0', 1e-60 * M**4, 0.0),
         ('prepared rhoVi=rho_Lambda, kappa=0', rhoL, 0.0),
         ('relaxed-to-zero total before EW: rhoVi=-(D_EW+D_QCD), kappa=0', -(D_EW + D_QCD[0]), 0.0),
         ('prepared rhoVi=rho_Lambda, kappa=1e-3 from T=1e18', rhoL, 1e-3))
for name, rVi, kap in cases:
    segs, reason = history(rVi, kap)
    Ns = np.concatenate([np.linspace(s.t[0], s.t[-1], 6001) for s in segs])
    Ys = np.concatenate([s.sol(np.linspace(s.t[0], s.t[-1], 6001)) for s in segs], axis=1)
    Tr = np.array([T_of_rr(x) for x in Ys[1]])
    rows = []
    for Tm in (1e16, 1e3, 160.0, 1.0, 1e-3, 1e-9, T_cmb):
        if Tr.min() > Tm * 1.5:
            rows.append({'T_target_GeV': Tm, 'not_reached': True})
            continue
        idx = int(np.argmin(abs(np.log(Tr / Tm))))
        yy = Ys[:, idx]
        vm = Vm(Tr[idx])
        rho_tot = yy[0] + yy[1] + yy[2] + vm
        H = np.sqrt(max(rho_tot, 1e-300) / (3 * M**2))
        rows.append({'T_rad_GeV': float(Tr[idx]), 'N': float(Ns[idx]), 'rho_V(q-sector)_GeV4': float(yy[0]),
                     'V_m(SM plateau)_GeV4': float(vm), 'Omega_vac_total': float((yy[0] + vm) / rho_tot),
                     'H_GeV': float(H), 'age_GeVinv': float(yy[3]), 'H*age': float(H * yy[3]),
                     'H/H_standard(T)': float(H / Hrad(Tr[idx], gstar(Tr[idx])))})
    yend = segs[-1].y[:, -1]
    hist[name] = {'rho_V_i_GeV4': rVi, 'kappa_hot': kap, 'end': reason, 'N_end': float(segs[-1].t[-1]),
                  'T_rad_end_GeV': float(T_of_rr(yend[1])), 'rows': rows,
                  'rho_V(q)_end_GeV4': float(yend[0]),
                  'total_vacuum_end_over_rho_Lambda': float((yend[0] + Vm(T_of_rr(yend[1]))) / rhoL)}
    print(name, reason, hist[name]['total_vacuum_end_over_rho_Lambda'], flush=True)
out['b_histories_from_T=1e18GeV'] = hist
out['b_note'] = ('A standard radiation era has H*age = 1/2 and H/H_standard(T) = 1. "formed?" cases: vacuum-dominated from '
                 'the start; rho_V falls only to ~9 M^4/(4 gamma N). Prepared case: the dissipation is negligible and the '
                 'final total vacuum energy is what was prepared. Relaxed-before-EW case: the EW and QCD plateaux are '
                 'released after the relaxation, total -> -(D_EW + D_QCD) < 0, turnaround and crunch. kappa = 1e-3 '
                 'from 1e18 GeV: drain gamma kappa^2 H_i^4/32 = 3.9e62 GeV^4, so the prepared value would also have to absorb it.')

# ------------------------------------------------------------------ (c) attractor numbers at t0 and at H = H0
tau0 = Estar * t0
rV_t0 = 9 * M**4 * (6 * gam * tau0)**(-2 / 3)
H_t0 = Estar * (6 * gam * tau0)**(-1 / 3)
N_t0 = (3 * tau0 / (4 * np.sqrt(gam)))**(2 / 3)
t_H0 = M**2 / (2 * gam * H0**3)
out['c_attractor_[D/C]'] = {
    'formulae': 'rho_V = 9 M^4 (6 gamma sqrt(3) M t)^(-2/3);  H = sqrt(3) M (6 gamma sqrt(3) M t)^(-1/3);  '
                'N = (3 sqrt(3) M t/(4 sqrt(gamma)))^(2/3);  Omega_M = 1/(4N);  t(H) = M^2/(2 gamma H^3)',
    'tau0=E_star*t0': tau0, 'rho_V(t0)_GeV4': rV_t0, 'rho_V(t0)/rho_Lambda': rV_t0 / rhoL,
    'H(t0)_GeV': H_t0, 'H(t0)/H0': H_t0 / H0, 'H(t0)*t0': H_t0 * t0, 'efolds_N(t0)': N_t0,
    'Omega_M(t0)=1/(4N)': 1 / (4 * N_t0),
    'age_when_H=H0_in_units_1/H0': t_H0 * H0, 'age_when_H=H0_Gyr': t_H0 * hbar / Gyr,
    'Omega_V_when_H=H0': 1 - 1 / (4 * (3 * Estar * t_H0 / (4 * np.sqrt(gam)))**(2 / 3)),
    'observed_[R]': {'H0*t0_LCDM': None, 'Omega_m': Om, 'Omega_Lambda': OL},
}
# LCDM H0 t0 for comparison
from scipy.integrate import quad
Or = Og_h2 / h_**2 * (1 + 0.2271 * 3.046)
H0t0 = quad(lambda a: 1 / (a * np.sqrt(Om / a**3 + Or / a**4 + (1 - Om - Or))), 1e-12, 1)[0]
out['c_attractor_[D/C]']['observed_[R]']['H0*t0_LCDM'] = H0t0
out['c_attractor_[D/C]']['H0*t0_from_data_[R]'] = H0 * t0

# ------------------------------------------------------------------ (d) attractor followed by unmodelled reheating
def N_end_for(target):
    return 9 * M**4 / (4 * gam * target)


Vinf = (1e16)**4
dinf = (4 * gam / 9) * Vinf**2 / M**4 * 60
out['d_attractor_then_reheating_[D/C]'] = {
    'requirement': 'rho_V,end = 9 M^4/(4 gamma N_end) = rho_Lambda + Delta_EW + Delta_QCD + D_trace + D_inflaton',
    'N_end_if_only_EW_QCD': N_end_for(rhoL + D_EW + D_QCD[0]),
    'fractional_tuning_of_N_end_needed': rhoL / (D_EW + D_QCD[0]),
    'H_at_end_GeV_if_rhoV_end=D_EW': np.sqrt(D_EW / (3 * M**2)),
    'radiation_fraction_at_end_1/(4N_end)': 1 / (4 * N_end_for(D_EW)),
    'radiation_density_at_end_GeV4': D_EW / (4 * N_end_for(D_EW)),
    'rho_r_needed_at_T=160GeV_GeV4': grho(160.0, 106.75),
    'inflaton_plateau_drain_[D]': {'formula': 'Delta rho_V = -(4 gamma/9)(V^2/M^4) N_inf  (V >> rho_V)',
                                   'V^(1/4)_GeV': 1e16, 'N_inf': 60, 'Delta_rho_V_GeV4': -dinf,
                                   'in_rho_Lambda': -dinf / rhoL},
}

# ------------------------------------------------------------------ (e) trace drain in the hot era
Dtr = {}
for TRH in (1e16, 1e9, 1e3, 160.0):
    for kap in (1e-3, 1e-2):
        HRH = Hrad(TRH, 106.75)
        D = gam * kap**2 * HRH**4 / 32
        Dtr['T_RH=%g,kappa=%g' % (TRH, kap)] = {'D_GeV4': D, 'D/rho_Lambda': D / rhoL}
H_qcd = Hrad(0.2, 61.75)
D_qcd_trace = gam * 0.3**2 * H_qcd**4 / 32
z_eq = 3400.0
H_eq = H0 * np.sqrt(2 * Om * (1 + z_eq)**3)
D_dust = gam / 24 * H_eq**4
out['e_trace_drains_[D/C]'] = {'formula_radiation_era': 'D = gamma kappa^2 H_RH^4 / 32 (constant kappa, H ~ a^-2)',
                               'formula_matter_era': 'D = gamma H_1^4 / 24 (from H_1 onwards, rho_V ~ 0)',
                               'radiation_era': Dtr,
                               'QCD_epoch_kappa=0.3,T=0.2GeV': {'D_GeV4': D_qcd_trace, 'D/rho_Lambda': D_qcd_trace / rhoL},
                               'matter_era_from_equality': {'H_eq_GeV': H_eq, 'D_GeV4': D_dust, 'D/rho_Lambda': D_dust / rhoL},
                               'EW_kick_scale_comparison': {'H(T=1.26TeV)^4/rho_Lambda': Hrad(1.26e3, 106.75)**4 / rhoL,
                                                            'note': 'the drain carries the KV/AHHKM H^4 scale times '
                                                                    'gamma kappa^2/32 and with negative sign'}}

# ------------------------------------------------------------------ (f) late-time observables
rho_g0 = Og_h2 / h_**2 * rho_crit
R0 = 3 * H0**2 * (Om + 4 * OL)
inj_today = gam / 36 * R0**2            # per e-fold
def inj_frac_z(z):
    Rz = 3 * H0**2 * Om * (1 + z)**3    # radiation traceless; Lambda negligible at high z
    return gam / 36 * Rz**2 / (rho_g0 * (1 + z)**4)
Hbbn = Hrad(1e-3, 10.75)
frac_bbn = gam * 1e-2**2 / 12 * (Hbbn / M)**2
frac_hot = gam * 1e-2**2 / 12 * (Hrad(1e16, 106.75) / M)**2
out['f_late_time_observables_[C]'] = {
    'produced_quanta': 'massless non-conformal fields (gravitons; KV 1601.00601 sec 5); photons are conformally '
                       'coupled at tree level and are not produced by a conformally flat FRW background [R: Parker 1969]',
    'energy_injected_per_efold_today_over_rho_gamma': inj_today / rho_g0,
    'energy_injected_per_efold_today_over_rho_Lambda': inj_today / rhoL,
    'w_plus_1_today_from_drain': (gam / 36) * R0**2 / (3 * rhoL),
    'energy_injected_per_efold_over_rho_gamma_at_z=2e6': inj_frac_z(2e6),
    'FIRAS_limits_[Q ledger row 15]': '|mu| < 9e-5, |y| < 1.5e-5',
    'graviton_fraction_per_efold_at_BBN_kappa=1e-2': frac_bbn,
    'graviton_fraction_per_efold_at_T=1e16_kappa=1e-2': frac_hot,
    'N_eff_limit_[Q ledger row 12]': 'N_nu = 2.898 +- 0.141 (BBN+CMB)',
    'G_variation': 'none at constant G (1601.00601, 2207.03453); with K\'(q) != 0 (QV1 structure) dK/K ~ H^2/m_q^2 '
                   '= 1.6e-121 today as in QV1 (../../equations/NOTES.md section 4.3)',
}

out['summary'] = {
    'efolds_needed': N_needed, 'efolds_standard': N_std,
    'rho_V(t0)_on_attractor_over_rhoL': rV_t0 / rhoL, 'H(t0)/H0_on_attractor': H_t0 / H0,
}
with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '../results/history.json'), 'w') as fh:
    json.dump(out, fh, indent=1, default=float)
print(json.dumps(out, indent=1, default=float)[:6000])
