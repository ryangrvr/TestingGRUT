#!/usr/bin/env python3
"""m1_critic quick checks (own code; algebra and 1D quadrature only).

Settles the contradictions between the four M1 tracks and their skeptics where a cheap, engine-free
calculation can do so, and fills the small gaps the tracks left (94f660d numbers nobody checked).
Everything here is a check of one supplied model; nothing validates or excludes it.

Items
  A  Option 1A holder: effective scale Lambda_eff = M_Y/xi and EFT ratios (holder claim 3 contradiction)
  B  94f660d numbers not covered by any track: H/m_X at the 20 MeV clock; -0.567645 % refit; FOLLOWUP -1.186 %/-1.74 %
  C  Holder late-charge exchange partner: Gauss-constraint arithmetic on the skeptic's saved window means
     (data quoted from tracks/m1_holder/skeptic/results/skeptic_evolve.json; the arithmetic is mine)
  D  Core interior for a slow X: Bogoliubov dispersion of the charged condensate at the Q = 2e5 centre
     (transport T13 mechanism contradiction)
  E  Formation-clock scaling of the free-streaming length (transport T11 contradiction), own background [E]
  F  Supplied target D = rho_DM/s versus the Planck-2018 value [R]
  G  Frozen mass test: coefficients, charge-conserving pulse rows, pt_m1tilt coherent inventories, and the
     H_i recalibration that would absorb each population (what 'failure' means)
  H  Checkpoint kurtosis <p^4>/<p^2>^2 of Codex's 96-bin gas (reproduction from the npz if present)
"""
import json, math, os, sys
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'critic_checks.json')
SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad'
R = {}

# ---------------------------------------------------------------- constants (GeV)
H_i = 134284273491.30815          # CANDIDATE_M1.json
M_Y = 4028528204739.244
mu = 961402368973.181
h_vev = 5371370939652.326
kappa = 0.0034586847211042523
lamX = 1e-3
xi = 30.746728738599167
alpha = 2.192144074179184
m_X = 0.09240896373753621
v = 0.924089637375362
g = 0.1
a_q = 0.1                         # lambda_X/g^2
M_red = 2.435e18
M_obj = 17324.832346421026
Q_ref = 200000
D = 4.4e-10

# ---------------------------------------------------------------- A. effective scale
Lam_eff = M_Y / xi
amps = {'X': 87.19556213287578, 'Y': 5.75424796328559, 'Sigma': 28.387181944361622}  # H_i units, f8186fa
R['A_effective_scale'] = {
    'Lambda_eff_over_H_i': Lam_eff / H_i, 'Lambda_eff_GeV': Lam_eff,
    'mu_over_H_i': mu / H_i, 'M_Y_over_H_i': M_Y / H_i,
    'amplitude_over_Lambda_eff': {k: val * H_i / Lam_eff for k, val in amps.items()},
    'amplitude_over_M_Y': {k: val * H_i / M_Y for k, val in amps.items()},
    'R_over_Lambda_eff2_at_R_12H2': 12.0 * H_i**2 / Lam_eff**2,
    'coefficient_relative_to_1e15_cutoff': xi * 1e15 / M_Y,
    'holder_term_over_mu_at_R12': abs(complex(math.cos(alpha), math.sin(alpha)) * xi / M_Y * 12 * H_i**2) / mu,
    'reading': 'the dimension-5 holder has its own scale M_Y/xi ~ H_i, below every field amplitude; the '
               '"margin 5" of m1_holder claim 3 is measured against M_Y, not against the operator scale; '
               'the skeptic\'s "not established" stands',
}

# ---------------------------------------------------------------- B. 94f660d numbers
gstar = 10.75
T20 = 0.020
H20 = math.sqrt(math.pi**2 * gstar / 90.0) * T20**2 / M_red
dM = 0.0085754931
R['B_94f660d'] = {
    'H_at_20MeV_GeV': H20, 'H_over_mX_at_20MeV': H20 / m_X, 'codex': 1.93010e-21,
    'refit_for_dM_0.0085754931': (1 + dM)**(-2.0/3.0) - 1, 'codex_refit': -0.00567645,
    'refit_fully_cooled_0.0084099365': (1 + 0.008409936530420273)**(-2.0/3.0) - 1, 'codex_fully_cooled': -0.00556762,
    'holder_refit_fixed_mass_per_charge': 1.0180645288398746**(-2.0/3.0) - 1, 'codex_followup': -0.0118646,
    'holder_plus_fully_cooled': (1.0180645288398746 * (1 + 0.008409936530420273))**(-2.0/3.0) - 1, 'codex_followup_combined': -0.0174,
    'delta_lna_over_20000_mX': 20000 * H20 / m_X, 'codex_delta_lna': 3.8602e-17,
}

# ---------------------------------------------------------------- C. holder exchange partner (Gauss arithmetic)
sk = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'm1_holder', 'skeptic', 'results', 'skeptic_evolve.json')
C = {'source': 'tracks/m1_holder/skeptic/results/skeptic_evolve.json window means [Q]; arithmetic own'}
try:
    d = json.load(open(sk))
    for tag in ('original', 'holder'):
        w = d[tag]
        dNX = w['NX_800_1000'] - w['NX_640_800']
        dNS6 = 6 * (w['NS_800_1000'] - w['NS_640_800'])
        dNY3 = 3 * (w['NY_800_1000'] - w['NY_640_800'])
        gauss_640 = w['NX_640_800'] + 3 * w['NY_640_800'] + 6 * w['NS_640_800']
        gauss_800 = w['NX_800_1000'] + 3 * w['NY_800_1000'] + 6 * w['NS_800_1000']
        C[tag] = {'dNX_640to800window': dNX, '6dN_Sigma': dNS6, '3dN_Y': dNY3,
                  'dNX_plus_6dNS': dNX + dNS6, 'share_carried_by_Sigma': -dNS6 / dNX if dNX else None,
                  'share_carried_by_Y': -dNY3 / dNX if dNX else None,
                  'gauss_sum_640_800': gauss_640, 'gauss_sum_800_1000': gauss_800}
    rr = d['comparison']['ratio_running_80_windows']
    C['ratio_running_windows_min_max'] = [rr['min'], rr['max']]
    C['ratio_half_spread_percent'] = 100 * (rr['max'] - rr['min']) / 2
    C['reading'] = ('total gauge charge N_X + 3N_Y + 6N_Sigma is zero in every window; the window-to-window change of N_X '
                    'is carried entirely by 6N_Sigma (Y carries < 1e-3 of it), so the slow late drift is X <-> Sigma '
                    'exchange, as the skeptic found, not light-Y quanta as m1_holder claims 11/20 say; the "+1.806 %" '
                    'is defined only to +-0.17 % at tau = 800')
except Exception as e:  # pragma: no cover
    C['error'] = repr(e)
R['C_holder_exchange_partner'] = C

# ---------------------------------------------------------------- D. core interior dispersion
S_c, F_c, Om = 0.1257717809745518, 2.6309605904405573, 0.8817204919524214
k_lab = 0.0684
nu = math.sqrt(1 + k_lab**2) - Om              # rotating-frame frequency of a slow free X
aF2 = a_q * F_c**2
# uniform-condensate Bogoliubov relation in the rotating frame (own derivation):
# (k^2 - nu^2)(k^2 - nu^2 + 2 a F^2) = 4 Omega^2 nu^2
z = (-2 * aF2 + math.sqrt(4 * aF2**2 + 16 * Om**2 * nu**2)) / 2      # z = k_in^2 - nu^2 (positive root)
k_in2 = z + nu**2
cs2 = aF2 / (2 * Om**2 + aF2)
diag_barrier = S_c**2 + 2 * aF2 - 1 + 0*k_lab**2
R['D_core_interior'] = {
    'centre_S': S_c, 'centre_F': F_c, 'Omega': Om, 'k_lab_over_mX': k_lab,
    'rotating_frame_nu': nu, 'aF2_centre': aF2,
    'uniform_condensate_dispersion': '(k^2 - nu^2)(k^2 - nu^2 + 2aF^2) = 4 Omega^2 nu^2; gapless branch nu -> 0 as k -> 0',
    'sound_speed_squared': cs2, 'skeptic': 0.31,
    'k_in_squared': k_in2, 'k_in': math.sqrt(k_in2) if k_in2 > 0 else None, 'skeptic_k_in': 0.22,
    'reviewer_diagonal_barrier_S2_plus_2aF2_minus_1': diag_barrier,
    'gapped_branch_nu2_at_k0': 4 * Om**2 + 2 * aF2,
    'reading': 'with the anomalous (off-diagonal) coupling kept, a slow X propagates in the core interior as a '
               'Bogoliubov phonon with k_in ~ 0.22 m_X (sound speed^2 ~ 0.31); the reviewer\'s "classically forbidden" '
               'premise (diagonal barrier +0.40) is wrong, the skeptic\'s mechanism is right; the geometric size of the '
               'cross section (160-265 hard sphere vs 332 coupled solve) is agreed within a factor 1.3-2 and the '
               'Gamma/H ~ 1e6 conclusion stands under both; the uniform relation is approximate at the centre (F\'\' != 0)',
}

# ---------------------------------------------------------------- E. clock scaling of the free-streaming length [E]
h0 = 0.674
Om_m_h2, Om_g_h2 = 0.143, 2.474e-5
Om_r_h2 = Om_g_h2 * (1 + 0.2271 * 3.044)
H0 = 100 * h0 * 1e3 / 3.0857e22          # s^-1
H0_GeV = H0 * 6.582119569e-25
T0 = 2.3486e-13                          # GeV (2.7255 K)
def gs(T):                                # crude g_*s table [E]
    if T < 5e-4: return 3.91
    if T < 0.1: return 10.75
    if T < 0.17: return 17.25
    if T < 1.0: return 61.75
    if T < 4.0: return 75.75
    return 86.25
def a_of_T(T):
    return T0 / T * (3.91 / gs(T))**(1.0/3.0)
def H_early(T):
    return math.sqrt(math.pi**2 * gs(T) / 90.0) * T**2 / M_red
a_eq = Om_r_h2 / Om_m_h2
def H_late(a):
    return H0_GeV * np.sqrt(Om_r_h2 / h0**2 * a**-4.0 + Om_m_h2 / h0**2 * a**-3.0)
GeVinv_to_pc = 1.97327e-16 / 3.0857e16   # GeV^-1 -> m -> pc
def lam_fs(Tf, pf, T_switch=2e-3, N=20000):
    """comoving distance to equality for p(a) = pf a_f/a, released at T_f [E]."""
    af = a_of_T(Tf)
    tot = 0.0
    # early part (T_f -> T_switch) in ln T with the g* table
    if Tf > T_switch:
        lnT = np.linspace(math.log(Tf), math.log(T_switch), N)
        T = np.exp(lnT)
        a = np.array([a_of_T(t) for t in T]); Hh = np.array([H_early(t) for t in T])
        p = pf * af / a; vv = p / np.sqrt(1 + p**2)
        # d = int v dt / a = int v da/(a^2 H); da = a dlna; dlna = -dlnT - (1/3) dln g_s (ignore g_s steps in da: use a directly)
        integrand = vv / (a * Hh)
        tot += -np.trapezoid(integrand, np.log(a))       # integrate in ln a (monotone), sign from decreasing a
        a_start = a_of_T(T_switch)
    else:
        a_start = af
    lna = np.linspace(math.log(a_start), math.log(a_eq), N)
    a = np.exp(lna); p = pf * af / a; vv = p / np.sqrt(1 + p**2)
    tot += np.trapezoid(vv / (a * H_late(a)), lna)
    return tot * GeVinv_to_pc * (1 + 0*af), af
pf = 0.684
rows = {}
for Tf in (2e-3, 2e-2, 2e-1, 2.0):
    lam, af = lam_fs(Tf, pf)
    p2 = pf * af / a_of_T(2e-3) if Tf > 2e-3 else pf
    rows['T_f_%g_MeV' % (Tf * 1e3)] = {'lambda_to_eq_pc': lam, 'p_at_2MeV_over_mX': p2,
                                       'v_at_2MeV': p2 / math.sqrt(1 + p2**2)}
base = rows['T_f_20_MeV']['lambda_to_eq_pc']
for k_ in rows: rows[k_]['ratio_to_20MeV'] = rows[k_]['lambda_to_eq_pc'] / base
R['E_clock_scaling'] = {
    'assumption': 'gas released at T_f with p_f = 0.684 m_X (gives p = 0.0684 at 2 MeV for T_f = 20 MeV); p redshifts as 1/a; '
                  'crude g_*s table above 2 MeV, standard late background below; z_eq = %.0f' % (1 / a_eq - 1),
    'rows': rows,
    'reviewer_claim': 'lambda_fs ~ T_f (200 MeV -> 366 pc, 2 MeV -> 3.7 pc)',
    'skeptic_claim': 'lambda_fs ~ 1/T_f (2 MeV -> 394 pc, 20 MeV -> 46 pc, 200 MeV -> 4.6 pc, 2 GeV -> 0.33 pc)',
    'reading': 'momentum redshifts, so an earlier (hotter) release gives a colder gas at 2 MeV and a SHORTER free-streaming '
               'length: lambda ~ 1/T_f up to logarithms; the skeptic is right, m1_transport T11/clock_scaling is CORRECTED '
               '(the conclusion "the scale is an input" survives; the direction does not). Absolute numbers here are [E] '
               '(crude background, no e+- annihilation detail); the 20 MeV row is within 15 % of the two tracks\' 36-46 pc',
}

# ---------------------------------------------------------------- F. target D versus Planck [R]
Oc_h2, sig = 0.1200, 0.0012
rho_c_h2 = 1.05368e-5    # GeV cm^-3
s0 = 2891.2              # cm^-3
D_planck = Oc_h2 * rho_c_h2 / s0
R['F_target_D'] = {'supplied_D_GeV': D, 'planck2018_D_GeV_recalled': D_planck,
                   'supplied_over_planck_minus_1': D / D_planck - 1, 'planck_rel_band_1sigma': sig / Oc_h2,
                   'lambda_shift_if_D_at_planck': (D_planck / D)**(2.0/3.0) - 1,
                   'reading': 'the +0.7009 % calibration step is comparable to the rounding of the supplied target (+0.6 %) '
                              'and inside its +-1 % band; the eight quoted digits of lambda are not meaningful (ledger skeptic N1/N9)'}

# ---------------------------------------------------------------- G. frozen mass test and what failure means
eps_Q = M_obj / Q_ref
c_eta = m_X / eps_Q - 1
c_p = 2 * m_X / eps_Q
def term(eta, p): return c_eta * (1 - eta) + c_p * p
def refit(eta, p): return (1 + term(eta, p))**(-2.0/3.0)
G = {'c_eta': c_eta, 'c_p': c_p, 'eta_min_p0': 1 - 0.01 / c_eta, 'p_max_eta1': 0.01 / c_p,
     'codex': {'c_eta': 0.06678046736338095, 'c_p': 2.133560934726762, 'eta_min': 0.8502556152297386, 'p_max': 0.004686999952631147}}
# charge-conserving pulse rows (skeptic N3): eta = 1 + (Nbar - Nplus)/Q0; phi quanta as future pairs
pulses = {'5%': (0.002088 * Q_ref, 17.0, 0.86), '10%': (1689.506, 89.356, 5.503), '20%': (0.03447 * Q_ref, 824.7, 53.7), '40%': (0.1459 * Q_ref, 10596.7, 589.5)}
prow = {}
for k_, (Nbar, Nplus, Nphi) in pulses.items():
    eta = 1 + (Nbar - Nplus) / Q_ref
    p = (Nbar + Nphi) / Q_ref
    prow[k_] = {'eta': eta, 'p': p, 'term': term(eta, p), 'ratio_to_0.01': term(eta, p) / 0.01, 'H_i_refit_factor': refit(eta, p)}
G['pulse_rows_charge_conserving'] = prow
# populations and the recalibration that would absorb them
pops = {
    'M1 matching hypothesis (eta 1, p 0)': (1.0, 0.0),
    'kin_solve gas-only R=5 (p = f/2.02855, f 0.0047584)': (1.0, 0.0047584 / 2.02855),
    'kin_solve gas-only R=7 (f 0.0288444)': (1.0, 0.0288444 / 2.02855),
    'pt_relax equilibrium Xbar/charge 0.07 (no annihilation)': (1.0, 0.07),
    'pt_relax equilibrium 0.11': (1.0, 0.11),
    'pt_m1tilt coherent tau-20 orbit, slow limit (1.49 Xbar + 0.18 phi -> 1.67 pairs/charge, no annihilation)': (1.0, 1.67),
    'pt_m1tilt coherent late orbit, slow limit (0.11-0.15 Xbar + 0.22-0.25 phi -> 0.35)': (1.0, 0.35),
    'pt_m1tilt adiabatic branch tau-20 (J_r = 1.76)': (1.0, 1.76),
    'ed73791 dilute masks 0.30-0.32 read as eta': (0.31, 0.0),
    'ed73791 dense masks 0.494 read as eta': (0.494, 0.0),
    '10% pulse kinetics tau=1000 (eta 0.94186, p 2.8e-6; frozen eps_Q)': (0.94186, 2.8e-6),
    'all charge in free X (eta 0, p 0)': (0.0, 0.0),
}
rows = {}
for name, (eta, p) in pops.items():
    t = term(eta, p)
    rows[name] = {'eta': eta, 'p': p, 'mass_over_D': 1 + t, 'excess': t, 'passes_1pct': t <= 0.01,
                  'H_i_refit_factor_to_absorb': refit(eta, p),
                  'pair_mass_fraction_of_total': (c_p * p) / (1 + t)}
G['populations'] = rows
G['survival_needed_for_pt_relax_0.07_to_pass'] = (0.01 / c_p) / 0.07
G['survival_needed_for_pt_m1tilt_1.67_to_pass'] = (0.01 / c_p) / 1.67
G['reading'] = ('the 1 % test is a benchmark on the calibration, not a falsification: any (eta, p) is absorbed by the H_i refit '
                'factor [1 + c_eta(1-eta) + c_p p]^(-2/3) because the dimensionless history is unchanged (Y ~ H_i^3/2). What a '
                'large p would change is the character of the dark matter: at p = 0.07 the symmetric pairs are 13 % of the mass, '
                'at p = 0.35 (pt_m1tilt late orbit) 43 %, at p = 1.67 (tau-20 slow limit) 78 %: the "asymmetric charged objects" '
                'picture is then gone unless annihilation at the actual heat removes them (kin_solve: possible at R <~ 5.8, '
                'e = 1e-6 band). eta is uncomputed; no mask fraction is eta.')
R['G_frozen_test'] = G

# ---------------------------------------------------------------- H. checkpoint kurtosis (reproduction from Codex data if present)
Hd = {'status': 'npz not found'}
for cand in (os.path.join(SP, 'review', 'm1_transport', 'codex', 'late_kin_h0025_96.npz'),
             os.path.join(SP, 'review', 'm1_transport_skeptic', 'codex', 'late_kin_h0025_96.npz')):
    if os.path.exists(cand):
        try:
            z = np.load(cand, allow_pickle=False)
            Hd = {'file': cand, 'keys': list(z.keys())[:40]}
            # try to find an energy grid and an occupation array for the positive species
            keys = list(z.keys())
            def pick(names):
                for n in names:
                    if n in keys: return z[n]
                return None
            E = pick(['E', 'energy', 'E_grid', 'energies', 'E_nodes', 'e_grid', 'eps'])
            fpos = pick(['Nplus', 'n_plus', 'f_plus', 'fplus', 'plus', 'N_plus', 'pos'])
            if E is not None and fpos is not None and E.shape == fpos.shape:
                E = np.asarray(E, float); fpos = np.asarray(fpos, float)
                # energies in m_X units; if E is kinetic use p^2 = E^2 - 1 when E >= 1 else p^2 = 2E
                if E.min() >= 1.0: p2 = E**2 - 1
                else: p2 = 2 * E
                w = fpos
                p2m = np.sum(w * p2) / np.sum(w); p4m = np.sum(w * p2**2) / np.sum(w)
                Hd.update({'kurtosis_p4_over_p2sq': p4m / p2m**2, 'codex_initial': 1.91713, 'v_rms': math.sqrt(p2m),
                           'codex_v_rms': 0.06838, 'energy_convention_guess': 'total' if E.min() >= 1 else 'kinetic'})
            else:
                Hd['note'] = 'array names not recognised; see keys'
                Hd['shapes'] = {k: list(z[k].shape) for k in keys[:40]}
        except Exception as e:
            Hd = {'file': cand, 'error': repr(e)}
        break
R['H_checkpoint_kurtosis_repro'] = Hd

json.dump(R, open(OUT, 'w'), indent=1, default=float)
print(json.dumps(R, indent=1, default=float))
