"""pt_relax skeptic (Claude, own code): independent recomputation of the scalar claims in pt_relax NOTES.md /
results/scales.json, by hand formulas (no import from pt_relax or Codex).  Writes results/recheck_scalars.json.
"""
import json, math, os
import numpy as np
from scipy.integrate import quad

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'recheck_scalars.json')
LX, G2, LP = 1e-3, 1e-2, 3e-2
AM, C = 8.869218738806, 1.9857962750634526
HI = 133349619398.33566; V = 0.924089637375362      # GeV; from the kin_source / hw spot point used by pt_relax
NQ, RHO = 8112.665383318339, 105481.24746711778
res = {}

# --- clock: a(s) = a_m + C s from tau = 20 + a_m s + C s^2/2 (d tau = a ds); check tau(34)
tau34 = 20 + AM * 34 + C * 34 ** 2 / 2
res['clock'] = dict(a_at_s34=AM + C * 34, tau_at_s34=tau34, note='H = C/a^2 in H_i units; Gamma_phys/H = Gamma_c a/C')

# --- kappa and onset
kappa = LP * V * V / HI ** 2
on = {}
for name, r2 in [('chi60', 3600.0), ('chi110', 12100.0), ('classical_s24_34', 10104.0), ('equilibrium', 3927.19)]:
    a = math.sqrt(G2 * r2 / kappa); s = (a - AM) / C
    on[name] = dict(a_onset=a, s_onset=s, Gamma_crit=1 / s, s_over_34=s / 34, mX_conf=math.sqrt(G2) * V * a / HI,
                    Gamma_crit_over_Teq=1 / s / 17.556)
res['kappa'] = kappa; res['onset'] = on
# g_* correction to the clock: a(eta) slope changes by (g_*s^{4/3}/g_*)^{1/2}-type factors between 10^14 GeV and 45 GeV
res['onset_note'] = 'a = a_m + C s assumes constant g_*; between g_* = 106.75 and ~86 (T ~ 45 GeV) the change in s_onset is O(10%)'

# --- temperatures at the onset
for name in ['equilibrium']:
    a = on[name]['a_onset']
    T_phys = 17.556 * HI / a; mX = math.sqrt(G2) * V
    res['T_at_onset'] = dict(a=a, T_hidden_GeV=T_phys, m_X_GeV=mX, T_over_mX=T_phys / mX)
a_ref = 7.047e12
res['T_at_a_ref'] = dict(T_hidden_GeV=17.556 * HI / a_ref, T_over_mX=17.556 * HI / a_ref / (math.sqrt(G2) * V),
                         T_hidden_GeV_19p4=19.40 * HI / a_ref)
# phi vacuum energy fraction at the onset: V0 = lambda v^4/4 in conformal units (a/H_i)^4
for name, a in [('a_onset_eq', on['equilibrium']['a_onset']), ('a_ref', a_ref)]:
    V0 = LP * V ** 4 / 4 * (a / HI) ** 4
    res[f'V0_over_rho_{name}'] = V0 / RHO

# --- free massless BE brackets (own): thermal energy = rho - (3/4) lambda^(1/3) nQ^(4/3)
Ec = 0.75 * LX ** (1 / 3) * NQ ** (4 / 3)
res['condensate_energy_fraction'] = Ec / RHO
z3 = 1.2020569
for dof in (2, 3):
    T = ((RHO - Ec) * 30 / (dof * math.pi ** 2)) ** 0.25
    res[f'free_BE_{dof}dof'] = dict(T=T, Xbar_per_net_mu0=z3 * T ** 3 / math.pi ** 2 / NQ)

# --- 2<->2 X Xbar -> phi phi rate, own route: Gamma_Xbar = <n_X sigma v_Mol>, massless Boltzmann/Bose thermal
#     |M|^2 = (2 g^2)^2, sigma = |M|^2/(32 pi s) (identical phi phi), thermal average over s with f1 f2 weights.
#     Gamma = (1/n_Xbar) * Int dPi1 dPi2 f1 f2 |M|^2 Phi_2 /(...)  -> use R_{2->2} = (1/2!)_{final} Int dPi1 dPi2 f1 f2 |M|^2 Phi2(s)
#     Phi2 = 1/(8 pi) massless; Int dPi1 dPi2 f1 f2 (indep. of s) = [Int dPi f]^2 = [n/ <2E>]... do it numerically with Bose f.
def dPi_moment(power, T=17.556, bose=True):
    f = (lambda E: 1 / math.expm1(E / T)) if bose else (lambda E: math.exp(-E / T))
    return quad(lambda p: p ** 2 / (2 * math.pi ** 2) / (2 * p) * f(p) * p ** power, 1e-9, 60 * T, limit=400)[0]
I0 = dPi_moment(0); nX = quad(lambda p: p ** 2 / (2 * math.pi ** 2) / math.expm1(p / 17.556), 1e-9, 2000, limit=400)[0]
R22 = 0.5 * I0 ** 2 * (2 * G2) ** 2 / (8 * math.pi)        # events per volume per time; each removes one Xbar
res['XXbar_to_phiphi'] = dict(n_X_one_dof=nX, R_per_vol=R22, Gamma_per_Xbar=R22 / nX, efolds=R22 / nX * 2.63e12,
                              note='massless Bose thermal, T = 17.556, no final-state Bose enhancement, no thermal masses; Gamma per Xbar = R/n_Xbar')

# --- cascade time to the thermal scale with p_max ~ s^(1/7) and s^(1/5)
for name, expo in [('1/7', 1 / 7), ('1/5', 1 / 5)]:
    res[f'cascade_s_to_3T_{name}'] = 200 * (53 / 8) ** (1 / expo)

# --- Rayleigh-Jeans sanity: energy per oscillator on N32 = (rho - E_cond) V / (3 N^3)
for N in (32, 48, 64):
    res[f'RJ_T_classical_est_N{N}'] = (RHO - Ec) * 12 ** 3 / (3 * N ** 3)

# --- FD4 maximum frequency: sqrt(3 * 16/3)/dx = 4/dx = N/3
res['FD4_omega_max'] = {N: N / 3 for N in (32, 40, 48, 64)}
# Bose spectrum peaks: number k^2 f peaks at k ~ 1.59 T, energy k^3 f at 2.82 T
from scipy.optimize import minimize_scalar
pk_n = minimize_scalar(lambda x: -(x ** 2 / math.expm1(x)), bounds=(0.1, 10), method='bounded').x
pk_e = minimize_scalar(lambda x: -(x ** 3 / math.expm1(x)), bounds=(0.1, 10), method='bounded').x
res['bose_peaks_over_T'] = dict(number=pk_n, energy=pk_e, at_T17p56=dict(number=pk_n * 17.556, energy=pk_e * 17.556))
# --- Hartree gap M_X - mu at the pt_relax end point (Hugenholtz-Pines violation of the Hartree treatment)
res['hartree_gap'] = dict(M_X=2.8355, mu=2.0394, gap=2.8355 - 2.0394, gap_over_T=(2.8355 - 2.0394) / 17.556)
# --- pt_relax's 'free BE' Xbar with the condensate mu applied to the thermal Xbar branch (mu/T suppression)
T2 = res['free_BE_2dof']['T']; mu = 2.0093
n_minus_mu = quad(lambda p: p ** 2 / (2 * math.pi ** 2) / math.expm1((p + mu) / T2), 1e-9, 2000, limit=400)[0]
res['free_BE_2dof']['Xbar_per_net_with_mu'] = n_minus_mu / NQ
json.dump(res, open(OUT, 'w'), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
