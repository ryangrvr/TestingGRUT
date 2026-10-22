"""v_bound: independent recomputation of the bound track's physical numbers (NOTES A.4-A.7, B.3 AB, B.5).
Reads ../results/v3_eotwash_digitize.json (own digitization); writes ../results/v3_numbers.json.  Pure numpy, < 1 s.
Inputs and provenance are listed in OUT['inputs'].
"""
import json
import os
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
dig = json.load(open(os.path.join(here, '..', 'results', 'v3_eotwash_digitize.json')))

# ---------------------------------------------------------------- inputs
c_kms = 299792.458
Mpc_km = 3.0856775814913673e19
hbar_eVs = 6.582119569e-16            # CODATA [R]
hbarc_eVm = 1.973269804e-7            # CODATA [R]
kB_GeV = 8.617333262e-14              # GeV/K [R]
M = 2.435e18                          # GeV, reduced Planck mass (BRIEF.md)
rhoL = 2.52e-47                       # GeV^4 (BRIEF.md)
H0 = 67.4 / Mpc_km * hbar_eVs         # eV (Planck 2018, ledger row 9 via bound)
Om, OL = 0.315, 0.685
Rsun = 6.957e8                        # m (IAU nominal) [R]
T0 = 2.7255                           # K [R]
gstar, gs_osc = 106.75, 106.75        # SM above the EW scale [R]
gs0, gstar0 = 3.909, 3.363            # entropy and energy d.o.f. today (3 massless nu, T_nu/T = (4/11)^(1/3)) [R]
Om_r = 9.1e-5                          # bound's value [R]; consistent with g*0 = 3.363, h = 0.674
OUT = {'inputs': {'H0 [eV]': H0, 'Omega_m': Om, 'Omega_L': OL, 'M [GeV]': M, 'rho_L [GeV^4]': rhoL, 'R_sun [m]': Rsun,
                  'g*(T_osc)': gstar, 'g*s0, g*0': [gs0, gstar0], 'Omega_r': Om_r}}

# ---------------------------------------------------------------- A.3/A.4 Eot-Wash
pur = dig['thick purple curves (Eot-Wash 2020 / 2007 / 2004), crossings in um']
l1, l13 = pur['alpha=1'][0], pur['alpha=1/3'][0]
l1_07, l13_07 = pur['alpha=1'][1], pur['alpha=1/3'][1]
scale = 38.6 / l1
lam13 = l13 * scale
dig_unc = max(abs(1 - scale), abs(1 - 56.0 / l1_07))
mmin = hbarc_eVm / (lam13 * 1e-6)                # eV
OUT['A4_EotWash'] = {'alpha=1 digitized [um] (paper 38.6)': l1, 'alpha=1/3 digitized [um]': l13,
                     'alpha=1/3 rescaled to 38.6 [um]': lam13, 'EW2007 alpha=1 [um] (Kapner abstract 56)': l1_07,
                     'EW2007 alpha=1/3 [um]': l13_07, 'digitization uncertainty (max of the two calibrations)': dig_unc,
                     'm_q lower bound [meV]': mmin * 1e3,
                     'm_q range [meV]': [mmin * 1e3 / (1 + dig_unc), mmin * 1e3 / (1 - dig_unc)],
                     'bound track': '54.0 +- 0.8 um, 3.65 (3.60-3.71) meV'}

# ---------------------------------------------------------------- Cassini
g, s = 2.1e-5, 2.3e-5
lo = g - 2 * s
# gamma - 1 = -(2/3) y/(1 + y/3), y = e^{-mb}:  solve -(2/3) y/(1+y/3) = lo
y = -lo / (2.0 / 3.0 + lo / 3.0)
xb = -np.log(y)
OUT['A4_Cassini'] = {'2sigma lower edge': lo, 'm b required': xb, 'm_q [eV] for b = R_sun': xb * hbarc_eVm / Rsun,
                     'm_q [eV] for b = 1.6 R_sun': xb * hbarc_eVm / (1.6 * Rsun),
                     'light q: (gamma-1)/sigma': 0.5 / s}

# ---------------------------------------------------------------- A.5 Omega_dyn
Od = lambda m_eV: 2.25 * Om**2 * (H0 / m_eV)**2
rho_c = rhoL / OL
m_nat = np.sqrt(2.0 / 3.0) * M * 1e9
m_cas = xb * hbarc_eVm / Rsun
OUT['A5_Omega_dyn'] = {'Omega_dyn max (Eot-Wash)': Od(mmin), 'rho_dyn max [GeV^4]': Od(mmin) * rho_c,
                       '1 + w max = 2 Omega_dyn/Omega_L': 2 * Od(mmin) / OL,
                       'Omega_dyn at Cassini-only bound': Od(m_cas), 'Omega_dyn at natural sqrt(2/3) M': Od(m_nat),
                       'm_q/H0 for Omega_dyn = Omega_L (formal)': 1.5 * Om / np.sqrt(OL),
                       'stage-1 m_q/H0 for rho_m^2/(12 K m^2) = rho_L': Om / (2 * np.sqrt(OL)),
                       'lab dK/K = rho_lab/(3 K m^2) at 1 g/cm^3, m = m_min': (1.0 / 2.3201e17) / (3 * M**2 * (mmin * 1e-9)**2)}

# ---------------------------------------------------------------- A.6 regime (ii)
zeq = 3400.0
Heq = np.sqrt(Om * (1 + zeq)**3 * 2)
OUT['A6_regime_ii'] = {'m_min/H0': mmin / H0, 'm_Cassini/H0': m_cas / H0, 'H(z_eq)/H0 (z_eq = 3400, rho_r = rho_m)': Heq,
                       'H(z_eq) [eV]': Heq * H0,
                       'G_cos/G_lab = 3/4 vs 1.011 +- 0.018 [sigma]': (1.011 - 0.75) / 0.018,
                       'same with ledger unrounded 1.01123 +- 0.01830': (1.01123136 - 0.75) / 0.01830192,
                       'dK/K to carry Omega_L: 2 sqrt(OL) H0/m': '2 sqrt(Omega_L) (H0/m_q) [V1h]'}

# ---------------------------------------------------------------- A.7 regime (iii) estimate (corrected growth factor)
mGeV = mmin * 1e-9
T_osc = (90.0 / (np.pi**2 * gstar))**0.25 * np.sqrt(mGeV * M)
alpha_s = 0.085
kappa = (49.0 / 6.0) * alpha_s**2 / (np.pi**2 * gstar / 30.0)          # V1h: (eps - 3p) = (49/6) alpha_s^2 T^4
T0GeV = T0 * kB_GeV
# relic (a^-3) over radiation (g* T^4) from T_osc to today, entropy conservation g_s a^3 T^3 = const:
#   [rho_relic/rho_rad]_0 / [..]_osc = (T_osc/T0) * (g*(osc)/g_s(osc)) * (g_s0/g*0)
growth_correct = (T_osc / T0GeV) * (gstar / gs_osc) * (gs0 / gstar0)
growth_bound = (T_osc / T0GeV) * (3.91 / gstar)**(1.0 / 3.0)            # what bound_numbers.py uses
cases = {}
for label, L2 in [('ln^2 = 1', 1.0), ('drive from T = 1e15 GeV', (2 * np.log(1e15 / T_osc))**2)]:
    f_on = 0.1 * kappa**2 * L2
    cases[label] = {'relic fraction at onset': f_on, 'Omega_relic (correct growth)': f_on * growth_correct * Om_r,
                    'Omega_relic (bound growth)': f_on * growth_bound * Om_r,
                    'overproduction vs Omega_c = 0.264 (correct)': f_on * growth_correct * Om_r / 0.2642}
tau1 = (192 * np.pi * M**2 * hbar_eVs * 1e-9)**(1.0 / 3.0)
OUT['A7_relic_estimate'] = {'T_osc [GeV]': T_osc, 'kappa_QCD': kappa, 'growth factor, correct': growth_correct,
                            'growth factor, bound code ((3.91/g*)^(1/3) T_osc/T0)': growth_bound,
                            'ratio correct/bound': growth_correct / growth_bound, 'cases': cases,
                            'm [GeV] with tau = 1 s for Gamma = m^3/(192 pi M^2)': tau1}

# ---------------------------------------------------------------- B.5 threshold
D_EW = 125.25**2 * 246.22**2 / 8.0         # m_H^2 v^2/8 (equations track 7) [R inputs]
K = M**2
R_dS = 4 * rhoL / K
OUT['B5_threshold'] = {'Delta_EW [GeV^4]': D_EW, 'dR = 2 Delta/K [GeV^2]': 2 * D_EW / K, 'R_dS [GeV^2]': R_dS,
                       'dR/R_dS': 2 * D_EW / K / R_dS, '= Delta_EW/(2 rho_L)': D_EW / (2 * rhoL)}

# ---------------------------------------------------------------- B.3 Appleby-Battye G_Mink/G = 1 + e^{2b}
R0_pl = 12 * (H0 / (M * 1e9))**2
b = 1.0
for _ in range(200):
    b = 0.5 * np.log(8 * b / (2 * R0_pl))     # their condition 8 b e^{-2b} << 2 R0 at equality (Planck units)
x = 1e-3 / hbarc_eVm                       # 1 mm in eV^-1
R0_eV2 = 12 * H0**2
b2 = 1.0
for _ in range(200):
    b2 = 0.5 * np.log(4 * b2 / (R0_eV2 * x**2))   # F''(R0) ~ (4b/R0) e^{-2b} = (1 mm)^2
# data-driven: |f_R0| = 1 - F'(R_today), F' = [1 + tanh(aR - b)]/2, a R_dS = 2b (a = 2b/R_dS for b >> 1),
# R_today/R_dS = (Om + 4 OL)/(4 OL)  =>  alpha_today = b (1 + Om/(2 OL)),  |f_R0| = 1/(1 + e^{2 alpha})
k_ = 1 + Om / (2 * OL)
dd = {}
for name, fb in [('|f_R0| <= 1.6e-5 (Cataneo+ 2015)', 1.62e-5), ('|f_R0| <= 2.6e-5 (DESI DR1 FS)', 2.57e-5),
                 ('|f_R0| <= 1.4e-8 (Desmond-Ferreira)', 1.4e-8)]:
    two_alpha = np.log(1 / fb - 1)
    dd[name] = {'b': two_alpha / (2 * k_), 'G_Mink/G = 1 + e^{2b} (>=)': 1 + np.exp(two_alpha / k_)}
OUT['B3_AB_GMink'] = {'AB as written (Planck units, F\'\'(R0) << 1): b, G_Mink/G': [b, np.exp(2 * b)],
                      'AB F\'\'(R0) << F\'(R0) (1 mm)^2: b, G_Mink/G': [b2, np.exp(2 * b2)],
                      'data-driven (|f_R0| bounds mapped to AB, [E])': dd,
                      'note': 'G_Mink/G ~ |f_R0|^(-1/(1 + Om/(2 OL))) = |f_R0|^-0.813 for the AB form'}

with open(os.path.join(here, '..', 'results', 'v3_numbers.json'), 'w') as fh:
    json.dump(OUT, fh, indent=1, default=float)
for k, v in OUT.items():
    print(k, json.dumps(v, indent=1, default=float))
