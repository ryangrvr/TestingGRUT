#!/usr/bin/env python3
"""(B) Vacuum-energy matching and the critical negative Lambda for the maser envelope.

Inputs copied from 2ee7f82 COMPARISON.json / PREDICTION.json (conservative H_floor,
cosmological_z_max corners, lower-distance corners). Reimplements the envelope with
an added constant negative vacuum term: H^2(z) >= Hf^2 (1+z)^3 - L, L = |Lambda| c^2/3 >= 0.
"""
import json, math
from scipy.integrate import quad
from scipy.optimize import brentq

c = 299792.458
Hf = 113.25181751823033           # conservative H_floor (PREDICTION.json)
rows = {  # name: (D_low, z_high_primary, z_high_stress)
    'UGC 3789': (39.5, 0.013614626766245808, 0.016153629314733964),
    'NGC 5765b': (96.9, 0.031021799463629485, 0.033604405226469014),
}
rho_c_over_h2_GeV_cm3 = 1.0536723728399925e-05   # EXCLUSION_AUDIT.json units
hbarc_GeV_cm = 1.973269804e-14
rho_c_over_h2_GeV4 = rho_c_over_h2_GeV_cm3*hbarc_GeV_cm**3

def DA(z, L):
    I = quad(lambda x: 1/math.sqrt(Hf**2*(1+x)**3 - L), 0, z, epsabs=0, epsrel=1e-13)[0]
    return c/(1+z)*I

def envelope(z, H):
    return 2*c/H/(1+z)*(1-1/math.sqrt(1+z))

out = {'rho_c_over_h2_GeV4': rho_c_over_h2_GeV4}
for name, (Dlow, zp, zs) in rows.items():
    for label, z in (('primary', zp), ('stress', zs)):
        assert abs(DA(z, 0) - envelope(z, Hf)) / envelope(z, Hf) < 1e-11
        if DA(z, 0) >= Dlow:
            out[f'{name} {label}'] = 'already compatible at Lambda=0'
            continue
        L = brentq(lambda L: DA(z, L) - Dlow, 0, Hf**2*0.999999, xtol=1e-12)
        H0 = math.sqrt(Hf**2 - L)
        OmL = -L/H0**2
        Om = Hf**2/H0**2
        omegaL = -L/1e4
        rho = omegaL*rho_c_over_h2_GeV4
        out[f'{name} {label}'] = {
            'critical_omega_Lambda=Omega_L*h^2': omegaL,
            'critical_Omega_Lambda (relative to model H0)': OmL,
            'model_H0_at_critical_km_s_Mpc': H0,
            'Omega_m_at_critical': Om,
            'q0_at_critical (flat, matter+Lambda)': Om/2 - OmL,
            'critical_rho_vac_GeV4': rho,
            'critical_rho_vac_as_-(x meV)^4, x=': (abs(rho))**0.25*1e12,
        }

# vacuum budget at v = 1 TeV, MSbar mu = 100 GeV (fixed-coefficient benchmark matching rho_Lambda(mu0)=0)
v = 1000.0
MP = 2.435e18   # reduced Planck mass, GeV
rho_new = -1.22652632365e-6*v**4
def cw(n, m, cst, mu=100.0):
    return n*m**4/(64*math.pi**2)*(math.log(m**2/mu**2) - cst)
sm = {}
for mt in (163.0, 172.5):
    for mu in (56.234, 100.0, 172.5):
        terms = {'top': cw(-12, mt, 1.5, mu), 'W': cw(6, 80.37, 5/6, mu), 'Z': cw(3, 91.19, 5/6, mu),
                 'h': cw(1, 125.1, 1.5, mu), 'b': cw(-12, 4.18, 1.5, mu)}
        sm[f'mt={mt},mu={mu}'] = {'total_GeV4': sum(terms.values()), **{k: round(t, 1) for k, t in terms.items()}}
L_AdS_GeVinv = math.sqrt(3*MP**2/abs(rho_new))
gstar = 106.75
T_turn = (30*abs(rho_new)/(math.pi**2*gstar))**0.25
rho_obs = 0.69*rho_c_over_h2_GeV4*0.70**2
out['vacuum_budget'] = {
    'rho_vac_new_scalar_one_loop_GeV4 (v=1TeV, rho_Lambda(mu0)=0, fixed coefficients)': rho_new,
    'L_AdS_m': L_AdS_GeVinv*hbarc_GeV_cm/100,
    'L_AdS_over_c_seconds': L_AdS_GeVinv*6.582119569e-25,
    'radiation_temperature_where_rho_rad=|rho_new|_GeV (turnaround scale)': T_turn,
    'SM_one_loop_vacuum_MSbar_Landau_GeV4 (tree V_SM(vac)=0)': sm,
    'observed_dark_energy_GeV4 (Omega_L=.69,h=.70)': rho_obs,
    'ratio |rho_new|/rho_obs': abs(rho_new)/rho_obs,
    'L_AdS_v at v/MP=1e-5 (quantum REPORT check)': math.sqrt(3/1.22652632365e-6)*1e5,
}
print(json.dumps(out, indent=1))
