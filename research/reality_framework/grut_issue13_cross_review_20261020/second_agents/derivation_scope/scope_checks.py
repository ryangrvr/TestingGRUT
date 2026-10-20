#!/usr/bin/env python3
"""(A) scope checks: nonminimal xi running and its effect on the averaged floor,
D-sector signs, omega_b <-> eta convention, empty-beam (clumped hidden matter) distance,
and gravitational-binding fraction needed to erase the stress-case mismatch."""
import json, math

gX, lp, lX = 0.1, 0.01, 0.001
c_, d, lD, alpha = gX**2, 0.01**2, 0.001, 0.0
mD0sq_GeV2 = (gX*1000.0)**2
v = 1000.0
mX = gX*v
mstar = v*(lp*lX)**0.25
loop = 16*math.pi**2
out = {}

# one-loop beta for xi_aa at xi=0:  16 pi^2 beta_{xi_ab} = lambda_abcd (xi_cd - delta_cd/6)
S = {'phi': 6*lp + 2*c_ + 2*c_ + 2*d,
     'X(each real comp.)': 6*lX + 2*lX + 2*c_ + 2*alpha,
     'D': 6*lD + 2*d + 2*alpha + 2*alpha}
beta = {k: -s/6/loop for k, s in S.items()}
out['beta_xi_at_mu0 (per e-fold)'] = beta
xiX_max = abs(beta['X(each real comp.)'])*abs(math.log(mstar/mX))
out['|xi_X| running mu0 -> m* (ln=%.3f)' % math.log(mstar/mX)] = xiX_max
xi_generous = 1e-3
# FRW-averaged improvement energy for L = -xi R |X|^2 (per real comp. -(xi/2) R q^2):
# rho_xi = 2 xi [3H^2 psi + 3H dpsi/dt - lap(psi)/a^2], psi=|X|^2. Average kills lap; dilution dpsi/dt=-3H psi
# => <rho_xi> = -12 xi H^2 <psi>.  With n_Q = 2 omega psi (rotating), ratio = 6 xi H^2/(omega m*).
H_GeV = 113.25/ (3.0856775814913673e19) * 6.582119569e-25       # (km/s/Mpc = 1/(3.0857e19 s))
out['H_floor_in_GeV'] = H_GeV
for label, omega in (('Q-ball omega=.8 mX', 0.8*mX), ('Q-matter omega=m*', mstar)):
    out[f'averaged xi-improvement / (m* n), xi=1e-3, {label}'] = 6*xi_generous*H_GeV**2/(omega*mstar)
MP = 2.435e18
psi_core = (math.sqrt(lp/lX)*v**2)/2     # |X|^2 = R*^2/2 in Q-matter interior
out['Jordan-frame 2 xi |X|^2/MP^2 inside Q-matter core, xi=1e-3'] = 2*xi_generous*psi_core/MP**2
out['pointwise wall estimate xi*mX^2/(m* omega) with xi=|beta|*1 (not averaged)'] = abs(beta['X(each real comp.)'])*mX**2/(mstar*0.8*mX)

# D sector sign conditions at benchmark (GRUT_CANDIDATE_1.md:211; PARAMETER_PROVENANCE.md:14)
out['D_sector'] = {'mD0^2_GeV2': mD0sq_GeV2, 'hD^2': d, 'lambda_D': lD,
                   'D_terms_nonnegative': mD0sq_GeV2 >= 0 and lD >= 0}
# generated alpha portal (quantum REPORT table): |alpha| <= 1.76e-8 for mu/mu0 in [1/4,4]
a_gen = 1.75753e-8
out['generated alpha: D-quadratic stays positive for R/v below'] = math.sqrt(0.01/a_gen)
out['generated alpha^2 vs lambda_X*lambda_D (quartic boundedness)'] = [a_gen**2, lX*lD]

# omega_b <-> eta convention: eta10 = 1e10 * omega_b * (rho_c/h^2)/(mb * n_gamma)
zeta3 = 1.2020569031595942
kB, hbar, cc = 1.380649e-23, 1.054571817e-34, 299792458.0
T0 = 2.72548
n_gamma_cm3 = 2*zeta3/math.pi**2*(kB*T0/(hbar*cc))**3/1e6
rhoc_h2 = 1.0536723728399925e-05
k_conv = 1e10*rhoc_h2/n_gamma_cm3          # eta10/omega_b = k_conv / mb[GeV]
out['n_gamma_cm3'] = n_gamma_cm3
out['eta10/omega_b for mb=.938'] = k_conv/0.938
out['mb implied by eta10=273.9 omega_b'] = k_conv/273.9
out['mb implied by eta10=273.279 omega_b'] = k_conv/273.279
out['fractional n_B difference vs mb=.938 (273.9, 273.279)'] = [0.938/(k_conv/273.9)-1, 0.938/(k_conv/273.279)-1]

# empty-beam (all hidden matter clumped; Dyer-Roeder alpha=0) distance with the same H floor
cc_km = 299792.458
Hf = 113.25181751823033
for name, z, Dlow in (('NGC 5765b stress', 0.033604405226469014, 96.9),
                      ('NGC 5765b primary', 0.031021799463629485, 96.9)):
    filled = 2*cc_km/Hf/(1+z)*(1-1/math.sqrt(1+z))
    empty = 2*cc_km/(5*Hf)*(1-(1+z)**-2.5)      # c * int dz/((1+z)^2 H), H = Hf (1+z)^1.5
    out[f'empty_vs_filled_beam {name}'] = {'filled_Mpc': filled, 'empty_Mpc': empty,
                                           'relative_increase': empty/filled-1,
                                           'D_low/empty': Dlow/empty}

# binding/collapse: fraction of the hidden floor energy that must be absent from the
# z~0 Friedmann source to erase the stress mismatch, at fixed charge
wb, mb = 0.02125, 0.94738
wh_floor = wb*mstar/mb
Hneed = Hf/1.154215714789063
wh_need = (Hneed/100)**2 - wb
out['hidden energy fraction that must be removed (stress)'] = 1 - wh_need/wh_floor
print(json.dumps(out, indent=1))
