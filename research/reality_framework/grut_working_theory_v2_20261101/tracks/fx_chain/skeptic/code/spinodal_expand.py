"""fx_chain_skeptic (Claude, own code): analytic small-k expansion of the M1 cold-fluid dispersion near the spinodal,
as an independent check of fx_chain frag_m1.py (which finds the roots of the full cubic numerically).
Linearised X (amplitude f, phase) + phi about the rotating homogeneous state gives, to O(k^4, z k^2) with z = nu^2,
   z (c + 4 Om^2 D) = -c k^2 - (D + 2 a F^2) k^4,   c = 2 F^2 (a D - 2 S^2),  D = U_phi''(S) + F^2,
so the spinodal is c = 0 (a D = 2 S^2), Gamma_peak = |c| / (4 Om sqrt(D (D + 2 a F^2))), k_peak^2 = |c| / (2 (D + 2aF^2)).
Also the a_ref-start entropy excess (network state functions) for the R = 4.97 heat.  [D] + [C]"""
import math, json, sys
from scipy.optimize import brentq, minimize_scalar
a, b, c3 = 0.1, 3.0, 0.2


def S_of_F(F):
    return brentq(lambda s: b * (s * s - 1) * s - 3 * c3 * (s - 1) ** 2 + F * F * s, 1e-14, 1.0, xtol=1e-16)


def st(F):
    S = S_of_F(F); Om = math.sqrt(S * S + a * F * F); D = b * (3 * S * S - 1) - 6 * c3 * (S - 1) + F * F
    return dict(F=F, S=S, Om=Om, D=D, n=100 * Om * F * F, c=2 * F * F * (a * D - 2 * S * S))


Fsp = brentq(lambda F: st(F)['c'], 1.7, 2.0, xtol=1e-15)
s0 = st(Fsp)
h = 1e-6
dc = (st(Fsp + h)['c'] - st(Fsp - h)['c']) / (2 * h); dn = (st(Fsp + h)['n'] - st(Fsp - h)['n']) / (2 * h)
c_over_delta = -dc * s0['n'] / dn          # c = (dc/dF) dF, delta = -dn/n_sp
G1 = abs(c_over_delta) / (4 * s0['Om'] * math.sqrt(s0['D'] * (s0['D'] + 2 * a * Fsp ** 2)))
K1 = math.sqrt(abs(c_over_delta) / (2 * (s0['D'] + 2 * a * Fsp ** 2)))
# minimum of Omega along the family (thermodynamic spinodal) must coincide with c = 0
Fmin = minimize_scalar(lambda F: st(F)['Om'], bracket=(1.7, 1.85, 2.0), tol=1e-13).x
out = dict(F_spinodal_c0=Fsp, n_spinodal=s0['n'], S=s0['S'], Omega=s0['Om'], D=s0['D'], F_min_Omega=Fmin,
           c_over_delta=c_over_delta, Gamma_over_delta=G1, k_over_sqrt_delta=K1,
           fx_chain=dict(Gamma_over_delta=0.774405344980941, k_over_sqrt_delta=0.9577685205840226, n_sp=233.22058697882017))
print(json.dumps(out, indent=1))
json.dump(out, open(sys.argv[1], 'w'), indent=1)
