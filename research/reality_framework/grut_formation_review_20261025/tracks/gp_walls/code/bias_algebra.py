#!/usr/bin/env python3
"""Exact algebra and small-number estimates for the proposed wall tilt (Claude review, track gp_walls).

V0 = lam/4 (phi^2-v^2)^2,  V_bias = Delta/4 (2 - 3x + x^3),  x = phi/v.
Checks: stationarity and energies at x=+/-1, factorisation, positivity (also with the Higgs portal),
curvatures, identity with Saikawa 1703.02576 eq. (2.23) plus a constant, first-order wall-tension
shift, false-vacuum decay by bubbles, and the size of the Z2-odd couplings it adds about +v.
"""
import json
import sys
from pathlib import Path
import numpy as np
import sympy as sp

x, lam, v, D, eps, phi, kap, lamH, b = sp.symbols('x lam v Delta epsilon phi kappa lambda_H b', real=True)
V0 = lam * v**4 / 4 * (x**2 - 1)**2
Vb = D / 4 * (2 - 3 * x + x**3)
Vt = V0 + Vb

out = {}
out['dV_at_plus1'] = str(sp.simplify(sp.diff(Vt, x).subs(x, 1)))
out['dV_at_minus1'] = str(sp.simplify(sp.diff(Vt, x).subs(x, -1)))
out['V_at_plus1'] = str(sp.simplify(Vt.subs(x, 1)))
out['V_at_minus1'] = str(sp.simplify(Vt.subs(x, -1)))
fact = (x - 1)**2 / 4 * (lam * v**4 * (x + 1)**2 + D * (x + 2))
out['factorisation_residual'] = str(sp.expand(Vt - fact))
br = lam * v**4 * (x + 1)**2 + D * (x + 2)
xm = sp.solve(sp.diff(br, x), x)[0]
out['bracket_min_location'] = str(sp.simplify(xm))
out['bracket_min_value'] = str(sp.simplify(br.subs(x, xm)))
# physical curvature m^2 = V''(phi) = V''(x)/v^2
out['m2_plus'] = str(sp.simplify(sp.diff(Vt, x, 2).subs(x, 1) / v**2))
out['m2_minus'] = str(sp.simplify(sp.diff(Vt, x, 2).subs(x, -1) / v**2))
# Saikawa (2.23): Delta V = eps v phi (phi^2/3 - v^2); V_bias = (4/3) eps v^4  (eq. 2.24)
saik = eps * v * (x * v) * ((x * v)**2 / 3 - v**2)
out['identity_with_Saikawa_2.23'] = str(sp.simplify(Vb - (saik.subs(eps, 3 * D / (4 * v**4)) + D / 2)))
# oddness: V_bias(x) + V_bias(-x) = Delta  => first-order tension shift vanishes for a kink centred at x=0
out['Vb(x)+Vb(-x)'] = str(sp.simplify(Vb + Vb.subs(x, -x)))

# Numbers at the stated point
LAM, VV, DD = 0.03, 0.924089637375362, 1e-22
MRED = 2.435e18
num = {lam: LAM, v: VV, D: DD}
sigma = 2 * np.sqrt(2 * LAM) * VV**3 / 3
w = np.sqrt(2 / LAM) / VV
nums = dict(sigma_GeV3=sigma, sigma_cuberoot_GeV=sigma**(1 / 3), wall_width_GeVinv=w,
            barrier_lam_v4_over_4=LAM * VV**4 / 4, Delta_over_barrier=DD / (LAM * VV**4 / 4),
            bracket_min_GeV4=float(sp.N(sp.sympify(out['bracket_min_value'], locals={'lam': lam, 'v': v, 'Delta': D}).subs(num), 30)),
            bracket_min_minus_Delta_GeV4=float(sp.N((sp.sympify(out['bracket_min_value'], locals={'lam': lam, 'v': v, 'Delta': D}) - D).subs(num), 30)),
            m2_plus_GeV2=float(sp.sympify(out['m2_plus'], locals={'lam': lam, 'v': v, 'Delta': D}).subs(num)),
            m2_minus_GeV2=float(sp.sympify(out['m2_minus'], locals={'lam': lam, 'v': v, 'Delta': D}).subs(num)),
            m2_shift_GeV2=float(3 * DD / (2 * VV**2)))
# first-order tension shift, numerically, for the unperturbed kink x = tanh(z/w)
z = np.linspace(-40 * w, 40 * w, 400001)
t = np.tanh(z / w)
f = DD / 4 * (2 - 3 * t + t**3) - DD * (z < 0)
nums['first_order_tension_shift_GeV3'] = float(np.trapezoid(f, z))
nums['first_order_tension_shift_if_centre_offset_by_w_GeV3'] = float(DD * w)    # gauge of wall position
nums['relative_second_order_tension_shift_scale'] = (DD / (LAM * VV**4))**2
# thin-wall bubble nucleation of the true vacuum inside a false domain (Coleman): B = 27 pi^2 sigma^4 / (2 Delta^3)
nums['thin_wall_bounce_exponent'] = 27 * np.pi**2 * sigma**4 / (2 * DD**3)
nums['thin_wall_critical_radius_GeVinv'] = 3 * sigma / DD
# wall acceleration from the pressure Delta (surface mass sigma)
nums['pressure_acceleration_GeV'] = DD / sigma

# Higgs portal (Codex debris packet): kappa/2 (phi^2-v^2)(|H|^2 - vH^2/2), Higgs quartic lamH (|H|^2-vH^2/2)^2
# Minimise over b = |H|^2 - vH^2/2 (unconstrained lower bound): lam_eff = lam - kappa^2/(4 lamH) multiplies the V0 term.
KAP, LAMH = 9.6e-10, 0.129
lam_eff = LAM - KAP**2 / (4 * LAMH)
nums['higgs_portal_lam_eff'] = lam_eff
nums['positivity_with_portal_requires_Delta_below_GeV4'] = 4 * lam_eff * VV**4
nums['positivity_margin'] = 4 * lam_eff * VV**4 / DD

# Z2-odd couplings about +v, s = phi - v
s = sp.symbols('s', real=True)
Vs = sp.expand(Vt.subs(x, (v + s) / v))
c3 = sp.Poly(Vs, s).coeff_monomial(s**3)
c1 = sp.Poly(Vs, s).coeff_monomial(s)
out['tadpole_coefficient_about_plus_v'] = str(sp.simplify(c1))
out['cubic_coefficient_about_plus_v'] = str(sp.simplify(c3))
c3_SSB = LAM * VV
c3_bias = DD / (4 * VV**3)
nums['cubic_s3_SSB_GeV'] = c3_SSB
nums['cubic_s3_from_tilt_GeV'] = c3_bias
nums['cubic_ratio'] = c3_bias / c3_SSB
# One-loop induced linear Higgs portal A phi |H|^2 from the tilt cubic x quadratic portal (order of magnitude):
# A ~ (3 Delta/(2 v^3)) * kappa / (16 pi^2)   [log ~ 1]
A_lin = 3 * DD / (2 * VV**3) * KAP / (16 * np.pi**2)
vH, mh = 246.22, 125.25
theta_quad = KAP * VV * vH / (mh**2)
theta_lin = A_lin * vH / mh**2
nums['induced_linear_portal_A_GeV'] = A_lin
nums['mixing_from_quadratic_portal'] = theta_quad
nums['mixing_from_induced_linear_portal'] = theta_lin
nums['mixing_ratio_lin_over_quad'] = theta_lin / theta_quad
# one-loop tadpole from the tilt cubic with an s loop ~ c3_feyn m_s^2/(32 pi^2): compare to tilt linear term 3Delta/(4v)
ms2 = 2 * LAM * VV**2
tad_loop = 6 * c3_bias * ms2 / (32 * np.pi**2)
nums['loop_tadpole_GeV3'] = tad_loop
nums['tilt_linear_term_GeV3'] = 3 * DD / (4 * VV)
nums['loop_tadpole_ratio'] = tad_loop / (3 * DD / (4 * VV))
# Planck-suppressed odd operator for comparison
nums['v5_over_Mred_GeV4'] = VV**5 / MRED
nums['coefficient_equivalent_to_Delta'] = DD / (VV**5 / MRED)
out['numbers'] = nums
txt = json.dumps(out, indent=1, default=float)
Path(sys.argv[1]).write_text(txt + '\n')
print(txt)
