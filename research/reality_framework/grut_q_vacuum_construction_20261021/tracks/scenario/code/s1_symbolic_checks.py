"""Symbolic checks of the closed FRW system used by the scenario track (sympy).

1. Constraint propagation (FRW content of the Bianchi identity): with
       C = 3 H^2 K + 3 H Kdot - rho_m - (eps - mu q),
       Hdot = R/6 - 2 H^2,  R = 2 (eps'(q) - mu + Gamma_q qdot)/K'(q)        [q-equation, with optional dissipation]
       Kddot = [S + 4 rho_V - K R]/3 - 3 H Kdot                               [trace of (g)]
       rho_m-dot = -3 H (rho_m + p_m) + Gamma_q qdot^2                         [matter, with the decay source]
   one has dC/dt = -4 H C identically, for generic eps(q), K(q) (K'' kept) and any Gamma_q(q).
   Without the matching source in the matter equation (Q = 0 but Gamma_q != 0) the identity fails: the source is
   fixed by the Bianchi identity, nabla^m T_mn = (mu(x) - mu) d_n q.
2. Reduction to the Y-variables of qv1_core.py with eps quartic about q0 and K linear:
       F(Y) = Y + 1.5 b3 Y^2 + (b3/2 + 2 b4) Y^3 + b4 Y^4,  R(Y), rho_V(Y) and the c-offset terms.
3. KV 0806.2805 shape: eps = (1/2chi)(-q^2/q0^2 + q^4/(3 q0^4)), K = q/3 gives b3 = 4/3, b4 = 1/3, beta = 1,
   m^2 = 2/(3 chi) (with K0 = q0/3 = 1).
4. Adiabatic reduction: substituting Y = S/(3 m^2) + O(m^-4) into the constraint gives
       3 H^2 = rho_m + [S^2/12 - H^2 S - H Sdot]/m^2 + O(m^-4)           (K0 = 1)
Run: python3 -I s1_symbolic_checks.py [outdir]
"""
import json
import os
import sys

import sympy as sp

out = {}
t = sp.symbols('t', real=True)
mu, Gq0 = sp.symbols('mu Gamma_q', real=True)
q = sp.Function('q')(t)
H = sp.Function('H')(t)
rho = sp.Function('rho')(t)
p = sp.Function('p')(t)
eps = sp.Function('eps')
K = sp.Function('K')
qd = sp.diff(q, t)
Kq, K1, K2 = K(q), sp.diff(K(q), q), sp.diff(K(q), q, 2)
e1 = sp.diff(eps(q), q)
rhoV = eps(q) - mu * q


def check(with_source, Gamma_q):
    R = 2 * (e1 - mu + Gamma_q * qd) / K1
    Hd = R / 6 - 2 * H ** 2
    S = rho - 3 * p
    # Kddot = K1 qdd + K2 qd^2 from the trace equation
    Kdd = (S + 4 * rhoV - Kq * R) / 3 - 3 * H * K1 * qd
    qdd = (Kdd - K2 * qd ** 2) / K1
    rhod = -3 * H * (rho + p) + (Gamma_q * qd ** 2 if with_source else 0)
    C = 3 * H ** 2 * Kq + 3 * H * K1 * qd - rho - rhoV
    dC = sp.diff(C, t)
    dC = dC.subs({sp.Derivative(H, t): Hd, sp.Derivative(rho, t): rhod})
    dC = dC.subs(sp.Derivative(q, (t, 2)), qdd)
    return sp.simplify(sp.expand(dC + 4 * H * C))


out['constraint_propagation_Gamma0'] = str(check(True, 0))
out['constraint_propagation_with_dissipation_and_source'] = str(check(True, Gq0))
out['constraint_propagation_with_dissipation_no_source'] = str(check(False, Gq0))

# 2. Y-variables
Y, V, m2, b3, b4, beta, c, G = sp.symbols('Y V m2 b3 b4 beta c Gamma', real=True)
K0, Kp, q0, e2 = sp.symbols('K0 Kp q0 e2', positive=True)
e3, e4 = sp.symbols('e3 e4', real=True)
y = K0 * Y / Kp                          # q - q0
mu0 = sp.Symbol('mu0')
dmu = -c / q0                            # c = -(mu - mu0) q0
epsinf_minus_mu0q = e2 * y ** 2 / 2 + e3 * y ** 3 / 6 + e4 * y ** 4 / 24   # eps_inf(q) - eps_inf(q0) - mu0 (q-q0), with eps_inf(q0) = mu0 q0
rhoV_Y = epsinf_minus_mu0q - dmu * (q0 + y)
epsp_minus_mu = e2 * y + e3 * y ** 2 / 2 + e4 * y ** 3 / 6 - dmu
Gq = sp.Rational(3, 2) * G * K0 * Kp ** 2 / (K0 * (1 + Y)) ** 2
qdot = K0 * V / Kp
R_Y = 2 * (epsp_minus_mu + Gq * qdot) / Kp
subs_shape = {e2: sp.Rational(3, 2) * m2 * Kp ** 2 / K0,                      # m^2 = (2/3) K0 e2/K'^2
              e3: 3 * b3 * e2 * Kp / K0, e4: 12 * b4 * e2 * Kp ** 2 / K0 ** 2,    # b3 = e3 K0/(3 e2 K'), b4 = e4 K0^2/(12 e2 K'^2)
              q0: K0 / (beta * Kp)}                                             # beta = K0/(K' q0)
rhoV_target = sp.Rational(3, 4) * m2 * K0 * Y ** 2 * (1 + b3 * Y + b4 * Y ** 2) + c * (1 + beta * Y)
R_target = 3 * m2 * Y * (1 + sp.Rational(3, 2) * b3 * Y + 2 * b4 * Y ** 2) + 2 * beta * c / K0 + 3 * G * V / (1 + Y) ** 2
out['rhoV_Y_matches'] = bool(sp.simplify((rhoV_Y - rhoV_target).subs(subs_shape).subs(subs_shape)) == 0)
out['R_Y_matches'] = bool(sp.simplify((R_Y - R_target).subs(subs_shape).subs(subs_shape)) == 0)
# trace equation in Y: K0 Vdot = [S + 4 rhoV - K0 (1+Y) R]/3 - 3 H K0 V
Ssym, Hs = sp.symbols('S H', real=True)
Vdot = ((Ssym + 4 * rhoV_target - K0 * (1 + Y) * R_target) / 3 - 3 * Hs * K0 * V) / K0
Fy = Y + sp.Rational(3, 2) * b3 * Y ** 2 + (b3 / 2 + 2 * b4) * Y ** 3 + b4 * Y ** 4
Vdot_target = (Ssym / 3 + (c / 3) * ((4 - 2 * beta) + 2 * beta * Y)) / K0 - m2 * Fy - 3 * Hs * V - G * V / (1 + Y)
out['Vdot_Y_matches'] = bool(sp.simplify(sp.expand(Vdot - Vdot_target)) == 0)

# 3. KV shape
chi, qq = sp.symbols('chi qq', positive=True)
q0v = 3                                  # K0 = q0/3 = 1
epsKV = (1 / (2 * chi)) * (-qq ** 2 / q0v ** 2 + qq ** 4 / (3 * q0v ** 4))
d = [sp.diff(epsKV, qq, k).subs(qq, q0v) for k in range(5)]
mu0KV = d[1]
Kp_KV = sp.Rational(1, 3)
m2KV = sp.simplify(sp.Rational(2, 3) * 1 * d[2] / Kp_KV ** 2)
out['KV'] = {'rhoV(q0) with mu0': str(sp.simplify(d[0] - mu0KV * q0v)),
             'm2': str(m2KV),
             'b3': str(sp.simplify(d[3] * 1 / (3 * d[2] * Kp_KV))),
             'b4': str(sp.simplify(d[4] * 1 / (12 * d[2] * Kp_KV ** 2))),
             'beta': str(sp.simplify(1 / (Kp_KV * q0v)))}

# 4. adiabatic reduction of the constraint
Sf = sp.Function('S')(t)
Hf = sp.Function('H')(t)
mm2 = sp.Symbol('m2', positive=True)
Y1 = Sf / (3 * mm2)
Cred = 3 * Hf ** 2 * (1 + Y1) + 3 * Hf * sp.diff(Y1, t) - sp.Rational(3, 4) * mm2 * Y1 ** 2
# 3H^2 = rho_m - [3 H^2 Y1 + 3 H Y1dot - rho_V(Y1)]  => q-correction:
corr = sp.expand(-(3 * Hf ** 2 * Y1 + 3 * Hf * sp.diff(Y1, t) - sp.Rational(3, 4) * mm2 * Y1 ** 2))
target = (Sf ** 2 / 12 - Hf ** 2 * Sf - Hf * sp.diff(Sf, t)) / mm2
out['adiabatic_correction_matches'] = bool(sp.simplify(corr - target) == 0)
out['adiabatic_correction'] = str(target)

for k_, v_ in out.items():
    print(k_, ':', v_)
assert out['constraint_propagation_Gamma0'] == '0'
assert out['constraint_propagation_with_dissipation_and_source'] == '0'
assert out['constraint_propagation_with_dissipation_no_source'] != '0'
assert out['rhoV_Y_matches'] and out['R_Y_matches'] and out['Vdot_Y_matches'] and out['adiabatic_correction_matches']
outdir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
with open(os.path.join(outdir, 'symbolic_checks.json'), 'w') as fh:
    json.dump(out, fh, indent=1)
