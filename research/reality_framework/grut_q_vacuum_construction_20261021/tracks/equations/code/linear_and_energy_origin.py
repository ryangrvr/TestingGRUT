"""Linearisation about Minkowski, scalar-tensor map, Newtonian limit, and the energy-origin
(eps -> eps + c) analysis of q-theory.  Generic functions K(q), eps(q); Taylor data at q0:
  K0=K(q0), K1=K'(q0), K2=K''(q0), e2=eps''(q0), e3=eps'''(q0), e4=eps''''(q0).
Results -> ../results/linear_and_energy_origin.json
"""
import json, os
import sympy as sp

out = {}
q, mu, c, R, phi = sp.symbols('q mu c R phi', real=True)
Kf, ef = sp.Function('K'), sp.Function('eps')
q0, K0, K1, K2, K3, e1, e2, e3, e4 = sp.symbols('q0 K0 K1 K2 K3 e1 e2 e3 e4', real=True)
T = sp.symbols('T', real=True)

# ---------------------------------------------------------------------------
# 2. Scalar-tensor map.  U(phi) = eps(q) - mu q with phi = K(q).
#    dU/dphi = (eps' - mu)/K',  d2U/dphi2 = [eps''/K' - (eps'-mu) K''/K'^2]/K'
# ---------------------------------------------------------------------------
Kq = Kf(q); eq_ = ef(q)
U = eq_ - mu * q
Up = sp.diff(U, q) / sp.diff(Kq, q)
Upp = sp.simplify(sp.diff(Up, q) / sp.diff(Kq, q))
out['U_prime_phi'] = str(sp.simplify(Up))
out['U_doubleprime_phi'] = str(Upp)
# BD omega=0 field eqs: R = 2 U'(phi) <=> eps' - K' R/2 = mu  (identical)
out['phi_eq_equals_q_eq'] = str(sp.simplify((R - 2 * Up) * sp.diff(Kq, q) / 2 - ((sp.diff(Kq, q) * R / 2) - (sp.diff(eq_, q) - mu))))

# trace equation: 3 Box K = T - 4 (eps - mu q) + K R, with R = 2(eps'-mu)/K'
# -> Box phi = dW/dphi + T/3,  dW/dphi = (2/3)(phi U' - 2U)
Rq = 2 * (sp.diff(eq_, q) - mu) / sp.diff(Kq, q)
boxphi = (T - 4 * U + Kq * Rq) / 3
out['Box_K_rhs'] = str(sp.simplify(boxphi))
out['check_BoxK_equals_(2/3)(phiU-2U)+T/3'] = str(sp.simplify(boxphi - (sp.Rational(2, 3) * (Kq * Up - 2 * U) + T / 3)))

# ---------------------------------------------------------------------------
# 4. Linearisation around Minkowski equilibrium: eps'(q0)=mu, eps(q0)-mu q0 = 0.
#    Box dK - m^2 dK = T/3 ;  m^2 = d/dK [ (2/3)(K U' - 2U) ] at q0
# ---------------------------------------------------------------------------
W1 = sp.Rational(2, 3) * (Kq * Up - 2 * U)
m2_expr = sp.diff(W1, q) / sp.diff(Kq, q)
# substitute Taylor data at q0 with equilibrium conditions
subsd = {}
m2_at = m2_expr.subs({sp.Derivative(ef(q), (q, 2)): e2, sp.Derivative(Kf(q), (q, 2)): K2})
m2_at = m2_at.subs({sp.Derivative(ef(q), q): mu, sp.Derivative(Kf(q), q): K1}).subs({ef(q): mu * q, Kf(q): K0})
m2_at = sp.simplify(m2_at)
out['m_q_squared'] = str(m2_at)
assert sp.simplify(m2_at - sp.Rational(2, 3) * K0 * e2 / K1**2) == 0
# Einstein frame: sigma with K = M^2 exp(sqrt(2/3) sigma/M), V = M^4 U/K^2 ; V''(sigma) at min
M = sp.symbols('M', positive=True)
sig = sp.symbols('sigma', real=True)
# use U as a generic function of phi with U(phi0)=U'(phi0)=0
Uf = sp.Function('Uphi')
phi_s = M**2 * sp.exp(sp.sqrt(sp.Rational(2, 3)) * sig / M)
Vsig = M**4 * Uf(phi_s) / phi_s**2
V2 = sp.diff(Vsig, sig, 2).subs(sig, 0)
u0, u1, u2 = sp.symbols('u0 u1 u2')
V2 = V2.subs({sp.Subs(sp.Derivative(Uf(sp.Symbol('_xi_1')), (sp.Symbol('_xi_1'), 2)), sp.Symbol('_xi_1'), M**2): u2})
V2 = sp.simplify(V2.doit())
# replace derivative objects generically
V2s = V2
for d, val in [(2, u2), (1, u1)]:
    V2s = V2s.replace(lambda e: isinstance(e, sp.Subs) and e.expr.is_Derivative and len(e.expr.variables) == d, lambda e: val)
V2s = V2s.replace(lambda e: isinstance(e, sp.Function) and e.func == Uf, lambda e: u0)
V2s = sp.simplify(V2s.subs({u0: 0, u1: 0}))
out['Einstein_frame_mass_squared (phi0=M^2, U=U\'=0)'] = str(V2s)
# compare with (2/3) phi0 U''(phi0) = (2/3) M^2 u2
assert sp.simplify(V2s - sp.Rational(2, 3) * M**2 * u2) == 0

# ---------------------------------------------------------------------------
# 4b. Newtonian limit. Static weak field ds^2 = -(1+2Psi)dt^2 + (1-2Phi)dx^2, dK = K - K0.
# Linearised equations (derived in NOTES 4.3):
#   00: 2 K0 lap(Phi) = rho + lap(dK)
#   ij: K0 (delta_ij lap - d_i d_j)(Psi - Phi) = -(delta_ij lap - d_i d_j) dK  ->  Psi - Phi = -dK/K0
#   trace: (lap - m^2) dK = -rho/3
# First verify the linearised Einstein tensor components symbolically.
# ---------------------------------------------------------------------------
t, x, y, z, ep = sp.symbols('t x y z epsilon', real=True)
X = [t, x, y, z]
Ps, Ph = sp.Function('Psi')(x, y, z), sp.Function('Phi')(x, y, z)
g = sp.diag(-(1 + 2 * ep * Ps), 1 - 2 * ep * Ph, 1 - 2 * ep * Ph, 1 - 2 * ep * Ph)
gi = sp.diag(*[1 / g[i, i] for i in range(4)])
Gam = [[[sum(gi[l, s] * (sp.diff(g[s, m], X[n]) + sp.diff(g[s, n], X[m]) - sp.diff(g[m, n], X[s])) for s in range(4)) / 2
         for n in range(4)] for m in range(4)] for l in range(4)]


def Riem(r, s, m, n):
    e = sp.diff(Gam[r][n][s], X[m]) - sp.diff(Gam[r][m][s], X[n])
    e += sum(Gam[r][m][l] * Gam[l][n][s] - Gam[r][n][l] * Gam[l][m][s] for l in range(4))
    return e


def lin(e):
    return sp.simplify(sp.diff(e, ep).subs(ep, 0))


Ric = [[lin(sum(Riem(r, s, r, n) for r in range(4))) for n in range(4)] for s in range(4)]
Rs1 = sp.simplify(-Ric[0][0] + Ric[1][1] + Ric[2][2] + Ric[3][3])
lap = lambda f: sp.diff(f, x, 2) + sp.diff(f, y, 2) + sp.diff(f, z, 2)
G00 = sp.simplify(Ric[0][0] + Rs1 / 2)
G11 = sp.simplify(Ric[1][1] - Rs1 / 2)
G12 = sp.simplify(Ric[1][2])
out['lin_G00_minus_2lapPhi'] = str(sp.simplify(G00 - 2 * lap(Ph)))
out['lin_G11_minus_(lap-d1d1)(Psi-Phi)'] = str(sp.simplify(G11 - (lap(Ps - Ph) - sp.diff(Ps - Ph, x, 2))))
out['lin_G12_minus_(-d1d2)(Psi-Phi)'] = str(sp.simplify(G12 + sp.diff(Ps - Ph, x, y)))
out['lin_R'] = str(sp.simplify(Rs1 - (2 * lap(Ps) - 4 * lap(Ph))))

# radial Yukawa solution for point mass Mp: check ODEs for r>0 and 1/r normalisation
r, m, Mp, Kz = sp.symbols('r m M_p K0', positive=True)
lapr = lambda f: sp.diff(r**2 * sp.diff(f, r), r) / r**2
dK = Mp * sp.exp(-m * r) / (12 * sp.pi * r)
Gs = 1 / (8 * sp.pi * Kz)
PhiN = -Gs * Mp / r * (1 - sp.exp(-m * r) / 3)
PsiN = -Gs * Mp / r * (1 + sp.exp(-m * r) / 3)
out['yukawa_trace_eq_r>0'] = str(sp.simplify(lapr(dK) - m**2 * dK))           # = -rho/3 = 0 for r>0
out['yukawa_dK_flux_check (should be Mp/3)'] = str(sp.limit(-4 * sp.pi * r**2 * sp.diff(dK, r), r, 0))  # flux = -int lap dK = rho_tot/3 -> Mp/3
out['yukawa_00_eq_r>0'] = str(sp.simplify(2 * Kz * lapr(PhiN) - lapr(dK)))
out['yukawa_Psi_minus_Phi_plus_dK/K0'] = str(sp.simplify(PsiN - PhiN + dK / Kz))
out['yukawa_00_flux_check (minus sign convention: should be -Mp)'] = str(sp.simplify(sp.limit(-4 * sp.pi * r**2 * sp.diff(2 * Kz * PhiN - dK, r), r, 0)))
out['gamma_PPN(r)'] = '(1 - exp(-m r)/3)/(1 + exp(-m r)/3)'
out['G_eff(r)'] = 'G_* (1 + exp(-m r)/3), G_* = 1/(8 pi K(q0))'

# ---------------------------------------------------------------------------
# 5a. Equilibrium with mu = mu0[eps + c]:  eps'(q0) = mu0, eps(q0) + c - mu0 q0 = 0
# ---------------------------------------------------------------------------
cc = sp.symbols('c', real=True)
qc = sp.Function('q0')(cc)
muc = sp.Function('mu0')(cc)
eqA = sp.diff(ef(qc), qc) - muc
eqB = ef(qc) + cc - muc * qc
d1 = sp.solve([sp.diff(eqA, cc), sp.diff(eqB, cc)], [sp.diff(qc, cc), sp.diff(muc, cc)], dict=True)[0]
dq_dc = sp.simplify(d1[sp.diff(qc, cc)].subs(muc, sp.diff(ef(qc), qc)))
dmu_dc = sp.simplify(d1[sp.diff(muc, cc)].subs(muc, sp.diff(ef(qc), qc)))
out['5a_dq0_dc'] = str(dq_dc)
out['5a_dmu0_dc'] = str(dmu_dc)
# second derivatives
d2q = sp.simplify(sp.diff(dq_dc, cc).subs(sp.diff(qc, cc), dq_dc))
d2mu = sp.simplify(sp.diff(dmu_dc, cc).subs(sp.diff(qc, cc), dq_dc))
out['5a_d2q0_dc2'] = str(d2q)
out['5a_d2mu0_dc2'] = str(d2mu)
out['5a_dK(q0)_dc'] = 'K\'(q0)/(q0 eps\'\'(q0))'

# ---------------------------------------------------------------------------
# 5b. Fixed mu.  Static (de Sitter / AdS) states of flat FRW or maximally symmetric space:
#   3 H^2 K = rhoV, R = 12 H^2 = 4 rhoV/K, eps'(q) - mu = K'(q) R/2 = 2 K' rhoV/K,  rhoV = eps + c - mu q
# Expand around the Minkowski point (c = c*, q = q0) in dc = c - c*.
# ---------------------------------------------------------------------------
dc = sp.symbols('delta_c', real=True)
s = sp.symbols('s', real=True)  # q - q0
a1, a2, a3 = sp.symbols('a1 a2 a3')
epsT = mu * q0 + mu * s + e2 * s**2 / 2 + e3 * s**3 / 6 + e4 * s**4 / 24   # eps(q0+s) with eps(q0)+c*-mu q0=0 -> eps(q0)=mu q0 - c*
cstar = sp.symbols('c_star')
epsT = epsT - cstar
KT = K0 + K1 * s + K2 * s**2 / 2 + K3 * s**3 / 6
rhoV = sp.expand(epsT + cstar + dc - mu * (q0 + s))
F = sp.diff(epsT, s) - mu - 2 * sp.diff(KT, s) * rhoV / KT
s_ser = a1 * dc + a2 * dc**2 + a3 * dc**3
Fs = sp.series(F.subs(s, s_ser), dc, 0, 4).removeO()
sol = sp.solve([sp.expand(Fs).coeff(dc, k) for k in (1, 2, 3)], [a1, a2, a3], dict=True)[0]
rhoV_ser = sp.series(rhoV.subs(s, s_ser).subs(sol), dc, 0, 4).removeO()
rhoV_ser = sp.collect(sp.simplify(sp.expand(rhoV_ser)), dc)
out['5b_q_shift_series'] = {k: str(sp.simplify(v)) for k, v in zip(['a1', 'a2', 'a3'], [sol[a1], sol[a2], sol[a3]])}
out['5b_rhoV_static_series_in_dc'] = str(rhoV_ser)
c2 = sp.simplify(sp.expand(rhoV_ser).coeff(dc, 2))
c3 = sp.simplify(sp.expand(rhoV_ser).coeff(dc, 3))
out['5b_rhoV_coeff_dc^1'] = str(sp.simplify(sp.expand(rhoV_ser).coeff(dc, 1)))
out['5b_rhoV_coeff_dc^2'] = str(c2)
out['5b_rhoV_coeff_dc^2_in_terms_of_m2'] = str(sp.simplify(c2.subs(e2, sp.Rational(3, 2) * sp.Symbol('m2') * K1**2 / K0)))
out['5b_rhoV_coeff_dc^3'] = str(c3)
# K at static point and H^2
H2_ser = sp.series((rhoV / (3 * KT)).subs(s, s_ser).subs(sol), dc, 0, 4).removeO()
out['5b_H2_static_series'] = str(sp.collect(sp.simplify(sp.expand(H2_ser)), dc))

# ---------------------------------------------------------------------------
# 5c. Oscillation energy deposited by a sudden drop Delta (post-transition mu):
#   displacement s0 = a1*Delta (pre-transition dS offset), Einstein-frame energy (1/2) m^2 dsigma^2,
#   dsigma = sqrt(3/(2 K0)) K1 s0
# ---------------------------------------------------------------------------
Dl, m2s = sp.symbols('Delta m2', positive=True)
s0 = sol[a1] * Dl
Eosc = sp.Rational(1, 2) * m2s * sp.Rational(3, 2) / K0 * (K1 * s0)**2
Eosc = sp.simplify(Eosc.subs(m2s, sp.Rational(2, 3) * K0 * e2 / K1**2))
out['5c_oscillation_energy_sudden_drop'] = str(Eosc)
out['5c_oscillation_energy_in_m2'] = str(sp.simplify(Eosc.subs(e2, sp.Rational(3, 2) * sp.Symbol('m2') * K1**2 / K0)))

for k_, v_ in out.items():
    print(k_, ':', v_)
with open(os.path.join(os.path.dirname(__file__), '..', 'results', 'linear_and_energy_origin.json'), 'w') as f:
    json.dump(out, f, indent=1)
