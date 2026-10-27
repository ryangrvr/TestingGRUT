"""Symbolic checks of the [D] statements in NOTES.md (sympy).
Usage: python3 -B derivations.py ../results/derivations.json
"""
import sys
import sympy as sp
from common import dump


def main(out):
    res = {}
    # 1. Salecker-Wigner optimum for a free-particle clock: s(T)^2 = s0^2 + (T/(2 M s0))^2 minimised over s0
    s0, T, M, v, w = sp.symbols('s0 T M v omega', positive=True)
    sT2 = s0**2 + (T / (2 * M * s0))**2
    sol = sp.solve(sp.diff(sT2, s0), s0)
    s0opt = [x for x in sol if x.is_positive][0]
    res['SW_s0opt_sq'] = str(sp.simplify(s0opt**2))
    res['SW_dT2_min'] = str(sp.simplify(sT2.subs(s0, s0opt) / v**2))
    # 2. GPP optimal clock: sigma(T) = t_P (t_P/T)^{1/3}, exponent int_0^T sigma w^2 dT'
    tP, Tp = sp.symbols('t_P Tp', positive=True)
    D = sp.integrate(tP * (tP / Tp)**sp.Rational(1, 3) * w**2, (Tp, 0, T))
    res['GPP_exponent'] = str(sp.simplify(D))
    # 3. harmonic orbit psi = A e^{iwt} + B e^{-iwt}: Q = w(A^2 - B^2), nbar = w B^2, eccentricity (A+B)/(A-B)
    A, B, e = sp.symbols('A B e', positive=True)
    nbarQ = B**2 / (A**2 - B**2)
    Bsol = sp.solve(sp.Eq((A + B) / (A - B), e), B)[0]
    res['harmonic_nbar_over_Q_of_ecc'] = str(sp.simplify(nbarQ.subs(B, Bsol)))
    # 4. quartic circular orbit at fixed Q: Omega^2 = lam r^2, Q = r^2 Omega, E = Omega^2 r^2/2 + lam r^4/4
    lam, r, Q = sp.symbols('lambda r Q', positive=True)
    Om = sp.sqrt(lam) * r
    rQ = sp.solve(sp.Eq(r**2 * Om, Q), r)[0]
    Ecirc = sp.simplify((Om**2 * r**2 / 2 + lam * r**4 / 4).subs(r, rQ))
    res['quartic_Emin_fixed_Q'] = str(Ecirc)
    # minimum over all orbits at fixed Q: E >= Q^2/(2 r^2) + lam r^4/4 minimised over r
    Eeff = Q**2 / (2 * r**2) + lam * r**4 / 4
    rmin = [x for x in sp.solve(sp.diff(Eeff, r), r) if x.is_positive][0]
    res['quartic_Emin_effective_potential'] = str(sp.simplify(Eeff.subs(r, rmin)))
    # 5. Ricci scalar in FLRW with a scalar inflaton (reduced Planck units): R = 6(2H^2 + Hdot), H^2 = rho/3, Hdot = -K
    K, V = sp.symbols('K V', positive=True)
    R = 6 * (2 * (K + V) / 3 - K)
    res['R_inflaton'] = str(sp.expand(R))
    res['xiR_over_H2_matter_avg'] = str(sp.Rational(-5, 6) * (4 * V - 2 * K).subs({K: sp.Symbol('rho') / 2, V: sp.Symbol('rho') / 2}) / (sp.Symbol('rho') / 3))
    # 6. CSL point-mass limit of the k-integral: int d^3k e^{-k^2 rC^2}(1 - cos k.d) = pi^{3/2}/rC^3 (1 - e^{-d^2/4rC^2})
    k, rC, d = sp.symbols('k r_C d', positive=True)
    I = sp.integrate(4 * sp.pi * k**2 * sp.exp(-k**2 * rC**2) * (1 - sp.sin(k * d) / (k * d)), (k, 0, sp.oo))
    res['CSL_point_integral_ratio'] = str(sp.simplify(I / (sp.pi**sp.Rational(3, 2) / rC**3)))
    # 7. DP: homogeneous sphere, int int mu mu /|x-y| = (6/5) m^2/R, cross term m^2/d for d >= 2R
    m, Rr, dd, Gs = sp.symbols('m R d G', positive=True)
    EG = Gs / 2 * (2 * sp.Rational(6, 5) * m**2 / Rr - 2 * m**2 / dd)
    res['DP_EG_sphere'] = str(sp.simplify(EG))
    # 8. torque = phase gradient of energy in the toy: T(delta) = -12 g 2 sqrt2 sin(delta) with g = 0.04
    dl = sp.symbols('delta', real=True)
    E = sp.Symbol('E0') + 2 * sp.Rational(4, 100) * 2 * sp.sqrt(2) * sp.cos(dl)
    res['toy_torque_6dEddelta_at_pi2'] = str(sp.nsimplify(6 * sp.diff(E, dl).subs(dl, sp.pi / 2)))
    res['toy_torque_numeric'] = float(6 * sp.diff(E, dl).subs(dl, sp.pi / 2))
    # 9. frustrated triangle: stationary canted point satisfies a sin A = b sin B = c sin(A+B) (law of sines)
    a, b, c, Aa, Bb = sp.symbols('a b c A B')
    Vt = -a * sp.cos(Aa) - b * sp.cos(Bb) + c * sp.cos(Aa + Bb)
    res['frustrated_gradient'] = [str(sp.diff(Vt, Aa)), str(sp.diff(Vt, Bb))]
    for kk, vv in res.items():
        print(kk, ':', vv)
    dump(res, out)


if __name__ == '__main__':
    main(sys.argv[1])
