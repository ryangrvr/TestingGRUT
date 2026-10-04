#!/usr/bin/env python3
"""CARD-01 pre-data model-only computations (NO DATA; no DESI product touched).

Card #1 v1 law (transcribed from scout-1 @ a2987fe, calc/wz_dark_energy.py, blob cd900011...):
    E_ref(z)^2 = 0.31 (1+z)^3 + 0.69,  R(z) = E_ref^2 / (1 + E_ref^2)   [H0*tau_2 = 1]
    w(z; eps)  = -1 + eps * R(z),      primary branch eps < 0, control eps > 0.

Registered CPL diagnostic (governance-frozen for Card #1 v1):
    w_CPL = w0 + wa f(z), f = z/(1+z);  J = int_0^2 [w_model - w_CPL]^2 dz  (unweighted, in z).
Because w_model + 1 = eps R and J is quadratic, (1+w0, wa) = eps * (c0, ca) with
(c0, ca) = argmin int_0^2 (R - c0 - ca f)^2 dz: the locus is EXACTLY linear and the
ratio wa/(1+w0) = ca/c0 is an exact eps-independent constant.  Computed here by
high-precision quadrature (mpmath, 30 digits) and cross-checked with scipy.
"""
import mpmath as mp
import numpy as np
from scipy.integrate import quad

mp.mp.dps = 30
OM, OL = mp.mpf("0.31"), mp.mpf("0.69")
ZMAX = mp.mpf(2)


def E2(z):
    return OM * (1 + z) ** 3 + OL


def R(z):
    e2 = E2(z)
    return e2 / (1 + e2)


def f(z):
    return z / (1 + z)


def gram():
    I = lambda g: mp.quad(g, [0, 1, ZMAX])
    A = mp.matrix([[I(lambda z: 1), I(f)], [I(f), I(lambda z: f(z) ** 2)]])
    b = mp.matrix([I(R), I(lambda z: R(z) * f(z))])
    c = mp.lu_solve(A, b)
    resid = I(lambda z: (R(z) - c[0] - c[1] * f(z)) ** 2)
    return c[0], c[1], resid, I(lambda z: (R(z) - I(R) / ZMAX) ** 2)


def scipy_check():
    Rn = lambda z: (0.31 * (1 + z) ** 3 + 0.69) / (1 + 0.31 * (1 + z) ** 3 + 0.69)
    fn = lambda z: z / (1 + z)
    A = np.array([[quad(lambda z: 1, 0, 2)[0], quad(fn, 0, 2)[0]],
                  [quad(fn, 0, 2)[0], quad(lambda z: fn(z) ** 2, 0, 2)[0]]])
    b = np.array([quad(Rn, 0, 2)[0], quad(lambda z: Rn(z) * fn(z), 0, 2)[0]])
    return np.linalg.solve(A, b)


def main():
    c0, ca, res, tot = gram()
    print("CARD-01 registered CPL diagnostic (unweighted LSQ in z on [0,2])")
    print(f"  c0 = (1+w0)/eps = {mp.nstr(c0, 15)}")
    print(f"  ca =  wa/eps    = {mp.nstr(ca, 15)}")
    print(f"  ratio wa/(1+w0) = ca/c0 = {mp.nstr(ca / c0, 15)}   (exact constant; locus exactly linear)")
    print(f"  fractional residual int(R-fit)^2 / int(R-mean)^2 = {mp.nstr(res / tot, 6)}")
    s = scipy_check()
    print(f"  scipy cross-check: c0={s[0]:.12f} ca={s[1]:.12f} ratio={s[1]/s[0]:.12f}")

    # UNREGISTERED HISTORICAL DIAGNOSTIC: value+slope at z=0 (scout-1 cpl_fit)
    R0 = R(mp.mpf(0)); dR0 = mp.diff(R, 0)
    print("\nUNREGISTERED HISTORICAL DIAGNOSTIC (local z=0 value+slope map, scout-1 cpl_fit):")
    print(f"  c0_loc = R(0) = {mp.nstr(R0, 15)}   ca_loc = R'(0) = {mp.nstr(dR0, 15)}   ratio = {mp.nstr(dR0 / R0, 15)}")
    print(f"  (analytic: R(0)=1/2, R'(0) = 3*OM/4 = {mp.nstr(3 * OM / 4, 15)}, ratio = 3*OM/2 = {mp.nstr(3 * OM / 2, 15)})")
    print("  note: the value 1.37 is NOT produced by either map; it is not found in the frozen scout-1 record.")

    # Sign facts on the projection domain and beyond
    print("\nSign / range facts (model only):")
    print(f"  R(0) = {mp.nstr(R0, 10)}, R(2) = {mp.nstr(R(ZMAX), 10)}, R(1100) = {mp.nstr(R(mp.mpf(1100)), 12)}, R(z->-1) = {mp.nstr(OL/(1+OL), 10)}")
    print("  R is strictly increasing in z (dR/dz = (dE2/dz)/(1+E2)^2 > 0) and 0 < R < 1 for all z > -1.")
    for eps in [mp.mpf("-0.1"), mp.mpf("-0.3")]:
        w0 = -1 + eps * c0; wa = eps * ca
        print(f"  eps={mp.nstr(eps,3)}: registered CPL (w0, wa) = ({mp.nstr(w0, 8)}, {mp.nstr(wa, 8)});"
              f" pointwise w on [0,2] in [{mp.nstr(-1 + eps * R(ZMAX), 8)}, {mp.nstr(-1 + eps * R0, 8)}]")
    print(f"  global-fit w0 > -1 possible on eps<0?  1+w0 = eps*c0 with c0 = {mp.nstr(c0, 10)} "
          f"{'> 0  => NO: w0 < -1 for every eps < 0' if c0 > 0 else '<= 0 => YES'}")
    print(f"  sign of wa on eps<0: wa = eps*ca, ca = {mp.nstr(ca, 10)} => wa {'< 0' if ca > 0 else '> 0'} for eps < 0")


if __name__ == "__main__":
    main()
