#!/usr/bin/env python3
"""kin_rates (Claude, own code): one-loop thermal masses of X and s and the phase of the radial field.

Finite-temperature one-loop effective potential (Dolan-Jackiw) for the background phi, with the X doublet and the
radial fluctuation in the loop:
    V(phi,T) = lP/4 (phi^2-v^2)^2 + (T^4/2pi^2) [ J_B(m_s^2(phi)/T^2) + 2 J_B(m_X^2(phi)/T^2) ],
    m_s^2(phi) = 3 lP phi^2 - lP v^2,   m_X^2(phi) = g^2 phi^2,
    J_B(y) = int_0^inf dx x^2 ln(1 - exp(-sqrt(x^2+y))).
The zero-temperature Coleman-Weinberg piece is absorbed in the renormalised (vacuum) parameters; its residual
field dependence is of relative order g^4/(16 pi^2) ~ 6e-6 and lX^2/(16 pi^2) and is dropped [A, size stated].
Negative m_s^2(phi) near phi = 0 is handled by the analytic continuation of J_B (real part), as usual.

Thermal masses are the curvatures of V at the thermal minimum, including the tadpole (J_B') and the
field-derivative-squared (J_B'') pieces:
    m_s^2(T)  = d^2V/dphi^2 |_{v_T}
    m_X^2(T)  = g^2 v_T^2 + (T^2/2pi^2) [ (6 lX + 2 lX) J_B'(y_X) + 2 g^2 J_B'(y_s) ]
      (x1-curvature: the x1 and x2 loops with lX and the s loop with g^2 s^2 |X|^2; high-T limit
       T^2 (lX/3 + g^2/12)).
High-T limits: v_T^2 = v^2 - T^2 (lP/4 + g^2/6)/lP,  T_c^2 = lP v^2 / (lP/4 + g^2/6) at this order.
The in-medium trilinear is lam3(T) = 2 g^2 v_T, the cubic 6 lP v_T; the in-medium width of s uses the in-medium
masses and the Bose enhancement 1 + f(E1) + f(E2).
"""
import numpy as np
from scipy import integrate, optimize
from scipy.special import kv
import model as M

TWO_PI2 = 2 * np.pi**2


EULER = 0.5772156649015329
LOG16PI2 = np.log(16 * np.pi**2)


def JB(y):
    """Bosonic thermal function J_B(y).  y >= 0: direct quadrature.  y < 0 (tachyonic background direction, only
    met for phi < v/sqrt3 in the phase search): real part of the Dolan-Jackiw high-T series
    -pi^4/45 + pi^2 y/12 - (pi/6) Re y^{3/2} - (y^2/32)[ln(|y|/16pi^2) + 2 gamma_E - 3/2] + O(y^3), adequate
    for |y| << 10 (the only negative y used are |y| <= lP v^2/T^2 < 0.03 at the temperatures of interest)."""
    y = float(y)
    if y >= 0:
        f = lambda x: x * x * np.log1p(-np.exp(-np.sqrt(x * x + y)))
        return integrate.quad(f, 0, 50 + np.sqrt(max(y, 0)), limit=200)[0]
    return -np.pi**4 / 45 + np.pi**2 * y / 12 - y * y / 32 * (np.log(-y) - LOG16PI2 + 2 * EULER - 1.5)


def JB1(y):
    """J_B'(y) = (1/2) int dx x^2 n_B(w)/w, w = sqrt(x^2+y)  (tadpole); J_B'(0) = pi^2/12.  y<0: series."""
    y = float(y)
    if y > 0:
        f = lambda x: 0.5 * x * x / (np.sqrt(x * x + y) * np.expm1(np.sqrt(x * x + y)))
        return integrate.quad(f, 0, 50 + np.sqrt(y), limit=200)[0]
    if y == 0:
        return np.pi**2 / 12
    return np.pi**2 / 12 - y / 16 * (np.log(-y) - LOG16PI2 + 2 * EULER - 1.0)


def JB2(y, h=None):
    """J_B''(y) by central difference for y > 0 (IR piece of the curvature); series for y <= 0."""
    y = float(y)
    if y <= 0:
        return -(np.log(max(-y, 1e-300)) - LOG16PI2 + 2 * EULER) / 16 if y < 0 else -1e3
    h = h or max(1e-3, 1e-2 * abs(y))
    h = min(h, 0.5 * y)
    return (JB1(y + h) - JB1(y - h)) / (2 * h)


def V_eff(phi, T, v=M.V, g=M.GX, lP=M.LP):
    ms2 = 3 * lP * phi**2 - lP * v**2
    mx2 = g**2 * phi**2
    V0 = lP / 4 * (phi**2 - v**2)**2
    if T <= 0:
        return V0
    return V0 + T**4 / TWO_PI2 * (JB(ms2 / T**2) + 2 * JB(mx2 / T**2))


def dV_dphi(phi, T, v=M.V, g=M.GX, lP=M.LP):
    ms2 = 3 * lP * phi**2 - lP * v**2
    mx2 = g**2 * phi**2
    d0 = lP * phi * (phi**2 - v**2)
    if T <= 0:
        return d0
    return d0 + T**2 / TWO_PI2 * (6 * lP * phi * JB1(ms2 / T**2) + 2 * 2 * g**2 * phi * JB1(mx2 / T**2))


def vev(T, v=M.V, g=M.GX, lP=M.LP):
    """Broken-phase minimum v_T (0 if none).  Bracket between the high-T estimate and v."""
    if T <= 0:
        return v
    hiT = v**2 - T**2 * (lP / 4 + g**2 / 6) / lP
    lo = 1e-4 * v
    hi = 1.05 * v
    f = lambda p: dV_dphi(p, T, v, g, lP)
    # scan for sign change
    grid = np.linspace(lo, hi, 60)
    vals = [f(p) for p in grid]
    for i in range(len(grid) - 1):
        if vals[i] < 0 and vals[i + 1] > 0:
            return optimize.brentq(f, grid[i], grid[i + 1], xtol=1e-14)
    return 0.0


def thermal_masses(T, v=M.V, g=M.GX, lX=M.LX, lP=M.LP):
    """Return dict with v_T, m_X(T), m_s(T), the pieces, and validity flags."""
    vT = vev(T, v, g, lP)
    if T <= 0:
        return dict(T=T, vT=vT, mX=g * vT, ms=np.sqrt(2 * lP) * vT, broken=True)
    T2 = T * T
    ms2_tree = 3 * lP * vT**2 - lP * v**2
    mx2_tree = g**2 * vT**2
    ys, yx = ms2_tree / T2, mx2_tree / T2
    j1s, j1x = JB1(ys), JB1(yx)
    # s curvature: d^2V/dphi^2 = 3 lP phi^2 - lP v^2 + (T^2/2pi^2)[6 lP J'(ys) + 4 g^2 J'(yx)]
    #                           + (1/2pi^2)[(6 lP phi)^2 J''(ys) + 2 (2 g^2 phi)^2 J''(yx)]
    tad_s = T2 / TWO_PI2 * (6 * lP * j1s + 4 * g**2 * j1x)
    j2s, j2x = JB2(ys), JB2(yx)
    ir_s = 1 / TWO_PI2 * ((6 * lP * vT)**2 * j2s + 2 * (2 * g**2 * vT)**2 * j2x)
    ms2 = ms2_tree + tad_s + ir_s
    # X curvature (x1 direction at x=0): g^2 vT^2 + (T^2/2pi^2)[(6+2) lX J'(yx) + 2 g^2 J'(ys)]; no J'' piece
    # (d m_i^2/d x1 = 0 at x1 = 0)
    pi_x = T2 / TWO_PI2 * (8 * lX * j1x + 2 * g**2 * j1s)
    mx2 = mx2_tree + pi_x
    return dict(T=T, vT=vT, broken=vT > 0, mX=np.sqrt(max(mx2, 0)), ms=np.sqrt(max(ms2, 0)),
                mX2_tree=mx2_tree, PiX=pi_x, ms2_tree=ms2_tree, tad_s=tad_s, ir_s=ir_s,
                lam3=2 * g**2 * vT, lam_s3=6 * lP * vT,
                hiT_vT=np.sqrt(max(v**2 - T2 * (lP / 4 + g**2 / 6) / lP, 0)),
                hiT_mX=np.sqrt(max(g**2 * (v**2 - T2 * (lP / 4 + g**2 / 6) / lP) + T2 * (lX / 3 + g**2 / 12), 0)))


def Tc_oneloop(v=M.V, g=M.GX, lP=M.LP):
    """Temperature where the broken minimum ceases to exist (one-loop, no daisy), and the high-T estimate."""
    Thi = np.sqrt(lP * v**2 / (lP / 4 + g**2 / 6))
    f = lambda T: vev(T, v, g, lP) - 1e-3 * v
    # the minimum disappears somewhat below Thi (cubic term); scan
    Ts = np.linspace(0.7 * Thi, 1.1 * Thi, 41)
    vs = [vev(T, v, g, lP) for T in Ts]
    Tdis = None
    for i in range(len(Ts) - 1):
        if vs[i] > 1e-3 * v and vs[i + 1] <= 1e-3 * v:
            Tdis = optimize.brentq(f, Ts[i], Ts[i + 1], xtol=1e-6)
            break
    # degeneracy temperature V(vT) = V(0)
    def dF(T):
        vT = vev(T, v, g, lP)
        return V_eff(vT, T, v, g, lP) - V_eff(0.0, T, v, g, lP)
    Tdeg = None
    if Tdis:
        Ts2 = np.linspace(0.7 * Thi, Tdis, 41)
        d = [dF(T) for T in Ts2]
        for i in range(len(Ts2) - 1):
            if d[i] < 0 and d[i + 1] > 0:
                Tdeg = optimize.brentq(dF, Ts2[i], Ts2[i + 1], xtol=1e-6)
                break
    return dict(T_hiT_estimate=Thi, T_minimum_disappears=Tdis, T_degenerate=Tdeg)


def bose(E, T, mu=0.0):
    x = (E - mu) / T
    return 1 / np.expm1(x) if x < 600 else 0.0


def width_s_vac(ms, mx, lam3):
    if ms <= 2 * mx:
        return 0.0
    return lam3**2 * np.sqrt(1 - 4 * mx**2 / ms**2) / (16 * np.pi * ms)


def width_s_medium(E, T, mu_p=0.0, mu_m=0.0, ms=M.MS, mx=M.MX, lam3=M.LAM3, n=64):
    """In-medium decay rate of an s of energy E (lab frame): (ms/E) Gamma_vac(ms,mx,lam3) <1 + f(E1) + f(E2)>,
    averaged over the decay angle; E1, E2 are the daughter lab energies, mu_p (X) and mu_m (Xbar) the chemical
    potentials.  The inverse-decay rate per unit volume follows by detailed balance (see thermal_avg.py)."""
    g0 = width_s_vac(ms, mx, lam3)
    if g0 == 0:
        return 0.0, 0.0
    p = np.sqrt(max(E * E - ms * ms, 0))
    pst = np.sqrt(ms * ms / 4 - mx * mx)
    c = np.linspace(-1, 1, n)
    E1 = (E * ms / 2 + p * pst * c) / ms
    E2 = E - E1
    fac = 1 + np.array([bose(e, T, mu_p) for e in E1]) + np.array([bose(e, T, mu_m) for e in E2])
    enh = np.trapezoid(fac, c) / 2
    return g0 * ms / E * enh, enh


def table(Ts_over_m):
    rows = []
    for t in Ts_over_m:
        T = t * M.MX
        d = thermal_masses(T)
        d['T_over_mX'] = t
        d['Gamma_s_medium_rest_over_vac'] = (width_s_medium(d['ms'], T, ms=d['ms'], mx=d['mX'], lam3=d['lam3'])[0]
                                             / M.GAMMA_S_VAC if d['ms'] > 2 * d['mX'] else 0.0)
        d['Gamma_s_vac_with_medium_masses'] = width_s_vac(d['ms'], d['mX'], d['lam3'])
        d['bose_enh_rest'] = width_s_medium(d['ms'], T, ms=d['ms'], mx=d['mX'], lam3=d['lam3'])[1] if d['ms'] > 2 * d['mX'] else 0.0
        rows.append({k: (float(v) if isinstance(v, (float, np.floating, int)) else bool(v)) for k, v in d.items()})
    return rows


if __name__ == '__main__':
    import json, sys
    print('J_B(0) =', JB(0), ' expected -pi^4/45 =', -np.pi**4 / 45)
    print("J_B'(0) =", JB1(0), ' pi^2/12 =', np.pi**2 / 12)
    print("J_B'(1e-6) =", JB1(1e-6))
    print('Tc:', Tc_oneloop())
    for t in (0.3, 1, 3, 6, 8, 12, 15, 17, 18, 20):
        d = thermal_masses(t * M.MX)
        print(f"T/m={t:5.1f}  vT/v={d['vT']/M.V:.4f} (hiT {d['hiT_vT']/M.V:.4f})  mX/mX0={d['mX']/M.MX:.4f} (hiT {d['hiT_mX']/M.MX:.4f})"
              f"  ms/ms0={d['ms']/M.MS:.4f}  PiX/mX0^2={d['PiX']/M.MX**2:.4f}  ir_s/ms0^2={d['ir_s']/M.MS**2:.2e}  broken={d['broken']}")
