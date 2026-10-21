"""Adversarial scheme test: a self-consistent (gap-equation) soft-phi daisy on top of the DR potential.
V_SC(phi,T) = V_DR(phi,T) - (T/12pi) Re(M^2 + 3 lp phi^2)^{3/2},  M^2 = c_DR(T) - (3 lp T/4pi) M,
so the curvature at phi=0 equals M^2 (self-consistent soft mass).  Tests whether 'second order in all
resummed schemes' is scheme-robust.  Writes results/gap_scheme_v.json."""
import json
import os
import numpy as np
from scipy.optimize import brentq
import vt_core as V

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
P = V.BENCH
MU = V.MU0
lp = P['lp']
PI = np.pi


def cDR(T):
    return 2 * float(V.dVds(np.array([0.0]), T, P, MU, 'DR')[0])


def Msc(T):
    b = 3 * lp * T / (4 * PI)
    c = cDR(T)
    if c <= 0:
        return 0.0
    return (-b + np.sqrt(b * b + 4 * c)) / 2


def dVds_SC(phi, T):
    M2 = Msc(T)**2
    g = V.dVds(phi, T, P, MU, 'DR')
    return g - T / (8 * PI) * 3 * lp * V.re_pow12(M2 + 3 * lp * np.asarray(phi)**2)


def V_SC(phi, T):
    M2 = Msc(T)**2
    return V.Veff(phi, T, P, MU, 'DR') - T / (12 * PI) * V.re_pow32(M2 + 3 * lp * np.asarray(phi)**2)


def mins(T):
    ph = np.geomspace(1e-7, 0.5, 20000)
    g = dVds_SC(ph, T)
    idx = np.where((g[:-1] < 0) & (g[1:] >= 0))[0]
    return [brentq(lambda x: float(dVds_SC(np.array([x]), T)[0]), ph[i], ph[i + 1], xtol=1e-16) for i in idx]


TDR = brentq(cDR, 1.5, 1.6, xtol=1e-14)
out = dict(T_DR=TDR)
Ts = TDR + np.concatenate([[1e-9, 1e-7, 1e-6], np.linspace(1e-5, 3e-3, 300)])
co = [(T, mins(T)) for T in Ts]
co = [(T, m) for T, m in co if m]
out['n_coexist'] = len(co)
if co:
    T1 = co[-1][0]

    def dV(T):
        m = max(mins(T))
        return float(V_SC(np.array([m]), T)[0] - V_SC(np.array([0.0]), T)[0]), m
    Tc = brentq(lambda T: dV(T)[0], TDR + 1e-12, T1, xtol=1e-14)
    out.update(order='first (spurious, perturbative)', Tc=Tc, phic_over_Tc=dV(Tc)[1] / Tc, T1_approx=T1,
               M_over_T_at_Tc=Msc(Tc) / Tc, window_over_Tc=(T1 - TDR) / Tc)
else:
    out.update(order='second', Tc=TDR)
print(out)
with open(os.path.join(OUT, 'gap_scheme_v.json'), 'w') as f:
    json.dump(out, f, indent=1, default=float)
