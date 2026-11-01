"""m1_vacuum SKEPTIC, part 4: cold-fluid dispersion of the homogeneous rotating condensate + phi, recomputed from the
6x6 first-order linear system (numpy eigvals of the companion matrix M(k)), not from the cubic's roots as Codex and the
track did.  Independent code (Claude skeptic, 1 Nov 2026).  sk_algebra.py has shown symbolically that det(M - nu) is
exactly Codex's cubic in z = nu^2; here the growth rate Gamma(k) = max Re nu is taken straight from the matrix.
Also: the k -> 0 limit against the thermodynamic sound speed, the band edge, a wide k scan (0 < k < 20) for the dense
cases to make sure 'no growing mode' is not a scanned-band artefact, and the growth rate of the fastest mode at the
track's spinodal density.
"""
import json, math, sys
import numpy as np
from scipy.optimize import brentq, minimize_scalar

a, b, c = 0.1, 3.0, 0.2
G2 = 0.01


def Uphi_p(S):  return b * (S * S - 1) * S - 3 * c * (S - 1)**2
def Uphi_pp(S): return b * (3 * S * S - 1) - 6 * c * (S - 1)
def S_of_F(F):  return brentq(lambda s: Uphi_p(s) + F * F * s, 1e-14, 1.0, xtol=1e-15)


def Mk(k, F, S, Om, D):
    M = np.zeros((6, 6))
    M[0, 3] = M[1, 4] = M[2, 5] = 1
    M[3, 0] = -(k * k + 2 * a * F * F); M[3, 4] = -2 * Om; M[3, 2] = -2 * F * S
    M[4, 1] = -k * k; M[4, 3] = 2 * Om
    M[5, 2] = -(k * k + D); M[5, 0] = -2 * F * S
    return M


def gamma(k, F, S, Om, D):
    return float(np.max(np.linalg.eigvals(Mk(k, F, S, Om, D)).real))


def analyse(F, kmax=3.0, nk=6001):
    S = S_of_F(F); Om = math.sqrt(S * S + a * F * F); D = Uphi_pp(S) + F * F
    ks = np.linspace(1e-4, kmax, nk)
    g = np.array([gamma(k, F, S, Om, D) for k in ks])
    ip = int(np.argmax(g))
    if g[ip] > 1e-12 and 0 < ip < nk - 1:
        r = minimize_scalar(lambda k: -gamma(k, F, S, Om, D), bounds=(ks[ip - 1], ks[ip + 1]), method='bounded', options={'xatol': 1e-12})
        kp, gp = float(r.x), float(-r.fun)
        iz = np.where(g > 1e-12)[0][-1]
        kedge = brentq(lambda k: gamma(k, F, S, Om, D) - 1e-12, ks[iz], ks[min(iz + 1, nk - 1)], xtol=1e-12) if iz < nk - 1 else float(ks[-1])
    else:
        kp, gp, kedge = 0.0, 0.0, 0.0
    # wide scan for residual growth up to k = 20 (dense cases)
    wide = max(gamma(k, F, S, Om, D) for k in np.linspace(kmax, 20.0, 2000))
    # k -> 0: z/k^2 from the matrix eigenvalue closest to zero
    k0 = 1e-3
    ev = np.linalg.eigvals(Mk(k0, F, S, Om, D))
    z0 = (ev[np.argmin(np.abs(ev))])**2
    # thermodynamic c_s^2 = n dOmega/dn / Omega along the family (central differences in F)
    dF = 1e-5
    def OmF(Fx): Sx = S_of_F(Fx); return math.sqrt(Sx * Sx + a * Fx * Fx)
    nF = lambda Fx: OmF(Fx) * Fx * Fx
    cs2 = nF(F) * (OmF(F + dF) - OmF(F - dF)) / (nF(F + dF) - nF(F - dF)) / Om
    n = Om * F * F / G2
    return dict(F=F, S=S, Omega=Om, D=D, n_over_mX3=n, peak_growth_over_mX=gp, peak_k_over_mX=kp, band_edge_k=kedge,
                max_growth_k_in_3_to_20=wide, z_over_k2_at_k_1em3=float(z0.real / k0**2), minus_cs2=-cs2,
                charge_per_fastest_cube=n * (2 * math.pi / kp)**3 if kp > 0 else None, growth_time=1 / gp if gp > 0 else None,
                codex_criterion_aD_lt_2S2=bool(a * D < 2 * S * S))


rows = []
for F in (0.1, 0.5, 1.0, 1.5, 1.8566256, 2.0, 2.406476756777611, 2.63, 2.6309606):
    r = analyse(F); rows.append(r)
    print('F=%.4f S=%.5f Om=%.6f n=%.2f | Gamma_peak=%.7f at k=%.6f edge=%.5f | wide(3..20)=%.1e | z/k2=%.5e -cs2=%.5e | crit=%s | 1/Gamma=%s cube=%s' % (
        r['F'], r['S'], r['Omega'], r['n_over_mX3'], r['peak_growth_over_mX'], r['peak_k_over_mX'], r['band_edge_k'], r['max_growth_k_in_3_to_20'],
        r['z_over_k2_at_k_1em3'], r['minus_cs2'], r['codex_criterion_aD_lt_2S2'], r['growth_time'], r['charge_per_fastest_cube']), flush=True)
out = dict(rows=rows,
           codex=dict(F1=dict(peak=0.1142789245915453, k=0.42), F0p5=dict(peak=0.02922175341675933, k=0.235), F1p5=dict(peak=0.18378951009285316, k=0.4475), F0p1=dict(peak=0.0011662639793450452, k=0.0475)),
           track=dict(F1=dict(peak=0.11427959138281993, k=0.4192293174108015, edge=0.6152586429760601), F0p5=dict(peak=0.02922187730020253, k=0.234653258138088), F1p5=dict(peak=0.18379067182072045, k=0.4465622028183507)))
json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)
