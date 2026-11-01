"""[fx_formation copy of tracks/m1_vacuum/code/coldfluid_m1.py; only change: main guard] Track m1_vacuum, item 3: linear stability of the homogeneous rotating condensate + phi in the M1 potential.
Independent code (Claude, 1 Nov 2026).  Derivation (NOTES.md): with X/v = (F + psi1 + i psi2) e^{-i Om t}/sqrt2,
S = S0 + s, perturbations ~ e^{nu t + i k x}, the linearised equations are
    (nu^2 + k^2 + 2 a F^2) psi1 + 2 Om nu psi2 + 2 F S s = 0
    (nu^2 + k^2) psi2 - 2 Om nu psi1 = 0
    (nu^2 + k^2 + D) s + 2 F S psi1 = 0,          D = U_phi''(S) + F^2,
whose determinant is a cubic in z = nu^2.  Growth rate Gamma = max Re nu over the roots.
"""
import json, math, sys
import numpy as np
from scipy.optimize import brentq

a, b, c = 0.1, 3.0, 0.2
G2 = 0.01


def Uphi_p(S):
    return b * (S * S - 1) * S - 3 * c * (S - 1)**2


def Uphi_pp(S):
    return b * (3 * S * S - 1) - 6 * c * (S - 1)


def S_of_F(F):
    return brentq(lambda s: Uphi_p(s) + F * F * s, 1e-14, 1.0, xtol=1e-15)


def cubic_coeffs(k2, F, S, Om, D):
    # det = (z+k2+2aF^2)(z+k2)(z+k2+D) + 4 Om^2 z (z+k2+D) - 4 F^2 S^2 (z+k2); build by polynomial multiplication
    p1 = np.polymul([1, k2 + 2 * a * F * F], [1, k2])
    p = np.polymul(p1, [1, k2 + D])
    p = np.polyadd(p, 4 * Om * Om * np.polymul([1, 0], [1, k2 + D]))
    p = np.polysub(p, 4 * F * F * S * S * np.array([1, k2]))
    return p


def growth(k, F, S, Om, D):
    z = np.roots(cubic_coeffs(k * k, F, S, Om, D))
    nu = np.sqrt(z.astype(complex))
    return float(max(np.max(nu.real), np.max((-nu).real), 0.0)), z


def analyse(F, ks):
    S = S_of_F(F)
    Om = math.sqrt(S * S + a * F * F)
    D = Uphi_pp(S) + F * F
    g = np.array([growth(k, F, S, Om, D)[0] for k in ks])
    ip = int(np.argmax(g))
    # refine the peak with a golden-section search between neighbours
    if 0 < ip < len(ks) - 1 and g[ip] > 0:
        lo, hi = ks[ip - 1], ks[ip + 1]
        for _ in range(60):
            m1 = lo + 0.381966 * (hi - lo); m2 = hi - 0.381966 * (hi - lo)
            if growth(m1, F, S, Om, D)[0] > growth(m2, F, S, Om, D)[0]:
                hi = m2
            else:
                lo = m1
        kp = 0.5 * (lo + hi); gp = growth(kp, F, S, Om, D)[0]
    else:
        kp, gp = float(ks[ip]), float(g[ip])
    # band edge: largest k with growth > 0
    kmax = float(ks[np.where(g > 0)[0][-1]]) if np.any(g > 0) else 0.0
    if kmax > 0 and kmax < ks[-1]:
        kmax = brentq(lambda k: growth(k, F, S, Om, D)[0] - 1e-12, kmax, kmax + (ks[1] - ks[0]), xtol=1e-10)
    # long-wave limit: put z = sigma k^2 in the determinant and keep O(k^2):
    #   D [2 a F^2 (sigma + 1) + 4 Om^2 sigma] = 4 F^2 S^2 (sigma + 1)
    #   sigma = F^2 (2 S^2 - a D) / [ D (a F^2 + 2 Om^2) - 2 F^2 S^2 ]   (sign set by 2S^2 - aD when the denominator > 0)
    k0 = 1e-3
    _, z0 = growth(k0, F, S, Om, D)
    zr = z0[np.argmin(np.abs(z0))]           # the small root
    lw_numeric = float(zr.real / k0**2)
    lw_denominator = D * (a * F * F + 2 * Om * Om) - 2 * F * F * S * S
    lw_formula = F * F * (2 * S * S - a * D) / lw_denominator
    # thermodynamic sound speed along the homogeneous family: c_s^2 = n dOm/dn / Om  (-> z/k^2 = -c_s^2)
    dF = 1e-5
    def OmF(Fx):
        Sx = S_of_F(Fx); return math.sqrt(Sx * Sx + a * Fx * Fx)
    nF = lambda Fx: OmF(Fx) * Fx * Fx
    dOm_dn = (OmF(F + dF) - OmF(F - dF)) / (nF(F + dF) - nF(F - dF))
    cs2 = nF(F) * dOm_dn / Om
    # determinant residual at the peak
    nu = complex(gp)
    M = np.array([[nu * nu + kp * kp + 2 * a * F * F, 2 * Om * nu, 2 * F * S], [-2 * Om * nu, nu * nu + kp * kp, 0], [2 * F * S, 0, nu * nu + kp * kp + D]])
    det_res = abs(np.linalg.det(M)) / max(1, np.linalg.norm(M)**3)
    n_over_mX3 = Om * F * F / G2
    return dict(F=F, S=S, Omega=Om, D=D, n_over_mX3=n_over_mX3, peak_growth_over_mX=gp, peak_k_over_mX=kp,
                unstable_band_kmax=kmax, longwave_z_over_k2_numeric=lw_numeric, longwave_z_over_k2_formula=lw_formula,
                minus_cs2_thermodynamic=-cs2, codex_criterion_aD_lt_2S2=bool(a * D < 2 * S * S),
                fastest_wavelength_charge=n_over_mX3 * (2 * math.pi / kp)**3 if kp > 0 and gp > 0 else None,
                det_residual_at_peak=float(det_res), growth_time_1_over_Gamma=1 / gp if gp > 0 else None)


# fx_formation copy: the module-level run below is wrapped in a main guard so the functions can be imported
if __name__ == '__main__':
    ks = np.linspace(1e-4, 3.0, 30001)
    rows = []
    for F in (0.1, 0.5, 1.0, 1.5, 1.8565, 2.0, 2.406476756777611, 2.63, 2.6309606):
        rw = analyse(F, ks)
        rows.append(rw)
        print('F=%.4f S=%.5f Om=%.6f n=%.2f  Gamma_peak=%.6f at k=%.5f  kmax=%.4f  z/k2 lw: num %.4e formula %.4e  -cs2 %.4e  crit=%s  n*lambda^3=%s' % (
            rw['F'], rw['S'], rw['Omega'], rw['n_over_mX3'], rw['peak_growth_over_mX'], rw['peak_k_over_mX'], rw['unstable_band_kmax'],
            rw['longwave_z_over_k2_numeric'], rw['longwave_z_over_k2_formula'], rw['minus_cs2_thermodynamic'], rw['codex_criterion_aD_lt_2S2'],
            rw['fastest_wavelength_charge']), flush=True)
    # the stability boundary in density: where the long-wave coefficient changes sign
    Fb = brentq(lambda F: (lambda S: 2 * S * S - a * (Uphi_pp(S) + F * F))(S_of_F(F)), 1.5, 2.3, xtol=1e-12)
    Sb = S_of_F(Fb)
    Omb = math.sqrt(Sb * Sb + a * Fb * Fb)
    out = dict(rows=rows, longwave_stability_boundary=dict(F=Fb, S=Sb, Omega=Omb, n_over_mX3=Omb * Fb * Fb / G2),
               codex_F1=dict(peak_growth_over_mX=0.1142789245915453, peak_k_over_mX=0.42, S=0.8212478761806628, Omega=0.8800273144233929),
               interpretation='growing modes exist only for n below the spinodal density where dOmega/dn < 0 (imaginary sound speed): '
                              'fragmentation of a dilute homogeneous condensate into dense objects (Kusenko-Shaposhnikov; Enqvist-McDonald), '
                              'not an instability of the bulk or of the finite objects; the bulk density and the object core are stable here')
    json.dump(out, open(sys.argv[1], 'w'), indent=1, default=float)
    print('longwave boundary', out['longwave_stability_boundary'])
