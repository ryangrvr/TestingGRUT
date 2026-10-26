#!/usr/bin/env python3
"""kin_rates: independent brute-force check of the X Xb -> X Xb pole treatment.
For an enlarged width Gamma' (so that a direct quadrature can resolve the Breit-Wigner), compute
  gamma_full  = (T/32pi^4) int ds lambda sigma_full(s) K1/sqrt s   with scipy.integrate.quad (breakpoints at the pole)
  gamma_RIS   = decay-rate density with the same Gamma'  (independent formula)
and compare gamma_full - gamma_RIS with thermal_rate(mode='sub') at the same width; also the eps -> 0 trend."""
import json, sys
import numpy as np
from scipy import integrate
from scipy.special import kve
import model as M
import xsec22 as X
from amp import kallen


def direct_full(T, P):
    ms2 = P.ms2; eps = P.eps_s
    smin = 4 * P.mX2
    smax = (np.sqrt(smin) + 40 * T)**2

    def f(s):
        s0, s1, s2 = X.sigma_parts(s, P)
        x = s - ms2
        D = 1 / (x + 1j * eps)
        sig = s0 + s1 * D.real + s2 * abs(D)**2
        return kallen(s, P.mX2, P.mX2) * sig * kve(1, np.sqrt(s) / T) * np.exp(-np.sqrt(s) / T) / np.sqrt(s)
    pts = [ms2 + k * eps for k in (-300, -30, -3, 0, 3, 30, 300)]
    tot = 0.0
    edges = [smin] + sorted(p for p in pts if smin < p < smax) + [ms2 + 1e4 * eps, smax]
    edges = sorted(set(e for e in edges if smin <= e <= smax))
    for a, b in zip(edges[:-1], edges[1:]):
        tot += integrate.quad(f, a, b, limit=400, epsrel=1e-10)[0]
    return T / (32 * np.pi**4) * tot


def main():
    out = {}
    # The identity on-shell BW = production x BR needs the physical width of the coupling used, so the width is
    # enlarged by scaling lam3 -> sqrt(fac) lam3 (Gamma_s' = fac Gamma_s), not by hand.
    for fac in (1000.0, 100.0, 10.0, 1.0):
        P = X.Params(lam3=np.sqrt(fac) * M.LAM3)
        rows = []
        for t in (0.3, 1.0, 3.0, 12.0):
            T = t * M.MX
            full = direct_full(T, P)
            ris = X.decay_rate_density(T, P)
            sub = X.thermal_rate('XXb_XXb', T, P, mode='sub')
            fullq = X.thermal_rate('XXb_XXb', T, P, mode='full')
            rows.append(dict(T_over_m=t, gamma_full_direct=full, gamma_full_quadmodule=fullq, gamma_RIS=ris,
                             full_minus_RIS=full - ris, gamma_sub_module=sub,
                             rel_diff_sub=(full - ris - sub) / (abs(sub) + 1e-300), rel_diff_full=(full - fullq) / full))
            print(f'Gamma={fac:g}xGamma_s T/m={t}: full {full:.6e} (module {fullq:.6e}) RIS {ris:.6e} full-RIS {full-ris:.4e} sub {sub:.4e} rel {rows[-1]["rel_diff_sub"]:.2e}')
        out[f'width_x{fac:g}'] = rows
    P = X.Params()
    out['physical_width_sub'] = {str(t): X.thermal_rate('XXb_XXb', t * M.MX, P, mode='sub') for t in (0.3, 1.0, 3.0, 12.0)}
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else '/dev/stdout', 'w'), indent=1)


if __name__ == '__main__':
    main()
