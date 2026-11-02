"""qm_micro (Claude, 2 Nov 2026): bulk Q-matter equation of state at T = 0 near Omega_b, phonon speed and its
density exponent u = dln c_s/dln n, and the Beliaev three-phonon decay estimate used for the trapping factor.  Own code.
P0(mu) = 100 max_S [(mu^2 - S^2)^2/(4a) - Uphi(S)] (F^2 = (mu^2 - S^2)/a eliminated), n = dP/dmu, e = mu n - P,
c_s^2 = dP/de = n/(mu dn/dmu).  Beliaev [E]: Gamma(k) = (1+u)^2 k^5 / (480 pi (e+P)), the universal superfluid
three-phonon result normalised to the dilute-BEC 3k^5/(640 pi m n) (u = 1/2); used only as an order of magnitude."""
import json, math, os, sys
import numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import qm_common as qc
import scatter1d as sc
a = qc.a_


def P0(mu):
    f = lambda S: -S * (mu * mu - S * S) / a - qc.Uphi_p(S)
    S = brentq(f, 0.01, 0.6, xtol=1e-15)
    return 100 * ((mu * mu - S * S) ** 2 / (4 * a) - qc.Uphi(S)), S


def main(out):
    B = qc.bulk(); Om = B['Omega']
    h = 1e-4
    mus = Om + h * np.arange(-3, 4)
    P = np.array([P0(m)[0] for m in mus])
    n = np.gradient(P, h); dn = np.gradient(n, h)
    i = 3
    cs2 = n[i] / (mus[i] * dn[i])
    c = np.sqrt(n / (mus * dn))
    u = float(np.gradient(np.log(c), np.log(n))[i])
    # phonon speed from the Bogoliubov dispersion at small w
    br = sc.bulk_branches(Om, B['S'], B['F'], 1e-3)
    c_disp = 1e-3 / br[0]['kx']
    eP = n[i] * Om
    rows = []
    for T in (0.005, 0.0075, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0):
        w = (1 - Om) + 1.5 * T
        b = sc.bulk_branches(Om, B['S'], B['F'], w)
        k = b[0]['kx']; vg = b[0]['vg']
        G = (1 + u) ** 2 * k ** 5 / (480 * math.pi * eP)
        rows.append(dict(T=T, w=w, k_phonon=k, v_g=vg, Gamma_Beliaev=G, mfp=vg / G))
    res = dict(P0_at_Omega_b=float(P[i]), n_B=float(n[i]), n_B_ref=B['n'], c_s_eos=float(math.sqrt(cs2)), c_s_dispersion=c_disp, u=u,
               e_plus_P=eP, rows=rows)
    json.dump(res, open(out, 'w'), indent=1, default=float)
    print(json.dumps({k: v for k, v in res.items() if k != 'rows'}, indent=1))
    for r in rows:
        print('T=%.4g w=%.4f k=%.4f v_g=%.3f Gamma_B=%.3e mfp=%.3e' % (r['T'], r['w'], r['k_phonon'], r['v_g'], r['Gamma_Beliaev'], r['mfp']))


if __name__ == '__main__':
    main(sys.argv[1])
