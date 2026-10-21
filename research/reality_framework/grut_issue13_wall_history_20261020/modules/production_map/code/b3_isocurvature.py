"""production_map / b3 [E]: isocurvature carried by the spectator relic (single-field spectral estimate).

Output: results/b3_isocurvature.json.  < 5 s.  Needs results/b1_spectator.json.

The spectator relic density in a Hubble patch is set by the local frozen field values phi_0 (b1: Y ~ amplitude^(3/2),
i.e. Y_i ~ |phi_0,i|^(3/2) per field for the condensate count, Y ~ (sum phi_0^2)^(3/4) for the floor).  It therefore
carries a superhorizon density perturbation uncorrelated with the inflaton curvature perturbation [D].

Spectral expansion of the Starobinsky-Yokoyama Fokker-Planck operator (per field, quartic V = lambda phi^4/4,
fields treated as independent):  Delta^2_f(k) = sum_n c_n^2 (2 Lambda_n) (k/(a_end H))^(2 Lambda_n),
c_n = <psi_0|f|psi_n> [R: Starobinsky-Yokoyama 1994; Markkanen-Rajantie-Stopyra-Tenkanen 2019, recalled, not fetched].
CMB scales: k/(a_end H) = exp(-N_*), N_* = 50-60 [P].
Limit used [R, retrospective, recalled not fetched]: Planck 2018 uncorrelated CDM isocurvature beta_iso <~ 0.04 with
P_zeta = 2.1e-9; for a relic fraction f of the CDM, f^2 P_S <~ 0.04 * 2.1e-9 / 0.96.
"""
import math

import numpy as np
from scipy.linalg import eigh_tridiagonal

from common import lphi, lX, lD, MP, Y_MAX_STRICT, load, dump

b1 = load("b1_spectator.json")
inst = b1["instant reheating, minimal couplings [C]"]
share = b1["conformal dynamics, instant reheating, minimal couplings (units H_e) [C]"]["with slow-roll velocity"][
    "share_J_by_field"]
Dc = 1 / (8 * math.pi ** 2)


def spectrum(lam, power, Nstar, nmodes=40):
    L = 14.0 / lam ** 0.25
    n = 6000
    x = np.linspace(-L, L, n)
    h = x[1] - x[0]
    U = lam * x ** 4 / 12            # v/3
    Up = lam * x ** 3 / 3
    Upp = lam * x ** 2
    W = Up ** 2 / (4 * Dc) - Upp / 2
    d = 2 * Dc / h ** 2 + W
    e = -Dc / h ** 2 * np.ones(n - 1)
    lam_n, psi = eigh_tridiagonal(d, e, select="i", select_range=(0, nmodes - 1))
    psi /= np.sqrt(h)
    psi0 = psi[:, 0] * np.sign(psi[n // 2, 0])
    f = np.abs(x) ** power
    mean = np.sum(psi0 ** 2 * f) * h
    c = np.array([np.sum(psi0 * f * psi[:, k]) * h for k in range(nmodes)])
    D2 = sum(c[k] ** 2 * 2 * lam_n[k] * math.exp(-2 * lam_n[k] * Nstar) for k in range(1, nmodes))
    return dict(Lambda1=float(lam_n[1]), Lambda2=float(lam_n[2]), Lambda1_over_sqrt_lambda=float(lam_n[1] / math.sqrt(lam)),
                mean_f=float(mean), P_S=float(D2 / mean ** 2))


out = {"module": "production_map", "part": "B3 isocurvature of the spectator relic [E]"}
rows = {}
for Ns in (50.0, 60.0):
    ph = spectrum(lphi, 1.5, Ns)
    xx = spectrum(lX, 1.5, Ns)          # X treated as one quartic real direction (signed radius) [E]
    dd = spectrum(lD, 1.5, Ns)
    PS = share[0] ** 2 * ph["P_S"] + share[1] ** 2 * xx["P_S"] + share[2] ** 2 * dd["P_S"]
    lim = 0.04 * 2.1e-9 / 0.96
    fmax = math.sqrt(lim / PS)
    Ymax_iso = fmax * Y_MAX_STRICT        # relic fraction f of the CDM <-> Y_N <= f * Y_max (m_X floor)
    YJ = inst["Y_J = coef (H/M_P)^(3/2) (condensate quanta conserved) [C]"]
    YF = inst["Y_floor = coef (H/M_P)^(3/2) (mu=0 number at the same energy) [C]"]
    rows[f"N*={Ns:g}"] = dict(phi=ph, X=xx, D=dd, P_S_total=PS, f_max=fmax, Y_N_max_iso=Ymax_iso,
                              H_thr_floor=MP * (Ymax_iso / YF) ** (2 / 3), H_thr_J=MP * (Ymax_iso / YJ) ** (2 / 3))
out["rows"] = rows
out["check: Lambda1/sqrt(lambda) for quartic = 0.08892 [R]"] = rows["N*=50"]["phi"]["Lambda1_over_sqrt_lambda"]
out["caveats"] = ("independent single-field FP per species (the phi-X coupling is ignored), f = |phi|^(3/2) per field, "
                  "equilibrium start, P_S evaluated at k = a_end H exp(-N_*); order-of-magnitude only")
dump("b3_isocurvature.json", out)
if __name__ == "__main__":
    import json
    print(json.dumps(out, indent=1, default=float))
