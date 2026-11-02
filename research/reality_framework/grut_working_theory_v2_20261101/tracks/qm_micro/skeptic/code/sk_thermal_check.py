"""qm_micro SKEPTIC (Claude, 2 Nov 2026): independent quadrature of the thermal-flux sticking average on qm_micro's
verified s(E, mu) table (stick_thermal.json; point values reproduced by sk_scatter.py to 1e-9).  Gauss-Laguerre in
the kinetic energy (weight e^{-eps/T}) with p^2 = 2 eps + eps^2, monotone (PCHIP) interpolation in ln eps, and the
table's own 16 Gauss-Legendre mu nodes, plus a check of the mu quadrature against a 64-point PCHIP re-integration."""
import json, math, os, sys
import numpy as np
from scipy.interpolate import PchipInterpolator
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, '..', '..', 'results', 'stick_thermal.json')))
Tk = np.array(d['T_kin']); mu = np.array(d['mu_nodes']); mw = np.array(d['mu_weights'])
sX = np.array(d['s_X']); dB = np.array(d['d_Xbar'])
ref = {('%g' % r['T_h']): r for r in d['thermal']}
xl, wl = np.polynomial.laguerre.laggauss(80)


def avg(tab, Th, fine_mu=False):
    eps = Th * xl
    p2 = 2 * eps + eps * eps
    # interpolate each mu column in ln eps (sqrt law below the table, constant above)
    S = np.empty((eps.size, mu.size))
    for j in range(mu.size):
        f = PchipInterpolator(np.log(Tk), tab[:, j])
        le = np.log(np.clip(eps, Tk[0], Tk[-1]))
        S[:, j] = f(le)
        lo = eps < Tk[0]
        S[lo, j] = tab[0, j] * np.sqrt(eps[lo] / Tk[0])
    if fine_mu:
        m2 = np.linspace(1e-4, 1, 4001)
        ang = np.array([np.trapezoid(PchipInterpolator(np.r_[0.0, mu], np.r_[0.0, S[i]])(m2) * m2, m2) / 0.5 for i in range(eps.size)])
    else:
        ang = (S * mu * mw).sum(1) / (mu * mw).sum()
    return float((wl * p2 * ang).sum() / (wl * p2).sum())


rows = []
for Th in (0.005, 0.01, 0.1, 0.3, 0.5, 1.0):
    r = dict(T_h=Th, s_X=avg(sX, Th), s_X_fine_mu=avg(sX, Th, True), d_Xbar=avg(dB, Th), d_Xbar_fine_mu=avg(dB, Th, True),
             qm_s_X=ref['%g' % Th]['s_X_boltz'], qm_d_Xbar=ref['%g' % Th]['d_Xbar_boltz'])
    rows.append(r)
    print('T_h=%.3g  s_X %.5f (fine-mu %.5f; qm_micro %.5f) | Xbar destroyed %.5f (fine-mu %.5f; qm_micro %.5f)' % (
        Th, r['s_X'], r['s_X_fine_mu'], r['qm_s_X'], r['d_Xbar'], r['d_Xbar_fine_mu'], r['qm_d_Xbar']))
json.dump(dict(notes=__doc__, rows=rows), open(sys.argv[1], 'w'), indent=1)
