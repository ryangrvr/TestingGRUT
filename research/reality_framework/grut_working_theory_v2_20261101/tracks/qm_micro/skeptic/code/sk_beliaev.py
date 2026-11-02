"""qm_micro SKEPTIC (Claude, 2 Nov 2026): EOS exponent u, phonon kinematics and the Beliaev mean free path from own
code (sk_scatter.P_of for P(mu); bulk dispersion roots of det M = 0).  Checks qm_micro eos_phonon.json.
Also: Beliaev 1 -> 2 decay needs a convex phonon branch (w/k rising); checked on the dispersion."""
import json, math, sys
import numpy as np
from scipy.optimize import brentq
import sk_scatter as ss
B = ss.bulk(); Om = B['Omega']
h = 2e-4
mus = Om + h * np.arange(-3, 4)
P = np.array([ss.P_of(m)[0] * 100 for m in mus])
n = np.gradient(P, h, edge_order=2); dn = np.gradient(n, h, edge_order=2)
c = np.sqrt(n / (mus * dn))
u = float(np.gradient(np.log(c), np.log(n))[3])
cs = float(c[3]); eP = float(Om * n[3])
F, S = B['F'], B['S']
a = ss.a
def det(w, k):
    M = np.array([[S * S + 3 * a * F * F - Om * Om - w * w + k * k, 2 * Om * w, 2 * S * F],
                  [2 * Om * w, S * S + a * F * F - Om * Om - w * w + k * k, 0.0],
                  [2 * S * F, 0.0, ss.d2Uphi(S) + F * F - w * w + k * k]])
    return np.linalg.det(M)
def w_ph(k):
    ws = np.linspace(1e-7, 6.0, 3001)
    dv = np.array([det(w, k) for w in ws])
    i = int(np.argmax(dv[:-1] * dv[1:] < 0))
    return brentq(det, ws[i], ws[i + 1], args=(k,), xtol=1e-15)
def k_of_w(w):
    return brentq(lambda k: w_ph(k) - w, 1e-3, 4.5, xtol=1e-12)
rows = []
for T in (0.005, 0.01, 0.1, 0.3, 1.0):
    for lab, eps in (('1.5T', 1.5 * T), ('2.5T', 2.5 * T)):
        w = (1 - Om) + eps
        k = k_of_w(w)
        vg = (w_ph(k + 1e-5) - w_ph(k - 1e-5)) / 2e-5
        G = (1 + u) ** 2 * k ** 5 / (480 * math.pi * eP)
        Grel = (1 + u - cs * cs) ** 2 * k ** 5 / (480 * math.pi * eP)
        rows.append(dict(T=T, eps=lab, w=w, k=k, w_over_k=w / k, v_g=vg, mfp=vg / G, mfp_if_vertex_1pu_minus_cs2=vg / Grel))
print(json.dumps(dict(u=u, c_s=cs, e_plus_P=eP), indent=0))
for r in rows:
    print('T=%.3g eps=%s w=%.4f k=%.4f w/k=%.4f vg=%.4f mfp=%.3e (vertex (1+u-cs^2): %.3e)' % (r['T'], r['eps'], r['w'], r['k'], r['w_over_k'], r['v_g'], r['mfp'], r['mfp_if_vertex_1pu_minus_cs2']))
json.dump(dict(notes=__doc__, u=u, c_s=cs, e_plus_P=eP, rows=rows), open(sys.argv[1], 'w'), indent=1)
