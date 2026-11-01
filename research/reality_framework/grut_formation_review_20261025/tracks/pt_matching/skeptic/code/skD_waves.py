"""Skeptic part D: spot check of the wave switch-on claim (pt_matching 3j) with scipy DOP853 (independent code).
Linear test modes chi_k'' = -(k^2 + m_w^2 + g^2 phi(s)^2) chi_k (actual phi) or g^2 M(s)^2 (smooth envelope), m_w^2 = 20,
driven by the homogeneous X+phi run; |beta_k|^2 = n_u + n_v - 1/2 from two real solutions started in the adiabatic vacuum,
averaged over 20 wave periods at M = 2 M_c,peri.  Also: the same modes with the coherent-orbit quartic modulation
2 lam_X r(s)^2 added (an isotropic stand-in for the omitted X-background term), to show that the reference is not adiabatic
before the transition either."""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sk_common import *
from scipy.integrate import solve_ivp

KS = np.array([1.0, 6.0]); MW2 = 20.0
t0 = time.time()
y4 = tau20_state(); info = orbit_info(y4, nsamp=1024)
Mperi = np.sqrt(G2 / LP) * info['rmin']; M0 = 0.5 * Mperi


def run(eps, ph, variant):
    kappa = eps * info['omega_r']
    st = phase_state(y4, info['T_r'], ph)
    nk = len(KS)
    def wmass(s, y):
        M = M0 * np.exp(kappa * s); r2 = y[0] ** 2 + y[1] ** 2
        base = KS ** 2 + MW2
        if variant == 'smooth': return base + G2 * M * M
        if variant == 'actual': return base + G2 * y[2] ** 2
        if variant == 'orbit_only': return base + 2 * LX * r2       # X-background modulation, phi frozen at 0
    def rhs(s, y):
        M = M0 * np.exp(kappa * s)
        r2 = y[0] ** 2 + y[1] ** 2; ph_ = y[2]
        if variant == 'orbit_only': ph_ = 0.0
        cx = LX * r2 + G2 * ph_ ** 2
        d = np.empty_like(y)
        d[0:3] = y[3:6]
        d[3] = -cx * y[0]; d[4] = -cx * y[1]
        d[5] = -(G2 * r2 + LP * ph_ ** 2 - LP * M * M) * ph_ if variant != 'orbit_only' else 0.0
        w2 = wmass(s, y)
        W = y[6:].reshape(nk, 4)          # chi_u, pi_u, chi_v, pi_v
        D = np.empty_like(W)
        D[:, 0] = W[:, 1]; D[:, 1] = -w2 * W[:, 0]; D[:, 2] = W[:, 3]; D[:, 3] = -w2 * W[:, 2]
        d[6:] = D.ravel()
        return d
    y0 = np.r_[st[0], st[1], 1e-6 * Mperi, st[2], st[3], 0.0, np.zeros(4 * nk)]
    w0 = np.sqrt(wmass(0.0, y0))
    W0 = np.zeros((nk, 4)); W0[:, 0] = 1 / np.sqrt(2 * w0); W0[:, 3] = np.sqrt(w0 / 2)
    y0[6:] = W0.ravel()
    s1 = np.log(2 * Mperi / M0) / kappa
    Tw = 20 * 2 * np.pi / np.sqrt(MW2 + G2 * (2 * Mperi) ** 2)
    ts = np.linspace(s1, s1 + Tw, 41)
    sol = solve_ivp(rhs, (0, s1 + Tw), y0, method='DOP853', rtol=1e-10, atol=1e-12, t_eval=ts)
    acc = np.zeros(nk)
    for j in range(len(ts)):
        y = sol.y[:, j]; w = np.sqrt(wmass(ts[j], y)); W = y[6:].reshape(nk, 4)
        nu = (W[:, 1] ** 2 + w ** 2 * W[:, 0] ** 2) / (2 * w); nv = (W[:, 3] ** 2 + w ** 2 * W[:, 2] ** 2) / (2 * w)
        acc += nu + nv - 0.5
    return acc / len(ts)


out = {}
for eps in [1e-3]:
    for variant in ['smooth', 'actual', 'orbit_only']:
        b = np.array([run(eps, ph, variant) for ph in np.arange(8) / 8])
        out[f'{variant}_eps{eps:g}'] = dict(k=KS, mean=b.mean(0), min=b.min(0), max=b.max(0))
        print(variant, eps, 'k', KS, 'mean |beta|^2', b.mean(0), 'min', b.min(0), 'max', b.max(0), round(time.time() - t0, 1))
save_json('skD_waves.json', out)
