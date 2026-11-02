"""fx_chain_skeptic (Claude, own code, written fresh): independent check of fx_chain item 3 (circular X orbit through
the M1 phi transition).  scipy DOP853 (not the copied C engine), instantaneous free-field antiparticle measure
n_- = [(p^2 + m^2 r^2)/(2 m) - q]/2 with m = g chi, averaged over the final window (M frozen at 40 M_c, 60 X periods).
Same model as pt_m1tilt mode 3: V = lX/4 r^4 + g^2/2 chi^2 r^2 + lP/4 (chi^2 - M^2)^2 - eps_t M (chi - M)^3,
M = M0 exp(kappa s), M0 = 0.5 M_c, kappa = eps_rate sqrt(6 lX) r0, q = 8117.8816 (tau-20 charge).
Also reports the eccentricity the start itself carries: the given r0 = (q/sqrt lX)^(1/3) is circular only for chi = 0,
and the tilted minimum chi0 at M0 shifts the radial equilibrium.  [C] homogeneous classical model."""
import sys, json, math, time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

LX, G2, LP = 1e-3, 1e-2, 3e-2
g = math.sqrt(G2)
Q = 8117.8816


def chi_min(M, r2, epst):
    f = lambda c: G2 * c * r2 + LP * c * (c * c - M * M) - 3 * epst * M * (c - M) ** 2
    return brentq(f, 0.0, M, xtol=1e-14)


def run(delta, eps_rate, epst=0.002, rtol=1e-11, ph=0.0):
    r0 = (Q / math.sqrt(LX)) ** (1 / 3)
    Mc = g * r0 / math.sqrt(LP); M0 = 0.5 * Mc; kap = eps_rate * math.sqrt(6 * LX) * r0
    Mend = 40 * Mc; s_end = math.log(Mend / M0) / kap
    ri = r0 * (1 + delta)
    c0 = chi_min(M0, ri * ri, epst)
    y0 = [ri, 0.0, 0.0, Q / ri, c0, 0.0]
    # radial equilibrium of the start with chi0 included: q^2/r^3 = lX r^3 + g^2 chi0^2 r
    req = brentq(lambda r: Q * Q / r ** 3 - LX * r ** 3 - G2 * c0 * c0 * r, 0.5 * r0, 1.5 * r0, xtol=1e-14)
    def Mof(s): return M0 * math.exp(kap * min(s, s_end))
    def rhs(s, y):
        x1, x2, p1, p2, c, pc = y
        r2 = x1 * x1 + x2 * x2; M = Mof(s)
        fX = LX * r2 + G2 * c * c
        fc = G2 * c * r2 + LP * c * (c * c - M * M) - 3 * epst * M * (c - M) ** 2
        return [p1, p2, -fX * x1, -fX * x2, pc, -fc]
    t0 = time.time()
    sol = solve_ivp(rhs, (0, s_end), y0, method='DOP853', rtol=rtol, atol=1e-10)
    y = sol.y[:, -1]
    mX = g * Mend; TX = 2 * math.pi / mX; Tw = 60 * TX
    ts = np.linspace(s_end, s_end + Tw, 6000)
    sol2 = solve_ivp(rhs, (s_end, s_end + Tw), y, method='DOP853', rtol=rtol, atol=1e-10, t_eval=ts)
    x1, x2, p1, p2, c, pc = sol2.y
    m = g * c
    q = x1 * p2 - x2 * p1
    nm = ((p1 ** 2 + p2 ** 2 + m * m * (x1 ** 2 + x2 ** 2)) / (2 * m) - q) / 2
    # phi quanta about the final minimum (free, m_phi^2 = 2 lP M^2 + g^2 r^2 ...): energy in oscillation / m_phi
    r2 = x1 ** 2 + x2 ** 2
    cbar = np.mean(c); mphi2 = G2 * np.mean(r2) + LP * (3 * cbar ** 2 - Mend ** 2) - 6 * epst * Mend * (cbar - Mend)
    Nphi = np.mean(0.5 * (pc ** 2 + mphi2 * (c - cbar) ** 2)) / math.sqrt(mphi2)
    return dict(delta=delta, eps_rate=eps_rate, epst=epst, rtol=rtol, r0=r0, M_c=Mc, chi0=c0, r_eq_start=req,
                delta_eff_start=(r0 * (1 + delta) - req) / req, Jr_over_q_start_est=1.2 * ((r0 * (1 + delta) - req) / req) ** 2,
                nminus_over_q=float(np.mean(nm) / Q), nminus_spread=float(np.std(nm) / Q), Nphi_over_q=float(Nphi / Q),
                q_drift=float(np.mean(q) / Q - 1), nfev=int(sol.nfev + sol2.nfev), sec=time.time() - t0)


if __name__ == '__main__':
    jobs = [(0.0, 1e-3), (0.0, 3e-4), (0.01, 1e-3), (1e-3, 1e-3)]
    if len(sys.argv) > 2:
        jobs = [(float(sys.argv[2]), float(sys.argv[3]))]
    out = []
    for d, e in jobs:
        r = run(d, e); out.append(r)
        print(json.dumps({k: (float('%.6g' % v) if isinstance(v, float) else v) for k, v in r.items()}), flush=True)
    json.dump(out, open(sys.argv[1], 'w'), indent=1)
