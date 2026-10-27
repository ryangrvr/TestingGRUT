"""en_rerun_skeptic -- own Floquet analysis of linear X fluctuations on the reduced homogeneous quartic orbit
(the orbit Codex's nonlinear stage starts from at tau 20). Different method from en_rerun's ir_tangent_scan
(which takes the max singular value of the tau 20->55 tangent map): here the orbit is autonomous in conformal
time s, so we find the radial period T_r, build the monodromy in the frame co-rotating with the orbit over one radial
period, and take the Floquet exponent mu(k) = ln|max multiplier|/T_r. Also the actual prepared-state growth of the
X variance from tau 20 to 55 per k for comparison.
usage: python3 -I floquet_reduced.py <MATCH_L12.npz> <out_json>"""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

LX = .001
m = np.load(sys.argv[1]); outp = sys.argv[2]
a_m = float(m['a_match']); C = float(m['C_match'])
x0 = np.array(m['light_mean_q'][:2], float); v0 = np.array(m['light_mean_p'][:2], float)


def s_of_tau(t):
    return 2*(t - 20)/(a_m + np.sqrt(a_m*a_m + 2*C*(t - 20)))


def orbit_rhs(s, y):
    x = y[:2]; return np.r_[y[2:], -LX*(x @ x)*x]


# radial period: successive minima of r^2 (event on d(r^2)/ds = 2 x.v crossing upward)
ev = lambda s, y: y[0]*y[2] + y[1]*y[3]
ev.direction = 1
sol = solve_ivp(orbit_rhs, (0, 6), np.r_[x0, v0], method='DOP853', rtol=1e-12, atol=1e-10, events=ev, dense_output=True)
tp = sol.t_events[0]
Tr = float(tp[2] - tp[1]); s_peri = float(tp[1])
xp = sol.sol(s_peri)
ang = lambda y: np.arctan2(y[1], y[0])
dtheta = float(np.unwrap([ang(sol.sol(t)) for t in np.linspace(s_peri, s_peri + Tr, 2001)])[-1] - ang(xp))
r_min2 = float(xp[:2] @ xp[:2]); r_max2 = float(max((sol.sol(t)[:2] @ sol.sol(t)[:2]) for t in np.linspace(s_peri, s_peri + Tr, 4001)))


def tangent_rhs(s, y, k):
    x = y[:2]; r2 = x @ x
    T = y[4:].reshape(4, 4)
    M = LX*r2*np.eye(2) + 2*LX*np.outer(x, x)
    dT = np.empty_like(T); dT[:2] = T[2:]; dT[2:] = -(k*k*T[:2] + M @ T[:2])
    return np.r_[y[2:4], -LX*r2*x, dT.ravel()]


def floquet(k):
    y0 = np.r_[xp, np.eye(4).ravel()]
    s = solve_ivp(tangent_rhs, (s_peri, s_peri + Tr), y0, args=(k,), method='DOP853', rtol=1e-11, atol=1e-12)
    T = s.y[4:, -1].reshape(4, 4)
    c, sn = np.cos(-dtheta), np.sin(-dtheta); Rm = np.array([[c, -sn], [sn, c]])
    Mono = np.block([[Rm, np.zeros((2, 2))], [np.zeros((2, 2)), Rm]]) @ T
    mult = np.linalg.eigvals(Mono)
    return float(np.log(np.max(abs(mult)))/Tr)


ks = np.unique(np.r_[np.linspace(.1, 8, 159), 2*np.pi/12*np.sqrt(np.arange(1, 13)), [1.683, 1.737, 1.8]])
mu = np.array([floquet(k) for k in ks])
s55 = s_of_tau(55.)
lat = {f'{2*np.pi/12*np.sqrt(n):.4f}': float(mu[np.argmin(abs(ks - 2*np.pi/12*np.sqrt(n)))]) for n in range(1, 13)}
res = dict(radial_period_s=Tr, orbit_rotation_per_radial_period=dtheta, r2_min=r_min2, r2_max=r_max2,
           s_tau55=float(s55), radial_periods_20_55=float(s55/Tr),
           k=ks.tolist(), mu=mu.tolist(),
           mu_max=float(mu.max()), k_at_mu_max=float(ks[mu.argmax()]),
           mu_at_1p737=float(mu[np.argmin(abs(ks - 1.737))]), mu_at_1p8=float(mu[np.argmin(abs(ks - 1.8))]),
           exp_mu_s55_at_1p737=float(np.exp(mu[np.argmin(abs(ks - 1.737))]*s55)),
           exp_mu_s55_max=float(np.exp(mu.max()*s55)),
           unstable_band=[float(ks[mu > 1e-3].min()), float(ks[mu > 1e-3].max())] if np.any(mu > 1e-3) else None,
           mu_L12_shells_below_1p8=lat)
json.dump(res, open(outp, 'w'), indent=1)
print({k: v for k, v in res.items() if k not in ('k', 'mu')})
for k, u in zip(ks, mu):
    if k < 3.2 and (abs(k - round(k*5)/5) < 1e-9 or abs(k - 1.737) < 1e-3 or abs(k-1.683) < 1e-3):
        print(f'k={k:.3f} mu={u:.4f}')
