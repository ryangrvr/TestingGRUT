"""Independent check of Codex 2335653: cold elliptic orbit (phi = 0) Floquet exponents, own code.
Units m_X = g v = 1: g = 0.1, v = 10, lambda_X = 1e-3, lambda_phi = 0.03.  rho_c^2 = lambda_phi v^2/g^2 = 300.
Invariants matched to the release: Q = n_Q(circle)/6.86259640973289, E_X = circle X-energy (excluding V0)."""
import numpy as np, json, sys
from scipy.integrate import solve_ivp
g, v, lX, lP = 0.1, 10.0, 1e-3, 0.03
rc = np.sqrt(lP * v * v / g**2)
om = np.sqrt(lX * rc**2); nQc = om * rc**2
Q = nQc / 6.86259640973289
EX = 0.5 * om**2 * rc**2 + lX * rc**4 / 4
p2 = Q / rc; p1 = -np.sqrt(2 * (EX - lX * rc**4 / 4) - p2**2)

def bg(t, y):
    q = y[:2]; r2 = q @ q
    return np.r_[y[2:], -lX * r2 * q]

def crossing(t, y):  # radius through rho_c, inward
    return np.hypot(y[0], y[1]) - rc
crossing.direction = -1
y0 = np.array([rc, 0.0, p1, p2])
s = solve_ivp(bg, [0, 200], y0, events=crossing, rtol=1e-13, atol=1e-13, dense_output=True)
ts = [t for t in s.t_events[0] if t > 1e-6]
Tr = ts[0]
yT = s.sol(Tr); dth = np.arctan2(yT[1], yT[0])

def var(t, Y, k, sector):
    y = Y[:4]; q = y[:2]; r2 = q @ q
    out = [bg(t, y)]
    if sector == 'X':
        M = lX * (r2 * np.eye(2) + 2 * np.outer(q, q)) + k * k * np.eye(2)
        P = Y[4:].reshape(4, 4)             # columns: perturbation states (dq, dp)
        dP = np.vstack([P[2:], -M @ P[:2]])
        out.append(dP.ravel())
    else:
        m2 = g * g * r2 - lP * v * v + k * k
        P = Y[4:].reshape(2, 2)
        dP = np.vstack([P[1:], -m2 * P[:1]])
        out.append(dP.ravel())
    return np.concatenate(out)

def floquet(k, sector):
    n = 4 if sector == 'X' else 2
    Y0 = np.r_[y0, np.eye(n).ravel()]
    r = solve_ivp(var, [0, Tr], Y0, args=(k, sector), rtol=1e-12, atol=1e-13)
    Mo = r.y[4:, -1].reshape(n, n)
    if sector == 'X':  # relative monodromy: undo the background rotation on coordinates and momenta
        c, sn = np.cos(-dth), np.sin(-dth); Rm = np.array([[c, -sn], [sn, c]])
        Mo = np.kron(np.eye(2), Rm) @ Mo
    ev = np.linalg.eigvals(Mo)
    return float(np.max(np.log(np.abs(ev))) / Tr), float(np.max(np.abs(ev)))

res = dict(rc=rc, omega=om, Q=Q, EX=EX, p1=p1, p2=p2, T_r=Tr, dtheta=dth)
res['X_k0.34148'] = floquet(0.34148460188515356, 'X')
res['phi_k0'] = floquet(0.0, 'phi')
res['phi_k0.34148'] = floquet(0.34148460188515356, 'phi')
ks = np.linspace(0.05, 0.8, 31)
scan = [(float(k), floquet(k, 'X')[0]) for k in ks]
kb, mb = max(scan, key=lambda t: t[1])
res['X_scan_peak'] = dict(k=kb, mu=mb)
res['codex'] = dict(k=0.34148460188515356, mu_X=0.1031462958281267, mu_phi_at_k0=0.8590822899833551, mu_phi_at_kX=0.8485739417008076)
json.dump(res, open(sys.argv[1], 'w'), indent=1)
for k_, v_ in res.items(): print(k_, v_)
