"""Skeptic part A: re-derive the orbit numbers, spectral identity, late-orbit matching, ramps and sudden range with
scipy (DOP853 + quad), i.e. with methods different from pt_matching (symplectic order-8 C map + mpmath)."""
import sys, os, time, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sk_common import *
from scipy.optimize import brentq, minimize_scalar
from scipy.integrate import solve_ivp

t0 = time.time()
out = {}
y4 = tau20_state(); q = charge(y4)
info = orbit_info(y4, nsamp=8192)
n, w, P = info.pop('n'), info.pop('w'), info.pop('P')
pos = n >= 1; neg2 = n <= -2
info['J_minus_A_per_q'] = info['J_per_q'] - info['A_minus_per_q']
info['diff_pos'] = np.sum((n * w * P)[pos]) / q
info['diff_neg2'] = np.sum(((np.abs(n) - 1) * np.abs(w) * P)[neg2]) / q
info['J_per_orig_net'] = info['J_per_q'] * q / Q_ORIG_NET
info['A_per_orig_net'] = info['A_minus_per_q'] * q / Q_ORIG_NET
info['spectral_identity_rel'] = info['J_spec_per_q'] / info['J_per_q'] - 1
out['tau20'] = info
print({k: v for k, v in info.items()})

# u = r^2 obeys u'' = 4E - 3 lam u^2 (independent of q) -> Lame structure; check numerically on the orbit
sol = integrate_orbit(y4, info['T_r'], t_eval=np.linspace(0, info['T_r'], 2001))
u = sol.y[0] ** 2 + sol.y[1] ** 2
v2 = sol.y[2] ** 2 + sol.y[3] ** 2
upp = 2 * v2 - 2 * LX * u * u   # d^2 u/ds^2 = 2|psi'|^2 + 2 Re(psi* psi'')
out['u_equation_residual'] = float(np.max(np.abs(upp - (4 * info['E'] - 3 * LX * u * u))) / (4 * info['E']))
print('u-equation residual', out['u_equation_residual'])

# late orbit: pure-quartic orbit with same q, apocentre/pericentre state, A_-/q = 0.4760441973431541/1.0090661119866706
A_target = 0.4760441973431541 / 1.0090661119866706
rc = (q * q / LX) ** (1 / 6)
f = lambda r0: orbit_info(np.array([r0, 0.0, 0.0, q / r0]), nsamp=2048)['A_minus_per_q'] - A_target
r0 = brentq(f, 10.0, 0.999 * rc, xtol=1e-11)
late = orbit_info(np.array([r0, 0.0, 0.0, q / r0]), nsamp=8192)
for k in ('n', 'w', 'P'): late.pop(k)
late['r0'] = r0; late['A_target'] = A_target
out['late_A_match'] = late
print('late', {k: late[k] for k in ('J_per_q', 'A_minus_per_q', 'rmin', 'rmax', 'T_r', 'omega_r')})

# ---------------------------------------------------------------- ramps with scipy (non-symplectic)
MOUT = 30.0
def ramp_rhs(T, lamramp, rtype):
    def fr(u):
        if u <= 0: return 0.0
        if u >= 1: return 1.0
        if rtype == 0: return 0.5 * (1 - np.cos(np.pi * u))
        a = np.exp(-1 / u); b = np.exp(-1 / (1 - u)); return a / (a + b)
    def rhs(t, y):
        fv = fr(t / T); m2 = (MOUT * fv) ** 2; lam = LX * (1 - fv) if lamramp else LX
        c = m2 + lam * (y[0] ** 2 + y[1] ** 2)
        return [y[2], y[3], -c * y[0], -c * y[1]]
    return rhs

def ramp(st, T, lamramp=1, rtype=0, rtol=1e-12):
    qq = charge(st)
    if T == 0:
        psi = st[0] + 1j * st[1]; pi_ = st[2] + 1j * st[3]
        return abs(MOUT * psi + 1j * pi_) ** 2 / (4 * MOUT) / qq      # n_- = m|B|^2, B e^{-ims} = (psi + i pi/m)/2
    sol = solve_ivp(ramp_rhs(T, lamramp, rtype), (0, T), st, method='DOP853', rtol=rtol, atol=1e-10)
    y = sol.y[:, -1]
    psi = y[0] + 1j * y[1]; pi_ = y[2] + 1j * y[3]
    if lamramp:
        return abs(MOUT * psi + 1j * pi_) ** 2 / (4 * MOUT) / qq
    Ef = energy4(y, LX, MOUT ** 2)
    return radial_quad(Ef, charge(y), LX, MOUT ** 2)['J_r'] / qq

Jq = info['J_per_q']; Tr = info['T_r']
ph8 = np.arange(8) / 8
st8 = [phase_state(y4, Tr, p) for p in ph8]
rows = {}
for T in [0.0, 0.2, 2.0, 20.0, 80.0]:
    a = np.array([ramp(s, T) for s in st8])
    rows[f'{T:g}'] = dict(values=a, min=a.min(), max=a.max(), mean=a.mean(), spread=np.ptp(a), max_dev=np.max(np.abs(a - Jq)))
    print('T', T, a.min(), a.max(), a.mean(), np.ptp(a), time.time() - t0)
out['ramps_8phase'] = rows
mo = np.array([ramp(s, 80.0, lamramp=0) for s in st8[::2]]) - Jq
out['mass_only_T80_dev'] = mo
print('mass-only devs', mo)
# sudden: continuous phase via dense output and optimisation
sol = integrate_orbit(y4, Tr)
def nsud(ph):
    y = sol.sol(ph * Tr); psi = y[0] + 1j * y[1]; pi_ = y[2] + 1j * y[3]
    return abs(MOUT * psi + 1j * pi_) ** 2 / (4 * MOUT) / q
grid = np.linspace(0, 1, 20001)[:-1]
vals = np.array([nsud(p) for p in grid])
i0, i1 = vals.argmin(), vals.argmax()
mn = minimize_scalar(nsud, bracket=(grid[i0] - 1e-4, grid[i0], grid[i0] + 1e-4), tol=1e-14)
mx = minimize_scalar(lambda p: -nsud(p), bracket=(grid[i1] - 1e-4, grid[i1], grid[i1] + 1e-4), tol=1e-14)
# closed form: n_- = |m psi + i pi|^2/(4 m q); its orbit average equals (<|pi|^2>/(2m) + m<r^2>/2 - q)/(2q)
out['sudden_continuous'] = dict(min=mn.fun, argmin=mn.x % 1, max=-mx.fun, argmax=mx.x % 1, mean=vals.mean(),
                                mean_analytic=None)
v2m = np.mean([abs(sol.sol(p * Tr)[2] + 1j * sol.sol(p * Tr)[3]) ** 2 for p in grid[::4]])
r2m = np.mean([abs(sol.sol(p * Tr)[0] + 1j * sol.sol(p * Tr)[1]) ** 2 for p in grid[::4]])
out['sudden_continuous']['mean_analytic'] = (v2m / (2 * MOUT) + MOUT * r2m / 2 - q) / (2 * q)
print('sudden', out['sudden_continuous'])
# T dependence (cosine and C-infinity), 16 phases
ph16 = np.arange(16) / 16; st16 = [phase_state(y4, Tr, p) for p in ph16]
td = []
for rt, Ts in [(0, [10, 20, 40, 80, 160]), (1, [10, 20, 40])]:
    for T in Ts:
        a = np.array([ramp(s, float(T), rtype=rt, rtol=1e-13) for s in st16])
        td.append(dict(ramp=['cosine', 'Cinf'][rt], T=T, max_abs_dev=np.max(np.abs(a - Jq)), mean_dev=np.mean(a - Jq)))
        print(td[-1], time.time() - t0)
cs = [x for x in td if x['ramp'] == 'cosine' and x['T'] >= 20]
out['T_dependence'] = td
out['cosine_slope'] = np.polyfit(np.log([x['T'] for x in cs]), np.log([x['max_abs_dev'] for x in cs]), 1)[0]
print('slope', out['cosine_slope'])
save_json('skA_orbit_ramps.json', out)
print('elapsed', time.time() - t0)
