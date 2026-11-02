#!/usr/bin/env python3
"""dq_heat track (Claude, own code; issue #13 decisive pass, 2 November 2026): COPY of tracks/fx_holder/code/fx_common.py.
Changes: scratch/binary paths point to this track's rptri_dq.c (the fixed copy of rptri.c), and run() also returns the
dense "D" state lines (y[0..18], H, R, Hdot at every sample).  Everything else is unchanged.
Original docstring:
fx_holder track (Claude, own code; issue #13 working-theory repair pass): shared model, preparation solvers and the
driver of rptri.c (our engine with mu_eff = mu + LOCK R + G phihat).

Model [Q CANDIDATE_M1 pins, as in m1_holder/code/holder_point.py]: H_i = 1.3428427349130815e11 GeV, M_Y = 30 H_i,
h = 40 H_i, mu = 7.159456 H_i, kappa = 3.4587e-3, lambda_X = 1e-3, lambda_Y = lambda_Sigma = 1 (or 0.5 for repair 1b),
xi_X = -5/6 (xi R = -10 H^2 in de Sitter); quadratic inflaton m_I = b H_i, Gamma_I = gam H_i (mode 0).

Copied (then modified) from our own earlier tracks, stated: the conditional-minimum preparation of
rp_reheating/code/rp.py::prepare (parameter dict instead of globals), the stationary-point and exact-Hessian method of
m1_holder/code/holder_point.py (sympy Hessian now with lambda_Y, lambda_Sigma symbolic), the yield and decay-ledger
conventions of rp_reheating/code/rp.py and analysis.py.
"""
import json, os, subprocess, time
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.optimize import minimize, root

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
SP = Path(os.environ.get('DQ_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/dq_heat'))
BIN = SP / 'build' / 'rptri_dq'
SRC = HERE / 'rptri_dq.c'

H_GEV = 134284273491.30815                 # M1 H_i [Q d4747cb]
MPL = 2.435e18
ETA2 = (H_GEV / MPL) ** 2
MU = 961402368973.181 / H_GEV               # 7.159456159...
BASE = dict(lx=1e-3, ly=1.0, ls=1.0, kap=0.0034586847211042523, M=4028528204739.244 / H_GEV, h=5371370939652.326 / H_GEV,
            mu=MU, xi=-5.0 / 6.0, cH=10.0)
LOCK_1A = (10j - MU) / 12.0                 # Option 1A: mu + 12 LOCK = 10 i  [Q f8186fa; reconstructed in m1_holder]
U_PREP = 10j                                # required mu_eff at the preparation slice
YQ_TARGET = 5.07941423e-9                   # [Q M1 calibration at E/Q = 0.9374 m_X]
GSTAR = 106.75
COLS = ['tau', 'a', 'H', 'R', 'NX', 'NY', 'NS', 'Gauss', 'EX', 'eps', 'NYH', 'NYL', 'Nsig', 'absX', 'absY', 'S_over_vev',
        'rho_r', 'rho_I', 'radfrac', 'Hcheck', 'Wcheck', 'w_I', 'm2curv', 'KV', 'Haux', 'rho', 'ure', 'uim', 'phihat', 'mL2']


def params(lam=1.0, **kw):
    p = dict(BASE); p['ly'] = lam; p['ls'] = lam; p.update(kw)
    return p


# ------------------------------------------------------------------------------------------------ potential
def potential(z, u, p):
    """hidden V without the curvature mass (as in the engine) and F = dV/dz*."""
    x, y, s = z
    x2, y2, s2 = abs(x) ** 2, abs(y) ** 2, abs(s) ** 2
    V = p['lx'] * x2 ** 2 + p['M'] ** 2 * y2 + p['ly'] * y2 ** 2 + p['ls'] * (s2 - p['h'] ** 2 / 2) ** 2 \
        + 2 * np.real(p['kap'] * np.conj(y) * x ** 3 + u * np.conj(s) * y * y)
    F = np.array([2 * p['lx'] * x2 * x + 3 * p['kap'] * y * np.conj(x) ** 2,
                  p['M'] ** 2 * y + 2 * p['ly'] * y2 * y + p['kap'] * x ** 3 + 2 * np.conj(u) * s * np.conj(y),
                  2 * p['ls'] * (s2 - p['h'] ** 2 / 2) * s + u * y * y])
    return float(V), F


def inflaton_I0(z, u, p, b):
    """I0/M_r from the Friedmann constraint at a = 1, H = H_i, zero velocities (mode 0 oscillator at rest)."""
    V, _ = potential(z, u, p)
    f = 1 - 2 * p['xi'] * ETA2 * abs(z[0]) ** 2
    need = 3 * f - ETA2 * V
    return float(np.sqrt(2 * need) / b)


def trilinear_G(z, p, b, u_target=U_PREP):
    """G = g_I e^{i alpha_I} such that mu + G phihat(0) = u_target, phihat(0) = I0/sqrt(eta2) (I0 self-consistent:
    V at the preparation already contains u_target)."""
    I0 = inflaton_I0(z, u_target, p, b)
    ph0 = I0 / np.sqrt(ETA2)
    return (u_target - p['mu']) / ph0, I0, ph0


# ------------------------------------------------------------------------------------------------ exact Hessian (sympy)
_r = sp.symbols('r1:7', real=True)
_ur, _ui, _R, _lx, _ly, _ls, _kap, _M, _h, _xi = sp.symbols('u_r u_i R lx ly ls kap M h xi', real=True)
_x = (_r[0] + sp.I * _r[1]) / sp.sqrt(2); _y = (_r[2] + sp.I * _r[3]) / sp.sqrt(2); _s = (_r[4] + sp.I * _r[5]) / sp.sqrt(2)
_u = _ur + sp.I * _ui
_cub = _kap * sp.conjugate(_y) * _x ** 3 + _u * sp.conjugate(_s) * _y ** 2
_V = sp.expand(_lx * (_x * sp.conjugate(_x)) ** 2 + _ly * (_y * sp.conjugate(_y)) ** 2 + _M ** 2 * _y * sp.conjugate(_y)
               + _ls * (_s * sp.conjugate(_s) - _h ** 2 / 2) ** 2 + _cub + sp.conjugate(_cub) + _xi * _R * _x * sp.conjugate(_x))
_args = (_r, _ur, _ui, _R, _lx, _ly, _ls, _kap, _M, _h, _xi)
_hess_f = sp.lambdify(_args, sp.hessian(_V, _r), 'numpy')
_grad_f = sp.lambdify(_args, [sp.diff(_V, v) for v in _r], 'numpy')
_V_f = sp.lambdify(_args, _V, 'numpy')
# derivative of V with respect to u (for the inflaton mixing): dV/du_r, dV/du_i as functions of r
_dVdu_f = sp.lambdify(_args, [[sp.diff(sp.diff(_V, _ur), v) for v in _r], [sp.diff(sp.diff(_V, _ui), v) for v in _r]], 'numpy')


def _pa(p):
    return (p['lx'], p['ly'], p['ls'], p['kap'], p['M'], p['h'], p['xi'])


def z_of_r(r):
    return np.array([(r[0] + 1j * r[1]), (r[2] + 1j * r[3]), (r[4] + 1j * r[5])]) / np.sqrt(2)


def r_of_z(z):
    return np.sqrt(2) * np.array([z[0].real, z[0].imag, z[1].real, z[1].imag, z[2].real, z[2].imag])


def Vr(r, u, p, R=12.0):
    return float(np.real(_V_f(r, u.real, u.imag, R, *_pa(p))))


def gradr(r, u, p, R=12.0):
    return np.array(_grad_f(r, u.real, u.imag, R, *_pa(p)), float)


def hessr(r, u, p, R=12.0):
    return np.array(_hess_f(r, u.real, u.imag, R, *_pa(p)), float)


def mixing_u(r, u, p, R=12.0):
    """d^2 V / (du_r dr), d^2 V / (du_i dr): 2 x 6."""
    return np.array(_dVdu_f(r, u.real, u.imag, R, *_pa(p)), float)


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def stationary(u, p, R=12.0, guess=None):
    """all-direction stationary point of V + xi R |X|^2 at fixed u (own guess as in holder_point.py), gauge-fixed X > 0."""
    if guess is None:
        x0 = np.sqrt(max(-p['xi'] * R, 1.0) / (2 * p['lx']))
        y0 = -p['kap'] * x0 ** 3 / p['M'] ** 2
        s0 = -(u / abs(u)) * p['h'] / np.sqrt(2)
        guess = np.array([x0, y0, s0])
    r0 = r_of_z(np.asarray(guess, complex))
    m = minimize(lambda r: Vr(r, u, p, R), r0, jac=lambda r: gradr(r, u, p, R), method='BFGS', options=dict(gtol=1e-9, maxiter=4000))
    sol = root(lambda r: gradr(r, u, p, R), m.x, method='lm', options=dict(xtol=1e-15, ftol=1e-15, maxiter=20000))
    z = z_of_r(sol.x); th = np.angle(z[0]); z = z * np.exp(-1j * th * np.array([1, 3, 6]))
    r = r_of_z(z)
    H = hessr(r, u, p, R); ev, U = np.linalg.eigh(H)
    return dict(z=z, r=r, grad=float(np.linalg.norm(gradr(r, u, p, R))), V=Vr(r, u, p, R), eig=ev, U=U,
                amps=np.abs(z), delta=float(wrap(np.angle(z[2]) - 6 * np.angle(z[0]))), delta_y=float(wrap(np.angle(z[1]) - 3 * np.angle(z[0]))))


# ------------------------------------------------------------------------------------------------ held-phase preparation
def ricci0(z, u, p, b, mode=0):
    V, F = potential(z, u, p)
    x2 = abs(z[0]) ** 2
    f = 1 - 2 * p['xi'] * ETA2 * x2
    need = 3 * f - ETA2 * V
    trI = 4 * need
    return (trI + ETA2 * (4 * V - 12 * p['xi'] * np.real(np.conj(z[0]) * F[0]))) / (f + 12 * p['xi'] ** 2 * ETA2 * x2)


def prepare_held(p, b, delta=0.0, u=None):
    """rp.prepare (own rp_reheating code) with a parameter dict: X = x > 0, Sigma = -i e^{i delta} s, Y free, curvature
    mass from the actual R(0); u = mu by default (the unheld/original preparation)."""
    u = p['mu'] if u is None else u
    direc = -1j * np.exp(1j * delta)

    def state(q):
        return np.array([q[0], q[1] + 1j * q[2], direc * q[3]], complex)

    def eq(q):
        z = state(q); V, F = potential(z, u, p); m = p['xi'] * ricci0(z, u, p, b)
        return [np.real(F[0] + m * z[0]), F[1].real, F[1].imag, np.real(F[2] * np.conj(direc))]
    x0 = np.sqrt(p['cH'] / (2 * p['lx']))
    q0 = [x0, -p['kap'] * x0 ** 3 / p['M'] ** 2, 0.0, p['h'] / np.sqrt(2)]
    sol = root(eq, q0, tol=1e-14)
    z = state(sol.x)
    return z, dict(q=sol.x.tolist(), residual=float(np.linalg.norm(eq(sol.x))), R0=float(ricci0(z, u, p, b)), success=bool(sol.success))


# ------------------------------------------------------------------------------------------------ engine driver
def build():
    BIN.parent.mkdir(parents=True, exist_ok=True)
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        subprocess.run(['gcc', '-O3', '-std=c11', '-Wall', '-o', str(BIN), str(SRC), '-lm'], check=True)


def run(tag, p, b, gam, z0=None, v0=None, state=None, lock=0j, G=0j, end=1000.0, rtol=1e-10, atol=1e-13, nsamp=5000, binw=20.0, save=True):
    build()
    g = lambda v: '%.17g' % float(v)
    args = [str(BIN), '0', g(b), g(gam), g(p['xi']), g(ETA2), g(p['cH']), g(p['lx']), g(p['ly']), g(p['ls']), g(p['kap']), g(p['M']),
            g(p['h']), g(p['mu']), g(lock.real), g(lock.imag), g(G.real), g(G.imag), g(end), g(rtol), g(atol), str(int(nsamp)), g(binw)]
    if state is None:
        v0 = np.zeros(3, complex) if v0 is None else np.asarray(v0, complex)
        args += ['F'] + [g(v) for c in z0 for v in (c.real, c.imag)] + [g(v) for c in v0 for v in (c.real, c.imag)]
    else:
        args += ['S'] + [g(v) for v in state]
    t0 = time.time()
    out = subprocess.run(args, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(tag + ': ' + out.stderr)
    S, Bn, W, Y, END, I0, Dl = [], [], {}, None, None, None, []
    for line in out.stdout.split('\n'):
        if line.startswith('S '):
            S.append([float(v) for v in line.split()[1:]])
        elif line.startswith('D '):
            Dl.append([float(v) for v in line.split()[1:]])
        elif line.startswith('B '):
            v = line.split(); Bn.append([float(x) for x in v[2:]])
        elif line.startswith('W') or line.startswith('A3040'):
            v = line.split(); W[v[0]] = dict(span=float(v[1]), **dict(zip(COLS, map(float, v[2:]))))
        elif line.startswith('Y '):
            Y = [float(v) for v in line.split()[1:]]
        elif line.startswith('I0 '):
            v = [float(x) for x in line.split()[1:]]; I0 = dict(I0=v[0], phihat0=v[1], u0=[v[2], v[3]], R0=v[4])
        elif line.startswith('END'):
            v = line.split(); END = dict(steps=int(v[1]), rejects=int(v[2]), max_gauss=float(v[3]), max_H_consistency=float(v[4]),
                                        max_hidden_work=float(v[5]), min_mL2=float(v[6]))
    S = np.array(S); Bn = np.array(Bn); Dl = np.array(Dl)
    if save:
        (SP / 'runs').mkdir(parents=True, exist_ok=True)
        np.save(SP / 'runs' / (tag + '.npy'), S); np.save(SP / 'runs' / (tag + '_bins.npy'), Bn); np.save(SP / 'runs' / (tag + '_D.npy'), Dl)
    return dict(tag=tag, S=S, bins=Bn, D=Dl, W=W, Y=Y, END=END, I0=I0, seconds=time.time() - t0, lam=p['ly'], b=b, gam=gam, end=end,
                lock=[lock.real, lock.imag], G=[G.real, G.imag], rtol=rtol)


def bin_col(rr, name):
    """bins array columns: t0, t1, then the NQ means in COLS order."""
    return rr['bins'][:, 2 + COLS.index(name)]


def window_stats(rr, name='NX', lo=400.0, hi=None, width=80.0):
    """mean over [lo, hi] and the spread of running `width`-wide windows built from the bins."""
    B = rr['bins']; hi = rr['end'] if hi is None else hi
    t0, t1 = B[:, 0], B[:, 1]; vals = bin_col(rr, name); dt = t1 - t0
    sel = (t0 >= lo - 1e-9) & (t1 <= hi + 1e-9)
    mean = float(np.sum(vals[sel] * dt[sel]) / np.sum(dt[sel]))
    bw = float(np.median(dt)); n = int(round(width / bw))
    idx = np.where(sel)[0]; run = []
    for k in range(len(idx) - n + 1):
        j = idx[k:k + n]; run.append((float(t0[j[0]]), float(np.sum(vals[j] * dt[j]) / np.sum(dt[j]))))
    rv = np.array([v for _, v in run])
    return dict(mean=mean, lo=lo, hi=hi, width=width, run_min=float(rv.min()), run_max=float(rv.max()), run_std=float(rv.std()),
                running=[[a, v] for a, v in run])


def window_mean(rr, name, lo, hi):
    B = rr['bins']; t0, t1 = B[:, 0], B[:, 1]; vals = bin_col(rr, name); dt = t1 - t0
    sel = (t0 >= lo - 1e-9) & (t1 <= hi + 1e-9)
    return float(np.sum(vals[sel] * dt[sel]) / np.sum(dt[sel]))


def yield_over_target(rr, NX, H=H_GEV):
    """Y_X = |N_X| H^3 / (a^3 s) at the end of the run (entropy a^3 s conserved after the inflaton decays), relative to
    the M1 target Y_Q = 5.07941423e-9; s = 2 pi^2 g* T^3/45, T from rho_r (rp.py convention)."""
    last = dict(zip(COLS, rr['S'][-1]))
    a, rr_ = last['a'], last['rho_r']
    rho_rad = rr_ * MPL ** 2 * H ** 2
    T = (30 * rho_rad / (np.pi ** 2 * GSTAR)) ** 0.25
    s = 2 * np.pi ** 2 / 45 * GSTAR * T ** 3
    Y = abs(NX) * H ** 3 / (a ** 3 * s)
    return float(Y), float(Y / YQ_TARGET), float(rr_ * a ** 4)


# decay ledger (rp_reheating/analysis.py, our hw_verify/hw_gates ledger with the (Rr/3)^(1/4) correction) [Q own, E]
def ledger(NYH, NYL, Nsig, EX, Rr, p, e=1e-7, H=H_GEV):
    M, h, mu, kap = p['M'] * H, p['h'] * H, p['mu'] * H, p['kap']
    mL = np.sqrt(M ** 2 - np.sqrt(2) * mu * h); mA = 6 * e * h
    GL = 3 * kap ** 2 * mL / (256 * np.pi ** 3); GA = e ** 2 * mA / (48 * np.pi)
    DEP_L = (mL / H) * np.sqrt(np.pi * H / (2 * GL)); DEP_A = (mA / H) * np.sqrt(np.pi * H / (2 * GA))
    corr = (Rr / 3) ** 0.25
    RY = (NYH + NYL) * DEP_L * corr / EX; RC = NYH * DEP_A * corr / EX; RS = 2 * Nsig * DEP_A * corr / EX
    return dict(R_Y=float(RY), R_cascade=float(RC), R_Sigma=float(RS), R_central=float(RY + RC + RS), reheat_factor=float(corr),
                deposit_per_light_Y=float(DEP_L), mL_over_H=float(mL / H), e=e)
