"""Common helpers for track pt_matching (independent code; no Codex module is imported).

Reads only numeric inputs: the tau=20 light mean field (MATCH_L12.npz, light_mean_q/p) from the reviewed
early packet, and Codex's published numbers for comparison.
"""
import ctypes, os, subprocess, json, math
import numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad'
F55 = SP + '/f55/research/reality_framework'
MATCH = F55 + '/grut_early_nonlinear_20261027/MATCH_L12.npz'
CODEX = F55 + '/grut_particle_transport_20261027'
RESULTS = os.path.join(os.path.dirname(HERE), 'results')
SCRATCH = SP + '/review/pt_matching'

LX, G2, LP = 1e-3, 1e-2, 3e-2   # lambda_X, g^2, lambda_phi (Codex/brief values)
Q_ORIG_NET = 8113.477086664833   # original full-state net charge, quoted from COHERENT_MATCHING.json [Q]


def tau20_state():
    d = np.load(MATCH)
    q = d['light_mean_q'][:2].astype(float); p = d['light_mean_p'][:2].astype(float)
    return np.array([q[0], q[1], p[0], p[1]])


def charge(y4):
    return y4[0] * y4[3] - y4[1] * y4[2]


# ---------------------------------------------------------------- C library
class Par(ctypes.Structure):
    _fields_ = [('mode', ctypes.c_int), ('lamramp', ctypes.c_int), ('ramp_type', ctypes.c_int), ('growth', ctypes.c_int),
                ('lamX0', ctypes.c_double), ('g2', ctypes.c_double), ('lamP', ctypes.c_double), ('T', ctypes.c_double),
                ('mout', ctypes.c_double), ('M0', ctypes.c_double), ('kappa', ctypes.c_double)]


def _lib():
    so = os.path.join(SCRATCH, 'libhom.so')
    src = os.path.join(HERE, 'hom.c')
    if not os.path.exists(so) or os.path.getmtime(so) < os.path.getmtime(src):
        os.makedirs(SCRATCH, exist_ok=True)
        subprocess.check_call(['gcc', '-O3', '-march=native', '-fPIC', '-shared', '-o', so, src, '-lm'])
    lib = ctypes.CDLL(so)
    dp = ctypes.POINTER(ctypes.c_double)
    lib.run.restype = ctypes.c_long
    lib.run.argtypes = [ctypes.POINTER(Par), dp, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
                        ctypes.c_int, ctypes.c_double, dp, ctypes.c_long]
    lib.run_fixed.restype = None
    lib.run_fixed.argtypes = [ctypes.POINTER(Par), dp, ctypes.c_double, ctypes.c_long, ctypes.c_int, dp]
    lib.energy.restype = ctypes.c_double
    lib.energy.argtypes = [ctypes.POINTER(Par), dp]
    lib.set_waves.restype = None
    lib.set_waves.argtypes = [ctypes.c_int, dp, ctypes.c_double, ctypes.c_int]
    lib.Mvalue.restype = ctypes.c_double
    lib.Mvalue.argtypes = [ctypes.POINTER(Par), ctypes.c_double]
    return lib


LIB = _lib()
DP = ctypes.POINTER(ctypes.c_double)


def make_par(mode=0, lamramp=1, ramp_type=0, growth=0, lamX0=LX, g2=G2, lamP=LP, T=1.0, mout=30.0, M0=0.0, kappa=0.0):
    return Par(mode, lamramp, ramp_type, growth, lamX0, g2, lamP, T, mout, M0, kappa)


def y8(x1, x2, phi, p1, p2, pphi, s=0.0, ps=0.0):
    return np.ascontiguousarray([x1, x2, phi, p1, p2, pphi, s, ps], dtype=float)


def c_run(par, y, s_end, h0, wb0=1.0, wbM=0.0, order=8, dsample=0.0, nmax=0):
    out = np.zeros((max(nmax, 1), 8)) if nmax else None
    n = LIB.run(ctypes.byref(par), y.ctypes.data_as(DP), s_end, h0, wb0, wbM, order, dsample,
                out.ctypes.data_as(DP) if nmax else None, nmax)
    return (out[:n] if nmax else n)


def c_fixed(par, y, h, n, order=8, store=True):
    out = np.zeros((n, 8)) if store else None
    LIB.run_fixed(ctypes.byref(par), y.ctypes.data_as(DP), h, n, order, out.ctypes.data_as(DP) if store else None)
    return out


def c_set_waves(k2, mw2, wslave=0):
    k2 = np.ascontiguousarray(k2, dtype=float)
    LIB.set_waves(len(k2), k2.ctypes.data_as(DP), float(mw2), int(wslave))


def c_energy(par, y):
    return LIB.energy(ctypes.byref(par), y.ctypes.data_as(DP))


# ---------------------------------------------------------------- central-orbit quadratures (mpmath)
def turning_u(E, q, lam, m2=0.0):
    """Roots of p_r^2 r^2 = 2E r^2 - q^2 - m2 r^4 - lam r^6/2 as a cubic in u=r^2 (u0<0<u1<u2)."""
    roots = mp.polyroots([mp.mpf(lam) / 2, mp.mpf(m2), -2 * mp.mpf(E), mp.mpf(q) ** 2], maxsteps=200, extraprec=200)
    rr = sorted([mp.re(x) for x in roots])
    return rr


def orbit_quadrature(E, q, lam, m2=0.0, dps=30):
    """J_r = (1/pi) int p_r dr, radial period, angle advance, by tanh-sinh quadrature in a smooth variable."""
    with mp.workdps(dps):
        E, q, lam, m2 = mp.mpf(E), mp.mpf(q), mp.mpf(lam), mp.mpf(m2)
        if lam > 0:
            u0, u1, u2 = turning_u(E, q, lam, m2)
            G = lambda u: (lam / 2) * (u - u0) / u      # p_r^2 = G(u) (u-u1)(u2-u)
        else:  # pure harmonic: m2 u^2 - 2E u + q^2 = 0
            disc = mp.sqrt(E * E - m2 * q * q)
            u1, u2 = (E - disc) / m2, (E + disc) / m2
            G = lambda u: m2 / u
        # u = u1 + (u2-u1) sin^2 t, t in [0, pi/2]; r = sqrt(u), dr = du/(2r); endpoint factors cancel analytically
        D = u2 - u1
        def U(t): return u1 + D * mp.sin(t) ** 2
        def fJ(t):
            u = U(t); return mp.sqrt(G(u)) * D * D * (mp.sin(t) * mp.cos(t)) ** 2 / mp.sqrt(u)
        def fT(t):
            u = U(t); return 1 / (mp.sqrt(u) * mp.sqrt(G(u)))
        def fA(t):
            u = U(t); return q / u / (mp.sqrt(u) * mp.sqrt(G(u)))
        J = mp.quad(fJ, [0, mp.pi / 4, mp.pi / 2]) / mp.pi
        Tr = 2 * mp.quad(fT, [0, mp.pi / 4, mp.pi / 2])
        Dth = 2 * mp.quad(fA, [0, mp.pi / 4, mp.pi / 2])
        return dict(J_r=float(J), T_r=float(Tr), dtheta=float(Dth), rmin=float(mp.sqrt(u1)), rmax=float(mp.sqrt(u2)))


def save_json(name, obj):
    os.makedirs(RESULTS, exist_ok=True)
    def conv(x):
        if isinstance(x, (np.floating,)): return float(x)
        if isinstance(x, (np.integer,)): return int(x)
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, mp.mpf): return float(x)
        raise TypeError(type(x))
    with open(os.path.join(RESULTS, name), 'w') as f:
        json.dump(obj, f, indent=1, default=conv)
