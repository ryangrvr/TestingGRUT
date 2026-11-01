"""Skeptic of track pt_m1tilt: common helpers (own code; nothing imported from the track's m1_*.py or m1_dop.c).
Orbits: tau-20 light mean field (MATCH_L12.npz light_mean_q/p) and the pure-quartic 'late' stand-in (r0 = 31.4946,
same q), as in pt_matching.  J_r, T_r, omega_r from a time-domain integration of the pure quartic orbit with scipy
DOP853 (J_r = (1/2pi) int_0^{T_r} p_r^2 ds), not from the track's quadrature.
Integrator: sk_sym.c (symplectic Kahan-Li order-8 composition in extended phase space, s-dependent step).
[C] homogeneous classical model; not a cosmological result."""
import ctypes, os, subprocess, json, time
import numpy as np
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad'
MATCH = SP + '/f55/research/reality_framework/grut_early_nonlinear_20261027/MATCH_L12.npz'
SCR = SP + '/review/pt_m1tilt_skeptic'
RUNS = SCR + '/runs'
RESULTS = os.path.join(os.path.dirname(HERE), 'results')
TRACK_RUNS = SP + '/review/pt_m1tilt/runs'

LX, G2, LP, EPST = 1e-3, 1e-2, 3e-2, 0.002
GG = np.sqrt(G2)
LATE_R0 = 31.494634874218747


class SkPar(ctypes.Structure):
    _fields_ = [('lamX', ctypes.c_double), ('g2', ctypes.c_double), ('lamP', ctypes.c_double), ('epst', ctypes.c_double),
                ('M0', ctypes.c_double), ('kappa', ctypes.c_double), ('om_fix', ctypes.c_double), ('fm', ctypes.c_double),
                ('growth', ctypes.c_int), ('method', ctypes.c_int)]


def _lib():
    os.makedirs(SCR, exist_ok=True)
    so = os.path.join(SCR, 'libsksym.so'); src = os.path.join(HERE, 'sk_sym.c')
    if not os.path.exists(so) or os.path.getmtime(so) < os.path.getmtime(src):
        subprocess.check_call(['gcc', '-O3', '-fPIC', '-shared', '-o', so, src, '-lm'])
    lib = ctypes.CDLL(so)
    dp = ctypes.POINTER(ctypes.c_double)
    lib.sk_run.restype = ctypes.c_long
    lib.sk_run.argtypes = [ctypes.POINTER(SkPar), dp, ctypes.c_double, ctypes.c_int]
    lib.sk_sample.restype = ctypes.c_long
    lib.sk_sample.argtypes = [ctypes.POINTER(SkPar), dp, ctypes.c_double, ctypes.c_int, ctypes.c_int, dp]
    lib.sk_E.restype = ctypes.c_double; lib.sk_E.argtypes = [ctypes.POINTER(SkPar), dp]
    lib.sk_V.restype = ctypes.c_double; lib.sk_V.argtypes = [ctypes.POINTER(SkPar), dp, ctypes.c_double]
    lib.sk_grad_py.restype = None; lib.sk_grad_py.argtypes = [ctypes.POINTER(SkPar), dp, ctypes.c_double, dp, dp]
    return lib


LIB = _lib()
DP = ctypes.POINTER(ctypes.c_double)


def arr(a):
    return np.ascontiguousarray(a, dtype=float)


# ------------------------------------------------------------------ orbits (own)
def charge(y4):
    return y4[0] * y4[3] - y4[1] * y4[2]


def quartic_rhs(t, y):
    c = LX * (y[0] ** 2 + y[1] ** 2)
    return [y[2], y[3], -c * y[0], -c * y[1]]


def tau20_state():
    d = np.load(MATCH)
    q = d['light_mean_q'][:2].astype(float); p = d['light_mean_p'][:2].astype(float)
    return np.array([q[0], q[1], p[0], p[1]])


def orbit_state(name):
    y4 = tau20_state()
    if name == 'late':
        q = charge(y4); y4 = np.array([LATE_R0, 0.0, 0.0, q / LATE_R0])
    return y4


def orbit_info(y4):
    """T_r from successive pericentre passages (r' = 0, r'' > 0) by dense output + root refinement; J_r from
    (1/2pi) int p_r^2 ds over one radial period; omega_r = 2 pi/T_r; Omega from the angle advance."""
    q = charge(y4)
    E = 0.5 * (y4[2] ** 2 + y4[3] ** 2) + 0.25 * LX * (y4[0] ** 2 + y4[1] ** 2) ** 2
    # integrate 3 rough periods with dense output
    rough_T = 2 * np.pi / np.sqrt(3 * LX * (2 * E / LX) ** 0.5)   # rough scale only
    sol = solve_ivp(quartic_rhs, (0, 6 * rough_T), y4, method='DOP853', rtol=1e-13, atol=1e-12, dense_output=True)
    ts = np.linspace(0, sol.t[-1], 200001); Y = sol.sol(ts)
    r = np.hypot(Y[0], Y[1]); rdot = (Y[0] * Y[2] + Y[1] * Y[3]) / r
    # pericentres: rdot crosses from - to +
    idx = np.where((rdot[:-1] < 0) & (rdot[1:] >= 0))[0]
    from scipy.optimize import brentq
    f = lambda t: (lambda y: (y[0] * y[2] + y[1] * y[3]))(sol.sol(t))
    tp = [brentq(f, ts[i], ts[i + 1], xtol=1e-14) for i in idx[:3]]
    Tr = tp[1] - tp[0]
    # J_r and angle advance over [tp0, tp1]
    from scipy.integrate import quad
    pr2 = lambda t: (lambda y: ((y[0] * y[2] + y[1] * y[3]) ** 2 / (y[0] ** 2 + y[1] ** 2)))(sol.sol(t))
    J = quad(pr2, tp[0], tp[1], epsabs=0, epsrel=1e-12, limit=400)[0] / (2 * np.pi)
    thdot = lambda t: (lambda y: ((y[0] * y[3] - y[1] * y[2]) / (y[0] ** 2 + y[1] ** 2)))(sol.sol(t))
    dth = quad(thdot, tp[0], tp[1], epsabs=0, epsrel=1e-12, limit=400)[0]
    y0 = sol.sol(tp[0]); rmin = np.hypot(y0[0], y0[1])
    # apocentre between tp0 and tp1
    ta = brentq(f, tp[0] + 0.3 * Tr, tp[0] + 0.7 * Tr, xtol=1e-14); ya = sol.sol(ta); rmax = np.hypot(ya[0], ya[1])
    return dict(q=q, E=E, T_r=Tr, omega_r=2 * np.pi / Tr, J_per_q=J / q, dtheta=dth, Omega=dth / Tr, rmin=rmin, rmax=rmax,
                Mperi=GG / np.sqrt(LP) * rmin, Mapo=GG / np.sqrt(LP) * rmax, t_peri0=tp[0])


def phase_state(y4, Tr, ph):
    if ph == 0:
        return np.array(y4, float)
    sol = solve_ivp(quartic_rhs, (0, ph * Tr), y4, method='DOP853', rtol=1e-13, atol=1e-12)
    return sol.y[:, -1].copy()


_ORB = {}


def orbit(name):
    if name not in _ORB:
        y4 = orbit_state(name); _ORB[name] = (y4, orbit_info(y4))
    return _ORB[name]


# ------------------------------------------------------------------ phi start: instantaneous minimum near 0
def chi_start(r2, M, epst=EPST):
    """Smallest root of dV/dchi = 0 at fixed r^2, M (Newton from the quadratic estimate), as the track does."""
    den = G2 * r2 - LP * M * M + 6 * epst * M * M
    p = 3 * epst * M ** 3 / den if den > 0 else 0.0
    for _ in range(60):
        d = p - M
        F = G2 * r2 * p + LP * (p * p - M * M) * p - 3 * epst * M * d * d
        dF = G2 * r2 + LP * (3 * p * p - M * M) - 6 * epst * M * d
        st = F / dF; p -= st
        if abs(st) < 1e-16 * max(1.0, abs(p)):
            break
    return p


# ------------------------------------------------------------------ measures (own formulas)
def free_measures(Y, s, par, q, epst=EPST):
    """Y[:, 8] samples.  n_-/q = [(E_X,free/m_X) - q]/(2q) with E_X,free = p_X^2/2 + m_X^2 r^2/2, m_X = g M;
    variant with m_X = g|chi|; N_phi/q with dchi = chi - M, dchi' = chi' - M', m_phi^2 = 2 lamP M^2.
    Also the energies of the budget."""
    M, Mp = Mof(par, s)
    r2 = Y[:, 0] ** 2 + Y[:, 1] ** 2; p2 = Y[:, 3] ** 2 + Y[:, 4] ** 2
    mX = GG * M
    EX = 0.5 * p2 + 0.5 * mX ** 2 * r2
    nminus = (EX / mX - q) / (2 * q)
    mXa = GG * np.abs(Y[:, 2])
    nminus_a = ((0.5 * p2 + 0.5 * mXa ** 2 * r2) / mXa - q) / (2 * q)
    mphi = np.sqrt(2 * LP) * M
    dchi = Y[:, 2] - M; dchid = Y[:, 5] - Mp
    Ephi = 0.5 * dchid ** 2 + 0.5 * mphi ** 2 * dchi ** 2
    Nphi = Ephi / mphi / q
    # exact phi potential energy about the vacuum and interaction energy
    Vphi = 0.25 * LP * (Y[:, 2] ** 2 - M ** 2) ** 2 - epst * M * dchi ** 3
    Eint = 0.5 * G2 * (Y[:, 2] ** 2 - M ** 2) * r2
    Equart = 0.25 * LX * r2 ** 2
    return dict(nminus=nminus, nminus_actual=nminus_a, Nphi=Nphi, EX=EX, Ephi_lin=Ephi, Ephi_exact=0.5 * dchid ** 2 + Vphi,
                Eint=Eint, Equart=Equart, mX=mX, mphi=mphi)


def Mof(par, s):
    s = np.asarray(s, float)
    if par.growth == 0:
        M = par.M0 * np.exp(par.kappa * s); Mp = par.kappa * M
    else:
        M = par.M0 * (1 + par.kappa * s); Mp = par.M0 * par.kappa * np.ones_like(s)
    return M, Mp


def s_of_M(par, M):
    return np.log(M / par.M0) / par.kappa if par.growth == 0 else (M / par.M0 - 1) / par.kappa


def tav(v, t):
    return float(np.trapezoid(v, t) / (t[-1] - t[0]))


class Sim:
    """Wrapper of the C integrator with the extended-phase-space state."""
    def __init__(self, par, q3, p3, nper=32, s0=0.0):
        self.par = par; self.nper = nper
        self.y = arr(np.r_[q3, p3, s0, 0.0])
        self.E0 = self.E(); self.steps = 0

    def E(self):
        return LIB.sk_E(ctypes.byref(self.par), self.y.ctypes.data_as(DP))

    @property
    def s(self):
        return self.y[6]

    def W(self):
        return -self.y[7]

    def to(self, s1):
        if s1 > self.y[6]:
            self.steps += LIB.sk_run(ctypes.byref(self.par), self.y.ctypes.data_as(DP), float(s1), self.nper)
        return self.y

    def sample(self, ds, ns):
        out = np.empty((ns, 8))
        self.steps += LIB.sk_sample(ctypes.byref(self.par), self.y.ctypes.data_as(DP), float(ds), int(ns), self.nper,
                                    out.ctypes.data_as(DP))
        return out

    def hext_err(self):
        """(E - E0 - W)/|E|: conservation of the extended Hamiltonian."""
        E = self.E(); return (E - self.E0 - self.W()) / abs(E), E

    def state(self):
        return dict(y=self.y.copy(), E0=self.E0, steps=self.steps)

    def restore(self, st):
        self.y = arr(st['y']); self.E0 = float(st['E0']); self.steps = int(st['steps'])


def make_par(info, eps, growth=0, epst=EPST, M_start_fac=0.5, method=8, om_fix=None, fm=1.5, rate_match_M=None):
    """kappa: exponential: eps*omega_r.  linear (growth 1): matched so that dlnM/ds at M = rate_match_M*Mperi equals
    eps*omega_r, i.e. kappa_lin = eps omega_r (rate_match_M Mperi)/M0."""
    M0 = M_start_fac * info['Mperi']
    if growth == 0:
        kappa = eps * info['omega_r']
    else:
        Mm = (rate_match_M if rate_match_M is not None else 1.0) * info['Mperi']
        kappa = eps * info['omega_r'] * Mm / M0
    if om_fix is None:
        # fastest early frequency: phi in the g^2 r^2 well at apocentre ~ g rmax; X ~ sqrt(3 lamX rmax^2); margin x1.5
        om_fix = 1.5 * max(GG * info['rmax'], np.sqrt(3 * LX) * info['rmax'], info['omega_r'])
    return SkPar(LX, G2, LP, epst, M0, kappa, om_fix, fm, growth, method)


def save_json(path, obj):
    def conv(x):
        if isinstance(x, (np.floating,)): return float(x)
        if isinstance(x, (np.integer,)): return int(x)
        if isinstance(x, np.ndarray): return x.tolist()
        raise TypeError(type(x))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w') as f:
        json.dump(obj, f, indent=1, default=conv)
