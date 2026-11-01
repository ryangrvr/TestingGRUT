"""Common helpers for track pt_m1tilt: the phi-condensation transition under Codex's GRUT-M1 potential (d4747cb).
Copied (minimal parts) from the pt_matching skeptic's sk_common.py; same orbits, same start (M0 = 0.5 M_c,peri), same
growth law M = M0 exp(kappa s), kappa = eps omega_r, same final window (60 X periods at M_end = 40 M_c,apo), same
work-integral and charge checks.  Only numeric input: the tau = 20 light mean field MATCH_L12.npz (light_mean_q/p).
Potential modes (see m1_dop.c): 1 = skeptic model (optional linear tilt beta), 3 = GRUT-M1 with tilt epst.
[C] homogeneous classical model; not a cosmological result."""
import ctypes, os, subprocess, json
import numpy as np
from scipy.integrate import solve_ivp, quad

HERE = os.path.dirname(os.path.abspath(__file__))
SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad'
F55 = SP + '/f55/research/reality_framework'
MATCH = F55 + '/grut_early_nonlinear_20261027/MATCH_L12.npz'
RESULTS = os.path.join(os.path.dirname(HERE), 'results')
SCR = SP + '/review/pt_m1tilt'

LX, G2, LP = 1e-3, 1e-2, 3e-2
GG = np.sqrt(G2)
EPST_M1 = 0.002          # Codex M1: epsilon = .002, lambda_phi = .03
LATE_R0 = 31.494634874218747


def tau20_state():
    d = np.load(MATCH)
    q = d['light_mean_q'][:2].astype(float); p = d['light_mean_p'][:2].astype(float)
    return np.array([q[0], q[1], p[0], p[1]])


def charge(y):
    return y[0] * y[3] - y[1] * y[2]


def energy4(y, lam=LX, m2=0.0):
    r2 = y[0] ** 2 + y[1] ** 2
    return 0.5 * (y[2] ** 2 + y[3] ** 2) + 0.5 * m2 * r2 + 0.25 * lam * r2 * r2


def orbit_rhs(lam, m2=0.0):
    def f(t, y):
        c = m2 + lam * (y[0] ** 2 + y[1] ** 2)
        return [y[2], y[3], -c * y[0], -c * y[1]]
    return f


def integrate_orbit(y4, T, lam=LX, m2=0.0, t_eval=None, rtol=1e-13):
    sol = solve_ivp(orbit_rhs(lam, m2), (0.0, T), y4, method='DOP853', rtol=rtol, atol=1e-12 * np.max(np.abs(y4)),
                    t_eval=t_eval, dense_output=t_eval is None)
    return sol


def radial_quad(E, q, lam=LX, m2=0.0):
    """J_r = (1/pi) int p_r dr, T_r, angle advance; u = r^2 roots of lam/2 u^3 + m2 u^2 - 2E u + q^2 = 0."""
    if lam > 0:
        rts = np.sort(np.roots([lam / 2, m2, -2 * E, q * q]).real)
        u0, u1, u2 = rts
        P = lambda u: lam / 2 * u ** 3 + m2 * u * u - 2 * E * u + q * q
        dP = lambda u: 1.5 * lam * u * u + 2 * m2 * u - 2 * E
        for _ in range(5):
            u1 -= P(u1) / dP(u1); u2 -= P(u2) / dP(u2)
        G = lambda u: (lam / 2) * (u - u0) / u
    else:
        disc = np.sqrt(E * E - m2 * q * q)
        u1, u2 = (E - disc) / m2, (E + disc) / m2
        G = lambda u: m2 / u
    D = u2 - u1
    U = lambda t: u1 + D * np.sin(t) ** 2
    fJ = lambda t: np.sqrt(G(U(t))) * D * D * (np.sin(t) * np.cos(t)) ** 2 / np.sqrt(U(t))
    fT = lambda t: 1.0 / (np.sqrt(U(t)) * np.sqrt(G(U(t))))
    fA = lambda t: q / U(t) / (np.sqrt(U(t)) * np.sqrt(G(U(t))))
    fR = lambda t: np.sqrt(U(t)) / np.sqrt(G(U(t)))     # for <r^2> over the radial period
    kw = dict(epsabs=0, epsrel=2e-14, limit=200)
    J = quad(fJ, 0, np.pi / 2, **kw)[0] / np.pi
    Tr = 2 * quad(fT, 0, np.pi / 2, **kw)[0]
    Dth = 2 * quad(fA, 0, np.pi / 2, **kw)[0]
    r2mean = 2 * quad(fR, 0, np.pi / 2, **kw)[0] / Tr
    return dict(J_r=J, T_r=Tr, dtheta=Dth, rmin=np.sqrt(u1), rmax=np.sqrt(u2), r2mean=r2mean)


def orbit_info(y4, lam=LX, m2=0.0):
    q = charge(y4); E = energy4(y4, lam, m2)
    o = radial_quad(E, q, lam, m2)
    return dict(q=q, E=E, J_per_q=o['J_r'] / q, T_r=o['T_r'], omega_r=2 * np.pi / o['T_r'], Omega=o['dtheta'] / o['T_r'],
                rmin=o['rmin'], rmax=o['rmax'], r2mean=o['r2mean'])


def phase_state(y4, Tr, ph):
    if ph == 0:
        return np.array(y4, float)
    sol = integrate_orbit(y4, ph * Tr, t_eval=[ph * Tr])
    return sol.y[:, -1].copy()


def orbit(name):
    y4 = tau20_state()
    if name == 'late':
        q = charge(y4); y4 = np.array([LATE_R0, 0.0, 0.0, q / LATE_R0])
    return y4, orbit_info(y4)


# ---------------------------------------------------------------- C library (own adaptive DOP853)
class SPar(ctypes.Structure):
    _fields_ = [('mode', ctypes.c_int), ('growth', ctypes.c_int), ('lamX', ctypes.c_double), ('g2', ctypes.c_double),
                ('lamP', ctypes.c_double), ('M0', ctypes.c_double), ('kappa', ctypes.c_double), ('beta', ctypes.c_double),
                ('epst', ctypes.c_double)]


def _lib():
    os.makedirs(SCR, exist_ok=True)
    so = os.path.join(SCR, 'libm1dop.so')
    src = os.path.join(HERE, 'm1_dop.c')
    if not os.path.exists(so) or os.path.getmtime(so) < max(os.path.getmtime(src), os.path.getmtime(os.path.join(HERE, 'dop853_coef.h'))):
        subprocess.check_call(['gcc', '-O2', '-fPIC', '-shared', '-o', so, src, '-lm'])
    lib = ctypes.CDLL(so)
    dp = ctypes.POINTER(ctypes.c_double)
    lib.m1_run.restype = ctypes.c_long
    lib.m1_run.argtypes = [ctypes.POINTER(SPar), dp, dp, ctypes.c_double, ctypes.c_double, ctypes.c_double, dp,
                           ctypes.POINTER(ctypes.c_long)]
    lib.m1_energy.restype = ctypes.c_double
    lib.m1_energy.argtypes = [ctypes.POINTER(SPar), ctypes.c_double, dp]
    return lib


LIB = _lib()
DP = ctypes.POINTER(ctypes.c_double)


class Runner:
    def __init__(self, par, y6, s0=0.0, rtol=1e-11, atol=1e-300, h=1e-3):
        self.par = par
        self.y = np.ascontiguousarray(np.r_[y6, 0.0], dtype=float)
        self.s = ctypes.c_double(s0); self.h = ctypes.c_double(h)
        self.rtol, self.atol = rtol, atol
        self.nrej = ctypes.c_long(0); self.nacc = 0
        self.E0 = self.energy()

    def to(self, s1):
        self.nacc += LIB.m1_run(ctypes.byref(self.par), self.y.ctypes.data_as(DP), ctypes.byref(self.s), s1, self.rtol,
                                self.atol, ctypes.byref(self.h), ctypes.byref(self.nrej))
        return self.y

    def energy(self):
        return LIB.m1_energy(ctypes.byref(self.par), self.s.value, self.y.ctypes.data_as(DP))


def Mof(par, s):
    s = np.asarray(s, float)
    if par.growth == 0:
        M = par.M0 * np.exp(par.kappa * s); Mp = par.kappa * M
    else:
        M = par.M0 * (1 + par.kappa * s); Mp = par.M0 * par.kappa * np.ones_like(s)
    return M, Mp


def vacuum_x(par):
    """Late vacuum chi = x M.  mode 3 (M1): x = 1 exactly (+v is the only stationary point).
    mode 1 with linear tilt beta: x^3 - x = beta, largest root (and the false vacuum on the negative side)."""
    if par.mode == 3:
        return (1.0, 1.0) if par.epst != 0 else (1.0, -1.0)   # epst = 0: Z2-symmetric model, vacuum on the side of chi
    rts = np.sort(np.roots([1.0, 0.0, -1.0, -par.beta]))
    rr = np.sort(rts[np.abs(rts.imag) < 1e-12].real)
    xp, xm = rr[-1], (rr[0] if (len(rr) == 3 and rr[0] < 0) else rr[-1])
    return xp, xm


def measures(par, s, Y, q):
    """n_-/q with the reference mass g x M (x = 1 for M1) and with the actual mass g|chi|; N_phi/q about chi = x M.
    m_phi^2 = V''(xM): M1 gives 2 lamP M^2 exactly; the linear tilt gives lamP (3x^2 - 1) M^2."""
    M, Mp = Mof(par, s)
    xp, xm = vacuum_x(par)
    x = np.where(Y[:, 2] >= 0, xp, xm)
    r2 = Y[:, 0] ** 2 + Y[:, 1] ** 2; pX2 = Y[:, 3] ** 2 + Y[:, 4] ** 2
    m = GG * np.abs(x) * M
    nM = ((0.5 * pX2 + 0.5 * m * m * r2) / m - q) / (2 * q)
    mf = GG * np.abs(Y[:, 2])
    nF = ((0.5 * pX2 + 0.5 * mf * mf * r2) / mf - q) / (2 * q)
    mphi = np.sqrt(LP * (3 * x * x - 1)) * M
    dphi = Y[:, 2] - x * M; dphid = Y[:, 5] - x * Mp
    Nphi = (0.5 * dphid ** 2 + 0.5 * mphi ** 2 * dphi ** 2) / mphi / q
    return nM, nF, Nphi, M


def phi_quasistatic(r2, M, epst, lamP=LP, g2=G2):
    """Instantaneous minimum of the M1 phi potential near chi = 0 at fixed r^2 and M (Newton from the quadratic
    estimate 3 epst M^3 / (g2 r^2 - lamP M^2 + 6 epst M^2)).  Returns (phi_newton, phi_quadratic_estimate)."""
    den = g2 * r2 - lamP * M * M + 6 * epst * M * M
    p = 3 * epst * M ** 3 / den if den > 0 else 0.0
    p0 = p
    for _ in range(50):
        d = p - M
        F = g2 * r2 * p + lamP * (p * p - M * M) * p - 3 * epst * M * d * d
        dF = g2 * r2 + lamP * (3 * p * p - M * M) - 6 * epst * M * d
        if dF <= 0: break
        step = F / dF; p -= step
        if abs(step) < 1e-15 * max(1.0, abs(p)): break
    return p, p0


def run_measure(y4, info, eps, phi0, M_end_fac=40.0, checkpoints=(5.0, 10.0, 20.0), growth=0, beta=0.0, epst=0.0,
                mode=3, rtol=1e-11, final_periods=60, nsamp=2000, M_start_fac=0.5, phi_dot0=0.0, phi0_rel=None):
    """phi0: initial chi in units of M_c,peri (used when phi0_rel is None); phi0_rel = 1 + delta multiplies phi0
    (seed-perturbation pairs).  phi0 = 'qs' puts chi at the instantaneous tilted minimum (mode 3)."""
    q = charge(y4)
    Mperi = np.sqrt(G2 / LP) * info['rmin']; Mapo = np.sqrt(G2 / LP) * info['rmax']
    M0 = M_start_fac * Mperi; kappa = eps * info['omega_r']
    par = SPar(mode, growth, LX, G2, LP, M0, kappa, beta, epst)
    r2_0 = y4[0] ** 2 + y4[1] ** 2
    if isinstance(phi0, str) and phi0 == 'qs':
        chi0, chi0_quad = phi_quasistatic(r2_0, M0, epst if mode == 3 else 0.0)
    else:
        chi0 = float(phi0) * Mperi; chi0_quad = None
    if phi0_rel is not None:
        chi0 *= phi0_rel
    y6 = [y4[0], y4[1], chi0 if mode != 2 else 0.0, y4[2], y4[3], phi_dot0 if mode != 2 else 0.0]
    R = Runner(par, y6, rtol=rtol)
    s_of = (lambda M: np.log(M / M0) / kappa) if growth == 0 else (lambda M: (M / M0 - 1) / kappa)
    res = dict(eps=eps, phi0=phi0 if isinstance(phi0, str) else float(phi0), chi0=float(chi0), chi0_over_Mperi=float(chi0 / Mperi),
               chi0_quadratic_estimate=chi0_quad, growth=growth, beta=beta, epst=epst, mode=mode, q=q, M0=M0,
               Mperi=Mperi, Mapo=Mapo, kappa=kappa, windows={})
    M_escape = None; M_dev = None; M_half = None
    stops = sorted(set(list(checkpoints) + [M_end_fac]))
    grid = np.geomspace(M0, M_end_fac * Mapo, 500)[1:]
    events = sorted([(Mg, None) for Mg in grid] + [(f * Mapo, f) for f in stops], key=lambda x: x[0])
    trace = []
    for Mg, f in events:
        if s_of(Mg) <= R.s.value and f is None:
            continue
        if s_of(Mg) > R.s.value:
            R.to(s_of(Mg))
        if mode != 2:
            if M_escape is None and abs(R.y[2]) > 0.1 * Mg:
                M_escape = Mg / Mperi
            if M_half is None and R.y[2] > 0.5 * Mg:
                M_half = Mg / Mperi
            if M_dev is None and mode == 3:
                # deviation from the quasi-static tilted minimum evaluated with the initial orbit's <r^2> (rough)
                pq, _ = phi_quasistatic(info['r2mean'], Mg, epst)
                if abs(R.y[2] - pq) > 0.1 * Mg:
                    M_dev = Mg / Mperi
            if f is None:
                trace.append((float(Mg / Mperi), float(R.y[2] / Mg)))
        if f is not None:
            s0 = R.s.value; Tw = final_periods * 2 * np.pi / (GG * Mg)
            ts = s0 + np.linspace(0, Tw, nsamp + 1)
            Y = np.empty((nsamp + 1, 7)); Y[0] = R.y
            for i, t in enumerate(ts[1:], 1):
                Y[i] = R.to(t)
            nM, nF, Nphi, M = measures(par, ts, Y, q)
            av = lambda v: float(np.trapezoid(v, ts) / (ts[-1] - ts[0]))
            qq = Y[:, 0] * Y[:, 4] - Y[:, 1] * Y[:, 3]
            E = R.energy()
            Edyn = abs(E) if mode == 2 else abs(0.5 * (Y[-1, 3] ** 2 + Y[-1, 4] ** 2 + Y[-1, 5] ** 2) + 0.5 * G2 * Y[-1, 2] ** 2 * (Y[-1, 0] ** 2 + Y[-1, 1] ** 2))
            res['windows'][f'{f:g}'] = dict(nminus=av(nM), nminus_actualmass=av(nF), Nphi=av(Nphi),
                                            nminus_std=float(np.std(nM)), qrel=float(np.max(np.abs(qq / q - 1))),
                                            work_rel=float(abs(E - R.E0 - R.y[6]) / abs(E)),
                                            work_abs=float(E - R.E0 - R.y[6]), E=float(E), E_dyn=float(Edyn))
    res['M_escape_over_Mperi'] = M_escape       # first grid M with |chi| > 0.1 M
    res['M_half_over_Mperi'] = M_half           # first grid M with chi > 0.5 M
    res['M_dev_over_Mperi'] = M_dev
    # settle: smallest grid M after which chi stays on the +M side (chi > 0.5 M at every later grid point)
    if trace:
        cm = np.array([t[1] for t in trace]); Ms = np.array([t[0] for t in trace])
        bad = np.where(cm <= 0.5)[0]
        res['M_settle_over_Mperi'] = float(Ms[bad[-1] + 1]) if len(bad) and bad[-1] + 1 < len(Ms) else (float(Ms[0]) if not len(bad) else None)
    res['trace_M_chiOverM'] = trace
    res['steps'] = R.nacc; res['rejected'] = R.nrej.value
    fin = res['windows'][f'{M_end_fac:g}']
    res.update(nminus=fin['nminus'], Nphi=fin['Nphi'], qrel=max(w['qrel'] for w in res['windows'].values()),
               work_rel=max(w['work_rel'] for w in res['windows'].values()),
               work_rel_dyn=max(abs(w['work_abs']) / w['E_dyn'] for w in res['windows'].values()))
    return res


def save_json(name, obj):
    os.makedirs(RESULTS, exist_ok=True)
    def conv(x):
        if isinstance(x, (np.floating,)): return float(x)
        if isinstance(x, (np.integer,)): return int(x)
        if isinstance(x, np.ndarray): return x.tolist()
        raise TypeError(type(x))
    with open(os.path.join(RESULTS, name), 'w') as f:
        json.dump(obj, f, indent=1, default=conv)
