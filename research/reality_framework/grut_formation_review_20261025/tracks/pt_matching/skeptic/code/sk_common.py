"""Common helpers for the pt_matching SKEPTIC track (independent of both Codex and pt_matching code).

Methods deliberately differ from pt_matching:
  * orbit integration: scipy DOP853 (adaptive explicit RK, rtol 1e-13) instead of an order-8 symplectic composition;
  * radial action / period: scipy.integrate.quad in the angle variable r^2 = u1 + (u2-u1) sin^2 t, roots by numpy;
  * transition: own C adaptive DOP853 (sk_dop.c) with a work integral as extra ODE component.
Only numeric input: the tau = 20 light mean field MATCH_L12.npz (light_mean_q/p) of the reviewed early packet."""
import ctypes, os, subprocess, json
import numpy as np
from scipy.integrate import solve_ivp, quad

HERE = os.path.dirname(os.path.abspath(__file__))
SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad'
F55 = SP + '/f55/research/reality_framework'
MATCH = F55 + '/grut_early_nonlinear_20261027/MATCH_L12.npz'
RESULTS = os.path.join(os.path.dirname(HERE), 'results')
SCR = SP + '/review/pt_matching_skeptic'
TRACK_SCAN = SP + '/review/pt_matching/scan'

LX, G2, LP = 1e-3, 1e-2, 3e-2
GG = np.sqrt(G2)
Q_ORIG_NET = 8113.477086664833  # [Q] COHERENT_MATCHING.json


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
        # refine the two positive roots by Newton
        P = lambda u: lam / 2 * u ** 3 + m2 * u * u - 2 * E * u + q * q
        dP = lambda u: 1.5 * lam * u * u + 2 * m2 * u - 2 * E
        for _ in range(5):
            u1 -= P(u1) / dP(u1); u2 -= P(u2) / dP(u2)
        G = lambda u: (lam / 2) * (u - u0) / u   # p_r^2 = G(u)(u-u1)(u2-u)
    else:
        disc = np.sqrt(E * E - m2 * q * q)
        u1, u2 = (E - disc) / m2, (E + disc) / m2
        G = lambda u: m2 / u
    D = u2 - u1
    U = lambda t: u1 + D * np.sin(t) ** 2
    # p_r dr with r = sqrt(u): dr = du/(2 sqrt u), du = 2 D sin cos dt, p_r = sqrt(G) D sin cos
    fJ = lambda t: np.sqrt(G(U(t))) * D * D * (np.sin(t) * np.cos(t)) ** 2 / np.sqrt(U(t))
    fT = lambda t: 1.0 / (np.sqrt(U(t)) * np.sqrt(G(U(t))))
    fA = lambda t: q / U(t) / (np.sqrt(U(t)) * np.sqrt(G(U(t))))
    kw = dict(epsabs=0, epsrel=2e-14, limit=200)
    J = quad(fJ, 0, np.pi / 2, **kw)[0] / np.pi
    Tr = 2 * quad(fT, 0, np.pi / 2, **kw)[0]
    Dth = 2 * quad(fA, 0, np.pi / 2, **kw)[0]
    return dict(J_r=J, T_r=Tr, dtheta=Dth, rmin=np.sqrt(u1), rmax=np.sqrt(u2))


def orbit_info(y4, lam=LX, m2=0.0, nsamp=4096):
    """Radial action by quadrature; lab negative-frequency statistic from the exact quasi-periodic Fourier series
    (psi e^{-i Omega t} periodic in T_r; sampled by scipy DOP853 dense output)."""
    q = charge(y4); E = energy4(y4, lam, m2)
    o = radial_quad(E, q, lam, m2)
    Tr, Dth = o['T_r'], o['dtheta']
    t = np.arange(nsamp) * Tr / nsamp
    sol = integrate_orbit(y4, Tr, lam, m2, t_eval=np.r_[t, Tr])
    psi = sol.y[0] + 1j * sol.y[1]
    closure = abs(psi[-1] - psi[0] * np.exp(1j * Dth)) / abs(psi[0])
    psi = psi[:-1]
    Om, wr = Dth / Tr, 2 * np.pi / Tr
    c = np.fft.fft(psi * np.exp(-1j * Om * t)) / nsamp
    n = np.fft.fftfreq(nsamp, 1.0 / nsamp).round().astype(int)
    w = Om + n * wr; P = np.abs(c) ** 2
    return dict(q=q, E=E, J_per_q=o['J_r'] / q, A_minus_per_q=-np.sum((w * P)[w < 0]) / q,
                A_plus_per_q=np.sum((w * P)[w > 0]) / q, q_spec_rel=np.sum(w * P) / q - 1,
                J_spec_per_q=np.sum(n * w * P) / q, T_r=Tr, omega_r=wr, Omega=Om, rmin=o['rmin'], rmax=o['rmax'],
                closure=closure, n=n, w=w, P=P)


def phase_state(y4, Tr, ph):
    if ph == 0:
        return np.array(y4, float)
    sol = integrate_orbit(y4, ph * Tr, t_eval=[ph * Tr])
    return sol.y[:, -1].copy()


# ---------------------------------------------------------------- C library (own adaptive DOP853)
class SPar(ctypes.Structure):
    _fields_ = [('mode', ctypes.c_int), ('growth', ctypes.c_int), ('lamX', ctypes.c_double), ('g2', ctypes.c_double),
                ('lamP', ctypes.c_double), ('M0', ctypes.c_double), ('kappa', ctypes.c_double), ('beta', ctypes.c_double)]


def _lib():
    os.makedirs(SCR, exist_ok=True)
    so = os.path.join(SCR, 'libskdop.so')
    src = os.path.join(HERE, 'sk_dop.c')
    if not os.path.exists(so) or os.path.getmtime(so) < max(os.path.getmtime(src), os.path.getmtime(os.path.join(HERE, 'dop853_coef.h'))):
        subprocess.check_call(['gcc', '-O2', '-fPIC', '-shared', '-o', so, src, '-lm'])
    lib = ctypes.CDLL(so)
    dp = ctypes.POINTER(ctypes.c_double)
    lib.sk_run.restype = ctypes.c_long
    lib.sk_run.argtypes = [ctypes.POINTER(SPar), dp, dp, ctypes.c_double, ctypes.c_double, ctypes.c_double, dp,
                           ctypes.POINTER(ctypes.c_long)]
    lib.sk_energy.restype = ctypes.c_double
    lib.sk_energy.argtypes = [ctypes.POINTER(SPar), ctypes.c_double, dp]
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
        self.nacc += LIB.sk_run(ctypes.byref(self.par), self.y.ctypes.data_as(DP), ctypes.byref(self.s), s1, self.rtol,
                                self.atol, ctypes.byref(self.h), ctypes.byref(self.nrej))
        return self.y

    def energy(self):
        return LIB.sk_energy(ctypes.byref(self.par), self.s.value, self.y.ctypes.data_as(DP))


def Mof(par, s):
    s = np.asarray(s, float)
    if par.growth == 0:
        M = par.M0 * np.exp(par.kappa * s); Mp = par.kappa * M
    else:
        M = par.M0 * (1 + par.kappa * s); Mp = par.M0 * par.kappa * np.ones_like(s)
    return M, Mp


def measures(par, s, Y, q):
    """n_-/q with the reference mass g M (track definition) and with the actual mass g|phi|; N_phi/q."""
    M, Mp = Mof(par, s)
    # late vacuum phi_min = x M with x^3 - x = beta (x = 1 without tilt); m_X = g x M, m_phi^2 = lamP (3x^2 - 1) M^2.
    # With a tilt the vacuum sits at phi = +x M; the measure below uses phi itself (no |phi|) when beta != 0.
    # minimum on the side of phi: x_+ = largest root of x^3 - x - beta; x_- = smallest real root (false vacuum) if present
    rts = np.sort(np.roots([1.0, 0.0, -1.0, -par.beta]))
    rr = np.sort(rts[np.abs(rts.imag) < 1e-12].real)
    xp, xm = rr[-1], (rr[0] if (len(rr) == 3 and rr[0] < 0) else rr[-1])
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


def run_measure(y4, info, eps, seed, M_end_fac=40.0, checkpoints=(5.0, 10.0, 20.0), growth=0, beta=0.0, mode=1,
                rtol=1e-11, final_periods=60, nsamp=2000, M_start_fac=0.5, phi_dot0=0.0):
    q = charge(y4)
    Mperi = np.sqrt(G2 / LP) * info['rmin']; Mapo = np.sqrt(G2 / LP) * info['rmax']
    M0 = M_start_fac * Mperi; kappa = eps * info['omega_r']
    par = SPar(mode, growth, LX, G2, LP, M0, kappa, beta)
    y6 = [y4[0], y4[1], seed * Mperi if mode == 1 else 0.0, y4[2], y4[3], phi_dot0 * seed * Mperi if mode == 1 else 0.0]
    R = Runner(par, y6, rtol=rtol)
    s_of = (lambda M: np.log(M / M0) / kappa) if growth == 0 else (lambda M: (M / M0 - 1) / kappa)
    res = dict(eps=eps, seed=seed, growth=growth, beta=beta, mode=mode, q=q, windows={})
    M_escape = None
    stops = sorted(set(list(checkpoints) + [M_end_fac]))
    grid = np.geomspace(M0, M_end_fac * Mapo, 500)[1:]
    # merge window starts into the grid
    events = sorted([(Mg, None) for Mg in grid] + [(f * Mapo, f) for f in stops], key=lambda x: x[0])
    for Mg, f in events:
        if s_of(Mg) <= R.s.value and f is None:
            continue
        if s_of(Mg) > R.s.value:
            R.to(s_of(Mg))
        if M_escape is None and mode == 1 and abs(R.y[2]) > 0.1 * Mg:
            M_escape = Mg / Mperi
        if f is not None:
            # window of final_periods X periods at this mass; sample uniformly; continue the SAME run afterwards
            s0 = R.s.value; Tw = final_periods * 2 * np.pi / (GG * Mg)
            ts = s0 + np.linspace(0, Tw, nsamp + 1)
            Y = np.empty((nsamp + 1, 7)); Y[0] = R.y
            for i, t in enumerate(ts[1:], 1):
                Y[i] = R.to(t)
            nM, nF, Nphi, M = measures(par, ts, Y, q)
            av = lambda v: float(np.trapezoid(v, ts) / (ts[-1] - ts[0]))
            qq = Y[:, 0] * Y[:, 4] - Y[:, 1] * Y[:, 3]
            E = R.energy()
            res['windows'][f'{f:g}'] = dict(nminus=av(nM), nminus_actualmass=av(nF), Nphi=av(Nphi),
                                            nminus_std=float(np.std(nM)), qrel=float(np.max(np.abs(qq / q - 1))),
                                            work_rel=float(abs(E - R.E0 - R.y[6]) / abs(E)))
    res['M_escape_over_Mperi'] = M_escape
    res['steps'] = R.nacc; res['rejected'] = R.nrej.value
    fin = res['windows'][f'{M_end_fac:g}']
    res.update(nminus=fin['nminus'], Nphi=fin['Nphi'], qrel=max(w['qrel'] for w in res['windows'].values()),
               work_rel=max(w['work_rel'] for w in res['windows'].values()))
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
