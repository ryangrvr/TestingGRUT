"""mf_gemini SKEPTIC: independent mild-waterfall calculation (linear modes + nonlinear delta N), own code.

Differences from the track (mf_gemini/code/waterfall_dn.py):
  * e-fold time N throughout (track: cosmic time); scipy RK45 for modes, LSODA for trajectories;
  * modes start in the ADIABATIC vacuum (omega^2 = k^2/a^2 + m^2, with the d omega/dN term) at k/aH = 50
    (track: WKB at k/aH = 100), from a deeper pre-history (chi heavy at the start: m^2/H^2 >= ~9);
  * chi_rms with two cuts: super-horizon only (q = 1, as the track) and including the sub-horizon tachyonic
    band k/a < |m_chi(t)| (q = |m|/H);
  * nonlinear trajectories integrate sigma directly until m_sigma = g chi > KDROP H, then DROP sigma
    (its energy fraction at the drop is recorded) instead of WKB-averaging it;
  * patches in which sigma leaves the tachyonic band before chi traps it are detected and counted as
    NO WATERFALL (N = infinity) instead of being extrapolated with the log law;
  * V_SR families: quad (convex), lin (constant slope), hill (concave hilltop, max above sigma_c, eta = 0.015),
    log (CW-like, V0 A/2 ln((sigma^2 + s0^2)/sigma_c^2), s0 = sigma_c: a model choice).
Code units: M_P = 1, H0 = sqrt(V0/3) = 1; hub = H0/M_P sets the quantum normalisation.  Delta = p_c/sigma_c with
p_c = V_SR'(sigma_c)/V0 (the track's definition).

Run:  python3 -B sk_waterfall.py case <name>   (one case -> scratch json)
      python3 -B sk_waterfall.py combine       (-> ../results/sk_waterfall.json)
      python3 -B sk_waterfall.py list
"""
import json
import os
import sys
import time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.special import erf

MP = 2.435e18
AS = 2.100e-9
H0GEV = 1.3e11
HUB = H0GEV / MP
V2 = 8.0 / 3.0
LAM = 27.0 / 16.0
P_STAR = np.sqrt(2 * HUB**2 / (8 * np.pi**2 * AS))
P_C_TRACK = P_STAR * np.exp(0.015 * 36.0)        # the track's default slope at sigma_c (3.18e-4)
SCR = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_gemini_skeptic/wf'
RES = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
KDROP = 100.0


class VSR:
    def __init__(self, kind, Delta, p_c=P_C_TRACK, eta_hill=0.015):
        self.kind, self.D = kind, Delta
        self.sc = p_c / Delta
        self.p = p_c
        if kind == 'hill':
            self.eta = eta_hill
            self.sm = self.sc + p_c / eta_hill            # maximum above sigma_c, slope p_c at sigma_c
        if kind == 'log':
            # V/V0 = A/2 ln((s^2 + s0^2)/sc^2), s0 = sc: slope at sc = A/(2 sc) = p_c -> A = 2 p_c sc
            self.A = 2 * p_c * self.sc
            self.s0 = self.sc

    def f(self, s):     # V_SR / V0
        k = self.kind
        if k == 'quad':
            return 0.5 * self.D * s**2
        if k == 'lin':
            return self.p * (s - self.sc)
        if k == 'hill':
            return -0.5 * self.eta * (s - self.sm)**2 + 0.5 * self.eta * self.sm**2
        return 0.5 * self.A * np.log((s**2 + self.s0**2) / self.sc**2)

    def df(self, s):
        k = self.kind
        if k == 'quad':
            return self.D * s
        if k == 'lin':
            return self.p + 0.0 * s
        if k == 'hill':
            return -self.eta * (s - self.sm)
        return self.A * s / (s**2 + self.s0**2)

    def monotone(self):
        return self.kind in ('lin', 'hill')


def model(vsr):
    g2 = LAM * V2 / vsr.sc**2
    V0 = 0.25 * LAM * V2**2   # = 3

    def V(s, c):
        return 0.25 * LAM * (c * c - V2)**2 + 0.5 * g2 * s * s * c * c + V0 * vsr.f(s)

    def Vs(s, c):
        return g2 * s * c * c + V0 * vsr.df(s)

    def Vc(s, c):
        return LAM * (c * c - V2) * c + g2 * s * s * c

    def m2chi(s):
        return g2 * s * s - LAM * V2
    return dict(g2=g2, V=V, Vs=Vs, Vc=Vc, m2=m2chi, V0=V0)


def sigma_bg(vsr, mdl, Npre, Nafter):
    """sigma-only background (chi = 0) in N from Npre e-folds (slow roll) before sigma_c."""
    if vsr.kind == 'quad':
        s0 = vsr.sc * min(np.exp(vsr.D * Npre), 8.0)
    elif vsr.kind == 'log':
        s0 = 2.0 * vsr.sc          # chi heavy (m^2 = 13.5 H^2); the log roll from 2 sc takes ~1.1/Delta e-folds
    else:
        s0 = min(vsr.sc + vsr.p * Npre, 8.0 * vsr.sc)
    V = lambda s: mdl['V'](s, 0.0)
    ds0 = -mdl['Vs'](s0, 0.0) / V(s0)

    def rhs(N, y):
        s, ds = y
        e = 0.5 * ds * ds
        return [ds, -(3 - e) * (ds + mdl['Vs'](s, 0.0) / V(s))]

    def evc(N, y):
        return y[0] - vsr.sc
    evc.direction = -1
    sol = solve_ivp(rhs, (0, 400), [s0, ds0], method='LSODA', rtol=1e-11, atol=1e-16, events=[evc],
                    dense_output=True, max_step=0.02)
    Nc = sol.t_events[0][0]
    Nend = min(Nc + Nafter, sol.t[-1])
    N = np.arange(0.0, Nend, 0.01)
    s, ds = sol.sol(N)
    eps = 0.5 * ds**2
    H2 = V(s) / (3 - eps)
    m2H2 = mdl['m2'](s) / H2
    return dict(Nc=Nc, N=N, s=s, ds=ds, eps=eps, H=np.sqrt(H2), m2H2=m2H2, sol=sol,
                spl_eps=CubicSpline(N, eps), spl_lnH=CubicSpline(N, np.log(np.sqrt(H2))),
                spl_m2=CubicSpline(N, mdl['m2'](s)), s0=s0)


def mode(bg, k, Ni, Nout):
    """u_k(N) in H0 units; adiabatic vacuum at Ni.  Integrated as v = a^{3/2} u (never tiny, so the explicit
    integrator stays inside its stability region during long heavy pre-histories):
        v'' - eps v' + (Omega - 9/4 + 3 eps/2) v = 0,  Omega = (k^2/a^2 + m^2)/H^2."""
    H = np.exp(bg['spl_lnH'](Ni))
    m2 = float(bg['spl_m2'](Ni))
    kp2 = (k * np.exp(-Ni))**2
    w = np.sqrt(kp2 + m2)
    dw2 = -2 * kp2 + float(bg['spl_m2'](Ni, 1))
    v0 = 1.0 / np.sqrt(2 * w)
    dv0 = v0 * (-1j * w / H - dw2 / (4 * w * w))

    def rhs(N, y):
        e = bg['spl_eps'](N)
        Hh = np.exp(bg['spl_lnH'](N))
        om = ((k * np.exp(-N))**2 + bg['spl_m2'](N)) / Hh**2 - 2.25 + 1.5 * e
        return [y[2], y[3], e * y[2] - om * y[0], e * y[3] - om * y[1]]
    sol = solve_ivp(rhs, (Ni, Nout[-1]), [v0.real, v0.imag, dv0.real, dv0.imag], method='RK45',
                    rtol=1e-8, atol=1e-12 * abs(v0), t_eval=Nout)
    a32 = np.exp(-1.5 * sol.t)
    v = sol.y[0] + 1j * sol.y[1]
    dv = sol.y[2] + 1j * sol.y[3]
    return a32 * v, a32 * (dv - 1.5 * v)


def trajectory(vsr, mdl, y0, N0):
    """full two-field from (s, ds, c, dc) at N0; drop sigma at g c > KDROP H; returns N_end or inf (no waterfall)."""
    g = np.sqrt(mdl['g2'])

    def H2of(y):
        s, ds, c, dc = y
        e = 0.5 * (ds * ds + dc * dc)
        return mdl['V'](s, c) / (3 - e), e

    def rhs(N, y):
        s, ds, c, dc = y
        H2, e = H2of(y)
        return [ds, -(3 - e) * ds - mdl['Vs'](s, c) / H2, dc, -(3 - e) * dc - mdl['Vc'](s, c) / H2]

    def ev_drop(N, y):
        return g * abs(y[2]) - KDROP * np.sqrt(H2of(y)[0])
    ev_drop.terminal = True
    ev_drop.direction = 1

    def ev_end(N, y):
        return H2of(y)[1] - 1.0
    ev_end.terminal = True
    ev_end.direction = 1

    def ev_fail(N, y):      # sigma beyond the far edge of the band while chi is not trapping it
        return (y[0] + 1.2 * vsr.sc) if mdl['m2'](y[0]) > 0 else 1.0
    ev_fail.terminal = True
    ev_fail.direction = -1
    s1 = solve_ivp(rhs, (N0, N0 + 400), y0, method='LSODA', rtol=1e-10, atol=1e-20,
                   events=[ev_drop, ev_end, ev_fail], max_step=0.05)
    if s1.t_events[2].size:
        return dict(N_end=np.inf, fail=True, N_fail=float(s1.t_events[2][0] - N0))
    if s1.t_events[1].size:
        return dict(N_end=float(s1.t_events[1][0]), fail=False, dropped=False)
    if not s1.t_events[0].size:
        return dict(N_end=np.nan, fail=False, note='no event')
    yd = s1.y_events[0][0]
    Nd = s1.t_events[0][0]
    H2d, ed = H2of(yd)
    s, ds, c, dc = yd
    rho_sig = 0.5 * H2d * ds * ds + 0.5 * mdl['g2'] * s * s * c * c + mdl['V0'] * (vsr.f(s) - vsr.f(0.0))
    frac_sig = rho_sig / (3 * H2d)
    Vc1 = lambda c: 0.25 * LAM * (c * c - V2)**2 + mdl['V0'] * vsr.f(0.0)

    def rhs2(N, y):
        c, dc = y
        e = 0.5 * dc * dc
        H2 = Vc1(c) / (3 - e)
        return [dc, -(3 - e) * dc - LAM * (c * c - V2) * c / H2]

    def ev_end2(N, y):
        return 0.5 * y[1]**2 - 1.0
    ev_end2.terminal = True
    ev_end2.direction = 1
    s2 = solve_ivp(rhs2, (Nd, Nd + 200), [c, dc], method='LSODA', rtol=1e-11, atol=1e-20, events=[ev_end2])
    return dict(N_end=float(s2.t_events[0][0]), fail=False, dropped=True, N_drop=float(Nd), frac_sigma_at_drop=float(frac_sig),
                H_end=float(np.sqrt(Vc1(s2.y_events[0][0][0]) / 2.0)))


CASES = {
    'quad_D1': ('quad', 1.0), 'quad_D0.3': ('quad', 0.3), 'quad_D0.1': ('quad', 0.1), 'quad_D0.03': ('quad', 0.03),
    'lin_D0.1': ('lin', 0.1), 'lin_D0.06': ('lin', 0.06), 'hill_D0.1': ('hill', 0.1), 'hill_D0.06': ('hill', 0.06),
    'log_D1': ('log', 1.0), 'log_D0.3': ('log', 0.3), 'log_D0.1': ('log', 0.1), 'log_D0.03': ('log', 0.03),
}


def run_case(name, fm=1e-2, dlnk=0.2):
    t0 = time.time()
    kind, D = CASES[name]
    vsr = VSR(kind, D)
    mdl = model(vsr)
    Npre = min(max(8.0, 1.5 / D), 45.0) if kind in ('quad',) else 12.0
    Nlin = 32.0
    bg = sigma_bg(vsr, mdl, Npre, Nlin + 2.0)
    Nc = bg['Nc']
    N = bg['N']
    Nout = np.arange(N[0] + 0.5, N[-1] - 0.05, 0.05)
    lnH = bg['spl_lnH'](Nout)
    aH = np.exp(Nout + lnH)
    m2H2 = bg['spl_m2'](Nout) / np.exp(2 * lnH)
    Nex = np.arange(N[0] - 2.0, Nout[-1] - 0.5, dlnk)
    ks = np.exp(Nex + bg['spl_lnH'](Nex))
    P = np.zeros((ks.size, Nout.size))
    DU = np.zeros((ks.size, Nout.size), complex)
    U = np.zeros((ks.size, Nout.size), complex)
    for i, k in enumerate(ks):
        Ni = brentq_N(lambda n: np.log(k) - n - float(bg['spl_lnH'](n)) - np.log(50.0), N[0], Nout[-1])
        sel = Nout > Ni + 1e-9
        u, du = mode(bg, k, Ni, Nout[sel])
        U[i, sel], DU[i, sel] = u, du
        P[i, sel] = k**3 * np.abs(u)**2 / (2 * np.pi**2)
    out = dict(name=name, kind=kind, Delta=D, sigma_c=vsr.sc, g_phys=float(np.sqrt(mdl['g2']) * HUB), N_pre=Npre,
               s_start_over_sc=float(bg['s0'] / vsr.sc), m2H2_at_start=float(bg['m2H2'][0]))
    res = {}
    for qlab in ['q1', 'qtach']:
        if qlab == 'q1':
            q = np.ones_like(Nout)
        else:
            q = np.maximum(1.0, np.sqrt(np.maximum(-m2H2, 0.0)))
        mask = ks[:, None] < (q * aH)[None, :]
        chi2 = (P * mask).sum(0) * dlnk
        chirms = np.sqrt(chi2) * HUB
        # mode-weighted velocity d chi_rms/dN
        dchi2 = (2 * np.real(np.conj(U) * DU) * ks[:, None]**3 / (2 * np.pi**2) * mask).sum(0) * dlnk
        dchirms = 0.5 * dchi2 / np.maximum(chi2, 1e-300) * chirms
        chi_trap = np.sqrt(2 * vsr.p * vsr.sc / 3.0) * np.sqrt(3.0) / np.sqrt(3.0)
        inb = np.abs(np.interp(Nout, N, bg['s'])) < vsr.sc
        ok = np.where((chirms >= fm * chi_trap) & inb & (Nout > Nc))[0]
        if ok.size == 0:
            res[qlab] = dict(ok=False, reason='chi_rms never reaches the matching value inside the band',
                             chirms_max_over_trap=float(chirms.max() / chi_trap))
            continue
        im = int(ok[0])
        Nm = Nout[im]
        s_m, ds_m = bg['sol'].sol(Nm)
        cs = np.exp(np.linspace(np.log(1e-6), np.log(5.0), 25)) if kind in ('lin', 'hill') else \
            np.exp(np.linspace(np.log(1e-4), np.log(5.0), 19))
        traj = []
        for c in cs:
            r = trajectory(vsr, mdl, [s_m, ds_m, c * chirms[im], c * dchirms[im]], Nm)
            traj.append(r)
        Ncur = np.array([r['N_end'] for r in traj])
        fails = np.array([r.get('fail', False) for r in traj])
        r1 = trajectory(vsr, mdl, [s_m, ds_m, chirms[im], dchirms[im]], Nm)
        ra = trajectory(vsr, mdl, [s_m, ds_m, 1.02 * chirms[im], 1.02 * dchirms[im]], Nm)
        rb = trajectory(vsr, mdl, [s_m, ds_m, 0.98 * chirms[im], 0.98 * dchirms[im]], Nm)
        dNdlnc = (ra['N_end'] - rb['N_end']) / np.log(1.02 / 0.98)
        # MC over chi_L = |Gaussian| chi_rms
        rng = np.random.default_rng(11)
        gs = np.abs(rng.standard_normal(1_000_000))
        lc = np.log(cs)
        good = np.isfinite(Ncur)
        cfail = cs[fails].max() if fails.any() else 0.0
        x = np.log(gs)
        Nx = np.interp(x, lc[good], Ncur[good])
        lo = x < lc[good][0]
        sl = (Ncur[good][1] - Ncur[good][0]) / (lc[good][1] - lc[good][0])
        Nx[lo] = Ncur[good][0] + (x[lo] - lc[good][0]) * sl
        hi = x > lc[good][-1]
        Nx[hi] = Ncur[good][-1] + (x[hi] - lc[good][-1]) * (Ncur[good][-1] - Ncur[good][-2]) / (lc[good][-1] - lc[good][-2])
        nofall = gs < cfail if fails.any() else np.zeros_like(gs, bool)
        fin = ~nofall
        dN = Nx[fin] - Nx[fin].mean()
        # linear delta N spectrum at matching
        dNdchi = dNdlnc / chirms[im]
        ex = ks < aH[im]
        Pz = dNdchi**2 * P[ex, im] * HUB**2
        ipk = int(np.argmax(Pz))
        Nb = r1['N_end'] - Nex[ex][ipk]
        res[qlab] = dict(
            ok=True, N_match_after_sc=float(Nm - Nc), chi_rms_match=float(chirms[im]), dlnchi_dN_match=float(dchirms[im] / chirms[im]),
            N_wf_classical=float(r1['N_end'] - Nc), N_wf_mean=float(Nx[fin].mean() - Nc), s_eff=float(-1.0 / dNdlnc),
            varN=float(dN.var()), varN_analytic=float(np.pi**2 / 8 * dNdlnc**2),
            P_tail={f'{z}': float((dN > z).mean()) for z in [0.5, 0.67, 1.0]},
            P_no_waterfall=float(nofall.mean()), c_fail_max=float(cfail),
            Pzeta_peak=float(Pz[ipk]), Pzeta_peak_Nexit_after_sc=float(Nex[ex][ipk] - Nc), N_before_end_peak=float(Nb),
            frac_sigma_at_drop=r1.get('frac_sigma_at_drop'), H_end=r1.get('H_end'),
            Ncurve=dict(c=[float(c) for c in cs], N_after_sc=[float(n - Nc) if np.isfinite(n) else None for n in Ncur]))
    out['results'] = res
    out['runtime_s'] = time.time() - t0
    return out


def brentq_N(f, a, b):
    """start time of a mode: first N with k/aH <= 50 (f decreasing); the background start if already below."""
    from scipy.optimize import brentq
    if f(a) <= 0:
        return a
    return brentq(f, a, b)


def regression():
    """exact de Sitter (H = 1, eps = 0): massless and m^2 = 2 spectra after exit, and m^2 = -4.5 growth rate."""
    out = {}
    for m2 in [0.0, 2.0, -4.5]:
        bg = dict(spl_eps=lambda N, d=0: 0.0 * N, spl_lnH=lambda N, d=0: 0.0 * N,
                  spl_m2=lambda N, d=0, m2=m2: (m2 + 0.0 * N) if d == 0 else 0.0 * N)
        k = np.exp(np.log(50.0))
        Nout = np.array([np.log(50.0) + 3.0, np.log(50.0) + 6.0])
        u, du = mode(bg, k, 0.0, Nout)
        Pk = k**3 * np.abs(u)**2 / (2 * np.pi**2)
        x = np.exp(-3.0)
        exact = {0.0: (1 / (2 * np.pi))**2 * (1 + x**2), 2.0: (1 / (2 * np.pi))**2 * x**2, -4.5: None}[m2]
        out[f'm2={m2}'] = dict(P_after_3_efolds=float(Pk[0]), exact=exact,
                               growth_rate_3_to_6=float(0.5 * np.log(Pk[1] / Pk[0]) / 3.0))
    out['s_exact'] = float((-3 + np.sqrt(9 + 18)) / 2)
    return out


def main():
    os.makedirs(SCR, exist_ok=True)
    if sys.argv[1] == 'list':
        print(list(CASES))
    elif sys.argv[1] == 'case':
        for name in sys.argv[2:]:
            fn = os.path.join(SCR, name + '.json')
            if os.path.exists(fn):
                continue
            r = run_case(name)
            with open(fn, 'w') as fh:
                json.dump(r, fh)
            print(name, json.dumps({q: {k: v for k, v in d.items() if k != 'Ncurve'} for q, d in r['results'].items()}),
                  'runtime', round(r['runtime_s'], 1), flush=True)
    elif sys.argv[1] == 'regress':
        print(json.dumps(regression(), indent=1))
    elif sys.argv[1] == 'combine':
        out = dict(note=__doc__.split('\n')[0], H0_GeV=H0GEV, p_c=P_C_TRACK, p_star=P_STAR, regression=regression(), cases=[])
        for name in CASES:
            fn = os.path.join(SCR, name + '.json')
            if os.path.exists(fn):
                with open(fn) as fh:
                    out['cases'].append(json.load(fh))
        with open(os.path.join(RES, 'sk_waterfall.json'), 'w') as fh:
            json.dump(out, fh, indent=1)


if __name__ == '__main__':
    main()
