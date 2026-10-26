#!/usr/bin/env python3
"""ls_threshold (Claude, own code): the case-L lump benchmark at the actual epsilon_net and yields of candidates
10 and 17 of b448f25, plus premise probes.  Thermodynamics from lsthermo.py (Bessel series).

Inputs [Q]: light target g_X=0.1, lambda_X=1e-3, lambda_phi=0.03, v=0.924089637 GeV; Y_Q=4.8294e-9;
DM 4.4e-10 GeV per visible entropy; lump proxy Q=200000, mu_c=0.861195 m_X, E_c/Q=0.923316 m_X;
threshold 4->2 coefficients A_a, A_b (gp_leftovers amplitudes.json); cold conversion sigma v
1.13039e-16 cm^3/s (bb211d1); coarse SM g* table (same nodes as the benchmark; monotone g*s below 1 MeV).

Modes
  proxy : the benchmark as specified (single chemical potential, isentropic hidden state, the X-bar fraction p
          relaxes to p_eq with Gamma = A_a n+^2 n-/4 + A_b n+^3/6 + A_b n+ n-^2/6 + (n_c/Q) sigma v).
          Own integrator: p_eq, n_c, n_Q, H tabulated on a log grid, interpolated linearly in ln a, and
          d ln p / d ln a integrated with scipy Radau (rtol 1e-9).
  ke    : premise probe.  Same rates, but the hidden state keeps the actual X-bar number (two chemical
          potentials, kinetic equilibrium) and every reaction conserves energy, so out-of-equilibrium
          annihilation produces entropy.  Rates carry the detailed-balance factor (1 - exp(-(mu+ + mu-)/T)).
          Lumps form when mu+ reaches mu_c; optionally anti-lumps when mu- reaches mu_c (charge conjugate).
          Operator splitting per step: adiabatic expansion at fixed numbers and entropy, then a backward-Euler
          reaction step at fixed energy density.
  s_channel (either mode): adds the on-shell light-scalar route X Xbar <-> s (inverse decay, m_s > 2 m_X)
          followed by s s -> X Xbar, with n_s = K(T) n+ n- (MB, relativistic) and the threshold
          sigma v(ss -> X Xbar) from our Berends-Giele amplitude (|M| = 0.026667).
"""
import time
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from scipy.special import kve
import lsthermo as TH

# ---------------- inputs [Q] ----------------
GX, LX, LP, V = 0.1, 1e-3, 0.03, 0.924089637375362
M = GX * V                      # m_X
MS = np.sqrt(2 * LP) * V        # m_s (radial phi quantum)
YQ, DM, MRED, Q = 4.8294009076988175e-9, 4.4e-10, 2.435e18, 200000
UC, EC = 0.8611951308217225, 0.9233161353726057
AA, AB = 0.0842146549366979, 0.002361976125302127
SV = 1.13039e-16 / 1.167e-17
SV_SS = 0.0001260281438027877   # GeV^-2, from amp42.amplitude (this track, run in NOTES)
XC = np.sqrt(LP * V * V / (2 * GX * GX))
RHO_C, NC, V0 = 3 * LX * XC**4, 2 * np.sqrt(2 * LX) * XC**3, LP * V**4 / 4

_TS = np.log([1e-6, 5e-4, 1e-3, .1, .15, .3, 1., 5., 100., 200., 1e19])
_G = [3.36, 3.36, 10.75, 10.75, 17.25, 61.75, 72.25, 86.25, 96.25, 106.75, 106.75]
_GS = [3.91, 3.91, 10.75, 10.75, 17.25, 61.75, 72.25, 86.25, 96.25, 106.75, 106.75]


def gstar(T): return float(np.interp(np.log(T), _TS, _G))
def gs(T): return float(np.interp(np.log(T), _TS, _GS))
def entropy(T): return 2 * np.pi**2 / 45 * gs(T) * T**3
def hubble(T): return np.sqrt(np.pi**2 * gstar(T) / 90) * T * T / MRED
def lna(T): return -np.log(gs(T)**(1 / 3) * T)


def K_s(t):
    """n_s^eq / (n_X^eq n_Xbar^eq) at zero chemical potentials, MB relativistic, in units m^-3."""
    r = MS / M
    return 2 * np.pi**2 * r * r * kve(2, r / t) / (t * kve(2, 1 / t)**2) * np.exp(-(r - 2) / t)


# ---------------- single chemical potential state (the benchmark's) ----------------
def state_single(nQ, sigma=None, energy=None, uc=None, ec=None):
    """Hidden state at net charge nQ with entropy sigma*nQ or energy density `energy`; lumps at uc."""
    uc = UC if uc is None else uc
    ec = EC if ec is None else ec
    n = nQ / M**3
    tgt = None if energy is None else energy / M**4

    def at(t):
        B = TH.bessel(t, uc)
        if TH.net_charge(t, uc, B) <= n:
            u = uc
        else:
            u = brentq(lambda x: TH.net_charge(t, x, B) - n, 0.0, uc, xtol=1e-15, rtol=1e-14)
        th = TH.pair(t, u, B)
        nc = max(n - th['nq'], 0.0) if u == uc else 0.0
        return u, nc, th

    def resid(lt):
        u, nc, th = at(np.exp(lt))
        return th['s'] - sigma * n if sigma is not None else th['rho'] + ec * nc - tgt

    lt = brentq(resid, np.log(1e-4), np.log(2e3), xtol=1e-13, rtol=1e-13)
    t = np.exp(lt); u, nc, th = at(t)
    return dict(t=t, u=u, nc=nc * M**3, npl=th['npl'] * M**3, nmi=th['nmi'] * M**3, s=th['s'] * M**3,
                rho=(th['rho'] + ec * nc) * M**4)


# ---------------- two chemical potentials (kinetic equilibrium) ----------------
def _species_fill(t, target, cap, B):
    """Solve n(t,u) = target for u <= cap; return (u, condensed excess) (dimensionless).
    ln n(u) is a log-sum of exponentials in u, hence convex and increasing: Newton started on the right of the
    root converges monotonically (no overshoot)."""
    k, k2, _ = B
    lt = np.log(target)
    lcap = TH.log_n(t, cap, B)
    if lcap <= lt:
        return cap, target - np.exp(lcap)
    u = min(1 + t * (lt - np.log(t * k2[0] / TH.TWO_PI2)), cap)   # MB k=1 estimate: n(u) >= target here
    if TH.log_n(t, u, B) < lt:
        u = cap
    for _ in range(200):
        e = np.exp(-(k - 1) * (1.0 - u) / t) * k2
        s1 = np.dot(e, 1 / k); s0 = e.sum()
        f = -(1.0 - u) / t + np.log(t / TH.TWO_PI2 * s1) - lt
        du = f * t * s1 / s0
        u -= du
        if abs(du) < 1e-15 * max(1.0, abs(u)):
            break
    return u, 0.0


def state_two(nQ, nm, S=None, E=None, antilumps=True, tguess=None):
    """Two-potential state: total X-bar number nm (gas + anti-lumps) is given, not equilibrated."""
    n, m_ = nQ / M**3, nm / M**3
    capm, ecm = (UC, EC) if antilumps else (0.995, 1.0)
    tgt = S / M**3 if S is not None else E / M**4

    def at(t):
        B = TH.bessel(t, max(UC, capm))
        um, cm = _species_fill(t, m_, capm, B)
        up, cp = _species_fill(t, n + m_, UC, B)
        a = TH.species(t, up, B); b = TH.species(t, um, B)
        return up, um, cp, cm, a, b

    def resid(lt):
        up, um, cp, cm, a, b = at(np.exp(lt))
        if S is not None:
            return np.log(max(a[3] + b[3], 1e-300) / tgt)
        return np.log((a[2] + b[2] + EC * cp + ecm * cm) / tgt)

    lo, hi = np.log(2e-5), np.log(2e3)
    if tguess is not None:
        a_, b_ = np.log(tguess) - 0.05, np.log(tguess) + 0.05
        ra, rb = resid(a_), resid(b_)
        while ra > 0 and a_ > lo:
            b_, rb = a_, ra; a_ -= 0.3; ra = resid(a_)
        while rb < 0 and b_ < hi:
            a_, ra = b_, rb; b_ += 0.3; rb = resid(b_)
        if ra <= 0 <= rb:
            lo, hi = a_, b_
    lt = brentq(resid, lo, hi, xtol=1e-13, rtol=1e-13)
    t = np.exp(lt); up, um, cp, cm, a, b = at(t)
    return dict(t=t, up=up, um=um, nc=cp * M**3, nac=cm * M**3, npl=a[0] * M**3, nmi=b[0] * M**3,
                s=(a[3] + b[3]) * M**3, rho=(a[2] + b[2] + EC * cp + ecm * cm) * M**4, P=(a[1] + b[1]) * M**4)


# ---------------- rates ----------------
def gamma_proxy(p, nQ, nc, t, s_channel=False):
    """Benchmark relaxation rate per X-bar (GeV), with optional on-shell s route (linearised)."""
    nmi = p * nQ; npl = nmi + nQ - nc
    g = AA * npl * npl * nmi / 4 + AB * npl**3 / 6 + AB * npl * nmi * nmi / 6 + nc / Q * SV
    if s_channel:
        K = K_s(t) / M**3
        g += SV_SS * K * K * npl * nmi * (npl + nmi) / (1 + K * (npl + nmi))
    return g


def rate_two(st, s_channel=False):
    """Net pair-destruction rate per volume (GeV^4) in the two-potential state."""
    npl, nmi, t = st['npl'], st['nmi'], st['t']
    x = (st['up'] + st['um']) / t
    db = -np.expm1(-x)
    r = npl * nmi * db * (AA * npl * nmi / 4 + AB * (npl**2 + nmi**2) / 6)
    r += SV * db * (st['nc'] / Q * nmi + st['nac'] / Q * npl)
    if s_channel:
        K = K_s(t) / M**3
        r += 0.5 * SV_SS * K * K * (npl * nmi)**2 * (-np.expm1(-2 * x)) / (1 + K * (npl + nmi))
    return r


# ---------------- runs ----------------
def onset(eps, R, ymult):
    nq0 = eps * NC; ye = YQ * ymult
    tv0 = brentq(lambda T: ye * entropy(T) - nq0, 1e-5, 1e4, xtol=1e-14, rtol=1e-14)
    E0 = RHO_C * (1 + R) + V0
    st = state_single(nq0, energy=E0)
    return nq0, ye, tv0, E0, st


def proxy_run(eps, R, ymult=1.0, endT=1e-8, N=1600, s_channel=False, keep=False):
    t0 = time.perf_counter()
    nq0, ye, tv0, E0, st0 = onset(eps, R, ymult)
    sigma = st0['s'] / nq0; p0 = st0['nmi'] / nq0
    Tv = np.geomspace(tv0, endT, N)
    L = np.array([lna(T) for T in Tv]); L -= L[0]
    nQ = ye * np.array([entropy(T) for T in Tv]); H = np.array([hubble(T) for T in Tv])
    rows = [st0] + [state_single(nq, sigma=sigma) for nq in nQ[1:]]
    lpeq = np.log([r['nmi'] / nq for r, nq in zip(rows, nQ)])
    cf = np.array([r['nc'] / nq for r, nq in zip(rows, nQ)])
    th = np.array([r['t'] for r in rows])
    lnq, lH = np.log(nQ), np.log(H)
    if np.any(np.diff(L) <= 0):
        raise ValueError('ln a not monotone')

    def rhs(x, y):
        lp = y[0]; p = np.exp(lp)
        nq = np.exp(np.interp(x, L, lnq)); h = np.exp(np.interp(x, L, lH))
        nc = np.interp(x, L, cf) * nq; t = np.exp(np.interp(x, L, np.log(th)))
        peq = np.exp(np.interp(x, L, lpeq))
        return [-(gamma_proxy(p, nq, nc, t, s_channel) / h) * (1 - peq / p)]

    sol = solve_ivp(rhs, (0, L[-1]), [np.log(p0)], method='Radau', rtol=1e-9, atol=1e-10, dense_output=keep)
    p = float(np.exp(sol.y[0, -1]))
    out = dict(eps=eps, R=R, ymult=ymult, endT=endT, N=N, s_channel=s_channel, Tvis_onset=tv0, sigma=sigma,
               Th0_over_m=st0['t'], cond0=st0['nc'] / nq0, p0=p0, p_final=p, f=2 * M * YQ * p / DM,
               peq_end=float(np.exp(lpeq[-1])), f_perfect_equilibrium_end=2 * M * YQ * np.exp(lpeq[-1]) / DM,
               max_cond=float(cf.max()), ode_status=int(sol.status), nfev=int(sol.nfev),
               runtime=time.perf_counter() - t0)
    if keep:
        sel = np.unique(np.r_[np.arange(0, N, max(1, N // 60)), N - 1])
        out['history'] = [dict(Tvis=Tv[i], Th_over_m=th[i], cond=cf[i], peq=float(np.exp(lpeq[i])),
                               p=float(np.exp(sol.sol(L[i])[0])),
                               G_over_H=gamma_proxy(float(np.exp(sol.sol(L[i])[0])), nQ[i], cf[i] * nQ[i], th[i], s_channel) / H[i])
                          for i in sel]
    return out


def ke_run(eps, R, ymult=1.0, endT=1e-8, N=1500, antilumps=True, s_channel=False, verbose=False):
    t0 = time.perf_counter()
    nq0, ye, tv0, E0, st0 = onset(eps, R, ymult)
    p = st0['nmi'] / nq0; sigma = st0['s'] / nq0
    Tv = np.geomspace(tv0, endT, N)
    hist = []; nroot = 0; nfall = 0; Lprev = lna(tv0); sig0 = sigma; tg = st0['t']
    prev = state_two(nq0, p * nq0, S=sigma * nq0, antilumps=antilumps, tguess=tg)
    a3prev = 1.0; ew_max = 0.0; ew_sum = 0.0
    for i in range(1, N):
        T = Tv[i]; nq = ye * entropy(T); dl = lna(T) - Lprev; Lprev = lna(T)
        h = hubble(np.sqrt(T * Tv[i - 1])); dt = dl / h
        st = state_two(nq, max(p, 1e-250) * nq, S=sigma * nq, antilumps=antilumps, tguess=tg)
        tg = st['t']
        # energy-work check of the adiabatic step: d(rho a^3) = -P d(a^3)
        a3 = a3prev * np.exp(3 * dl)
        res = (st['rho'] * a3 - prev['rho'] * a3prev + 0.5 * (st['P'] + prev['P']) * (a3 - a3prev)) / (st['rho'] * a3)
        ew_max = max(ew_max, abs(res)); ew_sum += res; a3prev = a3
        r = rate_two(st, s_channel)
        if abs(r * dt / nq) > 1e-10 * p:
            rho = st['rho']; nroot += 1
            peq = max(state_single(nq, energy=rho)['nmi'] / nq, 1e-250)

            def g(pp):
                s2 = state_two(nq, pp * nq, E=rho, antilumps=antilumps, tguess=tg)
                return pp - p + dt * rate_two(s2, s_channel) / nq
            lo, hi = (peq, p) if p > peq else (p, peq)
            if g(lo) * g(hi) < 0:
                try:
                    if lo > 0 and hi / lo > 10:      # many decades: solve in ln p'
                        p = np.exp(brentq(lambda x: g(np.exp(x)), np.log(lo), np.log(hi), xtol=1e-13, maxiter=400))
                    else:
                        p = brentq(g, lo, hi, xtol=1e-300, rtol=1e-11, maxiter=400)
                except RuntimeError:
                    nfall += 1
                    p = peq if abs(g(peq)) < abs(g(p)) else p
            else:
                nfall += 1
                p = peq if abs(g(peq)) < abs(g(p)) else p
            p = max(p, 1e-250)
            st = state_two(nq, p * nq, E=rho, antilumps=antilumps, tguess=tg)
            sigma = st['s'] / nq
        prev = st
        if i % max(1, N // 80) == 0 or i == N - 1:
            hist.append(dict(Tvis=T, Th_over_m=st['t'], up=st['up'], um=st['um'], cond=st['nc'] / nq,
                             anticond=st['nac'] / nq, p=p, sigma=sigma,
                             G_over_H=rate_two(st, s_channel) / max(p * nq, 1e-300) / h))
            if verbose:
                print({k: f'{v:.4g}' for k, v in hist[-1].items()}, flush=True)
    # optimistic bound on all later conversion, all charge in lumps, radiation era continuing (as benchmark)
    tail = ye * entropy(endT) / Q * SV / hubble(endT)
    return dict(eps=eps, R=R, ymult=ymult, endT=endT, N=N, antilumps=antilumps, s_channel=s_channel,
                sigma0=sig0, sigma_final=sigma, p_final=p, f=2 * M * YQ * p / DM, reaction_steps=nroot, root_fallbacks=nfall,
                energy_work_max_step_residual=ew_max, energy_work_cumulative_residual=ew_sum,
                tail_conversion_depth_bound=tail, history=hist, runtime=time.perf_counter() - t0)


def threshold(eps, ymult, lo=40.0, hi=200.0, iters=12, **kw):
    hist = []
    for _ in range(iters):
        mid = 0.5 * (lo + hi); r = proxy_run(eps, mid, ymult, **kw)
        hist.append(dict(R=mid, f=r['f'], sigma=r['sigma'], max_cond=r['max_cond']))
        if r['f'] > 0.01: hi = mid
        else: lo = mid
    return dict(eps=eps, ymult=ymult, bracket=[lo, hi], kw=kw, history=hist)
