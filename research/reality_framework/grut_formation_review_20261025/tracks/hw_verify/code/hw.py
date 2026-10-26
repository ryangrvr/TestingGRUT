#!/usr/bin/env python3
"""hw_verify track (Claude, own code): helpers for the independent check of ls_scaling's heavy-mediator window.

Release engine: hwrel.c (own conformal-time, conformal-field Yoshida-8 / DP5 integrator with run-time couplings).
Preparations (own minimisers, direct variables, Newton-polished, Hessian-checked):
  * 'codex'   : aligned radial minimum of V_infl = V - cH|X|^2 on X = x > 0, Y = -y, Sigma = -s, then
                Sigma -> -s e^{i delta} with Y left at its delta = 0 value (Codex 373d37e / b448f25 prescription);
  * 'relaxed' : conditional minimum over (x, Re Y, Im Y, s) at fixed arg X = 0, Sigma = -s e^{i delta}
                (ls_scaling's variant: what a heavy Y tracks during inflation when the phase is held).
Ledger [D, re-derived]: late-particle decay ledger of Codex b448f25 / ls_scaling,
  R_Y    = (N_YH + N_YL) (m_L/H) <a>_L / E_X,  <a> = Gamma(3/2) sqrt(2 H_i/Gamma) = sqrt(pi H_i/(2 Gamma)),
  R_Sig  = 2 N_sig (m_A/H) <a>_A / E_X,  R_casc = N_YH (m_A/H) <a>_A / E_X,  m_A = 6 e h,
  Gamma_L = 3 kap^2 m_L/(256 pi^3),  Gamma_A = e^2 m_A/(48 pi) (1 - 4 m_X^2/m_A^2)^{3/2}.
Threshold: case-L R_L(eps) at normalised yield from our sp_critic rows (own crossing + interpolation).
All quantities in H_i units unless the name says GeV.
"""
import json, os, subprocess, warnings
from pathlib import Path
import numpy as np
from scipy.optimize import minimize, fsolve

HERE = Path(__file__).resolve().parent
TRACKS = HERE.parents[1]
SCR = Path(os.environ.get('HWV_SCRATCH', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/'
                          'scratchpad/review/hw_verify/build'))
BIN = SCR / 'hwrel'
MRED, GSTAR = 2.435e18, 106.75
YQ = 4.8294009076988175e-9
MX_LIGHT = 0.1 * 0.924089637      # g_X v [Q]
H_CAND = 434867798206.6179         # Codex candidates' H_i (c_H = 10) [Q]
CH = 10.0
STATED = dict(lx=1e-3, ly=1.0, ls=1.0, kap=1e-3, M_GeV=1e12, h_GeV=1e11, mu_GeV=5e11, e=1e-7)
QN = ["NX", "NY", "NS", "EX", "eps", "NYH", "NYL", "Nsig", "Gauss", "H", "EYexc", "Esig", "absX", "absY", "absS",
      "rhophi_frac", "Ncirc", "EYfull"]


def build():
    src = HERE / 'hwrel.c'
    SCR.mkdir(parents=True, exist_ok=True)
    if not BIN.exists() or BIN.stat().st_mtime < src.stat().st_mtime:
        subprocess.run(['gcc', '-O3', '-o', str(BIN), str(src), '-lm'], check=True)


def heavy(m=30.0, r=0.45, hH=40.0, frac=0.5, ly=1.0, Hs=H_CAND, e=1e-7, cH=CH, lx=1e-3):
    """ls_scaling's heavy-direction family, re-derived: alpha6 = r kap^2/(2 M^2 (1-r^2)) [D],
    eta = alpha6 X0^2/lx with X0^2 = cH/(2 lx), eta_edge = r/12 [Q, gp_critic], frac = eta/eta_edge
    => kap = lx m sqrt(frac (1-r^2)/(3 cH)).  mu = r M^2/(sqrt2 h)."""
    kap = lx * m * np.sqrt(frac * (1 - r * r) / (3 * cH))
    M = m * Hs; h = hH * Hs; mu = r * M * M / (np.sqrt(2) * h)
    return dict(lx=lx, ly=ly, ls=1.0, kap=kap, M_GeV=M, h_GeV=h, mu_GeV=mu, e=e)


def pars(cp, H):
    return dict(lx=cp['lx'], ly=cp['ly'], ls=cp['ls'], kap=cp['kap'], M=cp['M_GeV'] / H, h=cp['h_GeV'] / H,
                mu=cp['mu_GeV'] / H)


# ---------------------------------------------------------------- potential and preparations (own)
def V(z, p, cH):
    x, y, s = z
    ax2, ay2, as2 = abs(x)**2, abs(y)**2, abs(s)**2
    return float(-cH * ax2 + p['lx'] * ax2**2 + p['ly'] * ay2**2 + p['ls'] * (as2 - p['h']**2 / 2)**2
                 + p['M']**2 * ay2 + 2 * np.real(p['kap'] * np.conj(y) * x**3 + p['mu'] * np.conj(s) * y * y))


def F(z, p, cH):
    """dV/dz* (complex force), inflationary potential."""
    x, y, s = z
    return np.array([2 * p['lx'] * abs(x)**2 * x + 3 * p['kap'] * y * np.conj(x)**2 - cH * x,
                     p['M']**2 * y + 2 * p['ly'] * abs(y)**2 * y + p['kap'] * x**3 + 2 * p['mu'] * s * np.conj(y),
                     2 * p['ls'] * (abs(s)**2 - p['h']**2 / 2) * s + p['mu'] * y * y])


def _grad4(q, p, cH, delta):
    """Gradient of V in the 4 real variables (x, Yr, Yi, s) with X = x, Sigma = -s e^{i delta}: dV/dx = 2 Re(F_X Xhat*),
    dV/dYr = 2 Re F_Y, dV/dYi = 2 Im F_Y, dV/ds = 2 Re(F_S Shat*)."""
    x, yr, yi, s = q
    z = np.array([x, yr + 1j * yi, -s * np.exp(1j * delta)])
    f = F(z, p, cH)
    return np.array([2 * np.real(f[0]), 2 * np.real(f[1]), 2 * np.imag(f[1]),
                     2 * np.real(f[2] * np.conj(-np.exp(1j * delta)))])


def _hess4(q, p, cH, delta, d=1e-6):
    Hm = np.zeros((4, 4))
    for j in range(4):
        e = np.zeros(4); e[j] = d * max(abs(q[j]), 1e-6)
        Hm[:, j] = (_grad4(q + e, p, cH, delta) - _grad4(q - e, p, cH, delta)) / (2 * e[j])
    return 0.5 * (Hm + Hm.T)


def _polish(q0, p, cH, delta):
    f4 = lambda q: V([q[0], q[1] + 1j * q[2], -q[3] * np.exp(1j * delta)], p, cH)
    r = minimize(f4, q0, jac=lambda q: _grad4(q, p, cH, delta), method='BFGS', options=dict(gtol=1e-12, maxiter=20000))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        q = fsolve(lambda q: _grad4(q, p, cH, delta), r.x, xtol=1e-15)
    if np.linalg.norm(_grad4(q, p, cH, delta)) > np.linalg.norm(_grad4(r.x, p, cH, delta)):
        q = r.x
    ev = np.linalg.eigvalsh(_hess4(q, p, cH, delta))
    return q, dict(V=f4(q), grad=float(np.linalg.norm(_grad4(q, p, cH, delta))), min_hess=float(ev.min()))


def aligned_min(cp, H, cH):
    """Codex-style aligned inflationary minimum (delta = 0): own multi-start in (x, Yr, Yi, s), kept if Y is real
    negative (Yi = 0 by symmetry at delta = 0)."""
    p = pars(cp, H)
    xr = np.sqrt(cH / (2 * p['lx']))
    best = None
    for fx in (0.7, 1.0, 1.4):
        for y0 in (-p['kap'] * (fx * xr)**3 / p['M']**2, -0.1, -1.0, -5.0):
            for s0 in (p['h'] / np.sqrt(2), 2 * p['h'] / np.sqrt(2)):
                q, info = _polish(np.array([fx * xr, y0, 0.0, s0]), p, cH, 0.0)
                if q[0] > 0 and q[3] > 0 and info['min_hess'] > 0 and (best is None or info['V'] < best[1]['V']):
                    best = (q, info)
    q, info = best
    return dict(x=float(q[0]), y=float(-q[1]), yi=float(q[2]), s=float(q[3]), **info)


def prep_codex(cp, H, cH, delta):
    m = aligned_min(cp, H, cH)
    z = np.array([m['x'], -m['y'] + 1j * m['yi'], -m['s'] * np.exp(1j * delta)])
    return z, dict(prep='codex', x=m['x'], Y=[-m['y'], m['yi']], s=m['s'], V_aligned=m['V'])


def prep_relaxed(cp, H, cH, delta):
    """Conditional minimum at the supplied phase; started from the Codex state and from Y-rotated variants, the
    lowest V kept.  Returns the state and the X-phase torque (dV/d arg X), which is not zero: the phase is held."""
    p = pars(cp, H)
    m = aligned_min(cp, H, cH)
    best = None
    for rot in (0.0, 0.5, -0.5, 1.0, -1.0, 2.0):
        y0 = -m['y'] * np.exp(1j * rot * delta / 2)
        q, info = _polish(np.array([m['x'], y0.real, y0.imag, m['s']]), p, cH, delta)
        if q[0] > 0 and q[3] > 0 and info['min_hess'] > 0 and (best is None or info['V'] < best[1]['V'] - 1e-9):
            best = (q, info)
    q, info = best
    z = np.array([q[0], q[1] + 1j * q[2], -q[3] * np.exp(1j * delta)])
    f = F(z, p, cH)
    torque = float(np.imag(np.conj(z[0]) * f[0]))      # dV/d(arg X) = -2 Im(X* F_X) ... sign: reported as Im(X* F_X)
    return z, dict(prep='relaxed', x=float(q[0]), Y=[float(q[1]), float(q[2])], s=float(q[3]), V=info['V'],
                   cond_grad=info['grad'], min_hess=info['min_hess'], torque_ImXF=torque)


# ---------------------------------------------------------------- release
def run(cp, H, z0, a_end=100.0, method=8, cpar=0.05, windows=(1, 10, 30, 100, 300, 1000), ns=60, bg=0, gam=0.0,
        cH=CH, v0=None):
    build()
    p = pars(cp, H)
    g = lambda v: '%.17g' % float(v)
    win = [w for w in windows if w <= a_end * (1 + 1e-12)]
    if win[-1] < a_end:
        win.append(a_end)
    args = [str(BIN), str(method), g(cpar), str(bg), g(gam), g(cH), g(p['lx']), g(p['ly']), g(p['ls']), g(p['kap']),
            g(p['M']), g(p['h']), g(p['mu']), g(a_end), str(ns), str(len(win))] + [g(w) for w in win]
    for zz in z0:
        args += [g(zz.real), g(zz.imag)]
    if v0 is not None:
        for zz in v0:
            args += [g(complex(zz).real), g(complex(zz).imag)]
    pr = subprocess.run(args, capture_output=True, text=True)
    if pr.returncode:
        raise RuntimeError(pr.stderr or pr.stdout[-500:])
    S, W, Fin, END = [], {}, None, None
    for line in pr.stdout.splitlines():
        f = line.split(); v = list(map(float, f[1:]))
        if f[0] == 'S':
            S.append(dict(a=v[0], **dict(zip(QN, v[1:]))))
        elif f[0] == 'W':
            W['%g-%g' % (v[0], v[1])] = dict(span=v[2], **dict(zip(QN, v[3:])))
        elif f[0] == 'F':
            Fin = v
        elif f[0] == 'END':
            END = dict(steps=int(v[0]), rejects=int(v[1]), max_gauss=v[2])
    return dict(samples=S, windows=W, final=Fin, end=END)


# ---------------------------------------------------------------- ledger, threshold, yield, survival
def masses(cp):
    M, h, mu, e = cp['M_GeV'], cp['h_GeV'], cp['mu_GeV'], cp['e']
    mL2 = M**2 - np.sqrt(2) * mu * h
    return dict(mYH=np.sqrt(M**2 + np.sqrt(2) * mu * h), mYL=np.sqrt(mL2) if mL2 > 0 else float('nan'),
                msig=np.sqrt(2 * cp['ls']) * h, mA=6 * e * h, r=np.sqrt(2) * mu * h / M**2)


def ledger(cp, H, NYH, NYL, Nsig, EX):
    m = masses(cp); kap, e = cp['kap'], cp['e']
    GL = 3 * kap**2 * m['mYL'] / (256 * np.pi**3)
    thr = (1 - 4 * MX_LIGHT**2 / m['mA']**2)**1.5
    GA = e**2 * m['mA'] / (48 * np.pi) * thr
    aL = np.sqrt(np.pi * H / (2 * GL)); aA = np.sqrt(np.pi * H / (2 * GA))
    dL = m['mYL'] / H * aL; dA = m['mA'] / H * aA
    RY = (NYH + NYL) * dL / EX; RS = 2 * Nsig * dA / EX; RC = NYH * dA / EX
    return dict(R_Y=RY, R_Sigma=RS, R_cascade=RC, R_lower=RY + RS, R_central=RY + RS + RC,
                deposit_per_Y=dL, deposit_per_A=dA, a_dec_L=aL, a_dec_A=aA, Gamma_L_GeV=GL, Gamma_A_GeV=GA,
                masses=m)


def _rl_table():
    """Own crossing of f = 0.01 (log-linear in f between bracketing R) from the sp_critic case-L rows."""
    rows = []
    for fn in ('rmax_relic.json', 'rmax_relic_L.json', 'rmax_relic_Llow.json'):
        for r in json.load(open(TRACKS / 'sp_critic' / 'results' / fn))['rows']:
            if r.get('case') == 'L' and 'f_extra' in r:
                rows.append((round(r['eps'], 6), r['R'], r['f_extra']))
    tab = []
    for eps in sorted({r[0] for r in rows}):
        pts = sorted({(R, f) for e, R, f in rows if e == eps})
        for (r0, f0), (r1, f1) in zip(pts, pts[1:]):
            if f0 <= 0.01 < f1 and r0 > 0:
                t = (np.log(0.01) - np.log(f0)) / (np.log(f1) - np.log(f0))
                tab.append((eps, float(np.exp(np.log(r0) + t * np.log(r1 / r0))), r0, r1))
                break
    return tab


RL_TAB = _rl_table()


def R_L(eps):
    """Case-L threshold at normalised yield, log-log interpolation in |eps| over our table; None outside it."""
    e = abs(eps)
    xs = [np.log(a) for a, *_ in RL_TAB]; ys = [np.log(b) for _, b, *_ in RL_TAB]
    if e < RL_TAB[0][0] or e > RL_TAB[-1][0]:
        return None
    return float(np.exp(np.interp(np.log(e), xs, ys)))


def yield_over_YQ(NX, H):
    T = (90 * MRED**2 * H**2 / (np.pi**2 * GSTAR))**0.25
    s = 2 * np.pi**2 / 45 * GSTAR * T**3
    return abs(NX) * H**3 / s / YQ


def alpha6_GeV(cp):
    M, h, mu, kap = cp['M_GeV'], cp['h_GeV'], cp['mu_GeV'], cp['kap']
    Sig = h / np.sqrt(2)
    return mu * kap**2 * Sig / (M**4 - 4 * (mu * Sig)**2)


def c6(cp):
    return alpha6_GeV(cp) * MRED**2


def kappa_max(cp):
    """Tree-level boundedness along the aligned direction [D, as ls_scaling]: kap < (lx/(1.5 2^{-1/3}))^{3/4} ly^{1/4}."""
    return (cp['lx'] / (1.5 * 2**(-1 / 3)))**0.75 * cp['ly']**0.25


def codex_envelope_bound(cp, H, cH, x0):
    """Codex 2908e91's conditional leading-EFT envelope bound evaluated at a given point [Q formula, C numbers]:
    R = r_i (1 + 2/cH)^{1/4}, |a^3 n_X| <= 12 |alpha6| R^6 / H_i (H_i units: alpha6_H = alpha6_GeV H^2)."""
    p = pars(cp, H)
    a6 = alpha6_GeV(cp) * H**2
    Renv = x0 * (1 + 2 / cH)**0.25
    NXmax = 12 * abs(a6) * Renv**6
    m = masses(cp)
    Sig = p['h'] / np.sqrt(2)
    gap2 = p['M']**2 - 2 * p['mu'] * Sig
    yforced = p['kap'] * Renv**3 / gap2
    return dict(alpha6_H=a6, R_env=Renv, NX_bound=NXmax, yield_bound_over_YQ=yield_over_YQ(NXmax, H),
                kappa_X2_over_M2=p['kap'] * Renv**2 / p['M']**2,
                X_freq_over_Y_gap=np.sqrt(2 * p['lx']) * Renv / np.sqrt(gap2),
                X_freq_over_Sigma_mass=np.sqrt(2 * p['lx']) * Renv / (np.sqrt(2 * p['ls']) * p['h']),
                Y_quartic_relative_force=2 * p['ly'] * yforced**2 / gap2,
                Sigma_relative_forced_shift=p['mu'] * yforced**2 / ((np.sqrt(2 * p['ls']) * p['h'])**2 * Sig))


def summarize(cp, H, w, prep_info=None, extra=None):
    L = ledger(cp, H, w['NYH'], w['NYL'], w['Nsig'], w['EX'])
    eps = w['eps']; RL = R_L(eps)
    out = dict(eps=eps, NX=w['NX'], EX=w['EX'], NY=w['NY'], NS=w['NS'], NYH=w['NYH'], NYL=w['NYL'], Nsig=w['Nsig'],
               EYexc=w['EYexc'], R_Y=L['R_Y'], R_Sigma=L['R_Sigma'], R_cascade=L['R_cascade'], R_lower=L['R_lower'],
               R_central=L['R_central'], R_L=RL, central_over_RL=(L['R_central'] / RL) if RL else None,
               yield_over_YQ=yield_over_YQ(w['NX'], H), c6=c6(cp), kappa=cp['kap'], kappa_max=kappa_max(cp),
               masses=L['masses'], deposit_per_Y=L['deposit_per_Y'], deposit_per_A=L['deposit_per_A'],
               couplings_over_H=dict(M=cp['M_GeV'] / H, h=cp['h_GeV'] / H, mu=cp['mu_GeV'] / H), H_GeV=H)
    if prep_info:
        out['prep_info'] = prep_info
    if extra:
        out.update(extra)
    return out
