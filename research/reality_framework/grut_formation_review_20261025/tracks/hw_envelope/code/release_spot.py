#!/usr/bin/env python3
"""hw_envelope (Claude, own code): independent three-field release at the ls_scaling heavy spot point, testing each
premise of Codex 2908e91's envelope bound directly on the full solution.

Own integrator (Python, scipy DOP853, cosmic time, physical fields, H_i = 1 units, instant radiation a = sqrt(2t)):
  V = lx|X|^4 + ly|Y|^4 + ls(|S|^2 - h^2/2)^2 + M^2|Y|^2 + 2Re(kap Y* X^3) + 2Re(mu S* Y^2),
  inflation only: -cH|X|^2 (sets the prepared state; removed instantly at a = 1; zero physical velocities).
Preparation (own minimiser): conditional minimum over |X|, complex Y, |Sigma| at arg X = 0, Sigma = -|Sigma|e^{i delta}
("Y-relaxed", as ls_scaling); option 'codex': Y kept at its delta = 0 aligned value (Codex b448f25 prescription).
Diagnostics along the solution:
  a^3 n_X (charge), E_X, eps = a^3 n_X / N_circ (Codex/ls convention), max a|X| versus the envelope R
  (premise test), the running bound integral 12 int a^3 |alpha6| |X|^6 dt, the actual source -6 Im(kap Y* X^3) against
  the forced-response EFT source -12 Im(alpha6 X^6) with instantaneous Sigma (matching test), the free Y excitation
  (Codex's excitation quadratic form in Y - Y_*, split into the two Sigma-frame modes -> N_YH, N_YL), the Sigma radial
  energy -> N_sigma, and the processed R_Y, R_Sigma, R_casc from the closed-form ledger of ls_scaling [D].
One-field EFT comparator: X alone with lam|X|^4 + beta6|X|^6 + (alpha6 X^6 + h.c.), same cH and phase.
Usage: python3 -B release_spot.py <out.json> [m=30] [a_end=30] [prep=relaxed|codex] [frac=0.5] [r=0.45] [hH=40]
"""
import json, math, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid
from scipy.optimize import minimize, fsolve

MRED, GSTAR, YQ, MXL = 2.435e18, 106.75, 4.8294009076988175e-9, 0.092408963737536
H_SPOT_GEV = 1.3335e11 * (434867798206.6179 * 5.889032093211876**(-2 / 3) / 1.3335e11)   # yield-normalised scale


def entropy(H):
    T = (90 * MRED**2 * H**2 / (math.pi**2 * GSTAR))**0.25
    return 2 * math.pi**2 / 45 * GSTAR * T**3


def couplings(m=30.0, r=0.45, hH=40.0, frac=0.5, cH=10.0, lx=1e-3, ly=1.0, ls=1.0, e=1e-7):
    kap = m * lx * math.sqrt(frac * 4 * (1 - r * r) / (12 * cH))     # eta = frac * r/12 with eta = alpha6 X0^2/lx
    mu = r * m * m / (math.sqrt(2) * hH)
    return dict(lx=lx, ly=ly, ls=ls, kap=kap, M=m, h=hH, mu=mu, cH=cH, e=e)


def V_total(x, Y, s, delta, p, cH):
    S = -s * np.exp(1j * delta)
    return (-cH * x * x + p['lx'] * x**4 + p['ly'] * abs(Y)**4 + p['ls'] * (s * s - p['h']**2 / 2)**2
            + p['M']**2 * abs(Y)**2 + 2 * p['kap'] * np.real(np.conj(Y)) * x**3 + 2 * p['mu'] * np.real(np.conj(S) * Y * Y))


def grad_total(q, delta, p, cH):
    x, Yr, Yi, s = q; Y = Yr + 1j * Yi; S = -s * np.exp(1j * delta)
    gx = -2 * cH * x + 4 * p['lx'] * x**3 + 6 * p['kap'] * x * x * Yr
    gY = p['M']**2 * Y + 2 * p['ly'] * abs(Y)**2 * Y + p['kap'] * x**3 + 2 * p['mu'] * S * np.conj(Y)
    gs = 4 * p['ls'] * s * (s * s - p['h']**2 / 2) - 2 * p['mu'] * np.real(np.exp(-1j * delta) * Y * Y)
    return np.array([gx, 2 * gY.real, 2 * gY.imag, gs])


def prepare(p, delta, mode='relaxed'):
    cH = p['cH']
    x_ref = math.sqrt(cH / (2 * p['lx']))

    def solve(d):
        f = lambda q: V_total(math.exp(q[0]), q[1] + 1j * q[2], math.exp(q[3]), d, p, cH)
        best = None
        for dx in (0.0, 0.1, -0.1):
            for y0 in (-0.5, -2.0, -5.0):
                q0 = [math.log(x_ref) + dx, y0, 0.0, math.log(p['h'] / math.sqrt(2))]
                rr = minimize(f, q0, method='Nelder-Mead', options=dict(xatol=1e-12, fatol=1e-16, maxiter=40000))
                rr = minimize(f, rr.x, method='BFGS', options=dict(gtol=1e-12))
                if best is None or rr.fun < best.fun:
                    best = rr
        q = np.array([math.exp(best.x[0]), best.x[1], best.x[2], math.exp(best.x[3])])
        q = fsolve(lambda v: grad_total(v, d, p, cH), q, xtol=1e-15)
        # Hessian positivity (numerical)
        Hm = np.zeros((4, 4))
        for j in range(4):
            e = np.zeros(4); e[j] = 1e-6 * max(abs(q[j]), 1e-6)
            Hm[:, j] = (grad_total(q + e, d, p, cH) - grad_total(q - e, d, p, cH)) / (2 * e[j])
        return q, float(np.linalg.norm(grad_total(q, d, p, cH))), float(np.linalg.eigvalsh(0.5 * (Hm + Hm.T)).min())
    if mode == 'relaxed':
        q, g, hmin = solve(delta)
        info = dict(mode=mode, x=q[0], Y=[q[1], q[2]], s=q[3], grad=g, hess_min=hmin)
        z0 = np.array([q[0], q[1] + 1j * q[2], -q[3] * np.exp(1j * delta)])
    else:   # Codex prescription: aligned (delta = 0) minimum, then rotate Sigma only
        q, g, hmin = solve(0.0)
        info = dict(mode=mode, x=q[0], Y=[q[1], q[2]], s=q[3], grad=g, hess_min=hmin)
        z0 = np.array([q[0], q[1] + 1j * q[2], -q[3] * np.exp(1j * delta)])
    info['x_over_lead'] = info['x'] / x_ref
    return z0, info


def run_full(p, z0, a_end=30.0, rtol=1e-10, atol=1e-13, npts=12000):
    lx, ly, ls, kap, M, h, mu = (p[k] for k in ('lx', 'ly', 'ls', 'kap', 'M', 'h', 'mu'))
    M2 = M * M; h22 = h * h / 2

    def rhs(t, w):
        X = complex(w[0], w[1]); Y = complex(w[2], w[3]); S = complex(w[4], w[5])
        Xd = complex(w[6], w[7]); Yd = complex(w[8], w[9]); Sd = complex(w[10], w[11])
        H3 = 1.5 / t
        ax2 = X.real * X.real + X.imag * X.imag; ay2 = Y.real * Y.real + Y.imag * Y.imag
        as2 = S.real * S.real + S.imag * S.imag
        Xc = X.conjugate()
        FX = 2 * lx * ax2 * X + 3 * kap * Y * Xc * Xc
        FY = M2 * Y + 2 * ly * ay2 * Y + kap * X * X * X + 2 * mu * S * Y.conjugate()
        FS = 2 * ls * (as2 - h22) * S + mu * Y * Y
        aX = -H3 * Xd - FX; aY = -H3 * Yd - FY; aS = -H3 * Sd - FS
        return [Xd.real, Xd.imag, Yd.real, Yd.imag, Sd.real, Sd.imag, aX.real, aX.imag, aY.real, aY.imag, aS.real, aS.imag]
    t0, t1 = 0.5, 0.5 * a_end**2
    w0 = [z0[0].real, z0[0].imag, z0[1].real, z0[1].imag, z0[2].real, z0[2].imag, 0, 0, 0, 0, 0, 0]
    tm = min(8.0, t1)
    t_eval = np.linspace(t0, tm, 4000) if t1 <= 8.0 else np.concatenate([np.linspace(t0, tm, 4000), np.geomspace(tm, t1, npts)[1:]])
    tic = time.time()
    sol = solve_ivp(rhs, (t0, t1), w0, method='DOP853', rtol=rtol, atol=atol, t_eval=t_eval)
    wall = time.time() - tic
    t = sol.t; a = np.sqrt(2 * t); a3 = a**3; a4 = a3 * a
    X = sol.y[0] + 1j * sol.y[1]; Y = sol.y[2] + 1j * sol.y[3]; S = sol.y[4] + 1j * sol.y[5]
    Xd = sol.y[6] + 1j * sol.y[7]; Yd = sol.y[8] + 1j * sol.y[9]; Sd = sol.y[10] + 1j * sol.y[11]
    nX = -2 * a3 * np.imag(np.conj(X) * Xd); nY = -2 * a3 * np.imag(np.conj(Y) * Yd); nS = -2 * a3 * np.imag(np.conj(S) * Sd)
    EX = a4 * (abs(Xd)**2 + lx * abs(X)**4)
    Ncirc = 2 * math.sqrt(2 * lx) * (EX / (3 * lx))**0.75
    eps = nX / Ncirc
    # forced response and its velocity (differentiate A Y + B Y* + J = 0)
    A = M2; B = 2 * mu * S; D = A * A - abs(B)**2
    J = kap * X**3; Jd = 3 * kap * X * X * Xd; Bd = 2 * mu * Sd
    Ys = (-A * J + B * np.conj(J)) / D
    J2 = Jd + Bd * np.conj(Ys)
    Ysd = (-A * J2 + B * np.conj(J2)) / D
    d = Y - Ys; dd = Yd - Ysd
    exc_codex = abs(dd)**2 + A * abs(d)**2 + 2 * np.real(mu * np.conj(S) * d * d)      # Codex's excitation form
    rot = np.exp(-0.5j * np.angle(S)); phid = np.imag(np.conj(S) * Sd) / abs(S)**2
    u = d * rot; ud = (dd - 0.5j * phid * d) * rot
    w1 = np.sqrt(M2 + 2 * mu * abs(S)); w2 = np.sqrt(M2 - 2 * mu * abs(S))
    E1 = ud.real**2 + w1**2 * u.real**2; E2 = ud.imag**2 + w2**2 * u.imag**2
    NYH = a3 * E1 / w1; NYL = a3 * E2 / w2
    sr = abs(S); srd = np.real(np.conj(S) * Sd) / sr
    Esig = srd**2 + ls * (sr**2 - h22)**2; msig = math.sqrt(2 * ls) * h; Nsig = a3 * Esig / msig
    # sources: actual vs EFT with instantaneous Sigma
    alpha6_t = mu * kap**2 * np.conj(S) / D
    src_act = -6 * np.imag(kap * np.conj(Y) * X**3) * a3
    src_eft = -12 * np.imag(alpha6_t * X**6) * a3
    src_D = -12 * np.imag(mu * np.conj(S) * Y * Y) * a3
    bound_integrand = 12 * abs(alpha6_t) * abs(X)**6 * a3
    cum_act = cumulative_trapezoid(src_act, t, initial=0); cum_eft = cumulative_trapezoid(src_eft, t, initial=0)
    cum_abs = cumulative_trapezoid(np.abs(src_act), t, initial=0); cum_bound = cumulative_trapezoid(bound_integrand, t, initial=0)
    aX = a * abs(X)
    late = t > 0.5 * (0.75 * a_end)**2           # window a in [0.75 a_end, a_end]
    def wmean(v):
        return float(np.trapezoid(v[late], a[late]) / (a[late][-1] - a[late][0]))
    res = dict(wall_s=wall, nfev=int(sol.nfev), status=int(sol.status), a_end=a_end,
               eps_end=float(eps[-1]), eps_window=wmean(eps), a3nX_end=float(nX[-1]), a3nX_window=wmean(nX),
               EX_end=float(EX[-1]), EX_window=wmean(EX), NYH_window=wmean(NYH), NYL_window=wmean(NYL),
               NY_window=wmean(NYH + NYL), Nsig_window=wmean(Nsig), exc_codex_a3_over_M_window=wmean(a3 * exc_codex / M),
               gauge_max=float(np.max(np.abs(nX + 3 * nY + 6 * nS)) / np.max(np.abs(nX))),
               nD_end=float((nX + 3 * nY)[-1]), nY_end=float(nY[-1]), nS_end=float(nS[-1]),
               max_aX=float(aX.max()), t_of_max_aX=float(t[aX.argmax()]), aX_end=float(aX[-1]),
               late_envelope_from_EX=float((EX[-1] / lx)**0.25),
               bound_integral_actual_X=float(cum_bound[-1]), integral_abs_source=float(cum_abs[-1]),
               charge_from_actual_source=float(cum_act[-1]), charge_from_eft_source=float(cum_eft[-1]),
               eft_over_actual_charge=float(cum_eft[-1] / cum_act[-1]),
               charge_by_a=[(float(aa), float(c)) for aa, c in zip((1.5, 2, 3, 5, 10, 20), np.interp([0.5 * x**2 for x in (1.5, 2, 3, 5, 10, 20)], t, nX))],
               bound_cum_by_a=[(float(aa), float(c)) for aa, c in zip((1.5, 2, 3, 5, 10, 20), np.interp([0.5 * x**2 for x in (1.5, 2, 3, 5, 10, 20)], t, cum_bound))],
               max_abs_source_ratio_eft_actual_early=float(np.max(np.abs(src_eft[t < 2]) ) / np.max(np.abs(src_act[t < 2]))),
               Ysd_check_initial=float(abs(Ysd[0])), Y_minus_Ystar_initial=float(abs(d[0])), Ystar_initial=float(abs(Ys[0])))
    hist = dict(a=a[::40].tolist(), nX=nX[::40].tolist(), eps=eps[::40].tolist(), aX=aX[::40].tolist(),
                NYH=NYH[::40].tolist(), NYL=NYL[::40].tolist(), src_act=src_act[::40].tolist(), src_eft=src_eft[::40].tolist(),
                cum_bound=cum_bound[::40].tolist())
    return res, hist


def run_eft(p, delta, x0=None, a_end=30.0):
    """One-field matched EFT: lam|X|^4 + beta6|X|^6 + alpha6 X^6 + h.c., fixed Sigma = -(h/sqrt2) e^{i delta}."""
    lx, kap, M, h, mu, cH = p['lx'], p['kap'], p['M'], p['h'], p['mu'], p['cH']
    Sig = -(h / math.sqrt(2)) * np.exp(1j * delta); D = M**4 - 4 * mu**2 * abs(Sig)**2
    a6 = mu * kap**2 * np.conj(Sig) / D; b6 = -M**2 * kap**2 / D
    if x0 is None:   # own Hubble-held radial minimum of -cH x^2 + lx x^4 + b6 x^6 + 2 Re(a6) x^6 along arg X = 0
        c = b6 + 2 * a6.real
        x0 = fsolve(lambda x: -2 * cH * x + 4 * lx * x**3 + 6 * c * x**5, math.sqrt(cH / (2 * lx)))[0]

    def rhs(t, w):
        X = complex(w[0], w[1]); Xd = complex(w[2], w[3])
        ax2 = abs(X)**2
        F = 2 * lx * ax2 * X + 3 * b6 * ax2**2 * X + 6 * np.conj(a6) * X.conjugate()**5
        acc = -1.5 / t * Xd - F
        return [Xd.real, Xd.imag, acc.real, acc.imag]
    t1 = 0.5 * a_end**2
    sol = solve_ivp(rhs, (0.5, t1), [x0, 0, 0, 0], method='DOP853', rtol=1e-11, atol=1e-13, dense_output=True)
    X = sol.y[0, -1] + 1j * sol.y[1, -1]; Xd = sol.y[2, -1] + 1j * sol.y[3, -1]
    a3 = a_end**3; nX = -2 * a3 * np.imag(np.conj(X) * Xd); EX = a3 * a_end * (abs(Xd)**2 + lx * abs(X)**4)
    return dict(x0=float(x0), alpha6=abs(a6), beta6=b6, eps=float(nX / (2 * math.sqrt(2 * lx) * (EX / (3 * lx))**0.75)),
                a3nX=float(nX), EX=float(EX))


def ledger(p, H_GeV, NYH, NYL, Nsig, EX):
    """Closed-form processed ledger of ls_scaling [D] (decay-averaged conformal deposit per quantum)."""
    M, h, mu = (p[k] * H_GeV for k in ('M', 'h', 'mu')); kap, e = p['kap'], p['e']
    mL = math.sqrt(M**2 - math.sqrt(2) * mu * h); mA = 6 * e * h
    GL = 3 * kap**2 * mL / (256 * math.pi**3); GA = e**2 * mA / (48 * math.pi) * (1 - 4 * MXL**2 / mA**2)**1.5
    depL = mL / H_GeV * math.sqrt(math.pi * H_GeV / (2 * GL)); depA = mA / H_GeV * math.sqrt(math.pi * H_GeV / (2 * GA))
    RY = (NYH + NYL) * depL / EX; RS = 2 * Nsig * depA / EX; RC = NYH * depA / EX
    return dict(R_Y=RY, R_Sigma=RS, R_casc=RC, R_central=RY + RS + RC, mL_over_H=mL / H_GeV)


if __name__ == '__main__':
    out = Path(sys.argv[1]); kw = dict(a.split('=') for a in sys.argv[2:])
    m = float(kw.get('m', 30)); a_end = float(kw.get('a_end', 30)); prep = kw.get('prep', 'relaxed')
    frac = float(kw.get('frac', 0.5)); r = float(kw.get('r', 0.45)); hH = float(kw.get('hH', 40))
    delta = math.pi / 2
    p = couplings(m=m, r=r, hH=hH, frac=frac)
    z0, info = prepare(p, delta, prep)
    print('prep', info, flush=True)
    res, hist = run_full(p, z0, a_end=a_end)
    print(json.dumps(res, indent=1), flush=True)
    H = 434867798206.6179 * 5.889032093211876**(-2 / 3)      # the yield-normalised scale (1.3335e11 GeV)
    yld = abs(res['a3nX_window']) * H**3 / entropy(H) / YQ
    led = ledger(p, H, res['NYH_window'], res['NYL_window'], res['Nsig_window'], res['EX_window'])
    # envelope bound with the leading relation and with the actual prepared amplitude
    lam = p['lx']; x_lead = math.sqrt(p['cH'] / (2 * lam)); a6 = r * p['kap']**2 / (2 * m**2 * (1 - r**2))
    bounds = {}
    for tag, x0 in (('lead', x_lead), ('actual', info['x'])):
        R = (x0**4 + x0**2 / lam)**0.25
        bounds[tag] = dict(R=R, a3nX_bound=12 * a6 * R**6, margin=12 * a6 * R**6 / abs(res['a3nX_window']),
                           max_aX_over_R=res['max_aX'] / R, yield_bound_over_YQ=12 * a6 * R**6 * H**3 / entropy(H) / YQ)
    eft = run_eft(p, delta, a_end=a_end)
    eft_same_x0 = run_eft(p, delta, x0=info['x'], a_end=a_end)
    summary = dict(couplings=p, H_GeV=H, prep=info, release=res, yield_over_YQ=yld, ledger=led, alpha6_Hunits=a6,
                   eta=a6 * x_lead**2 / lam, bounds=bounds, eft_one_field=eft, eft_one_field_same_x0=eft_same_x0,
                   eps_over_eta=abs(res['eps_window']) / (a6 * x_lead**2 / lam))
    out.write_text(json.dumps(dict(summary=summary, history=hist), indent=1, default=float))
    print(json.dumps({k: v for k, v in summary.items() if k != 'release'}, indent=1, default=float))
