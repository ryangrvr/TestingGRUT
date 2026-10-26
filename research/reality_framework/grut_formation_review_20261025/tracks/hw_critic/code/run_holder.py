#!/usr/bin/env python3
"""hw_critic (Claude, own code): independent check of hw_gates' D1 result (a Hubble-induced A-term holder that enters
the Y equation, switched off adiabatically, leaves little free Y).
D1 during inflation: mu -> mu e^{i delta} (so the Codex state, aligned moduli with Sigma rotated by delta, is the exact
inflationary minimum). After release the extra piece mu_P = mu (e^{i delta} - 1) decays as
  'exp':  mu_P(a) = mu_P exp(-(a - 1)/Dl);   'H': mu_P(a) = mu_P / a^2 (A-term tracking H in radiation);
  'instant': mu_P = 0 for a > 1 (= Codex preparation).
hw_gates (window 30-60): instant eps -0.159, R 1415; exp Dl = 0.1: eps -0.117, R 161; Dl = 0.3: eps -0.092, R 11.1;
Dl = 1: eps -0.065, R 2.4; H-tracking: eps -0.075, R 4.6.
Usage: python3 -B run_holder.py ../results/holder.json
"""
import cmath, json, math, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hwc

P = hwc.spot(); H = hwc.H_SPOT; D0 = math.pi / 2
MUP = P['mu'] * (cmath.exp(1j * D0) - 1.0)


def release_holder(z0, mode, Dl, a_end=60.0, win=(30.0, 60.0), nwin=121):
    """Radiation release with the time-dependent holder term V_P = 2 Re(mu_P(a) S* Y^2)."""
    p = P; lx, ly, ls, kap, M, h, mu = (p[k] for k in ('lx', 'ly', 'ls', 'kap', 'M', 'h', 'mu'))

    def muP(a):
        if mode == 'exp':
            return MUP * math.exp(-(a - 1.0) / Dl)
        if mode == 'H':
            return MUP / (a * a)
        return 0.0

    def rhs(N, u):
        a = math.exp(N); Hh = 1.0 / (a * a)
        X = complex(u[0], u[1]); Y = complex(u[2], u[3]); S = complex(u[4], u[5])
        vX = complex(u[6], u[7]); vY = complex(u[8], u[9]); vS = complex(u[10], u[11])
        fX, fY, fS = hwc.F((X, Y, S), p, 0.0)
        m = muP(a)
        fY = fY + 2 * m.conjugate() * S * Y.conjugate()
        fS = fS + m * Y * Y
        aX = (-3 * Hh * vX - fX) / Hh; aY = (-3 * Hh * vY - fY) / Hh; aS = (-3 * Hh * vS - fS) / Hh
        return [vX.real / Hh, vX.imag / Hh, vY.real / Hh, vY.imag / Hh, vS.real / Hh, vS.imag / Hh,
                aX.real, aX.imag, aY.real, aY.imag, aS.real, aS.imag]
    u0 = [c for zz in z0 for c in (zz.real, zz.imag)] + [0.0] * 6
    Nw = np.linspace(math.log(win[0]), math.log(win[1]), nwin)
    N_eval = np.concatenate([np.linspace(0.0, math.log(win[0]), 400, endpoint=False), Nw])
    sol = solve_ivp(rhs, (0.0, math.log(a_end)), u0, method='DOP853', rtol=1e-10, atol=1e-12, t_eval=N_eval)
    a = np.exp(sol.t); a3 = a ** 3
    X = sol.y[0] + 1j * sol.y[1]; Y = sol.y[2] + 1j * sol.y[3]; S = sol.y[4] + 1j * sol.y[5]
    Xd = sol.y[6] + 1j * sol.y[7]; Yd = sol.y[8] + 1j * sol.y[9]; Sd = sol.y[10] + 1j * sol.y[11]
    nX = -2 * a3 * np.imag(np.conj(X) * Xd); nY = -2 * a3 * np.imag(np.conj(Y) * Yd); nS = -2 * a3 * np.imag(np.conj(S) * Sd)
    EX = a3 * a * (np.abs(Xd) ** 2 + lx * np.abs(X) ** 4)
    eps = nX / (2 * math.sqrt(2 * lx) * (EX / (3 * lx)) ** 0.75)
    A = M * M; B = 2 * mu * S; D = A * A - np.abs(B) ** 2        # late window: holder gone (checked below)
    J = kap * X ** 3; Jd = 3 * kap * X * X * Xd; Bd = 2 * mu * Sd
    Ys = (-A * J + B * np.conj(J)) / D; J2 = Jd + Bd * np.conj(Ys); Ysd = (-A * J2 + B * np.conj(J2)) / D
    d = Y - Ys; dd = Yd - Ysd
    rot = np.exp(-0.5j * np.angle(S)); phid = np.imag(np.conj(S) * Sd) / np.abs(S) ** 2
    uu = d * rot; ud = (dd - 0.5j * phid * d) * rot
    w1 = np.sqrt(A + 2 * mu * np.abs(S)); w2 = np.sqrt(A - 2 * mu * np.abs(S))
    NYH = a3 * (ud.real ** 2 + w1 ** 2 * uu.real ** 2) / w1; NYL = a3 * (ud.imag ** 2 + w2 ** 2 * uu.imag ** 2) / w2
    sr = np.abs(S); srd = np.real(np.conj(S) * Sd) / sr
    Nsig = a3 * (srd ** 2 + ls * (sr ** 2 - h * h / 2) ** 2) / (math.sqrt(2 * ls) * h)
    W = slice(400, None)

    def wmean(v):
        return float(np.trapezoid(v[W], a[W]) / (a[W][-1] - a[W][0]))
    return dict(eps=wmean(eps), NX=wmean(nX), EX=wmean(EX), NYH=wmean(NYH), NYL=wmean(NYL), Nsig=wmean(Nsig),
                eps_end=float(eps[-1]), NX_end=float(nX[-1]), gauss=float(np.max(np.abs(nX + 3 * nY + 6 * nS)) / np.max(np.abs(nX))),
                nfev=int(sol.nfev), Rr_end=3.0, bg=0, gam=None, max_aX=float(np.max(a * np.abs(X))),
                holder_at_a30=abs(muP(30.0)) / abs(MUP), window=list(win), a_end=a_end)


def job(spec):
    mode, Dl = spec
    t0 = time.time()
    z0, info = hwc.prep('codex', D0, P)
    w = release_holder(z0, mode, Dl)
    o = hwc.summarize(P, w, H, tag='D1_%s_%s' % (mode, Dl), mode=mode, Dl=Dl, seconds=time.time() - t0,
                      holder_at_a30=w['holder_at_a30'])
    print('D1 %-8s Dl %-5s eps %+.5f NX %9.2f NYH %9.4g NYL %9.4g Nsig %8.2g R %9.4g RL %s Y/YQ %.3f %4.0fs'
          % (mode, Dl, o['eps'], o['NX'], o['NYH'], o['NYL'], o['Nsig'], o['R_central'],
             ('%.1f' % o['R_L_table']) if o['R_L_table'] else 'n/a', o['yield_over_YQ'], o['seconds']), flush=True)
    return o


if __name__ == '__main__':
    out = Path(sys.argv[1])
    jobs = [('instant', None), ('exp', 0.1), ('exp', 0.3), ('exp', 1.0), ('H', None)]
    rows = []
    with Pool(4) as pool:
        for o in pool.imap_unordered(job, jobs):
            rows.append(o)
            out.write_text(json.dumps(dict(method=__doc__, H_GeV=H, couplings=P, muP=[MUP.real, MUP.imag], rows=rows), indent=1, default=float))
    print('done')
