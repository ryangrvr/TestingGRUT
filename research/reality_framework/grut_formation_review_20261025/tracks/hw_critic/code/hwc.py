#!/usr/bin/env python3
"""hw_critic (Claude, own code): compact independent model of the ls_scaling heavy spot point, used to recompute the
numbers that the hw_verify / hw_envelope / hw_gates / ls_critic tracks dispute or that only one of them computed.

Units: H_i = 1 (the inflationary Hubble rate at release), fields physical. Potential (charges X:Y:Sigma = 1:3:6)
  V = -cX|X|^2 + lx|X|^4 + ly|Y|^4 + ls(|S|^2 - h^2/2)^2 + M^2|Y|^2 + 2Re(kap Y* X^3) + 2Re(mu S* Y^2),
cX = cH during inflation (and cH rho_phi/3 while an inflaton is present), 0 after an instant release.
Field equations z'' + 3H z' + dV/dz* = 0; comoving charges a^3 n_i = -2 a^3 Im(z_i* z_i') (Codex / ls convention).
Release backgrounds: (0) radiation from a = 1, H = 1/a^2; (1) inflaton matter rho_phi(1) = 3 decaying at width G
into radiation, H^2 = (rho_phi + rho_r)/3, with the X Hubble mass tied to rho_phi. Independent variable N = ln a.
Integrator: scipy DOP853 (rtol 1e-10, atol 1e-12), scalar complex arithmetic in the right-hand side.
Diagnostics in the late window: eps = a^3 n_X / N_circ, N_circ = 2 sqrt(2 lx) (E_X/3lx)^{3/4}, E_X = a^4(|X'|^2 + lx|X|^4);
free Y quanta about the forced response Y_* = (-AJ + B J*)/D split into the two Sigma-frame modes (N_YH, N_YL);
Sigma radial quanta N_sigma; Gauss residual. Ledger of ls_scaling / Codex b448f25 [D, re-derived]:
  R_Y = (N_YH + N_YL)(m_L/H)<a>_L/E_X, R_Sigma = 2 N_sigma (m_A/H)<a>_A/E_X, R_casc = N_YH (m_A/H)<a>_A/E_X,
  <a> = sqrt(pi H / 2 Gamma), Gamma_L = 3 kap^2 m_L/256 pi^3, Gamma_A = e^2 m_A (1 - 4 m_X^2/m_A^2)^{3/2}/48 pi.
Yield = |a^3 n_X| H^3 / s(T_R) with prompt radiation, g* = 106.75, M_red = 2.435e18 GeV.
"""
import cmath, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, minimize

MRED, GSTAR, YQ, MXL = 2.435e18, 106.75, 4.8294009076988175e-9, 0.1 * 0.924089637375362
H_CAND = 434867798206.6179                               # GeV, ls_scaling ladder scale [Q, own track]
H_SPOT = H_CAND * 5.889032057603274 ** (-2.0 / 3.0)     # yield-normalised scale (= 1.3335e11 GeV)
QCH = np.array([1.0, 3.0, 6.0])
# our sp_critic case-L table [Q, own track]; log-log interpolated (runs 7-9 % high, see CONSOLIDATED)
RL_TABLE = [(0.0677, 29.756595100551028), (0.1, 49.83550533984897), (0.15, 91.80895700634682),
            (0.19, 135.89463084582223), (0.240201575085664, 172.49014985200023)]


def spot(m=30.0, r=0.45, hH=40.0, frac=0.5, cH=10.0, lx=1e-3, ly=1.0, ls=1.0, e=1e-7):
    """ls_scaling parametrisation: eta = alpha6 X0^2/lx = frac * r/12 (companion cap), r = sqrt2 mu h / M^2."""
    kap = m * lx * math.sqrt(frac * 4.0 * (1.0 - r * r) / (12.0 * cH))
    mu = r * m * m / (math.sqrt(2.0) * hH)
    return dict(lx=lx, ly=ly, ls=ls, kap=kap, M=m, h=hH, mu=mu, cH=cH, e=e, r=r)


# ---------------------------------------------------------------- potential and complex gradient
def V(z, p, cX):
    X, Y, S = z
    aX, aY, aS = abs(X) ** 2, abs(Y) ** 2, abs(S) ** 2
    return (-cX * aX + p['lx'] * aX ** 2 + p['ly'] * aY ** 2 + p['ls'] * (aS - p['h'] ** 2 / 2) ** 2
            + p['M'] ** 2 * aY + 2 * (p['kap'] * Y.conjugate() * X ** 3 + p['mu'] * S.conjugate() * Y * Y).real)


def F(z, p, cX):
    """dV/dz_i^* (the force in z'' + 3Hz' + F = 0)."""
    X, Y, S = z
    Xc = X.conjugate()
    return (-cX * X + 2 * p['lx'] * (X * Xc).real * X + 3 * p['kap'] * Y * Xc * Xc,
            p['M'] ** 2 * Y + 2 * p['ly'] * (Y * Y.conjugate()).real * Y + p['kap'] * X ** 3 + 2 * p['mu'] * S * Y.conjugate(),
            2 * p['ls'] * ((S * S.conjugate()).real - p['h'] ** 2 / 2) * S + p['mu'] * Y * Y)


# ---------------------------------------------------------------- preparations (own minimisers)
def aligned(p, cX=None):
    """Aligned inflationary minimum X = x, Y = -y, S = -s (x, y, s > 0)."""
    cX = p['cH'] if cX is None else cX
    lx, ly, ls, kap, M, h, mu = (p[k] for k in ('lx', 'ly', 'ls', 'kap', 'M', 'h', 'mu'))

    def g(q):
        x, y, s = q
        return [-2 * cX * x + 4 * lx * x ** 3 - 6 * kap * y * x * x,
                4 * ly * y ** 3 + 2 * M * M * y - 2 * kap * x ** 3 - 4 * mu * s * y,
                4 * ls * s * (s * s - h * h / 2) - 2 * mu * y * y]
    best = None
    for x0 in (0.9, 1.0, 1.15):
        for y0 in (1.0, 3.0, 6.0):
            q, info, ier, _ = fsolve(g, [x0 * math.sqrt(cX / (2 * lx)), y0, h / math.sqrt(2)], xtol=1e-14, full_output=True)
            if ier == 1 and min(q) > 0:
                val = V((q[0], -q[1], -q[2]), p, cX)
                if best is None or val < best[0]:
                    best = (val, q)
    x, y, s = best[1]
    return np.array([x, -y, -s], complex), dict(x=x, y=y, s=s, V=best[0])


def prep(kind, delta, p, cX=None):
    """'codex': aligned minimum with Sigma rotated by delta (Y left at its delta = 0 value);
    'relaxed': conditional minimum over (|X|, Re Y, Im Y, |S|) at arg X = 0, S = -|S| e^{i delta} (ls_scaling)."""
    cX = p['cH'] if cX is None else cX
    z0, info = aligned(p, cX)
    ph = cmath.exp(1j * delta)
    if kind == 'codex':
        return np.array([z0[0], z0[1], -abs(z0[2]) * ph]), info
    if kind != 'relaxed':
        raise ValueError(kind)

    def resid(q):
        x, yr, yi, s = q
        z = (complex(x, 0), complex(yr, yi), -s * ph)
        fX, fY, fS = F(z, p, cX)
        return [2 * fX.real, 2 * fY.real, 2 * fY.imag, -2 * (fS / ph).real]   # dV/dx, dV/dYr, dV/dYi, dV/ds

    def hess4(q):
        Hm = np.zeros((4, 4))
        for j in range(4):
            e = np.zeros(4); e[j] = 1e-5 * max(1.0, abs(q[j]))
            Hm[:, j] = (np.array(resid(q + e)) - np.array(resid(q - e))) / (2 * e[j])
        return 0.5 * (Hm + Hm.T)
    best = None
    y0 = abs(z0[1])
    for ang in (0.0, delta / 2, -delta / 2, delta, math.pi / 2):
        for fx in (1.0, 0.95):
            yc = -y0 * cmath.exp(1j * ang)
            q0 = [fx * abs(z0[0]), yc.real, yc.imag, abs(z0[2])]
            q, info2, ier, _ = fsolve(resid, q0, xtol=1e-14, full_output=True)
            if ier != 1 or q[0] <= 0 or q[3] <= 0:
                continue
            hm = np.linalg.eigvalsh(hess4(q)).min()
            val = V((complex(q[0], 0), complex(q[1], q[2]), -q[3] * ph), p, cX)
            if hm > 0 and (best is None or val < best[0]):
                best = (val, q, hm)
    if best is None:
        raise RuntimeError('no relaxed conditional minimum at delta = %g' % delta)
    val, q, hm = best
    z = np.array([complex(q[0], 0), complex(q[1], q[2]), -q[3] * ph])
    return z, dict(x=q[0], Y=[q[1], q[2]], s=q[3], V=val, hess_min=hm, aligned=info)


# ---------------------------------------------------------------- canonical 6x6 Hessian, gauge mode projected
def canon(z):
    return math.sqrt(2) * np.ravel([[c.real, c.imag] for c in z])


def grad_canon(q, p, cX):
    z = (complex(q[0], q[1]) / math.sqrt(2), complex(q[2], q[3]) / math.sqrt(2), complex(q[4], q[5]) / math.sqrt(2))
    f = F(z, p, cX)
    return math.sqrt(2) * np.array([f[0].real, f[0].imag, f[1].real, f[1].imag, f[2].real, f[2].imag])


def hessian(z, p, cX, step=1e-4):
    q = canon(z)
    Hm = np.zeros((6, 6))
    for j in range(6):
        e = np.zeros(6); e[j] = step * max(1.0, abs(q[j]))
        Hm[:, j] = (grad_canon(q + e, p, cX) - grad_canon(q - e, p, cX)) / (2 * e[j])
    return 0.5 * (Hm + Hm.T)


def spectrum(z, p, cX):
    Hm = hessian(z, p, cX)
    g = canon(1j * QCH * z); g /= np.linalg.norm(g)
    P = np.eye(6) - np.outer(g, g)
    w, U = np.linalg.eigh(P @ Hm @ P)
    gi = int(np.argmax(np.abs(U.T @ g)))
    modes = []
    for k in range(6):
        if k == gi:
            continue
        v = U[:, k]
        dz = (v[0::2] + 1j * v[1::2]) / math.sqrt(2)
        dth = np.imag(dz / z); dlr = np.real(dz / z)
        wa = float(np.sum((np.abs(z) * dth) ** 2)); wr = float(np.sum((np.abs(z) * dlr) ** 2))
        modes.append(dict(m2=float(w[k]), ddelta=float(dth[2] - 6 * dth[0]), angular_weight=wa / (wa + wr),
                          dtheta=[float(t) for t in dth]))
    ang = [m for m in modes if m['angular_weight'] > 0.5]
    light = min(ang, key=lambda m: m['m2']) if ang else None
    out = dict(m2=[m['m2'] for m in modes], gauge_eig=float(w[gi]), modes=modes)
    if light:
        out['light_phase'] = dict(m2=light['m2'], mass_over_H=math.sqrt(max(light['m2'], 0.0)),
                                  f_eff=1.0 / abs(light['ddelta']), dtheta=light['dtheta'])
    return out


# ---------------------------------------------------------------- de Sitter stage (H = 1, inflationary potential)
def desitter(z0, N, p, held_delta=None, v0=None):
    """Evolve N e-folds at H = 1. held_delta = None: all six real components free. held_delta = delta: arg X = 0 and
    S = -s e^{i delta} are held (x, complex Y, s evolve). Returns (z, z') at the end."""
    cX = p['cH']
    if N <= 0:
        return np.array(z0, complex), (np.zeros(3, complex) if v0 is None else np.array(v0, complex))
    if held_delta is None:
        def rhs(t, u):
            z = (complex(u[0], u[1]), complex(u[2], u[3]), complex(u[4], u[5]))
            v = (complex(u[6], u[7]), complex(u[8], u[9]), complex(u[10], u[11]))
            f = F(z, p, cX)
            acc = [-3 * v[i] - f[i] for i in range(3)]
            return [v[0].real, v[0].imag, v[1].real, v[1].imag, v[2].real, v[2].imag,
                    acc[0].real, acc[0].imag, acc[1].real, acc[1].imag, acc[2].real, acc[2].imag]
        v0 = np.zeros(3, complex) if v0 is None else np.asarray(v0, complex)
        u0 = [c for zz in z0 for c in (zz.real, zz.imag)] + [c for vv in v0 for c in (vv.real, vv.imag)]
        sol = solve_ivp(rhs, (0.0, N), u0, method='DOP853', rtol=1e-12, atol=1e-14)
        u = sol.y[:, -1]
        return (np.array([complex(u[0], u[1]), complex(u[2], u[3]), complex(u[4], u[5])]),
                np.array([complex(u[6], u[7]), complex(u[8], u[9]), complex(u[10], u[11])]))
    ph = cmath.exp(1j * held_delta)

    def rhs(t, u):
        x, yr, yi, s, vx, vyr, vyi, vs = u
        fX, fY, fS = F((complex(x, 0), complex(yr, yi), -s * ph), p, cX)
        return [vx, vyr, vyi, vs, -3 * vx - fX.real, -3 * vyr - fY.real, -3 * vyi - fY.imag, -3 * vs + (fS / ph).real]
    u0 = [z0[0].real, z0[1].real, z0[1].imag, abs(z0[2]), 0, 0, 0, 0]
    if v0 is not None:
        u0[4:] = [v0[0].real, v0[1].real, v0[1].imag, (-v0[2] / ph).real]
    sol = solve_ivp(rhs, (0.0, N), u0, method='DOP853', rtol=1e-12, atol=1e-14)
    u = sol.y[:, -1]
    return (np.array([complex(u[0], 0), complex(u[1], u[2]), -u[3] * ph]),
            np.array([complex(u[4], 0), complex(u[5], u[6]), -u[7] * ph]))


# ---------------------------------------------------------------- release
def release(z0, p, v0=None, bg=0, gam=1.0, a_end=60.0, win=(30.0, 60.0), nwin=121, rtol=1e-10, atol=1e-12):
    lx, ly, ls, kap, M, h, mu, cH = (p[k] for k in ('lx', 'ly', 'ls', 'kap', 'M', 'h', 'mu', 'cH'))
    v0 = np.zeros(3, complex) if v0 is None else np.asarray(v0, complex)

    def rhs(N, u):
        z = (complex(u[0], u[1]), complex(u[2], u[3]), complex(u[4], u[5]))
        v = (complex(u[6], u[7]), complex(u[8], u[9]), complex(u[10], u[11]))
        rphi, rr = u[12], u[13]
        Hh = math.sqrt((rphi + rr) / 3.0)
        cX = cH * rphi / 3.0 if bg == 1 else 0.0
        f = F(z, p, cX)
        acc = [(-3 * Hh * v[i] - f[i]) / Hh for i in range(3)]
        vel = [v[i] / Hh for i in range(3)]
        dphi = -3 * rphi - (gam * rphi / Hh if bg == 1 else 0.0)
        drr = -4 * rr + (gam * rphi / Hh if bg == 1 else 0.0)
        return [vel[0].real, vel[0].imag, vel[1].real, vel[1].imag, vel[2].real, vel[2].imag,
                acc[0].real, acc[0].imag, acc[1].real, acc[1].imag, acc[2].real, acc[2].imag, dphi, drr, 1.0 / Hh]
    u0 = [c for zz in z0 for c in (zz.real, zz.imag)] + [c for vv in v0 for c in (vv.real, vv.imag)]
    u0 += [3.0, 0.0, 0.0] if bg == 1 else [0.0, 3.0, 0.0]
    Nw = np.linspace(math.log(win[0]), math.log(win[1]), nwin)
    Nearly = np.linspace(0.0, math.log(win[0]), 400, endpoint=False)
    N_eval = np.concatenate([Nearly, Nw])
    sol = solve_ivp(rhs, (0.0, math.log(a_end)), u0, method='DOP853', rtol=rtol, atol=atol, t_eval=N_eval)
    if sol.status != 0:
        raise RuntimeError(sol.message)
    a = np.exp(sol.t); a3 = a ** 3
    X = sol.y[0] + 1j * sol.y[1]; Y = sol.y[2] + 1j * sol.y[3]; S = sol.y[4] + 1j * sol.y[5]
    Xd = sol.y[6] + 1j * sol.y[7]; Yd = sol.y[8] + 1j * sol.y[9]; Sd = sol.y[10] + 1j * sol.y[11]
    rphi, rr = sol.y[12], sol.y[13]
    nX = -2 * a3 * np.imag(np.conj(X) * Xd); nY = -2 * a3 * np.imag(np.conj(Y) * Yd); nS = -2 * a3 * np.imag(np.conj(S) * Sd)
    EX = a3 * a * (np.abs(Xd) ** 2 + lx * np.abs(X) ** 4)
    eps = nX / (2 * math.sqrt(2 * lx) * (EX / (3 * lx)) ** 0.75)
    # free Y quanta about the forced response
    A = M * M; B = 2 * mu * S; D = A * A - np.abs(B) ** 2
    J = kap * X ** 3; Jd = 3 * kap * X * X * Xd; Bd = 2 * mu * Sd
    Ys = (-A * J + B * np.conj(J)) / D
    J2 = Jd + Bd * np.conj(Ys); Ysd = (-A * J2 + B * np.conj(J2)) / D
    d = Y - Ys; dd = Yd - Ysd
    rot = np.exp(-0.5j * np.angle(S)); phid = np.imag(np.conj(S) * Sd) / np.abs(S) ** 2
    u = d * rot; ud = (dd - 0.5j * phid * d) * rot
    w1 = np.sqrt(A + 2 * mu * np.abs(S)); w2 = np.sqrt(A - 2 * mu * np.abs(S))
    NYH = a3 * (ud.real ** 2 + w1 ** 2 * u.real ** 2) / w1
    NYL = a3 * (ud.imag ** 2 + w2 ** 2 * u.imag ** 2) / w2
    sr = np.abs(S); srd = np.real(np.conj(S) * Sd) / sr
    Nsig = a3 * (srd ** 2 + ls * (sr ** 2 - h * h / 2) ** 2) / (math.sqrt(2 * ls) * h)
    gauss = np.max(np.abs(nX + 3 * nY + 6 * nS)) / np.max(np.abs(nX))
    W = slice(len(Nearly), None)

    def wmean(v):
        return float(np.trapezoid(v[W], a[W]) / (a[W][-1] - a[W][0]))
    out = dict(eps=wmean(eps), NX=wmean(nX), EX=wmean(EX), NYH=wmean(NYH), NYL=wmean(NYL), Nsig=wmean(Nsig),
               eps_end=float(eps[-1]), NX_end=float(nX[-1]), gauss=float(gauss), nfev=int(sol.nfev),
               Rr_end=float(rr[-1] * a[-1] ** 4), rhophi_frac_end=float(rphi[-1] / (rphi[-1] + rr[-1])),
               max_aX=float(np.max(a * np.abs(X))), window=list(win), a_end=a_end, bg=bg, gam=gam if bg == 1 else None)
    out['hist'] = dict(a=a[::20].tolist(), eps=eps[::20].tolist(), nX=nX[::20].tolist(), NYH=NYH[::20].tolist(),
                       NYL=NYL[::20].tolist(), absX=np.abs(X)[::20].tolist())
    return out


# ---------------------------------------------------------------- ledger, yield, threshold table
def masses_GeV(p, H):
    M, h, mu = p['M'] * H, p['h'] * H, p['mu'] * H
    return dict(mYH=math.sqrt(M * M + math.sqrt(2) * mu * h), mYL=math.sqrt(M * M - math.sqrt(2) * mu * h),
                msig=math.sqrt(2 * p['ls']) * h, mA=6 * p['e'] * h)


def ledger(p, w, H, e=None):
    pp = dict(p); pp['e'] = p['e'] if e is None else e
    m = masses_GeV(pp, H); kap, ee = pp['kap'], pp['e']
    GL = 3 * kap ** 2 * m['mYL'] / (256 * math.pi ** 3)
    GA = ee ** 2 * m['mA'] / (48 * math.pi) * (1 - 4 * MXL ** 2 / m['mA'] ** 2) ** 1.5
    depL = m['mYL'] / H * math.sqrt(math.pi * H / (2 * GL)); depA = m['mA'] / H * math.sqrt(math.pi * H / (2 * GA))
    RY = (w['NYH'] + w['NYL']) * depL / w['EX']; RS = 2 * w['Nsig'] * depA / w['EX']; RC = w['NYH'] * depA / w['EX']
    return dict(R_Y=RY, R_Sigma=RS, R_casc=RC, R_central=RY + RS + RC, dep_A_in_EX=depA / w['EX'],
                dep_L_in_EX=depL / w['EX'], Gamma_A_GeV=GA, Gamma_L_GeV=GL, a_dec_A=math.sqrt(math.pi * H / (2 * GA)))


def entropy_RH(H):
    T = (90 * MRED ** 2 * H ** 2 / (math.pi ** 2 * GSTAR)) ** 0.25
    return 2 * math.pi ** 2 / 45 * GSTAR * T ** 3, T


def yield_over_YQ(NX, H):
    return abs(NX) * H ** 3 / entropy_RH(H)[0] / YQ


def alpha6_c6(p, H):
    M, h, mu, kap = p['M'] * H, p['h'] * H, p['mu'] * H, p['kap']
    Sig = h / math.sqrt(2); D = M ** 4 - 4 * (mu * Sig) ** 2
    a6 = mu * kap ** 2 * Sig / D
    return dict(alpha6_GeVm2=a6, c6=a6 * MRED ** 2, alpha6_Hm2=a6 * H * H, eta=a6 * H * H * (p['cH'] / (2 * p['lx'])) / p['lx'])


def R_L_table(eps):
    e = abs(eps)
    if e < RL_TABLE[0][0] or e > RL_TABLE[-1][0]:
        return None
    xs = [math.log(a) for a, _ in RL_TABLE]; ys = [math.log(b) for _, b in RL_TABLE]
    return float(math.exp(np.interp(math.log(e), xs, ys)))


def summarize(p, w, H, tag='', **extra):
    L = ledger(p, w, H)
    o = dict(tag=tag, eps=w['eps'], NX=w['NX'], EX=w['EX'], NYH=w['NYH'], NYL=w['NYL'], Nsig=w['Nsig'], NY=w['NYH'] + w['NYL'],
             gauss=w['gauss'], nfev=w['nfev'], eps_end=w['eps_end'], max_aX=w['max_aX'], **L)
    o['yield_over_YQ'] = yield_over_YQ(w['NX'], H)
    o['R_L_table'] = R_L_table(w['eps'])
    o['ratio_table'] = (o['R_central'] / o['R_L_table']) if o['R_L_table'] else None
    if w['bg'] == 1:
        dil = (w['Rr_end'] / 3.0) ** 0.75
        o.update(gamma_over_H=w['gam'], Rr_end=w['Rr_end'], entropy_dilution=dil, yield_diluted_over_YQ=o['yield_over_YQ'] / dil,
                 R_central_reheat=o['R_central'] * (w['Rr_end'] / 3.0) ** 0.25, rhophi_frac_end=w['rhophi_frac_end'])
    o.update(extra)
    return o
