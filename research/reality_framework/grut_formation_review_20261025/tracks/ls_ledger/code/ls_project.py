#!/usr/bin/env python3
"""ls_ledger (Claude, own code): late-state projection and Sigma radial decomposition in PHYSICAL fields, cosmic time.

Separate from the C diagnostics in sp_yledger/code/yrelease.c (used for the releases) and from Codex's
modes.py: the forced Y response is solved with the complex inverse, the Y modes are taken in the frame
co-rotating with arg(Sigma)/2, and the Sigma radial count is examined three ways:
  (1) Codex/sp_yledger formula   N_sig = a^3 [rdot^2 + lambda_S (r^2 - h^2/2)^2] / m_sig;
  (2) a continuation of the undamped release over ~20 Sigma periods (own DOP853 call via scipy), with a
      least-squares split of delta r = r - h/sqrt2 into a static offset, the free m_sig oscillation and the
      Y-forced lines (2 w_H, 2 w_L, w_H + w_L, w_H - w_L);
  (3) the predicted static (forced) offsets: centrifugal r thetadot^2/m^2 and Y-induced -mu Re(Y^2 e^{-i theta})/m^2.
Units: H_i = 1; radiation a = sqrt(2 t), H = 1/(2 t); comoving numbers in H_i^3 units (as Codex).
"""
import numpy as np
from scipy.integrate import solve_ivp

MY, HV, MU = 1e12, 1e11, 5e11
LX, LY, LS, KAP = 1e-3, 1.0, 1.0, 1e-3


def pars(H):
    return MY / H, HV / H, MU / H


def force(z, H):
    M, h, mu = pars(H)
    X, Y, S = z
    return np.array([2 * LX * abs(X)**2 * X + 3 * KAP * Y * np.conj(X)**2,
                     M * M * Y + 2 * LY * abs(Y)**2 * Y + KAP * X**3 + 2 * mu * S * np.conj(Y),
                     2 * LS * (abs(S)**2 - h * h / 2) * S + mu * Y * Y])


def project(t, z, zd, H):
    """Snapshot projection at cosmic time t (a = sqrt(2t)) of physical fields z and dz/dt."""
    M, h, mu = pars(H)
    a = np.sqrt(2 * t); a3 = a**3
    X, Y, S = z; Xd, Yd, Sd = zd
    NX = -2 * a3 * np.imag(np.conj(X) * Xd)
    NY = -2 * a3 * np.imag(np.conj(Y) * Yd)
    NS = -2 * a3 * np.imag(np.conj(S) * Sd)
    EX = a**4 * (abs(Xd)**2 + LX * abs(X)**4)
    Ncirc = 2 * np.sqrt(2 * LX) * (EX / (3 * LX))**0.75
    # forced response: M^2 y + 2 mu S y* = -J (complex inverse), and its exact time derivative
    M2 = M * M; D = M2 * M2 - 4 * mu * mu * abs(S)**2
    J = KAP * X**3
    yp = (-M2 * J + 2 * mu * S * np.conj(J)) / D
    Jd = 3 * KAP * X**2 * Xd + 2 * mu * Sd * np.conj(yp)
    ypd = (-M2 * Jd + 2 * mu * S * np.conj(Jd)) / D
    d, dd = Y - yp, Yd - ypd
    r = abs(S); th = np.angle(S)
    rd = np.real(np.conj(S) * Sd) / r
    thd = np.imag(np.conj(S) * Sd) / r**2
    rot = np.exp(-0.5j * th)
    u = d * rot; ud = (dd - 0.5j * thd * d) * rot
    wH = np.sqrt(M2 + 2 * mu * r); wL = np.sqrt(M2 - 2 * mu * r)
    NYH = a3 * (ud.real**2 + wH**2 * u.real**2) / wH
    NYL = a3 * (ud.imag**2 + wL**2 * u.imag**2) / wL
    msig = np.sqrt(2 * LS) * h; r0 = h / np.sqrt(2)
    Erad = rd**2 + LS * (r * r - h * h / 2)**2
    Nsig = a3 * Erad / msig
    Nsig_quad = a3 * (rd**2 + msig**2 * (r - r0)**2) / msig
    dr_centrifugal = r * thd**2 / msig**2
    dr_Yforced = -mu * np.real(Y * Y * np.exp(-1j * th)) / msig**2
    return dict(a=a, NX=NX, NY=NY, NS=NS, Gauss=NX + 3 * NY + 6 * NS, EX=EX, Ncirc=Ncirc, eps=NX / Ncirc,
                NYH=NYH, NYL=NYL, Nsig=Nsig, Nsig_quadratic=Nsig_quad,
                ESig_angular_a4=a**4 * r * r * thd**2, theta_dot=thd,
                dr_centrifugal=dr_centrifugal, dr_Yforced=dr_Yforced,
                Nsig_from_static_centrifugal=a3 * msig * dr_centrifugal**2,
                yp_over_Y=abs(yp) / max(abs(Y), 1e-300))


def rhs(t, v, H):
    z = v[0:6:2] + 1j * v[1:6:2]; zd = v[6:12:2] + 1j * v[7:12:2]
    F = force(z, H); acc = -3 / (2 * t) * zd - F
    out = np.empty(12)
    out[0:6:2], out[1:6:2] = zd.real, zd.imag
    out[6:12:2], out[7:12:2] = acc.real, acc.imag
    return out


def sigma_decomposition(t0, z, zd, H, periods=20, nper=64, rtol=1e-11, atol=1e-15):
    """Continue the undamped release from (t0, z, zd) for `periods` Sigma periods; split delta r."""
    M, h, mu = pars(H)
    msig = np.sqrt(2 * LS) * h; r0 = h / np.sqrt(2)
    T = 2 * np.pi / msig * periods
    wH = np.sqrt(M * M + np.sqrt(2) * mu * h); wL = np.sqrt(M * M - np.sqrt(2) * mu * h)
    # sample finely enough for the fastest forced line 2 wH
    n = int(max(periods * nper, 40 * T * 2 * wH / (2 * np.pi)))
    ts = np.linspace(t0, t0 + T, n)
    v0 = np.concatenate([np.ravel(np.column_stack([z.real, z.imag])), np.ravel(np.column_stack([zd.real, zd.imag]))])
    sol = solve_ivp(rhs, [t0, t0 + T], v0, method='DOP853', t_eval=ts, rtol=rtol, atol=atol, args=(H,))
    assert sol.success
    S = sol.y[4] + 1j * sol.y[5]; Sd = sol.y[10] + 1j * sol.y[11]
    Y = sol.y[2] + 1j * sol.y[3]
    r = abs(S); dr = r - r0
    rd = np.real(np.conj(S) * Sd) / r
    thd = np.imag(np.conj(S) * Sd) / r**2
    a = np.sqrt(2 * ts); a3 = a**3
    Nsig_formula = a3 * (rd**2 + LS * (r * r - h * h / 2)**2) / msig
    # least squares: delta r = c0 + sum_k [A_k cos w_k t + B_k sin w_k t]  (frequencies fixed; Delta a/a ~ 4e-4)
    freqs = dict(free=msig, twoH=2 * wH, twoL=2 * wL, HplusL=wH + wL, HminusL=wH - wL)
    cols = [np.ones_like(ts)]
    for w in freqs.values():
        cols += [np.cos(w * (ts - t0)), np.sin(w * (ts - t0))]
    A = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(A, dr, rcond=None)
    resid = dr - A @ coef
    am = np.mean(a3)
    out = dict(t0=t0, T=T, a_range=[float(a[0]), float(a[-1])], samples=n,
               Nsig_formula_mean=float(np.mean(Nsig_formula)),
               Nsig_formula_relative_ripple=float(np.std(Nsig_formula) / np.mean(Nsig_formula)),
               dr_mean=float(np.mean(dr)), dr_static_fit=float(coef[0]),
               dr_centrifugal_pred=float(np.mean(r * thd**2) / msig**2),
               dr_Yforced_pred=float(-mu * np.mean(np.real(Y * Y * np.exp(-1j * np.angle(S)))) / msig**2),
               N_static=float(am * msig * coef[0]**2),
               fit_residual_rms_over_dr_rms=float(np.std(resid) / np.std(dr)))
    lines = {}
    for i, (k, w) in enumerate(freqs.items()):
        amp2 = coef[1 + 2 * i]**2 + coef[2 + 2 * i]**2
        # oscillator energy of a line at frequency w: E = rdot^2 + m^2 dr^2 averaged -> (w^2 + m^2) amp^2 / 2;
        # counted with the Codex formula (divided by m_sig).
        lines[k] = dict(omega=float(w), amplitude=float(np.sqrt(amp2)),
                        N_counted_by_formula=float(am * (w * w + msig * msig) * amp2 / 2 / msig))
    out['lines'] = lines
    out['N_free'] = lines['free']['N_counted_by_formula']
    out['N_forced_lines_total'] = float(sum(v['N_counted_by_formula'] for k, v in lines.items() if k != 'free'))
    return out
