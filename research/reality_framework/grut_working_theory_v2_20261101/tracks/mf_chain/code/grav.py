#!/usr/bin/env python3
"""mf_chain: gravitational production of hidden quanta at the abrupt end of the hybrid (finite-k, beyond the homogeneous
engine) [C for the mode integration, E for the post-switch tail].

Free scalar of mass m and curvature coupling xi on the benchmark's own a(t), R(t) (raw R through the homogeneous stage):
  u = a phi,  u'' + [k^2 + a^2 (m^2 + (xi - 1/6) R)] u = 0  (conformal time; integrated in cosmic time, vectorised RK4).
Start: adiabatic vacuum in the de Sitter stage before the slice (low k from t = -2.5, high k from t = -0.3: deep
sub-horizon or heavy); read-out of the adiabatic number at T_READ (after the switch to the cycle-mean fluid R, so that
the read-out is clean).  Fields (M1r, model choices): X-like (xi_X = -5/6, m^2 = 20.2 H_i^2 so that the inflationary
m^2 + xi R = 10.2 matches the prepared state's lightest mass^2; bracket m^2 = 15, 30), Y_L (22.25 H_i, xi = 0, complex),
Y_H (36.1 H_i), sigma_Sigma (40 H_i, real).
Tail after the switch (fluid part; coherent or not, chi chi -> graviton* -> pairs has this rate for non-relativistic chi):
  per real d.o.f. rho_dot = A^2 omega/(64 pi), A = |xi - 1/6| R_osc, omega = 2 m_chi (R_osc_chi) and 2 m_sig (R_osc_sig);
  comoving energy a^4 rho_dot, quanta rate 2 rho_dot/omega.
Outputs (results/grav_<tag>.json): comoving energy E (conformal, a = 1 at the slice) and number N per field, E/E_X of the
history's own release, Delta N_YL, Delta N_YH -> ledger shifts.
Usage: python3 -B grav.py BENCH_JSON TAG [--tread 1.0]
"""
import argparse, json, math
import numpy as np
import mfc_common as M

FIELDS = {'X_m2_20.2': dict(m2=20.2, xi=-5.0 / 6.0, dof=2), 'X_m2_15': dict(m2=15.0, xi=-5.0 / 6.0, dof=2),
          'X_m2_30': dict(m2=30.0, xi=-5.0 / 6.0, dof=2), 'Y_L': dict(m2=22.25 ** 2, xi=0.0, dof=2),
          'Y_H': dict(m2=36.1 ** 2, xi=0.0, dof=2), 'sigma_Sigma': dict(m2=40.0 ** 2, xi=0.0, dof=1)}


def col(h, k):
    return np.array([np.nan if v is None else float(v) for v in h[k]])


def background(d):
    h = d['history']
    t, N, R = col(h, 't'), col(h, 'N'), col(h, 'R')
    N0 = np.interp(0.0, t, N)
    return t, N - N0, R, col(h, 'H'), col(h, 'R_osc_chi'), col(h, 'R_osc_sig')


def modes(tb, Nb, Rb, ks, m2, xi, t0, t1, dt):
    a_of = lambda t: np.exp(np.interp(t, tb, Nb)) if t >= tb[0] else math.exp(t)
    R_of = lambda t: np.interp(t, tb, Rb) if t >= tb[0] else 12.0

    def om2(t):
        a = a_of(t)
        return ks ** 2 + a * a * (m2 + (xi - 1.0 / 6.0) * R_of(t))
    w0 = np.sqrt(np.maximum(om2(t0), 1e-30))
    u = 1.0 / np.sqrt(2 * w0) + 0j
    w = -1j * w0 * u
    t = t0
    nst = int(math.ceil((t1 - t0) / dt))
    h = (t1 - t0) / nst

    def f(t, u, w):
        a = a_of(t)
        return w / a, -om2(t) * u / a
    for _ in range(nst):
        k1u, k1w = f(t, u, w)
        k2u, k2w = f(t + h / 2, u + h / 2 * k1u, w + h / 2 * k1w)
        k3u, k3w = f(t + h / 2, u + h / 2 * k2u, w + h / 2 * k2w)
        k4u, k4w = f(t + h, u + h * k3u, w + h * k3w)
        u = u + h / 6 * (k1u + 2 * k2u + 2 * k3u + k4u)
        w = w + h / 6 * (k1w + 2 * k2w + 2 * k3w + k4w)
        t += h
    W = np.sqrt(np.maximum(om2(t1), 1e-30))
    n = (np.abs(w) ** 2 + W ** 2 * np.abs(u) ** 2) / (2 * W) - 0.5
    wr = np.abs(np.imag(np.conj(u) * w) * 2 + 1)         # Wronskian check: 2 Im(u* w) = -1
    return n, W, float(np.max(wr))


def run(path, tag, tread=1.0):
    d = json.load(open(path))
    tb, Nb, Rb, Hb, Roc, Ros = background(d)
    rel = json.load(open(M.ENG / 'release' / ('release_%s.json' % tag)))
    EX = rel['result']['E_X']
    mchi, msig = d['params']['m_chi_over_Hi'], d['params']['m_sig_over_Hi']
    out = dict(tag=tag, bench=str(path), E_X=EX, t_read=tread, fields={})
    hom_end = float(tb[np.isnan(Roc) & (tb >= 0)].max())
    for name, F in FIELDS.items():
        klo = np.geomspace(0.05, 40.0, 36)
        khi = np.geomspace(40.0, 1500.0, 48)[1:]
        nlo, Wlo, e1 = modes(tb, Nb, Rb, klo, F['m2'], F['xi'], -2.5, tread, 1.0e-4)
        nhi, Whi, e2 = modes(tb, Nb, Rb, khi, F['m2'], F['xi'], -0.3, tread, 2.0e-5)
        # control: identical initial data and integrator on pure de Sitter continued to T_READ (no end of inflation);
        # its 'production' is the adiabatic-start ambiguity plus integrator error, subtracted mode by mode
        tds = np.array([-3.0, tread + 1.0]); Nds = tds.copy(); Rds = np.array([12.0, 12.0])
        clo, _, _ = modes(tds, Nds, Rds, klo, F['m2'], F['xi'], -2.5, tread, 1.0e-4)
        chi_, _, _ = modes(tds, Nds, Rds, khi, F['m2'], F['xi'], -0.3, tread, 2.0e-5)
        ks = np.r_[klo, khi]; n_raw = np.r_[nlo, nhi]; n_ctl = np.r_[clo, chi_]; W = np.r_[Wlo, Whi]
        n = n_raw - n_ctl
        lk = np.log(ks)
        dN = F['dof'] * ks ** 3 * np.maximum(n, 0) / (2 * math.pi ** 2)
        Nk = float(np.trapezoid(dN, lk)); Ek = float(np.trapezoid(dN * W, lk))
        # tail after the homogeneous stage (fluid part), perturbative
        s = (tb > max(hom_end, tread)) & np.isfinite(Roc)
        tt, aa = tb[s], np.exp(Nb[s])
        A_c = abs(F['xi'] - 1.0 / 6.0) * Roc[s]; A_s = abs(F['xi'] - 1.0 / 6.0) * Ros[s]
        m = math.sqrt(F['m2'])
        bc = np.sqrt(np.maximum(1 - (2 * m / (2 * mchi)) ** 2, 0)); bs = np.sqrt(np.maximum(1 - (2 * m / (2 * msig)) ** 2, 0))
        rdot_c = F['dof'] * A_c ** 2 * (2 * mchi) / (64 * math.pi) * bc
        rdot_s = F['dof'] * A_s ** 2 * (2 * msig) / (64 * math.pi) * bs
        Et = float(np.trapezoid((rdot_c + rdot_s) * aa ** 4, tt)) if len(tt) > 1 else 0.0
        Nt = float(np.trapezoid((2 * rdot_c / (2 * mchi) + 2 * rdot_s / (2 * msig)) * aa ** 3, tt)) if len(tt) > 1 else 0.0
        out['fields'][name] = dict(m2=F['m2'], xi=F['xi'], dof=F['dof'], N_modes=Nk, E_modes=Ek, N_tail=Nt, E_tail=Et,
                                   E_total_over_EX=(Ek + Et) / EX, N_total=Nk + Nt, wronskian_err=max(e1, e2),
                                   N_modes_raw=float(np.trapezoid(F['dof'] * ks ** 3 * np.maximum(n_raw, 0) / (2 * math.pi ** 2), lk)),
                                   spectrum=[dict(k=float(k), n=float(x), n_raw=float(y), n_control=float(z))
                                             for k, x, y, z in zip(ks[::4], n[::4], n_raw[::4], n_ctl[::4])])
        print('%s %-12s N_modes %.3e E_modes %.3e | tail N %.3e E %.3e | E/E_X %.3e | Wr %.1e' % (
            tag, name, Nk, Ek, Nt, Et, (Ek + Et) / EX, max(e1, e2)), flush=True)
    M.dump(M.RES / 'grav' / ('grav_%s.json' % tag), out)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('bench'); ap.add_argument('tag'); ap.add_argument('--tread', type=float, default=1.0)
    a = ap.parse_args()
    run(a.bench, a.tag, a.tread)
