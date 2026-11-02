#!/usr/bin/env python3
"""inf_potential (Claude, own code): linear fluctuations on the post-inflation background of each benchmark.

(1) Floquet chart (Minkowski, no expansion): for the homogeneous oscillation of amplitude A (turning point on the
    inflationary side, phi_I = A), integrate two real solutions of  dphi_k'' + (k^2 + V''(phi(t))) dphi_k = 0  over one
    period; the Floquet exponent is mu_k = ln|lambda_max| / T.  Reported: mu_max/m_I and mu_max/H_i at the amplitudes the
    real background passes through (cf. Amin et al. 1106.3335; Lozanov-Amin 1608.01213, 1710.06851) [C].
(2) Expanding background, inflaton self-resonance: complex modes dphi_k (Bunch-Davies-like start at eps_H = 1,
    dphi = e^{-i w t}/(a^{3/2} sqrt(2w)), w = k/a for the sub-horizon band) on the full background with Gamma_I = 0 or 0.5 H_i
    (fixed-step RK4, background integrated in the same state; the phenomenological width damps the modes as it damps the
    background).  Occupations n_k against w^2 = k^2/a^2 + m_I^2; the excess variance int n_k/(w a^3) d^3k/(2pi)^3 against
    the current amplitude^2 marks backreaction/fragmentation when it reaches ~0.1-1 [C on the linear theory; E for the
    fragmentation criterion].  Metric perturbations are dropped (sub-horizon band, k/aH >~ 10).
(3) Control of the Born Y-pair estimate (yy_harmonics.py): light-Y modes  y_k'' + 3H y_k' + (k^2/a^2 + m_L^2 + 2 c_L phi_I(t)) y_k = 0
    from vacuum, counted with the adiabatic number at the end, against the Born count of the same channel [C].
Units: M_P = 1, time 1/sqrt(V0) internally; reported in H_i (H at S1) and m_I.
Usage: python3 modes.py <benchmark> [floquet|resonance|ycontrol] ...   Output: ../results/modes_<benchmark>.json (merged).
"""
import json, math, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0, str(Path(__file__).resolve().parent))
import inf_common as ic

RES = Path(__file__).resolve().parent.parent / 'results'
ALPHA_I = math.atan2(10.0, -7.159456159)
SIG = 40 / math.sqrt(2)
M_L = 22.249


def load(name, Hi_index=0):
    d = json.load(open(RES / ('bench_%s.json' % name)))
    rec = d['by_Hi'][Hi_index]
    kind = {'T1': 'T', 'E1': 'E', 'HT6t': 'HT', 'HT6': 'HT', 'HT4': 'HT', 'HT8': 'HT', 'PL2': 'PL'}[name]
    m = ic.Model(kind, **rec['params'])
    return m, rec


def floquet(m, rec, amps):
    """Minkowski Floquet exponents at turning-point amplitudes phi_I = A (on the inflationary side)."""
    pmin = m.phi_min()
    sgn = 1.0 if rec['phiI_S1_MP'] > 0 else -1.0
    mI = math.sqrt(m.m2_min())
    hs = rec['H_i_GeV'] / ic.MPL_GEV / math.sqrt(rec['V0_MP4'])
    ks = np.linspace(0.02, 2.5, 125) * mI
    out = []
    for A in amps:
        f0 = pmin + sgn * A
        E = m.v(f0)
        # period: integrate the background until the velocity returns to zero on the same side twice
        def rhs(t, y):
            return [y[1], -m.dv(y[0])]
        ev = lambda t, y: y[1]
        ev.direction = 0
        sol = solve_ivp(rhs, (0, 200 * 2 * math.pi / mI), [f0, 0.0], method='DOP853', rtol=1e-11, atol=1e-14, events=ev, dense_output=True)
        te = sol.t_events[0]
        te = te[te > 1e-9]
        T = te[1]
        # monodromy for all k simultaneously
        def rhs2(t, y):
            f = sol.sol(t)[0]
            d2 = m.d2v(f)
            n = len(ks)
            u1, v1, u2, v2 = y[:n], y[n:2 * n], y[2 * n:3 * n], y[3 * n:]
            w2 = ks ** 2 + d2
            return np.concatenate([v1, -w2 * u1, v2, -w2 * u2])
        n = len(ks)
        y0 = np.concatenate([np.ones(n), np.zeros(n), np.zeros(n), np.ones(n)])
        s2 = solve_ivp(rhs2, (0, T), y0, method='DOP853', rtol=1e-10, atol=1e-12)
        Y = s2.y[:, -1]
        u1, v1, u2, v2 = Y[:n], Y[n:2 * n], Y[2 * n:3 * n], Y[3 * n:]
        tr = u1 + v2
        lam = np.abs(tr / 2) + np.sqrt(np.maximum((tr / 2) ** 2 - 1, 0))
        mu = np.log(np.maximum(lam, 1.0)) / T
        j = int(np.argmax(mu))
        out.append(dict(A_MP=A, V_over_V0=E, omega_over_mI=2 * math.pi / T / mI, mu_max_over_mI=float(mu[j] / mI),
                        k_at_max_over_mI=float(ks[j] / mI), mu_max_over_Hi=float(mu[j] / hs),
                        band_width_over_mI=float(np.sum(mu > 0.5 * mu[j]) * (ks[1] - ks[0]) / mI) if mu[j] > 0 else 0.0))
        print('  floquet A=%.4e MP V/V0=%.3f w/mI=%.3f mu_max/mI=%.4f at k/mI=%.3f  mu_max/H_i=%.2f' % (
            A, E, out[-1]['omega_over_mI'], out[-1]['mu_max_over_mI'], out[-1]['k_at_max_over_mI'], out[-1]['mu_max_over_Hi']), flush=True)
    return out


def rk4_run(m, rec, gam_Hi, Nafter, ks, extra=None, steps_per_osc=48):
    """fixed-step RK4 of background + complex inflaton modes from eps_H = 1 (a = 1) for Nafter e-folds."""
    bg = ic.background(m, 75.0)
    V0 = rec['V0_MP4']
    hs = rec['H_i_GeV'] / ic.MPL_GEV / math.sqrt(V0)
    f, fN = ic.bg_at(bg, 0.0)
    h, eps = ic.h_of(m, f, fN)
    fd = h * fN
    mI = math.sqrt(m.m2_min())
    dt = 2 * math.pi / (mI * steps_per_osc)
    g = gam_Hi * hs
    nk = len(ks)
    w0 = ks.astype(float)
    z0 = 1 / np.sqrt(2 * w0)
    zr, zi = z0.copy(), np.zeros(nk)
    zdr = -1.5 * h * z0
    zdi = -w0 * z0
    st = np.concatenate([[0.0, f, fd, 0.0], zr, zi, zdr, zdi])

    def deriv(s):
        lna, f, fd, rr = s[:4]
        a = math.exp(lna)
        rho = 0.5 * fd * fd + m.v(f) + rr
        H = math.sqrt(rho / 3)
        zr, zi, zdr, zdi = s[4:4 + nk], s[4 + nk:4 + 2 * nk], s[4 + 2 * nk:4 + 3 * nk], s[4 + 3 * nk:]
        w2 = (ks / a) ** 2 + m.d2v(f)
        return np.concatenate([[H, fd, -(3 * H + g) * fd - m.dv(f), g * fd * fd - 4 * H * rr], zdr, zdi, -(3 * H + g) * zdr - w2 * zr, -(3 * H + g) * zdi - w2 * zi])
    t = 0.0; rows = []; last_turn = None
    pmin = m.phi_min()
    amp = abs(f - pmin)
    nextrec = 0.0
    prev_fd = fd
    while st[0] < Nafter:
        k1 = deriv(st); k2 = deriv(st + 0.5 * dt * k1); k3 = deriv(st + 0.5 * dt * k2); k4 = deriv(st + dt * k3)
        st = st + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt
        fdn = st[2]
        if np.sign(fdn) != np.sign(prev_fd):
            amp = abs(st[1] - pmin)
        prev_fd = fdn
        if st[0] >= nextrec:
            a = math.exp(st[0])
            zr, zi = st[4:4 + nk], st[4 + nk:4 + 2 * nk]
            zdr, zdi = st[4 + 2 * nk:4 + 3 * nk], st[4 + 3 * nk:]
            # occupation of chi = a^{3/2} dphi against w^2 = k^2/a^2 + m_I^2 (the minimum's mass; adiabatic once the
            # amplitude is small, a diagnostic while it is large)
            a32 = a ** 1.5
            Hn = math.sqrt((0.5 * st[2] ** 2 + m.v(st[1]) + st[3]) / 3)
            w = np.sqrt((ks / a) ** 2 + mI ** 2)
            cr, ci = a32 * zr, a32 * zi
            cdr, cdi = a32 * (zdr + 1.5 * Hn * zr), a32 * (zdi + 1.5 * Hn * zi)
            nkk = np.maximum((cdr ** 2 + cdi ** 2 + w ** 2 * (cr ** 2 + ci ** 2)) / (2 * w) - 0.5, 0.0)
            # excess (non-vacuum) variance: int k^2 dk/(2 pi^2) n_k/(w a^3), comoving k -> per ln k: k^3 n_k/(2 pi^2 w a^3)
            P = ks ** 3 * nkk / (2 * math.pi ** 2 * w * a ** 3)
            var = float(np.trapezoid(P, np.log(ks))) * V0                          # times V0 -> M_P^2
            jmax = int(np.argmax(nkk))
            rows.append(dict(N_after=float(st[0]), t_Hi=t * hs, a=a, amp_MP=amp, var_excess_MP2=var, ratio=var / amp ** 2,
                             n_k_max=float(nkk.max()), k_peak_phys_over_mI=float(ks[jmax] / a / mI)))
            nextrec += 0.05
    return rows


def resonance(m, rec, gam_Hi, Nafter):
    mI = math.sqrt(m.m2_min())
    ks = mI * np.logspace(-1.3, 0.6, 160)          # comoving k (a = 1 at the end) in sqrt(V0) units
    rows = rk4_run(m, rec, gam_Hi, Nafter, ks)
    frag = None
    for r in rows:
        if r['ratio'] > 0.1 and frag is None:
            frag = r
    print('  resonance Gamma=%.1f: ratio at N=0.5,1,1.5,2: %s ; n_k,max at N=1,2: %s ; ratio>0.1 at %s' % (
        gam_Hi, ['%.2e' % np.interp(x, [r['N_after'] for r in rows], [r['ratio'] for r in rows]) for x in (0.5, 1, 1.5, 2)],
        ['%.2e' % np.interp(x, [r['N_after'] for r in rows], [r['n_k_max'] for r in rows]) for x in (1, 2)],
        ('N=%.2f t_Hi=%.2f k/a/mI=%.2f' % (frag['N_after'], frag['t_Hi'], frag['k_peak_phys_over_mI'])) if frag else 'never'), flush=True)
    return dict(gamma_over_Hi=gam_Hi, N_after_max=Nafter, rows=rows[::2], first_ratio_0p1=frag)


def ycontrol(m, rec, gam_Hi, tmax_Hi):
    """light-Y modes on the background, LL channel only, compared with the Born count."""
    bg = ic.background(m, 75.0)
    V0 = rec['V0_MP4']; Hi = rec['H_i_GeV']
    hs = Hi / ic.MPL_GEV / math.sqrt(V0)
    Hi_MP = Hi / ic.MPL_GEV
    gI = abs(10j - 7.159456159) * Hi_MP / abs(rec['phiI_S1_MP'])
    cL = gI * abs(math.cos(ALPHA_I)) * SIG                  # H_i units, phi_I in H_i
    pmin = m.phi_min()
    sol = ic.post_inflation(m, bg, gam_Hi * hs, tmax_Hi / hs, rtol=1e-10, atol=1e-14, dense=True)
    # work in H_i units: time tau = t_code * hs; phi_I/H_i = (phi - pmin)/Hi_MP
    ks = np.linspace(0.25, 250.0, 500)                      # comoving k in H_i (a = 1 at the end)
    nk = len(ks)
    def bgH(tau):
        y = sol.sol(tau / hs)
        lna, f, fd, rr = y
        rho = 0.5 * fd * fd + m.v(f) + rr
        return math.exp(lna), math.sqrt(rho / 3) / hs, (f - pmin) / Hi_MP
    def deriv(tau, s):
        a, H, ph = bgH(tau)
        yr, yi, ydr, ydi = s[:nk], s[nk:2 * nk], s[2 * nk:3 * nk], s[3 * nk:]
        w2 = (ks / a) ** 2 + M_L ** 2 + 2 * cL * ph
        return np.concatenate([ydr, ydi, -3 * H * ydr - w2 * yr, -3 * H * ydi - w2 * yi])
    a0, H0, ph0 = bgH(0.0)
    w0 = np.sqrt(ks ** 2 + M_L ** 2 + 2 * cL * ph0)
    y0 = np.concatenate([1 / np.sqrt(2 * w0), np.zeros(nk), -1.5 * H0 / np.sqrt(2 * w0), -w0 / np.sqrt(2 * w0)])
    s = solve_ivp(deriv, (0, tmax_Hi), y0, method='DOP853', rtol=1e-9, atol=1e-12, t_eval=[tmax_Hi])
    a, H, ph = bgH(tmax_Hi)
    yr, yi, ydr, ydi = [s.y[i * nk:(i + 1) * nk, -1] for i in range(4)]
    w = np.sqrt((ks / a) ** 2 + M_L ** 2 + 2 * cL * ph)
    # adiabatic number of chi = a^{3/2} y: (|chi'|^2 + w^2 |chi|^2)/(2 w) - 1/2 with chi' = a^{3/2}(y' + 1.5 H y)
    a32 = a ** 1.5
    cr, ci = a32 * yr, a32 * yi
    cdr, cdi = a32 * (ydr + 1.5 * H * yr), a32 * (ydi + 1.5 * H * yi)
    nkk = np.maximum((cdr ** 2 + cdi ** 2 + w ** 2 * (cr ** 2 + ci ** 2)) / (2 * w) - 0.5, 0.0)
    N_LL = float(np.trapezoid(ks ** 2 * nkk, ks) / (2 * math.pi ** 2))   # comoving number (a = 1 at the end), H_i^3
    print('  ycontrol Gamma=%.1f tmax=%g: linear-mode N_YL(LL) = %.4e (peak n_k %.3e at k=%.2f)' % (gam_Hi, tmax_Hi, N_LL, nkk.max(), ks[np.argmax(nkk)]), flush=True)
    return dict(gamma_over_Hi=gam_Hi, t_max_Hi=tmax_Hi, N_YL_LL_linear=N_LL, n_k_peak=float(nkk.max()), k_peak=float(ks[np.argmax(nkk)]),
                spectrum=[dict(k=float(k), n=float(n)) for k, n in zip(ks[::10], nkk[::10])])


def main(name, tasks):
    m, rec = load(name)
    try:
        out = json.load(open(RES / ('modes_%s.json' % name)))
    except Exception:
        out = dict(benchmark=name)
    out['method'] = __doc__
    t0 = time.time()
    pmin = m.phi_min()
    if 'floquet' in tasks:
        A1 = max(abs(r['phiI_turn']) for r in rec['osc_gamma0']['turning'][:2])
        amps = [A1 * x for x in (1.0, 0.7, 0.5, 0.35, 0.2, 0.1, 0.05)]
        out['floquet'] = floquet(m, rec, amps)
    if 'resonance' in tasks:
        out['resonance'] = [resonance(m, rec, 0.0, 2.5), resonance(m, rec, 0.5, 2.5)]
    if 'ycontrol' in tasks:
        out['ycontrol'] = [ycontrol(m, rec, 0.5, 12.0)]
    json.dump(out, open(RES / ('modes_%s.json' % name), 'w'), indent=1, default=float)
    print(name, 'done %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
