#!/usr/bin/env python3
"""hw_verify (Claude, own code): is the Y-relaxed state the physical preparation for a heavy Y during inflation?

Spot point (m = 30, r = 0.45, h = 40 H, frac = 0.5, c_H = 10), delta = pi/2.
(a) Linearised relaxation at the relaxed conditional minimum: canonical mode frequencies omega/H of the held-phase
    directions (x, Re Y, Im Y, s); with Hubble friction 3H each underdamped mode decays as e^{-3N/2} in amplitude,
    e^{-3N} in number.  Also the free (unheld) X-phase and Sigma-phase directions.
(b) De Sitter relaxation (H = 1, Hubble mass on, full potential, cosmic time, DOP853) from the Codex-style state
    with the two supplied phases HELD (arg X = 0, arg Sigma = delta + pi fixed; x, complex Y, s evolve), for N
    e-folds; then an instant release carrying the acquired velocities.  Gives the residual Y displacement, the
    residual free-mode number and the debris R as functions of N.
(c) The same from the relaxed state (control: it should sit still) and with the phases FREE (shows the X-phase
    torque unwinding the supplied phase).
(d) Linear interpolation between the Codex and relaxed states, zero velocity (intermediate preparations).
(e) Estimate [E]: de Sitter excitation of a field with m >> H is Boltzmann suppressed (e^{-2 pi m/H}); a crude
    massless upper bound for the Y number from fluctuations is also given.
Usage: python3 -B prep.py ../results/prep.json     (~3 min on 4 cores)
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hw

H = hw.H_CAND
CP = hw.heavy()
P = hw.pars(CP, H)
DELTA = np.pi / 2


def desitter(z0, N, held=True, v0=None):
    """Evolve N e-folds in de Sitter (H = 1) with 3H friction and the inflationary potential; cosmic time."""
    ph = -np.exp(1j * DELTA)
    if held:
        def rhs(t, u):
            x, yr, yi, s, vx, vyr, vyi, vs = u
            z = np.array([x, yr + 1j * yi, s * ph])
            f = hw.F(z, P, hw.CH)
            ax = -3 * vx - np.real(f[0]); ay = -3 * (vyr + 1j * vyi) - f[1]
            a_s = -3 * vs - np.real(f[2] * np.conj(ph))
            return [vx, vyr, vyi, vs, ax, ay.real, ay.imag, a_s]
        u0 = [z0[0].real, z0[1].real, z0[1].imag, abs(z0[2]), 0, 0, 0, 0]
        if v0 is not None:
            u0[4:] = [v0[0].real, v0[1].real, v0[1].imag, np.real(v0[2] * np.conj(ph))]
        if N == 0:
            return z0, np.zeros(3, complex)
        sol = solve_ivp(rhs, [0, N], u0, method='DOP853', rtol=1e-12, atol=1e-14)
        u = sol.y[:, -1]
        return np.array([u[0], u[1] + 1j * u[2], u[3] * ph]), np.array([u[4], u[5] + 1j * u[6], u[7] * ph])
    def rhs(t, u):
        z = u[0:6:2] + 1j * u[1:6:2]; v = u[6::2] + 1j * u[7::2]
        a = -3 * v - hw.F(z, P, hw.CH)
        out = np.empty(12)
        out[0:6:2], out[1:6:2], out[6::2], out[7::2] = v.real, v.imag, a.real, a.imag
        return out
    u0 = np.zeros(12); u0[0:6:2], u0[1:6:2] = z0.real, z0.imag
    if N == 0:
        return z0, np.zeros(3, complex)
    sol = solve_ivp(rhs, [0, N], u0, method='DOP853', rtol=1e-12, atol=1e-14)
    u = sol.y[:, -1]
    return u[0:6:2] + 1j * u[1:6:2], u[6::2] + 1j * u[7::2]


def release(args):
    tag, z, v = args
    r = hw.run(CP, H, np.asarray(z, complex), a_end=100.0, windows=(1, 10, 30, 100), v0=np.asarray(v, complex))
    w = r['windows']['30-100']
    o = hw.summarize(CP, H, w, extra=dict(tag=tag, controls=r['end']))
    o['state'] = dict(z=[[c.real, c.imag] for c in np.asarray(z, complex)], v=[[c.real, c.imag] for c in np.asarray(v, complex)])
    return o


def ymode_numbers(z, v):
    """Free-mode numbers of the Y displacement about the forced response at a = 1 (same split as the release code)."""
    x, y, s = z
    A = P['M']**2; B = 2 * P['mu'] * s; J = P['kap'] * x**3; D = A * A - abs(B)**2
    yp = (-A * J + B * np.conj(J)) / D
    Jp = 3 * P['kap'] * x * x * v[0] + 2 * P['mu'] * v[2] * np.conj(yp)
    ypp = (-A * Jp + B * np.conj(Jp)) / D
    d = y - yp; dd = v[1] - ypp
    rot = np.exp(-1j * np.angle(B) / 2)
    u = rot * d; w = rot * dd
    wH, wL = np.sqrt(A + abs(B)), np.sqrt(A - abs(B))
    return dict(abs_dY=abs(d), NYH=(w.real**2 + wH**2 * u.real**2) / wH, NYL=(w.imag**2 + wL**2 * u.imag**2) / wL)


def main(out):
    t0 = time.time()
    zc, ic = hw.prep_codex(CP, H, hw.CH, DELTA)
    zr, ir = hw.prep_relaxed(CP, H, hw.CH, DELTA)
    res = dict(codex_state=ic, relaxed_state=ir)
    # (a) linearised modes about the relaxed minimum
    q = np.array([ir['x'], ir['Y'][0], ir['Y'][1], ir['s']])
    H4 = hw._hess4(q, P, hw.CH, DELTA)
    w2, U = np.linalg.eigh(0.5 * H4)     # canonical: L = xdot^2 + |Ydot|^2 + sdot^2 - V
    modes = []
    for k in range(4):
        om2 = w2[k]; vec = U[:, k]
        modes.append(dict(omega2_over_H2=float(om2), omega_over_H=float(np.sqrt(max(om2, 0))),
                          underdamped=bool(om2 > 2.25), amplitude_decay_rate_over_H=1.5 if om2 > 2.25 else float(1.5 - np.sqrt(2.25 - om2)),
                          composition=dict(x=float(vec[0]), ReY=float(vec[1]), ImY=float(vec[2]), s=float(vec[3]))))
    # phase directions (free): second derivative of V w.r.t. the two phases at the relaxed point, canonical
    zr_ = np.array(zr)
    def Vph(th, ps):
        z = np.array([zr_[0] * np.exp(1j * th), zr_[1], zr_[2] * np.exp(1j * ps)])
        return hw.V(z, P, hw.CH)
    e = 1e-4
    Vtt = (Vph(e, 0) - 2 * Vph(0, 0) + Vph(-e, 0)) / e**2
    Vpp = (Vph(0, e) - 2 * Vph(0, 0) + Vph(0, -e)) / e**2
    Vtp = (Vph(e, e) - Vph(e, -e) - Vph(-e, e) + Vph(-e, -e)) / (4 * e * e)
    gX = (Vph(e, 0) - Vph(-e, 0)) / (2 * e); gS = (Vph(0, e) - Vph(0, -e)) / (2 * e)
    K = np.diag([abs(zr_[0])**2, abs(zr_[2])**2])          # kinetic: |X|^2 thdot^2 + |S|^2 psdot^2
    Wm = 0.5 * np.array([[Vtt, Vtp], [Vtp, Vpp]])
    ev = np.linalg.eigvals(np.linalg.solve(K, Wm)).real
    res['linearised'] = dict(held_modes=modes, phase_mass2_over_H2=sorted(ev.tolist()),
                             phase_gradient=dict(dV_dthetaX=gX, dV_dthetaS=gS,
                                                 slow_roll_dthetaX_per_efold=-gX / (3 * 2 * abs(zr_[0])**2),
                                                 slow_roll_dthetaS_per_efold=-gS / (3 * 2 * abs(zr_[2])**2)),
                             note='canonical frequencies; number of an underdamped mode decays as e^{-3N}')
    res['displacement_codex_vs_relaxed'] = dict(dX=float(abs(zc[0] - zr[0])), dY=float(abs(zc[1] - zr[1])),
                                                dS=float(abs(abs(zc[2]) - abs(zr[2]))), ymodes_codex=ymode_numbers(zc, np.zeros(3)),
                                                ymodes_relaxed=ymode_numbers(zr, np.zeros(3)))
    # (b), (c): de Sitter stages
    Ns = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
    stages = [('codex_held_N%g' % N, *desitter(zc, N, True)) for N in Ns]
    stages += [('relaxed_held_N%g' % N, *desitter(zr, N, True)) for N in (0.5, 2.0)]
    stages += [('relaxed_free_N%g' % N, *desitter(zr, N, False)) for N in (0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0)]
    stages += [('codex_free_N%g' % N, *desitter(zc, N, False)) for N in (0.1, 0.5, 2.0)]
    fs = [0.0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
    stages += [('interp_f%g' % f, zc + f * (zr - zc), np.zeros(3, complex)) for f in fs]
    print('stages prepared', len(stages), '%.0fs' % (time.time() - t0), flush=True)
    with Pool(4) as pool:
        rel = pool.map(release, stages)
    for (tag, z, v), r in zip(stages, rel):
        r['ymodes_at_release'] = ymode_numbers(np.asarray(z, complex), np.asarray(v, complex))
        r['phases_at_release'] = dict(argX=float(np.angle(z[0])), argS_minus_pi=float(np.angle(-z[2])),
                                      argY=float(np.angle(z[1])), absX=float(abs(z[0])), absY=float(abs(z[1])), absS=float(abs(z[2])))
        print(tag, 'eps %.5f NYH %.4g NYL %.4g R %.4g RL %s |dY| %.3g' % (r['eps'], r['NYH'], r['NYL'], r['R_central'], r['R_L'],
                                                                           r['ymodes_at_release']['abs_dY']), flush=True)
    res['releases'] = rel
    # fit of the decay rate of the free Y number with N (held, Codex start), N = 0.2..2
    pts = [(N, r['NYH'] + r['NYL']) for N, r in zip(Ns, rel[:len(Ns)]) if 0.2 <= N <= 2.0]
    A = np.vstack([np.ones(len(pts)), [p[0] for p in pts]]).T
    coef = np.linalg.lstsq(A, np.log([p[1] for p in pts]), rcond=None)[0]
    res['NY_decay_fit'] = dict(N_range=[0.2, 2.0], d_lnNY_dN=float(coef[1]), predicted=-3.0)
    # (e) fluctuation estimates
    mY = P['M']
    res['fluctuation_estimates'] = dict(
        boltzmann_factor_exp_minus_2pi_m_over_H=float(np.exp(-2 * np.pi * mY)),
        massless_upper_bound_NY_per_H3=float(np.sqrt(mY**2) * (1 / (2 * np.pi))**2),
        heavy_variance_estimate_H2=float((1 / (2 * np.pi))**2 * (1 / mY)**2),
        note='[E] a field with m = 30 H is not excited by de Sitter fluctuations: particle production ~ e^{-2 pi m/H}; '
             'even the massless variance (H/2pi)^2 would give only ~ m (H/2pi)^2 = %.3g quanta per H^-3' % (mY / (4 * np.pi**2)))
    Path(out).write_text(json.dumps(dict(method=__doc__, H_GeV=H, delta=DELTA, results=res), indent=1, default=float))
    print('done %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main(sys.argv[1])
