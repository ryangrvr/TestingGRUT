#!/usr/bin/env python3
"""Full nonlinear radial time-domain check of first- vs second-order charge drain (track xcharge).

Evolves the complete classical equations (no linearization) with
   chi_tt = lap chi - (S^2 + 2A|chi|^2) chi - d eta conj(chi)^(d-1),    X = v chi
   S_tt   = lap S   - B (S^2-1) S - 2|chi|^2 S,                         phi = v S
eta = lam v^(d-2)/m_X^2 (real), switched on smoothly. Two runs (eta and 0) from the same
ground state; their core-charge difference isolates the breaking effect, cancelling the
common discretization breathing. Radial variables u = rho chi, w = rho (S-1); 4th-order
finite differences, RK4, polynomial sponge for outgoing waves.

Outputs: q_core(tau) = g^2 Q_core inside rho < Rc, for later fitting of oscillation (O(eta))
and secular drift (O(eta^2)).
"""
import numpy as np
from profile import A as AQ, B as BQ


def run(bg, d, eta0, h=0.04, dt=0.01, Rmax=110.0, Rs=50.0, smax=3.0, Rc=25.0, T=320.0, ton=40.0, every=5):
    Om = bg.Om
    r = np.arange(0, Rmax + h / 2, h)
    n = r.size
    F, S = bg.FS(r)
    u = (r * F / np.sqrt(2)).astype(complex)
    ut = -1j * Om * u
    w = r * (S - 1)
    wt = np.zeros_like(w)
    sig = np.where(r > Rs, smax * ((r - Rs) / (Rmax - Rs))**2, 0.0)
    rin = np.where(r > 0, r, 1.0)
    ic = int(round(Rc / h))
    # Simpson-like weights for core integral up to Rc
    wq = np.ones(ic + 1); wq[1:-1:2] = 4; wq[2:-1:2] = 2; wq *= h / 3

    def lap(f, odd=True):
        g = np.empty(n + 4, dtype=f.dtype)
        g[2:-2] = f
        g[1] = -f[1]; g[0] = -f[2]          # odd extension (u, w vanish at origin)
        g[-2] = 0; g[-1] = 0
        return (-g[4:] + 16 * g[3:-1] - 30 * g[2:-2] + 16 * g[1:-3] - g[:-4]) / (12 * h * h)

    def acc(u, ut, w, wt, eta):
        chi = u / rin; chi[0] = 0
        Sf = 1 + w / rin; Sf[0] = 1 + (w[1] / r[1])  # w~rho*(S-1); only used with factor r below
        a2 = np.abs(chi)**2
        au = lap(u) - (Sf * Sf + 2 * AQ * a2) * u - sig * ut
        if eta != 0.0:
            au -= d * eta * r * np.conj(chi)**(d - 1)
        aw = lap(w) - r * (BQ * (Sf * Sf - 1) * Sf + 2 * a2 * Sf) - sig * wt
        au[0] = 0; aw[0] = 0; au[-1] = 0; aw[-1] = 0
        return au, aw

    def env(t):
        return 1.0 if t >= ton else np.sin(0.5 * np.pi * t / ton)**2

    ts, qs = [], []
    t = 0.0
    nsteps = int(round(T / dt))
    for k in range(nsteps + 1):
        if k % every == 0:
            q = 4 * np.pi * np.sum(wq * (-2) * np.imag(np.conj(u[:ic + 1]) * ut[:ic + 1]))
            ts.append(t); qs.append(q)
        if k == nsteps:
            break
        e1 = eta0 * env(t); e2 = eta0 * env(t + dt / 2); e3 = eta0 * env(t + dt)
        a1u, a1w = acc(u, ut, w, wt, e1)
        u2 = u + dt / 2 * ut; ut2 = ut + dt / 2 * a1u; w2 = w + dt / 2 * wt; wt2 = wt + dt / 2 * a1w
        a2u, a2w = acc(u2, ut2, w2, wt2, e2)
        u3 = u + dt / 2 * ut2; ut3 = ut + dt / 2 * a2u; w3 = w + dt / 2 * wt2; wt3 = wt + dt / 2 * a2w
        a3u, a3w = acc(u3, ut3, w3, wt3, e2)
        u4 = u + dt * ut3; ut4 = ut + dt * a3u; w4 = w + dt * wt3; wt4 = wt + dt * a3w
        a4u, a4w = acc(u4, ut4, w4, wt4, e3)
        u = u + dt / 6 * (ut + 2 * ut2 + 2 * ut3 + ut4)
        ut = ut + dt / 6 * (a1u + 2 * a2u + 2 * a3u + a4u)
        w = w + dt / 6 * (wt + 2 * wt2 + 2 * wt3 + wt4)
        wt = wt + dt / 6 * (a1w + 2 * a2w + 2 * a3w + a4w)
        t += dt
    return np.array(ts), np.array(qs)


def analyse(ts, dq, d, Om, t0=120.0):
    """Least-squares fit dq = a + b t + sum_{k=1..3} [c_k cos(k d Om t) + s_k sin(k d Om t)] on t>=t0."""
    m = ts >= t0
    t = ts[m]; y = dq[m]
    cols = [np.ones_like(t), t]
    for k in (1, 2, 3):
        cols += [np.cos(k * d * Om * t), np.sin(k * d * Om * t)]
    Amat = np.vstack(cols).T
    coef, *_ = np.linalg.lstsq(Amat, y, rcond=None)
    resid = y - Amat @ coef
    osc1 = float(np.hypot(coef[2], coef[3]))
    return dict(drift=float(coef[1]), osc_amp_harm1=osc1, resid_rms=float(np.sqrt(np.mean(resid**2))))


if __name__ == '__main__':
    import json, sys, time
    from profile import ground_state
    from forced import Bg
    from scipy.integrate import simpson
    d = int(sys.argv[1]); etas = [float(x) for x in sys.argv[2].split(',')]
    out_path = sys.argv[3]
    h = float(sys.argv[4]) if len(sys.argv) > 4 else 0.04
    dt = float(sys.argv[5]) if len(sys.argv) > 5 else 0.01
    T = float(sys.argv[6]) if len(sys.argv) > 6 else 420.0
    Om, sol = ground_state(200000.)
    bg = Bg(sol, Om)
    x = np.linspace(0, bg.R, 40001); Fx, _ = bg.FS(x)
    I_d = 4 * np.pi * simpson(x * x * (Fx / np.sqrt(2))**d, x=x)
    t0 = time.time()
    ts, q0 = run(bg, d, 0.0, h=h, dt=dt, T=T)
    res = dict(d=d, h=h, dt=dt, T=T, Omega=Om, I_d=float(I_d), q0_start=float(q0[0]), q0_end=float(q0[-1]),
               q0_ptp=float(np.ptp(q0)), runs=[])
    for eta in etas:
        _, qp = run(bg, d, eta, h=h, dt=dt, T=T)
        _, qm = run(bg, d, -eta, h=h, dt=dt, T=T)
        even = 0.5 * ((qp - q0) + (qm - q0)); odd = 0.5 * ((qp - q0) - (qm - q0))
        ae = analyse(ts, even, d, Om, t0=160.0); ao = analyse(ts, odd, d, Om, t0=160.0)
        row = dict(eta=eta, even_drift=ae['drift'], even_resid=ae['resid_rms'], odd_osc_amp=ao['osc_amp_harm1'],
                   odd_drift=ao['drift'], odd_resid=ao['resid_rms'],
                   predicted_first_order_amp=float(2 * eta * I_d / Om),
                   plus_drift=analyse(ts, qp - q0, d, Om, t0=160.0)['drift'],
                   minus_drift=analyse(ts, qm - q0, d, Om, t0=160.0)['drift'])
        res['runs'].append(row)
        np.save(out_path.replace('.json', f'_eta{eta:g}.npy'), np.vstack([ts, qp - q0, qm - q0]))
        print(d, eta, row, f'{time.time()-t0:.0f}s', flush=True)
        json.dump(res, open(out_path, 'w'), indent=1)
