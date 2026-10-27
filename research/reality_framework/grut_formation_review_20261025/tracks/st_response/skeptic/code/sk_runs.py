#!/usr/bin/env python3
"""st_response skeptic (Claude, own code): chunked, checkpointed runs of sk_tangent.py on a 0.002 tau grid (plus
Codex's archived sample times, so values can be compared at identical tau), for ALL 18 of Codex's momenta.

Usage: python3 -B sk_runs.py <job> [max_seconds]     (re-run until it prints DONE; resumes from its checkpoint)
Jobs:
  q3      (3, 0.5), 18 Codex momenta + the old-archive momentum 2.122989207841491, tau 0-800
  q5      (5, 1),   18 Codex momenta, tau 0-800
  q3tight (3, 0.5), k = 2.0, 2.13, 2.4 at rtol 1e-12 / atol 1e-14, tau 0-160   (own tolerance control)
  q3full  (3, 0.5), k = 2.13 with the full nonminimal geometry (f and xi eta2 terms) is NOT implemented; instead
          q3fd below tests the linearisation itself.
  q3fd    (3, 0.5), k = 2.13 and 2.4: NONLINEAR two-point antisymmetric finite difference (no Hessian, no dual
          numbers): points x = 0 and x = pi/k carry z + d0 and z + d1, d0(0) = +eps dz, d1(0) = -eps dz, the spectral
          Laplacian of the antisymmetric mode is -(k^2/a^2)(d0 - d1)/2, and F(z + d) - F(z) is evaluated directly.
          The response (d0 - d1)/(2 eps) is odd in eps (error O(eps^2)).  eps = 1e-6, tau 0-200.
Outputs (scratch): runs/<job>/chunk_###.npz with float64 A, dN, N_X, a on the union grid; checkpoint.npz.
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import sk_tangent as S

SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')
CODEX = SP / 'b73/research/reality_framework/grut_source_timing_20261026'
OUT = SP / 'review/st_response_skeptic/runs'
KS = [.02, .05, .1, .2, .4, .8, 1.2, 1.6, 2., 2.13, 2.4, 2.8, 3.2, 4., 5.6, 5.8, 6., 8.]   # [Q] Codex PROTOCOL.md
KOLD = 2.122989207841491
JOBS = {'q3': dict(b=3.0, gam=0.5, ks=KS + [KOLD], tmax=800.0),
        'q5': dict(b=5.0, gam=1.0, ks=KS, tmax=800.0),
        'q3tight': dict(b=3.0, gam=0.5, ks=[2.0, 2.13, 2.4], tmax=160.0, rtol=1e-12, atol=1e-14),
        'q3fd': dict(b=3.0, gam=0.5, ks=[2.13, 2.4], tmax=200.0, fd=1e-6)}
CHUNK = 25.0
FINE = 0.002


class FDModel(S.Model):
    """Background + two lattice points' deviations per k (d0, d1: 6 + 6 real each)."""

    def __init__(self, b, gam, ks, eps):
        super().__init__(b, gam, ks); self.eps = eps

    def rhs(self, t, y):
        z, v, H, R = self.geo(y)
        out = np.empty_like(y)
        out[0] = H; out[1] = y[2]; out[2] = -(3 * H + self.gam) * y[2] - self.b ** 2 * y[1]
        out[3] = self.gam * y[2] ** 2 - 4 * H * y[3]
        F0 = np.array(S.force(*z))
        acc = -3 * H * v - F0; acc[0] -= S.XI * R * z[0]
        out[4:10:2], out[5:10:2] = v.real, v.imag
        out[10:16:2], out[11:16:2] = acc.real, acc.imag
        a2 = np.exp(2 * y[0])
        T = y[16:].reshape(self.nk, 2, 12)
        d = T[:, :, 0:6:2] + 1j * T[:, :, 1:6:2]; dv = T[:, :, 6:12:2] + 1j * T[:, :, 7:12:2]
        zz = z[None, None, :] + d
        Fp = np.stack(S.force(zz[..., 0], zz[..., 1], zz[..., 2]), axis=-1)
        lap = (d[:, 0, :] - d[:, 1, :]) / 2                     # antisymmetric amplitude
        dacc = -3 * H * dv - (Fp - F0[None, None, :])
        dacc[:, 0, :] -= (self.k2 / a2)[:, None] * lap
        dacc[:, 1, :] += (self.k2 / a2)[:, None] * lap
        dacc[:, :, 0] -= S.XI * R * d[:, :, 0]
        o = out[16:].reshape(self.nk, 2, 12)
        o[:, :, 0:6:2], o[:, :, 1:6:2] = dv.real, dv.imag
        o[:, :, 6:12:2], o[:, :, 7:12:2] = dacc.real, dacc.imag
        return out


def fd_initial(model):
    y0, z0, dz, info = S.initial(S.Model(model.b, model.gam, model.ks))
    y = np.zeros(16 + 24 * model.nk); y[:16] = y0[:16]
    T = y[16:].reshape(model.nk, 2, 12)
    for p, sg in ((0, 1.0), (1, -1.0)):
        T[:, p, 0:6:2], T[:, p, 1:6:2] = sg * model.eps * dz.real, sg * model.eps * dz.imag
    return y


def fd_to_tangent(model, Y):
    """Map FD state to the tangent layout of sk_tangent: response = (d0 - d1)/(2 eps)."""
    T = Y[:, 16:].reshape(len(Y), model.nk, 2, 12)
    resp = (T[:, :, 0, :] - T[:, :, 1, :]) / (2 * model.eps)
    return np.concatenate([Y[:, :16], resp.reshape(len(Y), -1)], axis=1)


def main(job, budget):
    spec = JOBS[job]; tmax = spec['tmax']
    out = OUT / job; out.mkdir(parents=True, exist_ok=True)
    ck = out / 'checkpoint.npz'
    if spec.get('fd'):
        model = FDModel(spec['b'], spec['gam'], spec['ks'], spec['fd']); tang = S.Model(spec['b'], spec['gam'], spec['ks'])
    else:
        model = S.Model(spec['b'], spec['gam'], spec['ks']); tang = model
    if ck.exists():
        c = np.load(ck); t0, y = float(c['t']), c['y']
    else:
        t0 = 0.0
        y = fd_initial(model) if spec.get('fd') else S.initial(model)[0]
    grid = np.load(CODEX / ('quiet3.npz' if spec['b'] == 3.0 else 'quiet5.npz'))['times']
    oldgrid = np.load(CODEX.parent / 'grut_reheating_modes_20261026' / 'PHASE_PROFILES.npz')['times']
    start = time.time()
    while t0 < tmax - 1e-9:
        t1 = min(tmax, t0 + CHUNK)
        fine = np.arange(round(t0 / FINE), round(t1 / FINE) + 1) * FINE
        te = np.unique(np.round(np.r_[fine, grid[(grid >= t0) & (grid <= t1)], oldgrid[(oldgrid >= t0) & (oldgrid <= t1)]], 12))
        te = te[(te >= t0) & (te <= t1)]; te[0] = t0; te[-1] = t1
        sol = solve_ivp(model.rhs, (t0, t1), y, method='DOP853', t_eval=te, rtol=spec.get('rtol', 1e-11),
                        atol=spec.get('atol', 1e-12))
        if not sol.success:
            raise RuntimeError(sol.message)
        Y = sol.y.T
        if spec.get('fd'):
            Y = fd_to_tangent(model, Y)
        ob = S.observables(tang, Y)
        idx = int(round(t0 / CHUNK))
        np.savez(out / ('chunk_%03d.npz' % idx), t=sol.t, a=ob['a'], H=ob['H'], NX=ob['NX'], absX=ob['absX'],
                 A=ob['A'], dN=ob['dN'], ks=np.asarray(spec['ks']), nfev=sol.nfev)
        y = sol.y[:, -1].copy(); t0 = t1
        np.savez(ck, t=t0, y=y)
        print('chunk', idx, 't', t1, 'nfev', sol.nfev, 'elapsed', round(time.time() - start, 1), flush=True)
        if time.time() - start > budget:
            print('PAUSE', flush=True); return
    print('DONE', job, flush=True)


if __name__ == '__main__':
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 500.0)
