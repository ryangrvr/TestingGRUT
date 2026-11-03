"""Skeptic S1: independent exact mode-function solver for the homogeneous control of lf_fermion.
Claude, lf_fermion skeptic (BRIEF_LF.md; issue #13 context).  Model choices, not validations.

Independent of lf1d / lfline in method and code:
  * per comoving momentum k, the 2x2 problem i d psi/dt = (1/a) [[k, c], [c*, -k]] psi, c = a M e^{i theta(t)},
    integrated in COSMIC time with the 4th-order two-point Gauss Magnus integrator (exact SU(2) exponential of the
    Magnus generator, including the commutator term), not Strang splitting in conformal time;
  * no FFT, no line grid: the homogeneous problem is diagonal in k;
  * background: line mean of the recorded lines.bin (read-only, same file lf_fermion used); theta = v X / f between
    records by cubic Hermite in t (values and cosmic-time derivatives) [model choice; option 'spline' uses a cubic
    spline through the values only, to test the interpolation];
  * measurement: positive-energy projection at the instantaneous phase, eps_3D = (1/2 pi^2) sum k^2 w n dk / a^4
    (same definition as lf1d.eps3d), with and without |k| > 3 a M; baseline = mean over t < 0.15 (as closure_cmp).
Usage: python3 -B sk_hom_exact.py RUN LINE SPECIES(chi|sig) [nsub_factor] [interp] [kmax_over_dk]
"""
import json, math, sys, time
from pathlib import Path
import numpy as np

SPEC = {'chi': dict(col=2, f=2000.0, Ms=[400.0, 200.0, 120.0], names=['N2', 'N3', 'N4']),
        'sig': dict(col=0, f=1200.0, Ms=[250.0], names=['N1'])}
SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2')
RUNS = SCR / 'lf_fermion' / 'runs'
OUT = Path(__file__).resolve().parents[1] / 'results'


def read_lines(rd):
    # independent reader of the documented lines.bin layout (int32 N, nl; int32 ax[nl], i1[nl], i2[nl]; then records
    # of float64 t, a, H, F[nl, N, 4])
    with open(rd / 'lines.bin', 'rb') as f:
        N, nl = np.fromfile(f, dtype=np.int32, count=2)
        np.fromfile(f, dtype=np.int32, count=3 * nl)
        raw = np.fromfile(f, dtype=np.float64)
    rec = 3 + nl * N * 4
    raw = raw[:(raw.size // rec) * rec].reshape(-1, rec)
    t = raw[:, 0]
    keep = np.concatenate([[True], np.diff(t) > 0])       # drop duplicated records at restarts
    raw = raw[keep]
    return raw[:, 0], raw[:, 1], raw[:, 2], raw[:, 3:].reshape(-1, nl, N, 4)


def su2_exp(g):
    """exp(-i g.sigma) acting as 2x2 matrix, g real (3, ...)."""
    n = np.sqrt(g[0] ** 2 + g[1] ** 2 + g[2] ** 2)
    cs, sn = np.cos(n), np.sin(n)
    s = np.where(n > 0, sn / np.where(n > 0, n, 1), 1.0)
    # -i g.sigma = -i (gx sx + gy sy + gz sz)
    U00 = cs - 1j * s * g[2]
    U11 = cs + 1j * s * g[2]
    U01 = -1j * s * (g[0] - 1j * g[1])
    U10 = -1j * s * (g[0] + 1j * g[1])
    return U00, U01, U10, U11


def hvec(k, a, M, th):
    """H_t = (1/a)[[k, c],[c*, -k]] = h.sigma with c = a M e^{i th}: h = (M cos th, -M sin th, k / a)."""
    return np.array([np.broadcast_to(M * np.cos(th), k.shape), np.broadcast_to(-M * np.sin(th), k.shape), k / a])


def pos_spinor(k, c):
    w = np.sqrt(k * k + np.abs(c) ** 2)
    nu = np.sqrt(np.abs(c) ** 2 + (w - k) ** 2)
    u0, u1 = c / nu, (w - k) / nu
    bad = nu < 1e-300
    u0 = np.where(bad, 1.0, u0)
    u1 = np.where(bad, 0.0, u1)
    return u0, u1, w


def neg_spinor(k, c):
    w = np.sqrt(k * k + np.abs(c) ** 2)
    nv = np.sqrt(np.abs(c) ** 2 + (w + k) ** 2)
    v0, v1 = c / nv, -(w + k) / nv
    bad = nv < 1e-300
    return np.where(bad, 1.0, v0), np.where(bad, 0.0, v1)


def main():
    run, line, sp = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    fac = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
    interp = sys.argv[5] if len(sys.argv) > 5 else 'hermite'
    kmul = int(sys.argv[6]) if len(sys.argv) > 6 else 4096
    S = SPEC[sp]
    rd = RUNS / run
    info = json.loads((rd / 'info.json').read_text())
    v, L = info['P']['v'], info['L']
    t, a, H, F = read_lines(rd)
    X = F[:, line, :, S['col']].mean(axis=1)
    Xd = F[:, line, :, S['col'] + 1].mean(axis=1)
    vf = v / S['f']
    th, thd = vf * X, vf * Xd
    if interp == 'spline':
        from scipy.interpolate import CubicSpline
        csp = CubicSpline(t, th)
    dk = 2 * math.pi / L
    k = dk * np.arange(-kmul, kmul)
    Ms = np.array(S['Ms'])[:, None]
    kk = np.broadcast_to(k, (len(S['Ms']), len(k)))
    c0 = a[0] * Ms * np.exp(1j * th[0])
    p0, p1 = neg_spinor(kk, c0)
    p0 = p0.astype(complex); p1 = p1.astype(complex)
    meas = []
    t_end = float(info.get('t_end', 1.5))
    t_stop = min(1.1 + 1e-9, t[-1])
    t0w = time.time()
    g1, g2 = 0.5 - math.sqrt(3) / 6, 0.5 + math.sqrt(3) / 6
    next_m = t[0]
    for n in range(len(t) - 1):
        if t[n] > t_stop:
            break
        h = t[n + 1] - t[n]
        # substeps so that the per-step rotation angle stays <= 0.25 rad
        rate = max(abs(thd[n]), abs(thd[n + 1]), np.abs(k).max() / min(a[n], a[n + 1]), Ms.max())
        nsub = max(1, int(math.ceil(fac * rate * h / 0.25)))
        hs = h / nsub
        for q in range(nsub):
            taus = [(q + g1) / nsub, (q + g2) / nsub]
            hv = []
            for tau in taus:
                if interp == 'spline':
                    thm = float(csp(t[n] + tau * h))
                else:
                    h00, h10 = 2 * tau ** 3 - 3 * tau ** 2 + 1, tau ** 3 - 2 * tau ** 2 + tau
                    h01, h11 = -2 * tau ** 3 + 3 * tau ** 2, tau ** 3 - tau ** 2
                    thm = h00 * th[n] + h10 * h * thd[n] + h01 * th[n + 1] + h11 * h * thd[n + 1]
                am = a[n] * (a[n + 1] / a[n]) ** tau          # a interpolated geometrically
                hv.append(hvec(kk, am, Ms, thm))
            h1, h2 = hv
            gen = 0.5 * hs * (h1 + h2) + (math.sqrt(3) / 6) * hs * hs * np.cross(h2, h1, axis=0)
            U00, U01, U10, U11 = su2_exp(gen)
            p0, p1 = U00 * p0 + U01 * p1, U10 * p0 + U11 * p1
        if t[n + 1] >= next_m - 1e-9:
            an = a[n + 1]
            c = an * Ms * np.exp(1j * th[n + 1])
            u0, u1, w = pos_spinor(kk, c)
            nk = np.abs(np.conj(u0) * p0 + np.conj(u1) * p1) ** 2
            wt = kk ** 2 * w * dk / (2 * math.pi ** 2) / an ** 4
            cut = np.abs(kk) > 3 * an * Ms
            row = dict(t=float(t[n + 1]), a=float(an))
            for j, nm in enumerate(S['names']):
                row[nm] = dict(e=float(np.sum(wt[j] * nk[j])), e_cut=float(np.sum((wt[j] * nk[j])[cut[j]])),
                               nmax=float(nk[j].max()), k_at_half=float(np.abs(k)[nk[j] > 0.5].max()) if np.any(nk[j] > 0.5) else 0.0)
            meas.append(row)
            next_m += 0.01
    tm = np.array([r['t'] for r in meas])
    res = dict(run=run, line=line, sp=sp, fac=fac, interp=interp, kmax=float(k.max()), v=v, meas=meas,
               wall_s=time.time() - t0w, label='C', note=__doc__)
    for nm in S['names']:
        e = np.array([r[nm]['e_cut'] for r in meas])
        b = e[tm < 0.15].mean()
        res.setdefault('Q_cut_minus_base_over_v2', {})[nm] = list((e - b) / v ** 2)
    res['t'] = list(tm)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / ('sk_hom_%s_l%d_%s_f%g_%s.json' % (run, line, sp, fac, interp))).write_text(json.dumps(res, indent=0))
    print('done', run, line, sp, fac, interp, 'wall %.0f s' % res['wall_s'])


if __name__ == '__main__':
    main()
