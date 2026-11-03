"""lfline: evolve the axial singlets on one recorded lattice line (non-back-reacting), with Fermi statistics.
Claude, lf_fermion track (issue #13 context).  Uses lf1d.Line (regressed in regress_homog.py).

Background: sigma, chi and their cosmic-time derivatives on a line of the 64^3 lf3d run, at every lattice step
(lines.bin).  theta = v phi / f.  Between records: cubic Hermite in t (values and derivatives), Fourier (band-limited)
interpolation from the 64 lattice sites to the fermion grid.  Comoving coordinates, conformal time eta, mass a(t) M.
Grid: N_c sites until t_ref, then N_f (zero-padded in k: the extra modes are filled with fresh vacuum, which is exact
while the band stays below the coarse Nyquist; checked by the occupation near the old Nyquist).
Time step: chosen per record interval so that deta (2 k_Nyq + max|a theta_dot|) <= cfl * 2 pi (no aliased
resonance of the discrete-time map; cfl = 0.5 default).
Measurements (every dt_meas): per species the 3D-weighted energy density (physical, H^4 units) with
  (i) 'proj': positive-energy projection at the line-mean phase (exact for a homogeneous background);
  (ii) 'chir': chiral occupation (R for k > 0, L for k < 0) minus its vacuum value |c|^2/(|c|^2 + (w + |k|)^2), which
       does not depend on theta (exact for |k| >> aM, used for the inhomogeneous comparison);
  both also with the cut |k| > 3 a M.  Plus n(k) binned, and the occupation near the coarse Nyquist.
Usage: python3 -B lfline.py RUN LINE SPECIES(chi|sig) NR SEED T_END [N_c N_f t_ref cfl tag wall]
A tag starting with 'hom' replaces the line by its line mean (homogeneous control).
Writes scratch/lfq/<RUN>/<tag>_l<LINE>_<SPECIES>.npz (state checkpoint) and .json (measurements); resumable.
"""
import json, math, sys, time
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
import lf1d, lines_io, lfcommon as LC

V = 367046.6168416639           # v in design H units (P['v'] of the run; the 8.4e10 run reads its own)
SPEC = {'chi': dict(col=2, f=2000.0, Ms=[400.0, 200.0, 120.0], names=['N2', 'N3', 'N4']),
        'sig': dict(col=0, f=1200.0, Ms=[250.0], names=['N1'])}


def finterp(c, Nf):
    """band-limited interpolation of real periodic samples c[..., n] to Nf points."""
    n = c.shape[-1]
    C = np.fft.rfft(c, axis=-1)
    out = np.zeros(c.shape[:-1] + (Nf // 2 + 1,), dtype=complex)
    out[..., :n // 2] = C[..., :n // 2]
    out[..., n // 2] = 0.5 * C[..., n // 2]
    return np.fft.irfft(out, n=Nf, axis=-1) * (Nf / n)


def measure(ln, a, theta_bar):
    res = {}
    n_proj = ln.occupation(a, theta_bar)
    k = ln.k
    for j, M in enumerate(ln.Ms):
        c2 = (a * M) ** 2
        w = np.sqrt(k * k + c2)
        vac = c2 / (c2 + (w + np.abs(k)) ** 2)
        R = np.mean(np.abs(ln.psik[j, :, 0, :]) ** 2, axis=0)
        L = np.mean(np.abs(ln.psik[j, :, 1, :]) ** 2, axis=0)
        n_ch = np.where(k > 0, R, L) - vac
        wt = k * k * w * ln.dk / (2 * np.pi ** 2) / a ** 4
        cut = np.abs(k) > 3 * a * M
        res[j] = dict(e_proj=float(np.sum(wt * n_proj[j])), e_proj_cut=float(np.sum((wt * n_proj[j])[cut])),
                      e_chir=float(np.sum(wt * n_ch)), e_chir_cut=float(np.sum((wt * n_ch)[cut])),
                      n_ch=n_ch)
    return res


def main():
    run, line, sp = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    nr, seed, t_end = int(sys.argv[4]), int(sys.argv[5]), float(sys.argv[6])
    Nc = int(sys.argv[7]) if len(sys.argv) > 7 else 4096
    Nf = int(sys.argv[8]) if len(sys.argv) > 8 else 8192
    t_ref = float(sys.argv[9]) if len(sys.argv) > 9 else 0.34
    cfl = float(sys.argv[10]) if len(sys.argv) > 10 else 0.5
    tag = sys.argv[11] if len(sys.argv) > 11 else 'base'
    wall = float(sys.argv[12]) if len(sys.argv) > 12 else 540.0
    dt_meas = 0.01
    S = SPEC[sp]
    rd = LC.RUNS / run
    info = json.loads((rd / 'info.json').read_text())
    v = info['P']['v']
    D = lines_io.read_lines(rd)
    t, a_r = D['t'], D['a']
    X = D['F'][:, line, :, S['col']]
    Xd = D['F'][:, line, :, S['col'] + 1]
    vf = v / S['f']
    if tag.startswith('hom'):
        # homogeneous control: the line-mean field and its derivative (same history, no inhomogeneity)
        X = np.repeat(X.mean(axis=1, keepdims=True), X.shape[1], axis=1)
        Xd = np.repeat(Xd.mean(axis=1, keepdims=True), Xd.shape[1], axis=1)
    out_dir = LC.SP / 'lfq' / run
    out_dir.mkdir(parents=True, exist_ok=True)
    base = out_dir / ('%s_l%d_%s' % (tag, line, sp))
    L = info['L']
    t0w = time.time()
    if base.with_suffix('.npz').exists():
        z = np.load(base.with_suffix('.npz'), allow_pickle=True)
        irec, N = int(z['irec']), int(z['N'])
        ln = lf1d.Line(N, L, S['Ms'], n_r=1, phases='ones')
        ln.nr = nr
        ln.psik = z['psik']
        meas = json.loads(base.with_suffix('.json').read_text())['meas']
        nsteps = int(z['nsteps'])
    else:
        irec = 0
        N = Nc
        th0 = float(np.mean(vf * X[0]))
        ln = lf1d.Line(N, L, S['Ms'], n_r=nr, seed=seed, theta0=th0, a0=a_r[0], dtype=np.complex128)
        meas = []
        nsteps = 0
    th_c = finterp(vf * X[irec], N)
    thd_c = finterp(vf * Xd[irec], N)
    next_meas = (math.floor(t[irec] / dt_meas + 1e-9) + 1) * dt_meas if meas else t[irec]
    while irec < len(t) - 1 and t[irec] < t_end - 1e-12:
        if N == Nc and t[irec] >= t_ref and Nf > Nc:
            # refine: zero-pad in k, new modes filled with vacuum
            kn = 2 * np.pi * np.fft.fftfreq(Nf, d=L / Nf)
            new = np.zeros(ln.psik.shape[:-1] + (Nf,), dtype=ln.psik.dtype)
            ko = ln.k
            idx = np.round(ko / ln.dk).astype(int) % Nf
            new[..., idx] = ln.psik
            rng = np.random.default_rng(seed + 1000)
            fresh = np.ones(Nf, bool)
            fresh[idx] = False
            for j, M in enumerate(ln.Ms):
                vv, _, _ = lf1d.spinors(kn, a_r[irec] * M * np.exp(1j * float(np.mean(th_c))))
                for r in range(ln.psik.shape[1]):
                    xi = np.exp(2j * np.pi * rng.random(Nf))
                    new[j, r][:, fresh] = (vv * xi)[:, fresh]
            # amplitude convention: |psi_k|^2 = 1 per filled mode on either grid (the ifft/fft pair is unnormalised
            # consistently), so no rescaling is needed
            ln.N, ln.k, ln.psik = Nf, kn, new
            N = Nf
            th_c = finterp(vf * X[irec], N)
            thd_c = finterp(vf * Xd[irec], N)
        th_n = finterp(vf * X[irec + 1], N)
        thd_n = finterp(vf * Xd[irec + 1], N)
        h = t[irec + 1] - t[irec]
        a0_, a1_ = a_r[irec], a_r[irec + 1]
        kny = math.pi * N / L
        thp = max(a0_ * np.abs(thd_c).max(), a1_ * np.abs(thd_n).max())
        amid = 0.5 * (a0_ + a1_)
        deta_rec = h / amid
        nsub = max(1, int(math.ceil(deta_rec * (2 * kny + thp) / (cfl * 2 * math.pi))))
        for q in range(nsub):
            tau = (q + 0.5) / nsub
            h00, h10 = 2 * tau ** 3 - 3 * tau ** 2 + 1, tau ** 3 - 2 * tau ** 2 + tau
            h01, h11 = -2 * tau ** 3 + 3 * tau ** 2, tau ** 3 - tau ** 2
            thm = h00 * th_c + h10 * h * thd_c + h01 * th_n + h11 * h * thd_n
            am = a0_ + (a1_ - a0_) * tau
            ln.step(0.0, deta_rec / nsub, am, thm)
        nsteps += nsub
        irec += 1
        th_c, thd_c = th_n, thd_n
        if t[irec] >= next_meas - 1e-9 or irec == len(t) - 1:
            r = measure(ln, a_r[irec], float(np.mean(th_c)))
            k = ln.k
            edges = np.array([0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5]) * kny
            row = dict(t=float(t[irec]), a=float(a_r[irec]), N=N, nsteps=nsteps)
            for j, nm in enumerate(S['names']):
                d = r[j]
                row[nm] = dict(e_proj=d['e_proj'], e_proj_cut=d['e_proj_cut'], e_chir=d['e_chir'], e_chir_cut=d['e_chir_cut'],
                               n_kbin=[float(np.mean(d['n_ch'][(np.abs(k) >= lo) & (np.abs(k) < hi)])) if np.any((np.abs(k) >= lo) & (np.abs(k) < hi)) else 0.0
                                       for lo, hi in zip(edges[:-1], edges[1:])],
                               n_edge=float(np.mean(d['n_ch'][np.abs(k) > 0.9 * kny])))
            meas.append(row)
            next_meas += dt_meas
            if time.time() - t0w > wall:
                break
    np.savez(base.with_suffix('.npz'), psik=ln.psik, irec=irec, N=N, nsteps=nsteps)
    base.with_suffix('.json').write_text(json.dumps(dict(run=run, line=line, sp=sp, nr=nr, seed=seed, Nc=Nc, Nf=Nf, t_ref=t_ref,
                                                         cfl=cfl, L=L, meas=meas), indent=0))
    print('%s line %d %s t=%.4f N=%d steps=%d wall=%.0fs' % (run, line, sp, t[irec], N, nsteps, time.time() - t0w), flush=True)
    return 0 if t[irec] >= t_end - 1e-9 or irec == len(t) - 1 else 3


if __name__ == '__main__':
    sys.exit(main())
