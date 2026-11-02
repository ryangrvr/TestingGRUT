"""qm_micro item 1b driver (Claude, 2 Nov 2026): wave packets of X / Xbar on the planar M1 Q-matter slab, classical
lattice (lattice.c), full X-phi potential with the M1 tilt.  Own code.

Setup (1D: Ny = 1; quasi-1D: Ny > 1 periodic in y with one transverse mode k_perp = 2 pi / (Ny dy)):
  * the static interface for the SAME second-order stencil (mirror ends), by sparse Newton with Omega free and the
    pin F(0) = F_B/2 (Omega comes out at Omega_b);
  * background run: psi = F e^{-i Omega t}, S static;  packet run: background + pure positive-frequency (X) or
    negative-frequency (Xbar) Gaussian packet of peak amplitude eps in the vacuum, moving toward the slab;
  * analysis on the difference (packet - background) at the end: in the vacuum region x > x_m the free-field
    decomposition into X, Xbar and phi quanta with the lattice dispersion; in x < x_m the excess charge and the
    rotating-frame energy (E - Omega Q) deposited in the slab;
  * prediction: the frequency-domain solver scatter1d.solve(scheme='fd2', h = dx) on the same discrete profile
    (exact image of the lattice in space), averaged over the packet's incident number spectrum.
Usage: python3 lat_driver.py <out.json> <case-set>
"""
import json, math, os, subprocess, sys, time
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import qm_common as qc
import scatter1d as sc

SP = os.environ.get('QM_SCR', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/qm_micro')
LAT = os.path.join(SP, 'lattice')
a_ = qc.a_


def discrete_profile(Lb, Lv, dx):
    B = qc.bulk()
    x = -Lb + dx * np.arange(int(round((Lb + Lv) / dx)) + 1)
    n = x.size
    P = qc.planar_profile(Lb=60, Lv=60, h=0.01)
    F = np.interp(x, P['x'], P['F'], left=B['F'], right=0.0)
    S = np.interp(x, P['x'], P['S'], left=B['S'], right=1.0)
    Om = B['Omega']
    ipin = int(np.argmin(np.abs(x)))
    i2 = 1 / dx ** 2
    for it in range(60):
        Fm = np.r_[F[1], F[:-1]]; Fp = np.r_[F[1:], F[-2]]
        Sm = np.r_[S[1], S[:-1]]; Spp = np.r_[S[1:], S[-2]]
        RF = (Fp + Fm - 2 * F) * i2 - (S * S + a_ * F * F - Om * Om) * F
        RS = (Spp + Sm - 2 * S) * i2 - qc.Uphi_p(S) - F * F * S
        R = np.r_[RF, RS, F[ipin] - 0.5 * B['F']]
        nr = np.max(np.abs(R))
        if nr < 1e-13:
            break
        idx = np.arange(n)
        rows, cols, vals = [], [], []
        def add(r, c, v):
            rows.append(r); cols.append(c); vals.append(v)
        # Laplacian with mirror ends
        for off in (0, n):
            add(idx + off, idx + off, np.full(n, -2 * i2))
            up = idx[:-1]; add(up + off, up + 1 + off, np.full(n - 1, i2))
            dn = idx[1:]; add(dn + off, dn - 1 + off, np.full(n - 1, i2))
            add(np.array([off]), np.array([1 + off]), np.array([i2]))           # mirror at 0: F_{-1} = F_1
            add(np.array([n - 1 + off]), np.array([n - 2 + off]), np.array([i2]))
        add(idx, idx, -(S * S + 3 * a_ * F * F - Om * Om))
        add(idx, idx + n, -2 * S * F)
        add(idx + n, idx + n, -(qc.Uphi_pp(S) + F * F))
        add(idx + n, idx, -2 * F * S)
        add(idx, np.full(n, 2 * n), 2 * Om * F)
        add(np.array([2 * n]), np.array([ipin]), np.array([1.0]))
        J = sp.csc_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(2 * n + 1, 2 * n + 1))
        d = spsolve(J, -R)
        F, S, Om = F + d[:n], S + d[n:2 * n], Om + d[2 * n]
    return dict(x=x, F=F, S=S, Omega=Om, Omega_b=Om, resid=nr, its=it, S_B=B['S'], F_B=B['F'], n_B=B['n'])


def lat_k2(k, d):
    return (2 / d * np.sin(k * d / 2)) ** 2


def free_decomp(dpsi, dvel, dS, dvS, dx, kp2):
    """X, Xbar, phi quanta (rescaled units, per unit transverse length) of free fields on a region (1D arrays, the
    transverse mode already projected out)."""
    n = dpsi.size
    k = 2 * np.pi * np.fft.fftfreq(n, d=dx)
    w = np.sqrt(1 + lat_k2(k, dx) + kp2)
    wf = np.sqrt(6 + lat_k2(k, dx) + kp2)
    P = np.fft.fft(dpsi) * dx; V = np.fft.fft(dvel) * dx
    L = n * dx
    NX = np.sum(np.abs(w * P + 1j * V) ** 2 / (4 * w)) / L
    NXb = np.sum(np.abs(w * P - 1j * V) ** 2 / (4 * w)) / L
    Sk = np.fft.fft(dS) * dx; Vk = np.fft.fft(dvS) * dx
    Nphi = np.sum((np.abs(Vk) ** 2 + wf ** 2 * np.abs(Sk) ** 2) / (2 * wf)) / L
    specX = np.abs(w * P + 1j * V) ** 2 / (4 * w) / L
    specXb = np.abs(w * P - 1j * V) ** 2 / (4 * w) / L
    specphi = np.abs(wf * Sk + 1j * Vk) ** 2 / (2 * wf) / L
    # right-movers: X with k > 0; Xbar (psi ~ e^{+i w t + i q x}) with q < 0; phi with k > 0
    dirn = dict(X_right=float(specX[k > 0].sum()), X_left=float(specX[k < 0].sum()),
                Xbar_right=float(specXb[k < 0].sum()), Xbar_left=float(specXb[k > 0].sum()),
                phi_right=float(specphi[k > 0].sum()), phi_left=float(specphi[k < 0].sum()))
    return NX, NXb, Nphi, k, specX, specXb, w, dirn


def packet_geometry(T, theta_deg):
    E = 1 + T; p = math.sqrt(E * E - 1); kx = p * math.cos(math.radians(theta_deg))
    sigma = max(28.0, 6.0 / kx)
    x0 = max(220.0, 5 * sigma + 70.0)
    vmin = packet_vmin(kx, sigma, E)
    vmax = (kx + 3 / (2 * sigma)) / E
    t_end = (x0 + 5 * sigma + 60.0) / vmin
    Lb = round(0.5 * 0.95 * t_end + 120.0, -1)
    Lv = round(max(x0 + 5 * sigma, vmax * t_end - x0 + 5 * sigma) + 80.0, -1)
    return sigma, x0, t_end, Lb, Lv


def packet_vmin(kx, sigma, E):
    kmin = max(kx - 3 / (2 * sigma), 0.5 * kx)
    return kmin / math.sqrt(1 + kmin ** 2 + (E * E - 1 - kx * kx))


def run_case(prof, T, inc, eps, Ny=1, theta_deg=0.0, sigma=28.0, x0=220.0, dx=0.1, dt=0.05, label=''):
    t0 = time.time()
    x = prof['x']; Nx = x.size; Om = prof['Omega']
    E = 1 + T; p = math.sqrt(E * E - 1)
    th = math.radians(theta_deg)
    if Ny > 1:
        kperp = p * math.sin(th)
        dy = 2 * math.pi / (Ny * kperp)
        kp2 = float(lat_k2(kperp, dy))
    else:
        kperp, dy, kp2 = 0.0, 1.0, 0.0
    # carrier k_x on the lattice so that omega = E exactly
    kx2_hat = E * E - 1 - kp2
    kx = 2 / dx * math.asin(math.sqrt(kx2_hat) * dx / 2)
    vgx = (math.sin(kx * dx) / dx) / E                    # d omega / d k_x on the lattice
    # time to bring the whole reflected packet (down to its k - 3 sigma_k component) past x_m = 40
    xm = 40.0
    vmin = packet_vmin(kx, sigma, E)
    t_end = (x0 + 5 * sigma + xm + 20) / vmin
    nsteps = int(math.ceil(t_end / dt))
    # fields
    yv = dy * np.arange(Ny)
    F2 = np.repeat(prof['F'][:, None], Ny, 1); S2 = np.repeat(prof['S'][:, None], Ny, 1)
    pr0, pi0 = F2.copy(), np.zeros_like(F2)
    vr0, vi0 = np.zeros_like(F2), -Om * F2
    vS0 = np.zeros_like(F2)
    # packet: X moving left psi ~ e^{-i E t - i kx x + i kperp y}; Xbar moving left psi ~ e^{+i E t + i kx x + i kperp y}
    env = eps * np.exp(-(x - x0) ** 2 / (4 * sigma ** 2))
    sgn = -1 if inc == 'X' else +1
    # pure frequency sign via FFT along x with the lattice dispersion
    ph1 = np.exp(1j * sgn * kx * x)[:, None] * np.exp(1j * kperp * yv)[None, :]
    dpsi = env[:, None] * ph1
    k = 2 * np.pi * np.fft.fftfreq(Nx, d=dx)
    wk = np.sqrt(1 + lat_k2(k, dx) + kp2)
    Pk = np.fft.fft(dpsi, axis=0)
    dvel = np.fft.ifft(-1j * (1 if inc == 'X' else -1) * wk[:, None] * Pk, axis=0)
    def write(fn, pr, pi, S, vr, vi, vS):
        np.concatenate([a.ravel() for a in (pr, pi, S, vr, vi, vS)]).astype(np.float64).tofile(fn)
    tag = '%s_%s_T%g_e%g_Ny%d_th%g' % (label, inc, T, eps, Ny, theta_deg)
    fin_b, fout_b = os.path.join(SP, tag + '_bg.in'), os.path.join(SP, tag + '_bg.out')
    fin_p, fout_p = os.path.join(SP, tag + '_pk.in'), os.path.join(SP, tag + '_pk.out')
    write(fin_b, pr0, pi0, S2, vr0, vi0, vS0)
    write(fin_p, pr0 + dpsi.real, pi0 + dpsi.imag, S2, vr0 + dvel.real, vi0 + dvel.imag, vS0)
    outs = []
    for fi, fo in ((fin_b, fout_b), (fin_p, fout_p)):
        r = subprocess.run([LAT, fi, fo, str(Nx), str(Ny), repr(dx), repr(dy), repr(dt), str(nsteps)], capture_output=True, text=True, check=True)
        outs.append([float(v) for v in r.stdout.split()[1::2]])
    def load(fn):
        d = np.fromfile(fn, dtype=np.float64).reshape(6, Nx, Ny)
        return d
    B = load(fout_b); Pq = load(fout_p)
    D = Pq - B
    # project on the transverse mode
    proj = np.exp(-1j * kperp * yv)[None, :] / Ny
    dpsi_e = ((D[0] + 1j * D[1]) * proj).sum(1)
    dvel_e = ((D[3] + 1j * D[4]) * proj).sum(1)
    dS_e = (D[2] * proj).sum(1); dvS_e = (D[5] * proj).sum(1)
    mv = x > xm
    # initial incident numbers (same decomposition on the initial packet, whole vacuum side)
    ini = (x > xm)
    dpsi0 = (dpsi * proj).sum(1); dvel0 = (dvel * proj).sum(1)
    NX0, NXb0, _, _, specin, specinb, wk_in, _ = free_decomp(dpsi0[ini], dvel0[ini], np.zeros(ini.sum()), np.zeros(ini.sum()), dx, kp2)
    Nin = NX0 if inc == 'X' else NXb0
    NX, NXb, Nphi, kk, specX, specXb, _, dirn = free_decomp(dpsi_e[mv], dvel_e[mv], dS_e[mv].real, dvS_e[mv].real, dx, kp2)
    # charge deposited in the slab (x < x_m): the leapfrog-conserved staggered charge sum(pi vr' - pr vi'), v' = v + dt/2 a
    def accel(Fd):
        pr, pi, S = Fd[0], Fd[1], Fd[2]
        def lap(f):
            fm = np.r_[f[1:2], f[:-1]]; fp = np.r_[f[1:], f[-2:-1]]
            L = (fp + fm - 2 * f) / dx ** 2
            if Ny > 1:
                L = L + (np.roll(f, 1, 1) + np.roll(f, -1, 1) - 2 * f) / dy ** 2
            return L
        r2 = pr * pr + pi * pi
        m = S * S + a_ * r2
        return lap(pr) - m * pr, lap(pi) - m * pi
    wtr = np.ones((Nx, 1)); wtr[0] = wtr[-1] = 0.5        # trapezoid end weights: exact conservation with mirror ends
    def Qstag(Fd, sel):
        ar, ai = accel(Fd)
        vr = Fd[3] + 0.5 * dt * ar; vi = Fd[4] + 0.5 * dt * ai
        return np.sum((wtr * (Fd[1] * vr - Fd[0] * vi))[sel]) * dx * dy
    Qin_slab = Qstag(Pq, ~mv) - Qstag(B, ~mv)
    Qvac_end = Qstag(Pq, mv) - Qstag(B, mv)
    I0b = np.stack([pr0, pi0, S2, vr0, vi0, vS0]); I0p = np.stack([pr0 + dpsi.real, pi0 + dpsi.imag, S2, vr0 + dvel.real, vi0 + dvel.imag, vS0])
    allm = np.ones(Nx, bool)
    Qdiff0 = Qstag(I0p, allm) - Qstag(I0b, allm)
    Qdiff1 = Qstag(Pq, allm) - Qstag(B, allm)
    # BVP prediction on the same discrete profile and stencil, averaged over the incident spectrum
    sel = (prof['x'] >= -60) & (prof['x'] <= 60)
    Ip = sc.Interface(profile=dict(x=prof['x'][sel], F=prof['F'][sel], S=prof['S'][sel], Omega_b=Om), h=dx)
    Ip.h = dx
    kin = 2 * np.pi * np.fft.fftfreq(ini.sum(), d=dx)
    win = np.sqrt(1 + lat_k2(kin, dx) + kp2)
    spin = specin if inc == 'X' else specinb
    wts = spin / spin.sum()
    keep = (wts > 1e-6 * wts.max()) & (np.abs(kin) > 1e-9)
    Rp, Xp, Xbp, php, trp = 0.0, 0.0, 0.0, 0.0, 0.0
    for wv, wt in zip(win[keep], wts[keep]):
        r = sc.solve(Ip, (wv - Om) if inc == 'X' else (wv + Om), math.sqrt(kp2), inc, scheme='fd2')
        if r.get('closed'):
            continue
        Xp += wt * r['out'].get('out_X', 0.0); Xbp += wt * r['out'].get('out_Xbar', 0.0); php += wt * r['out'].get('out_phi', 0.0)
        trp += wt * r['trans_total']
    rc = sc.solve(Ip, (E - Om) if inc == 'X' else (E + Om), math.sqrt(kp2), inc, scheme='fd2')
    res = dict(label=label, inc=inc, T=T, E=E, eps=eps, Ny=Ny, theta_deg=theta_deg, kperp=kperp, kp_lattice=math.sqrt(kp2),
               kx=kx, vgx=vgx, sigma=sigma, x0=x0, dx=dx, dy=dy, dt=dt, nsteps=nsteps, t_end=nsteps * dt, Nx=Nx,
               N_in=Nin, frac_out_X=NX / Nin, frac_out_Xbar=NXb / Nin, frac_out_phi=Nphi / Nin,
               frac_dir={k_: v_ / Nin for k_, v_ in dirn.items()},
               slab_charge_over_Nin=Qin_slab / (Nin * (dy * Ny if Ny > 1 else 1.0)),
               vac_charge_over_Nin=Qvac_end / (Nin * (dy * Ny if Ny > 1 else 1.0)),
               staggered_charge_diff_start_end=[Qdiff0 / (Nin * (dy * Ny if Ny > 1 else 1.0)), Qdiff1 / (Nin * (dy * Ny if Ny > 1 else 1.0))],
               pred_packet=dict(out_X=Xp, out_Xbar=Xbp, out_phi=php, trans=trp),
               pred_carrier=dict(out=rc['out'], trans=rc['trans_total']),
               energy_drift_bg=(outs[0][2] - outs[0][0]) / abs(outs[0][0]), charge_drift_bg=(outs[0][3] - outs[0][1]) / abs(outs[0][1]),
               energy_drift_pk=(outs[1][2] - outs[1][0]) / abs(outs[1][0]), charge_drift_pk=(outs[1][3] - outs[1][1]) / abs(outs[1][1]),
               dQ_packet_run=(outs[1][3] - outs[1][1]), runtime_s=time.time() - t0)
    for f in (fin_b, fout_b, fin_p, fout_p):
        os.remove(f)
    return res


CASES = {
    'lin1d': [('X', T, 1e-3, 1, 0.0) for T in (0.005, 0.02, 0.1, 0.3, 1.0)] + [('Xbar', T, 1e-3, 1, 0.0) for T in (0.02, 0.1, 0.5)],
    'nl1d': [('X', 0.1, e, 1, 0.0) for e in (0.03, 0.1, 0.3, 0.6)] + [('X', 0.02, e, 1, 0.0) for e in (0.1, 0.3)] + [('Xbar', 0.1, e, 1, 0.0) for e in (0.1, 0.3)],
    'q1d': [('X', 0.1, 1e-3, 8, 30.0), ('X', 0.3, 1e-3, 8, 60.0), ('Xbar', 0.1, 1e-3, 8, 45.0)],
    'q1d_b': [('X', 0.1, 1e-3, 8, 60.0), ('X', 0.3, 0.1, 16, 60.0)],
}


def main(out, which):
    t0 = time.time()
    res = dict(notes='qm_micro lattice wave-packet runs; fractions are of the incident quanta number', cases=[])
    prev = None
    if prev:
        res = prev
    for inc, T, eps, Ny, th in CASES[which]:
        sigma, x0, t_end, Lb, Lv = packet_geometry(T, th)     # phonons reflected at the far end cannot return in time
        prof = discrete_profile(Lb, Lv, 0.1)
        r = run_case(prof, T, inc, eps, Ny=Ny, theta_deg=th, sigma=sigma, x0=x0, label=which)
        r['Lb'], r['Lv'], r['profile_Omega_minus_bulk'] = Lb, Lv, prof['Omega'] - qc.bulk()['Omega']
        res['cases'].append(r)
        print('%s %s T=%g eps=%g Ny=%d th=%g: out X %.6f Xbar %.6f phi %.2e (left-movers left %.1e) | pred X %.6f Xbar %.6f phi %.2e | slab Q %.6f | drift E %.1e Q %.1e | %.0fs' % (
            which, inc, T, eps, Ny, th, r['frac_dir']['X_right'], r['frac_dir']['Xbar_right'], r['frac_dir']['phi_right'], r['frac_dir']['X_left'] + r['frac_dir']['Xbar_left'], r['pred_packet']['out_X'],
            r['pred_packet']['out_Xbar'], r['pred_packet']['out_phi'], r['slab_charge_over_Nin'], r['energy_drift_pk'], r['charge_drift_pk'], r['runtime_s']), flush=True)
        json.dump(res, open(out, 'w'), indent=1, default=float)
    print('total %.0fs' % (time.time() - t0))


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
