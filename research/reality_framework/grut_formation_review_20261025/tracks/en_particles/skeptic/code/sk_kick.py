"""en_particles SKEPTIC: kicked twin runs at a DIFFERENT start time with an own driver and own pole analysis.

Driver (reproduction-type use of Codex's unmodified early3d.cpp, compiled into scratch): start from a regenerated
snapshot, add +-eps*eta to pi_psi (eta: complex white noise, own seed), evolve with the relabelled clock
(tau0=0, a0=1, C=0, so the engine's tau equals the elapsed conformal time; the Hamiltonian is autonomous),
snapshot every ds, reduce each snapshot to psi_k immediately.
Analysis (own code): shell-averaged retarded function G_n(D) = sum dpsi_k(D) eta_k^* / sum |eta_k|^2,
dpsi = (psi(+eps) - psi(-eps))/2. Poles by a matrix-pencil fit G(D) ~ sum_j c_j exp(z_j D) (no taper, no window);
branch weights Z_j = Re(i w_j c_j) for Im z_j = w_j (free field: 1/2 per branch). Also the en_particles transform
(cos^2 taper, rho = -2 Im G~) for a like-for-like comparison.
Usage: python3 -I sk_kick.py run <engine> <state.bin> <outdir> <eps> [T=5] [ds=0.05] [seed=11] [threads=2]
       python3 -I sk_kick.py analyse <outdir_plus> <outdir_minus> <out.json>
"""
import sys, os, json, struct, subprocess, glob
import numpy as np

G2, LX, LP = 1e-2, 1e-3, 3e-2


def read_state(path):
    with open(path, 'rb') as f:
        N, L, s, bits = struct.unpack('<iddi', f.read(24))
    return N, L, s, np.fromfile(path, dtype='<f8' if bits == 64 else '<f4', offset=24).reshape(6, N, N, N).astype(float)


def write_init(path, q, p, L):
    with open(path, 'wb') as f:
        f.write(b'GRUT3D1\n'); f.write(struct.pack('<id4d', q.shape[-1], L, 0., G2, LX, LP))
        np.asarray(q, '<f8').tofile(f); np.asarray(p, '<f8').tofile(f)


def run(engine, start, outdir, eps, T=5., ds=.05, seed=11, threads=2):
    eps = float(eps); T = float(T); ds = float(ds); seed = int(seed)
    os.makedirs(outdir, exist_ok=True)
    N, L, s0, a = read_state(start)
    rng = np.random.default_rng(seed)
    eta = (rng.normal(size=(N, N, N))+1j*rng.normal(size=(N, N, N)))/np.sqrt(2); eta -= eta.mean()
    a[3] += eps*eta.real; a[4] += eps*eta.imag
    init = os.path.join(outdir, 'init.bin'); write_init(init, a[:3], a[3:], L)
    pref = os.path.join(outdir, 'st')
    subprocess.run([engine, init, pref, '0.004', str(T), str(ds), str(ds), '4', str(threads), '64', '0', '1', '0'],
                   check=True, stdout=subprocess.DEVNULL)
    files = sorted(glob.glob(pref+'_s*.bin')); ss = []; psi = []
    for fpath in files:
        n_, L_, s_, b = read_state(fpath)
        if ss and abs(s_-ss[-1]) < 1e-9: os.unlink(fpath); continue
        ss.append(s_); psi.append(np.fft.fftn(b[0]+1j*b[1], norm='forward')); os.unlink(fpath)
    os.unlink(init)
    np.savez(os.path.join(outdir, 'resp.npz'), s=np.array(ss), psi=np.array(psi), eta_k=np.fft.fftn(eta, norm='forward'),
             eps=eps, N=N, L=L, s0=s0, start=start)
    print('run done', outdir, len(ss), 'samples, s0', s0)


def matrix_pencil(y, dt, K):
    """Fit y(t_n) = sum_j c_j exp(z_j t_n), n=0..M-1, with K poles (Hua-Sarkar)."""
    M = len(y); Lp = M//2
    Y = np.array([y[i:i+Lp+1] for i in range(M-Lp)])
    U, S, Vh = np.linalg.svd(Y, full_matrices=False)
    V = Vh.conj().T[:, :K]
    V1, V2 = V[:-1, :], V[1:, :]
    lam = np.linalg.eigvals(np.linalg.pinv(V1)@V2)
    z = np.log(lam)/dt
    A = np.exp(np.outer(np.arange(M)*dt, z)); c = np.linalg.lstsq(A, y, rcond=None)[0]
    resid = float(np.linalg.norm(A@c-y)/np.linalg.norm(y))
    return z, c, resid


def analyse(dp, dm, out):
    P = np.load(os.path.join(dp, 'resp.npz')); Mn = np.load(os.path.join(dm, 'resp.npz'))
    eps = float(P['eps']); s = P['s']; D = s-s[0]; ds = D[1]-D[0]; N = int(P['N']); L = float(P['L'])
    eta = P['eta_k'].ravel()
    d = ((P['psi']-Mn['psi'])/2).reshape(len(s), -1)/eps
    n = np.fft.fftfreq(N)*N; sh = np.rint(np.sqrt(n[:, None, None]**2+n[None, :, None]**2+n[None, None, :]**2)).astype(int).ravel()
    th = 2*np.pi*np.fft.fftfreq(N); dx = L/N; one = ((4/3)*(2-2*np.cos(th))-(2-2*np.cos(2*th))/12)/dx**2
    k2 = (one[:, None, None]+one[None, :, None]+one[None, None, :]).ravel()
    T = D[-1]; taper = np.cos(np.pi*D/(2*T))**2; trap = np.full(len(D), ds); trap[[0, -1]] = ds/2
    wg = np.linspace(-20, 20, 4001); dw = wg[1]-wg[0]; E = np.exp(-1j*np.outer(wg, D))
    rows = []
    for shell in range(1, 21):
        sel = sh == shell
        G = (d[:, sel]*np.conj(eta[sel])[None, :]).sum(axis=1)/(np.abs(eta[sel])**2).sum()
        kf = float(np.sqrt(k2[sel].mean()))
        # like-for-like transform
        rho = -2*(E@(G*taper*trap)).imag
        ip = int(np.argmax(np.where(wg > 0, rho, -np.inf))); im = int(np.argmax(np.where(wg < 0, -rho, -np.inf)))
        Zp = float((wg*rho)[wg > 0].sum()*dw/(2*np.pi)); Zm = float((wg*rho)[wg < 0].sum()*dw/(2*np.pi))
        # matrix pencil, K=4 and K=6
        mp = {}
        for K in (2, 4, 6):
            z, c, res = matrix_pencil(G, ds, K)
            Z = np.real(1j*z.imag*c)
            order = np.argsort(-np.abs(Z))
            mp[K] = dict(resid=res, poles=[dict(omega=float(z[i].imag), gamma=float(-z[i].real), Z=float(Z[i])) for i in order])
        main = mp[4]['poles']
        pos = [p for p in main if p['omega'] > 0]; neg = [p for p in main if p['omega'] < 0]
        wp = pos[0]['omega'] if pos else None; wm = neg[0]['omega'] if neg else None
        rows.append(dict(shell=shell, k=2*np.pi/L*shell, k_fd4=kf, modes=int(sel.sum()),
                         G0=float(abs(G[0])), Gprime0=float(((G[1]-G[0])/ds).real),
                         rho_peak_pos=float(wg[ip]), rho_peak_neg=float(wg[im]), Z_plus_rho=Zp, Z_minus_rho=Zm,
                         m2_rho_peak_kfd4=float(wg[ip]**2-kf**2),
                         pencil=mp, pencil_omega_pos=wp, pencil_omega_neg=wm,
                         m2_pencil_pos=(wp**2-kf**2) if wp else None, m2_pencil_neg=(wm**2-kf**2) if wm else None,
                         Z_pencil_sum_K4=float(sum(p['Z'] for p in main))))
        r = rows[-1]
        print('%2d k=%.2f kf=%.2f rho+ %.2f rho- %.2f Z+ %.3f Z- %.3f | pencil w+ %s w- %s m2+ %s Zsum %.3f res4 %.1e res6 %.1e' % (
            shell, r['k'], kf, r['rho_peak_pos'], r['rho_peak_neg'], Zp, Zm,
            '%.3f' % wp if wp else '-', '%.3f' % wm if wm else '-', '%.1f' % r['m2_pencil_pos'] if wp else '-',
            r['Z_pencil_sum_K4'], mp[4]['resid'], mp[6]['resid']), flush=True)
    json.dump(dict(source='en_particles_skeptic own driver/analysis; reproduction-type engine runs', eps=eps,
                   s0=float(P['s0']), T=float(T), start=str(P['start']), shells=rows), open(out, 'w'), indent=1)


if __name__ == '__main__':
    if sys.argv[1] == 'run': run(*sys.argv[2:])
    else: analyse(*sys.argv[2:])
