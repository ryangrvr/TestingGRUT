"""en_particles SKEPTIC: where in k does the late-time nonzero-mode X-bar sit, and does that depend on the box?

Own driver around Codex's unmodified engine (compiled into scratch; relabelled clock tau0=s_start, a0=1, C=0).
From a start state, evolve to s_a without output, then sample every ds up to s_b in chunks; for every sample
accumulate per-mode time averages <|psi_k|^2>, <|pi_k|^2>, <Im psi_k^* pi_k> (uniform weights) and keep the
zero-mode time series. Output: per-shell virial counts n_- = 1/2[sqrt(<|pi|^2><|psi|^2>) - <q>] per net charge
(per-mode, then summed per continuum shell), the zero-mode Hann spectral split, and the box-scale fractions.
Usage: python3 -I sk_box.py run <engine> <start.bin> <workdir> <s_a> <s_b> <ds> <out.json> [threads=2] [chunk=1.0]
       python3 -I sk_box.py window <window.npz> <i0> <out.json>      (same reduction on an existing L12 window)
"""
import sys, os, json, struct, subprocess, glob
import numpy as np

G2, LX, LP = 1e-2, 1e-3, 3e-2


def read_state(path):
    with open(path, 'rb') as f:
        N, L, s, bits = struct.unpack('<iddi', f.read(24))
    return N, L, s, np.fromfile(path, dtype='<f8' if bits == 64 else '<f4', offset=24).reshape(6, N, N, N).astype(float)


def write_init(path, a, L):
    with open(path, 'wb') as f:
        f.write(b'GRUT3D1\n'); f.write(struct.pack('<id4d', a.shape[-1], L, 0., G2, LX, LP)); np.asarray(a, '<f8').tofile(f)


class Acc:
    def __init__(self, N, L):
        self.N, self.L = N, L; self.A = 0.; self.B = 0.; self.q = 0.; self.n = 0; self.zero = []; self.s = []; self.Q = []

    def add(self, s, a):
        psi = np.fft.fftn(a[0]+1j*a[1], norm='forward'); pi = np.fft.fftn(a[3]+1j*a[4], norm='forward')
        self.A = self.A+np.abs(psi)**2; self.B = self.B+np.abs(pi)**2; qk = np.imag(np.conj(psi)*pi); self.q = self.q+qk
        self.n += 1; self.zero.append(psi[0, 0, 0]); self.s.append(s); self.Q.append(float(qk.sum()))

    def result(self, tag):
        N, L = self.N, self.L; A = self.A/self.n; B = self.B/self.n; q = self.q/self.n; Q = float(np.mean(self.Q))
        nv = .5*(np.sqrt(A*B)-q)
        n = np.fft.fftfreq(N)*N; sh = np.rint(np.sqrt(n[:, None, None]**2+n[None, :, None]**2+n[None, None, :]**2)).astype(int)
        kk = 2*np.pi/L*sh; nz = sh > 0
        z = np.array(self.zero); ss = np.array(self.s); M = len(z); h = np.hanning(M); ds = ss[1]-ss[0]
        P = np.abs(np.fft.fft(h*z))**2/(M*(h**2).sum()); w = 2*np.pi*np.fft.fftfreq(M, ds)
        shells = []
        for m in range(1, int(sh.max())+1):
            sel = sh == m
            if sel.any(): shells.append(dict(shell=m, k=float(2*np.pi/L*m), modes=int(sel.sum()), nminus=float(nv[sel].sum()/Q)))
        tot = float(nv[nz].sum()/Q)
        return dict(tag=tag, N=N, L=L, s=[float(ss[0]), float(ss[-1])], samples=M, Q=Q,
                    nminus_nonzero_virial=tot, nminus_zero_virial=float(nv[~nz].sum()/Q),
                    nminus_zero_hann=float((np.abs(w)*P)[w < 0].sum()/Q),
                    frac_k_le_0p6=float(nv[nz & (kk <= .6)].sum()/Q)/tot, frac_k_le_1p1=float(nv[nz & (kk <= 1.1)].sum()/Q)/tot,
                    frac_k_le_1p6=float(nv[nz & (kk <= 1.6)].sum()/Q)/tot, frac_lowest_shell=shells[0]['nminus']/tot,
                    lowest_shell_k=shells[0]['k'], shells=shells[:30])


def run(engine, start, work, sa, sb, ds, out, threads=2, chunk=1.0):
    sa, sb, ds, chunk = float(sa), float(sb), float(ds), float(chunk); os.makedirs(work, exist_ok=True)
    N, L, s0, a = read_state(start)
    init = os.path.join(work, 'init.bin'); pref = os.path.join(work, 'st')
    if sa > s0+1e-9:  # pre-evolution, only the final snapshot is kept
        write_init(init, a, L)
        subprocess.run([engine, init, pref, '0.004', str(sa-s0), str(sa-s0), str(sa-s0), '4', str(threads), '64', str(s0), '1', '0'],
                       check=True, stdout=subprocess.DEVNULL)
        files = sorted(glob.glob(pref+'_s*.bin')); _, _, t, a = read_state(files[-1])
        for fpath in files: os.unlink(fpath)
        print('pre-evolved to s=%.5f' % (s0+t), flush=True)
    acc = Acc(N, L); acc.add(sa, a); cur = sa
    while cur < sb-1e-9:
        D = min(chunk, sb-cur); write_init(init, a, L)
        subprocess.run([engine, init, pref, '0.004', str(D), str(ds), str(ds), '4', str(threads), '64', str(cur), '1', '0'],
                       check=True, stdout=subprocess.DEVNULL)
        files = sorted(glob.glob(pref+'_s*.bin')); last = 0.
        for fpath in files[1:]:
            _, _, t, b = read_state(fpath)
            if t > last+1e-9: acc.add(cur+t, b); last = t; a = b
            os.unlink(fpath)
        os.unlink(files[0]); cur += D
        print('chunk to s=%.4f samples=%d' % (cur, acc.n), flush=True)
    os.unlink(init)
    res = acc.result('engine run from '+os.path.basename(start)); json.dump(res, open(out, 'w'), indent=1); report(res)


def window(path, i0, out):
    d = np.load(path); s = d['s']; i0 = int(i0); N = int(d['N']); L = float(d['L'])
    acc = Acc(N, L)
    for j in range(i0, len(s)):
        psi = d['psi'][j].astype(complex); pi = d['ppsi'][j].astype(complex)
        acc.A = acc.A+np.abs(psi)**2; acc.B = acc.B+np.abs(pi)**2; qk = np.imag(np.conj(psi)*pi); acc.q = acc.q+qk
        acc.n += 1; acc.zero.append(psi[0, 0, 0]); acc.s.append(float(s[j])); acc.Q.append(float(qk.sum()))
    res = acc.result('window '+path.split('/')[-2]); json.dump(res, open(out, 'w'), indent=1); report(res)


def report(r):
    print('%s N%d L%g s %.2f-%.2f: n- nonzero %.4f zero(vir) %.4f zero(hann) %.4f | lowest shell k=%.3f frac %.3f, k<=0.6 %.3f, k<=1.1 %.3f, k<=1.6 %.3f' % (
        r['tag'], r['N'], r['L'], r['s'][0], r['s'][1], r['nminus_nonzero_virial'], r['nminus_zero_virial'], r['nminus_zero_hann'],
        r['lowest_shell_k'], r['frac_lowest_shell'], r['frac_k_le_0p6'], r['frac_k_le_1p1'], r['frac_k_le_1p6']))
    print('   shells:', ' '.join('%.2f:%.3f' % (x['k'], x['nminus']) for x in r['shells'][:12]))


if __name__ == '__main__':
    if sys.argv[1] == 'run': run(*sys.argv[2:])
    else: window(*sys.argv[2:])
