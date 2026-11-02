#!/usr/bin/env python3
"""dq_heat (Claude, own code), item 2: finite-k mediator production during the M1r release.

Linear fluctuations of the hidden fields on the homogeneous M1r background, integrated in rpmodes.c together with the
background (rptri_dq.c dynamics, unchanged):
  delta r'' + 3 H delta r' + (k^2/a^2) delta r + Hess V(t) delta r = 0,  canonical real components of (Y, Sigma) [subsystem
  'YS', X frozen as background] or of (X, Y, Sigma) ['XYS'].  The Y mass matrix is time dependent through mu_eff(t) =
  mu + G phihat(t) (the holder), the Y and Sigma backgrounds, and (XYS) the 3 kap X^2 mixing.
Initial state: exact de Sitter Bunch-Davies for every eigenmode of the inflationary Hessian at the prepared point (the
background is static before the slice, H = H_i, so the eigenmodes decouple there):
  delta r_alpha(0) = (sqrt(pi)/2) e^{i(pi nu/2 + pi/4)} H^(1)_nu(k),  nu = (9/4 - m_alpha^2)^{1/2} (imaginary for heavy modes),
  delta r_alpha'(0) = -(sqrt(pi)/2) e^{i(...)} [3/2 H_nu(k) + k H_nu'(k)]   (a = 1, H = 1, -eta = 1 at the slice).
Read-out: instantaneous eigenmodes of k^2/a^2 + Hess V(t); n_beta(k) = sum_alpha a^3 (|e_beta . v'_alpha|^2 + w^2 |e_beta . v_alpha|^2)/(2 w) - 1/2;
N_beta = int dk k^2 n_beta / (2 pi^2) per comoving H_i^3 (the frame of N_X, N_YL).  Late eigenmodes are identified as
Y_L (M^2 - sqrt2 mu h = 495), Y_H (1305), sigma (2 lambda_Sigma h^2 = 1600), the Sigma phase (Goldstone, massless) and,
in XYS, the X modes.
Usage: python3 -B ymodes.py YS 100 [tag]  |  python3 -B ymodes.py XYS 40 [tag]
"""
import json, subprocess, sys, time
from pathlib import Path
import numpy as np
import mpmath as mp
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dq_model as fc

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
SRC = HERE / 'rpmodes.c'
BINM = fc.SP / 'build' / 'rpmodes'
B, GAM, LAM = 3.0, 0.5, 0.5
mp.mp.dps = 50
RTOL = sys.argv[4] if len(sys.argv) > 4 else '1e-11'


def build():
    BINM.parent.mkdir(parents=True, exist_ok=True)
    if not BINM.exists() or BINM.stat().st_mtime < SRC.stat().st_mtime:
        subprocess.run(['gcc', '-O3', '-std=c11', '-o', str(BINM), str(SRC), '-lm'], check=True)


def hessian6(z, u, m2c, p):
    """numpy port of rpmodes.c hessian6 (canonical real comps x1 x2 y1 y2 s1 s2)."""
    x, y, s = z
    A = np.zeros((3, 3), complex); Bm = np.zeros((3, 3), complex)
    A[0, 0] = 4 * p['lx'] * abs(x) ** 2 + m2c; A[0, 1] = 3 * p['kap'] * np.conj(x) ** 2; A[1, 0] = 3 * p['kap'] * x * x
    A[1, 1] = p['M'] ** 2 + 4 * p['ly'] * abs(y) ** 2; A[1, 2] = 2 * np.conj(u) * np.conj(y); A[2, 1] = 2 * u * y
    A[2, 2] = 2 * p['ls'] * (2 * abs(s) ** 2 - p['h'] ** 2 / 2)
    Bm[0, 0] = 2 * p['lx'] * x * x + 6 * p['kap'] * y * np.conj(x); Bm[1, 1] = 2 * p['ly'] * y * y + 2 * np.conj(u) * s
    Bm[2, 2] = 2 * p['ls'] * s * s
    Hs = np.zeros((6, 6))
    for a in range(3):
        for c in range(3):
            Hs[2 * a, 2 * c] = A[a, c].real + Bm[a, c].real; Hs[2 * a + 1, 2 * c + 1] = A[a, c].real - Bm[a, c].real
            Hs[2 * a, 2 * c + 1] = -A[a, c].imag + Bm[a, c].imag; Hs[2 * a + 1, 2 * c] = A[a, c].imag + Bm[a, c].imag
    return Hs


def bd_desitter(k, m2):
    """physical-field BD mode of mass^2 m2 in de Sitter (H = 1) at a = 1 (-eta = 1): (delta r, delta r_dot)."""
    nu = mp.sqrt(mp.mpf(9) / 4 - mp.mpf(m2))
    pref = mp.sqrt(mp.pi) / 2 * mp.exp(1j * (mp.pi * nu / 2 + mp.pi / 4))
    Hn = mp.hankel1(nu, k); Hm = mp.hankel1(nu - 1, k)
    Hp = Hm - nu / k * Hn
    dr = pref * Hn
    drd = -pref * (mp.mpf(3) / 2 * Hn + k * Hp)
    return complex(dr), complex(drd)


def main(sub, end, tag=None):
    t0 = time.time()
    build()
    p = fc.params(LAM)
    st = fc.stationary(fc.U_PREP, p)
    zh = st['z']
    G, I0, ph0 = fc.trilinear_G(zh, p, B)
    NF, OFF = (4, 2) if sub == 'YS' else (6, 0)
    m2c0 = p['xi'] * 12.0
    H6 = hessian6(zh, fc.U_PREP, m2c0, p)
    ev6 = np.linalg.eigvalsh(H6)
    Hsub = H6[OFF:OFF + NF, OFF:OFF + NF]
    w2, E = np.linalg.eigh(Hsub)
    KG = np.geomspace(float(sys.argv[6]) if len(sys.argv) > 6 else 0.03, float(sys.argv[5]) if len(sys.argv) > 5 else 80.0, int(sys.argv[7]) if len(sys.argv) > 7 else 40)
    rows = []
    wr = []
    for k in KG:
        for al in range(NF):
            dr, drd = bd_desitter(float(k), max(float(w2[al]), 0.0))
            wr.append((dr * np.conj(drd) - np.conj(dr) * drd).imag)
            col = []
            for i in range(NF):
                v, vd = E[i, al] * dr, E[i, al] * drd
                col += [v.real, v.imag, vd.real, vd.imag]
            rows.append(col)
    mf = fc.SP / 'runs' / f'modes_{sub}.txt'
    with open(mf, 'w') as f:
        f.write('%d %d %d\n' % (len(KG), NF, OFF))
        f.write(' '.join('%.17g' % k for k in KG) + '\n')
        for col in rows:
            f.write(' '.join('%.17g' % v for v in col) + '\n')
    g = lambda v: '%.17g' % float(v)
    nsamp = int(round(end / 0.25))
    args = [str(BINM), '0', g(B), g(GAM), g(p['xi']), g(fc.ETA2), g(p['cH']), g(p['lx']), g(p['ly']), g(p['ls']), g(p['kap']), g(p['M']),
            g(p['h']), g(p['mu']), '0', '0', g(G.real), g(G.imag), g(end), RTOL, '1e-16', str(nsamp), '0', str(mf), '2', 'F']
    args += [g(v) for c in zh for v in (c.real, c.imag)] + ['0'] * 6
    out = subprocess.run(args, capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(out.stderr)
    Dl, Ml, END = {}, [], None
    for line in out.stdout.split('\n'):
        if line.startswith('D '):
            v = np.array(line.split()[1:], float); Dl[round(v[0], 9)] = v
        elif line.startswith('M '):
            v = np.array(line.split()[1:], float); Ml.append(v)
        elif line.startswith('END'):
            END = line.split()[1:]
    print('engine done %.1fs, %d M lines, END %s' % (time.time() - t0, len(Ml), END), flush=True)
    Q = np.array([1.0, 3.0, 6.0])
    res = dict(method=__doc__, subsystem=sub, NF=NF, OFF=OFF, kgrid=KG.tolist(), end=end, engine_END=END,
               prep_hessian6=ev6.tolist(), prep_sub_eigs=w2.tolist(), wronskian_max_dev=float(max(abs(np.array(wr) - 1))), times=[])
    lk = np.log(KG)
    for v in Ml:
        t = v[0]; data = v[1:].reshape(len(KG), NF, NF, 4)
        D = Dl[round(t, 9)]
        y = D[1:20]; Hh, Rr, Hd = D[20], D[21], D[22]
        a = np.exp(y[0]); z = np.array([y[3] + 1j * y[4], y[5] + 1j * y[6], y[7] + 1j * y[8]])
        u = p['mu'] + G * (y[1] / np.sqrt(fc.ETA2))
        Hs = hessian6(z, u, p['xi'] * Rr, p)[OFF:OFF + NF, OFF:OFF + NF]
        wm2, Em = np.linalg.eigh(Hs)
        # label eigenvectors by field content (Y share, Sigma share, X share)
        comp = []
        for b in range(NF):
            e = Em[:, b]; full = np.zeros(6); full[OFF:OFF + NF] = e ** 2
            comp.append(dict(m2=float(wm2[b]), X=float(full[0] + full[1]), Y=float(full[2] + full[3]), S=float(full[4] + full[5])))
        nk = np.zeros((len(KG), NF)); ek = np.zeros((len(KG), NF))
        for ik, k in enumerate(KG):
            # rescaled v = a^{3/2} delta r: v'' + Omega^2 v = 0, Omega^2 = k^2/a^2 + M^2 - 9H^2/4 - 3Hdot/2 (cosmic time)
            om2 = (k / a) ** 2 + wm2 - 2.25 * Hh ** 2 - 1.5 * Hd
            om = np.sqrt(np.maximum(om2, 1e-30))
            V = data[ik, :, :, 0] + 1j * data[ik, :, :, 1]       # [alpha, i]
            Vd = data[ik, :, :, 2] + 1j * data[ik, :, :, 3]
            q = a ** 1.5 * (V @ Em); pq = a ** 1.5 * ((Vd + 1.5 * Hh * V) @ Em)    # [alpha, beta]
            nk[ik] = np.sum(np.abs(pq) ** 2 + om2 * np.abs(q) ** 2, axis=0) / (2 * om) - 0.5
            ek[ik] = np.sqrt((k / a) ** 2 + np.maximum(wm2, 0))
        Nb = [float(np.trapezoid(KG ** 3 * nk[:, b] / (2 * np.pi ** 2), lk)) for b in range(NF)]
        Eb = [float(np.trapezoid(KG ** 3 * nk[:, b] * ek[:, b] * a / (2 * np.pi ** 2), lk)) for b in range(NF)]   # conformal energy
        res['times'].append(dict(tau=float(t), a=float(a), H=float(Hh), eig=comp, N=Nb, E_conformal=Eb,
                                 n_k=nk.T.tolist(), absY=float(abs(z[1])), absS=float(abs(z[2])), u=[float(u.real), float(u.imag)]))
        if abs(t - round(t)) < 1e-9 and (round(t) in (0, 1, 2, 3, 5, 10, 20, 30, 40, 50, 75, 100, 150, 200) or t == Ml[-1][0]):
            print('tau=%6.2f a=%7.3f ' % (t, a) + ' '.join('[m2=%.1f X%.2f Y%.2f S%.2f N=%.3e]' % (c['m2'], c['X'], c['Y'], c['S'], n)
                                                    for c, n in zip(comp, Nb)), flush=True)
    res['seconds'] = time.time() - t0
    out_file = RES / ('ymodes_%s.json' % (tag or sub))
    json.dump(res, open(out_file, 'w'), default=float)
    print('wronskian dev', res['wronskian_max_dev'], 'prep sub eigs', w2, 'seconds', res['seconds'])


if __name__ == '__main__':
    main(sys.argv[1], float(sys.argv[2]), sys.argv[3] if len(sys.argv) > 3 else None)
