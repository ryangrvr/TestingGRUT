"""en_rerun_skeptic -- OWN 3D evolution code (no Codex module imported, no Codex binary called).

Purpose: continue the main N48 history from a Codex-engine snapshot with an independently written integrator and
measure the reference projections DENSELY in tau (en_rerun sampled them only every 20 tau, i.e. about twice per
late orbit of the coherent mean).

Discretisation (same model and discretisation as MATHEMATICS.md, own implementation):
  H = 1/2 pi^2 + 1/2 |grad chi|^2_FD4 + lX/4 r^4 + g2/2 phi^2 r^2 + lp/4 phi^4, r^2 = chi1^2+chi2^2
  FD4 Laplacian applied SPECTRALLY with its exact symbol (Codex uses the real-space 13-point stencil);
  fourth-order Yoshida triple jump (same coefficients), same output schedule (tau grid .05 accumulated exactly as
  a running double sum from 20.05, steps h = min(hmax, target - s)), radiation clock tau = 20 + a_m s + C s^2/2.
Diagnostics and projections: own code (proj_own.py).

usage: python3 -I own_engine.py <snapshot.bin> <tau_end> <diag_every> <proj_every> <out_json> [save_final.bin]
"""
import sys, json, struct, time
import numpy as np
from scipy import fft as sfft
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import proj_own as P

LX, G2, LP = .001, .01, .03
A_M, CC, TAU0 = 8.869218738806, 1.9857962750634526, 20.0
HMAX, SAMPLE = .004, .05
W1 = 1/(2 - 2**(1/3)); W0 = 1 - 2*W1
WORKERS = 2


def clock(s):
    return TAU0 + A_M*s + CC*s*s/2


def eta_of_tau(t):
    return 2*(t - TAU0)/(A_M + np.sqrt(A_M*A_M + 2*CC*(t - TAU0)))


def read_snapshot(path):
    with open(path, 'rb') as f:
        N, L, s, bits = struct.unpack('<iddi', f.read(24))
        a = np.fromfile(f, dtype='<f8' if bits == 64 else '<f4').reshape(6, N, N, N).astype(float)
    return N, L, s, a[:3].copy(), a[3:].copy()


def write_snapshot(path, N, L, s, q, p):
    with open(path, 'wb') as f:
        f.write(struct.pack('<iddi', N, L, s, 64))
        np.asarray(q, '<f8').tofile(f); np.asarray(p, '<f8').tofile(f)


def main():
    snap, tau_end, diag_every, proj_every, outp = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), sys.argv[5]
    save_final = sys.argv[6] if len(sys.argv) > 6 else None
    N, L, s, q, p = read_snapshot(snap)
    sym = P.fd4_symbol_r(N, L)

    def force(q):
        qh = sfft.rfftn(q, axes=(1, 2, 3), workers=WORKERS)
        lap = sfft.irfftn(-sym*qh, s=(N, N, N), axes=(1, 2, 3), workers=WORKERS)
        r2 = q[0]*q[0] + q[1]*q[1]; ph2 = q[2]*q[2]
        fX = G2*ph2 + LX*r2
        lap[0] -= fX*q[0]; lap[1] -= fX*q[1]; lap[2] -= (LP*ph2 + G2*r2)*q[2]
        return lap

    # reproduce the running output-time sum exactly (Codex: next_tau = tau0+sample, then next_tau += sample)
    nt = TAU0 + SAMPLE
    targets_done = 0
    while True:
        tg = eta_of_tau(nt)
        if tg > s + 1e-12:
            break
        nt += SAMPLE; targets_done += 1
    end_s = eta_of_tau(tau_end)
    rows = []; k = 0; t0 = time.time()
    rows.append(P.diagnostics(q, p, L, sym, s, clock(s), proj=True))
    while s < end_s - 1e-13:
        target = min(end_s, eta_of_tau(nt))
        while s < target - 1e-13:
            h = min(HMAX, target - s)
            p += (W1*h/2)*force(q); q += (W1*h)*p
            p += ((W1+W0)*h/2)*force(q); q += (W0*h)*p
            p += ((W0+W1)*h/2)*force(q); q += (W1*h)*p
            p += (W1*h/2)*force(q)
            s += h
        s = target; nt += SAMPLE; k += 1
        if k % diag_every == 0 or s >= end_s - 1e-13:
            rows.append(P.diagnostics(q, p, L, sym, s, clock(s), proj=(k % proj_every == 0 or s >= end_s - 1e-13)))
    if save_final:
        write_snapshot(save_final, N, L, s, q, p)
    json.dump(dict(start_snapshot=snap, N=N, L=L, tau_end=tau_end, seconds=time.time()-t0, rows=rows), open(outp, 'w'), default=float)
    print('done', N, 'tau', clock(s), 'rows', len(rows), 'sec', round(time.time()-t0, 1), flush=True)


if __name__ == '__main__':
    main()
