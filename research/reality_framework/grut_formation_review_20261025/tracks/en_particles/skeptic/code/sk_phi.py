"""en_particles SKEPTIC (own code): neutral phi quanta per net charge from a window, per-mode virial definition
n_k = sqrt(<|phi_k|^2><|pi_k|^2>) (real field, rfft half space with weights 2/1), uniform time weights.
Usage: python3 -I sk_phi.py <window.npz> <i0> <out.json>"""
import sys, json
import numpy as np
def main(path, i0, out):
    d = np.load(path); s = d['s']; i0 = int(i0); N = int(d['N']); L = float(d['L'])
    a0 = np.arange(N//2+1); wr = np.where((a0 == 0) | (a0 == N//2), 1., 2.)[None, None, :]
    A = 0.; B = 0.; Q = []
    for j in range(i0, len(s)):
        A = A+np.abs(d['phi'][j].astype(complex))**2; B = B+np.abs(d['pphi'][j].astype(complex))**2
        p = d['psi'][j].astype(complex); q = d['ppsi'][j].astype(complex); Q.append(float(np.imag(np.conj(p)*q).sum()))
    M = len(s)-i0; A /= M; B /= M; Q = float(np.mean(Q))
    nk = wr*np.sqrt(A*B); nk[0, 0, 0] = 0.
    res = dict(window=[float(s[i0]), float(s[-1])], phi_quanta_per_net_virial=float(nk.sum()/Q))
    json.dump(res, open(out, 'w'), indent=1); print(res)
if __name__ == '__main__': main(*sys.argv[1:])
