#!/usr/bin/env python3
"""dq_heat (Claude, own code), item 1, piece B (gate1b_pieceB.py adapted to M1r): the release-induced Goldstone angle.

In each superhorizon patch the release converts a perturbation of the prepared state into a late gauge-orbit offset
alpha_G = Delta theta_Sigma(inf)/6 (late F is Sigma dominated), which adds to the inflationary Goldstone angle of piece A.
Separate-universe responses d alpha_G / d q_alpha along every physical eigenmode alpha of the inflationary 6x6 Hessian at the
M1r prepared point (canonical real components; the gauge zero mode is piece A itself and is excluded), by symmetric
differences of rptri_dq runs (amplitudes 1e-3 and 1e-2 H_i), theta_Sigma window mean over tau = 300-400 plus the analytic
radiation tail -N_Sigma / (|Sigma|^2 (2K)^{3/2} t^{1/2}).
Superhorizon power of each heavy mode at the end of inflation: the exact de Sitter BD mode (ymodes.bd_desitter),
P_alpha(k) = k^3 |delta r_alpha(k)|^2 / (2 pi^2) (gate1b used the order-of-magnitude e^{-3N/2}).
Piece-A Goldstone power: P_theta = (H/2pi F_inf)^2 (massless BD, F_inf = 275.93 H_i).
N_A,B = int dlnk [dN_A/dlnk]_A(k) * P_B(k)/P_theta, restricted to k <= 0.3 (superhorizon during the release), with the
piece-A spectra of vector_modes.json.  Usage: python3 -B sep_universe.py
"""
import json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dq_model as fc
from ymodes import hessian6, bd_desitter

RES = Path(__file__).resolve().parent.parent / 'results'
B, GAM, LAM = 3.0, 0.5, 0.5
END = 400.0


def theta_inf(rr, K):
    D = rr['D']; t = D[:, 0]; y = D[:, 1:20]
    S = y[:, 7] + 1j * y[:, 8]; vS = y[:, 13] + 1j * y[:, 14]; a = np.exp(y[:, 0])
    th = np.unwrap(np.angle(S))
    NS = -2 * np.imag(np.conj(S) * vS) * a ** 3
    tail = -NS / (np.abs(S) ** 2 * (2 * K) ** 1.5 * np.sqrt(np.maximum(t, 1e-9)))
    sel = t >= 300.0
    return float(np.mean(th[sel] + tail[sel])), float(np.mean(NS[sel])), float(th[0])


def main():
    t0 = time.time()
    p = fc.params(LAM)
    st = fc.stationary(fc.U_PREP, p)
    zh = st['z']
    G, I0, ph0 = fc.trilinear_G(zh, p, B)
    H6 = hessian6(zh, fc.U_PREP, p['xi'] * 12.0, p)
    w2, E = np.linalg.eigh(H6)
    r0 = fc.r_of_z(zh)
    K = json.load(open(RES / 'background.json'))['late']['K_end']
    base = fc.run('su_base', p, B, GAM, z0=zh, G=G, end=END, nsamp=2000, binw=0, save=False)
    th0, NS0, _ = theta_inf(base, K)
    resp = []
    for al in range(6):
        if w2[al] < 1e-6:
            continue                        # gauge zero mode = piece A
        row = dict(mode=al, m2=float(w2[al]), vec=E[:, al].tolist())
        for amp in (1e-3, 1e-2):
            th = {}
            for sg in (+1, -1):
                z = fc.z_of_r(r0 + sg * amp * E[:, al])
                # the prepared slice is at rest: G (holder) unchanged; the inflaton follows from the constraint in the engine
                rr = fc.run('su_%d_%+d' % (al, sg), p, B, GAM, z0=z, G=G, end=END, nsamp=2000, binw=0, save=False)
                th[sg] = theta_inf(rr, K)
            row['dthetaS_dq_amp%.0e' % amp] = (th[1][0] - th[-1][0]) / (2 * amp)
            row['dalphaG_dq_amp%.0e' % amp] = row['dthetaS_dq_amp%.0e' % amp] / 6.0
            row['dNS_dq_amp%.0e' % amp] = (th[1][1] - th[-1][1]) / (2 * amp)
        resp.append(row)
        print('mode %d m2=%.3f  dalphaG/dq = %.5e (1e-3) %.5e (1e-2)' % (al, w2[al], row['dalphaG_dq_amp1e-03'], row['dalphaG_dq_amp1e-02']), flush=True)
    F_inf = float(np.sqrt(2 * np.sum(np.array([1, 9, 36]) * np.abs(zh) ** 2)))
    P_theta = (1 / (2 * np.pi * F_inf)) ** 2
    VM = json.load(open(RES / 'vector_modes.json'))
    out = dict(method=__doc__, base=dict(theta_S_inf=th0, N_Sigma=NS0, winding_total=th0 - float(np.angle(zh[2]))), K=K, F_inf=F_inf,
               P_theta_pieceA=P_theta, hessian_eigs=w2.tolist(), responses=resp, per_e=[])
    for case in VM['cases']:
        if case['variant'] != 'full':
            continue
        e = case['e']; ks = np.array([s['k'] for s in case['spectrum']]); dN = np.array([s['dN_dlnk'] for s in case['spectrum']])
        sel = ks <= 0.3
        PB = np.zeros_like(ks)
        for r in resp:
            g = r['dalphaG_dq_amp1e-03']
            for i, k in enumerate(ks):
                if sel[i]:
                    dr, _ = bd_desitter(float(k), r['m2'])
                    PB[i] += g * g * k ** 3 * abs(dr) ** 2 / (2 * np.pi ** 2)
        NB = float(np.trapezoid((dN * PB / P_theta)[sel], np.log(ks[sel])))
        ik = int(np.argmin(abs(ks - case['k_star'])))
        out['per_e'].append(dict(e=e, N_A_pieceA=case['N_A'], N_A_pieceB_SU=NB, ratio_B_over_A=NB / case['N_A'],
                                 PB_over_Ptheta_at_kstar=float(PB[ik] / P_theta), k_star=case['k_star'],
                                 R_vec_pieceB_SU=NB * case['deposit_per_vector_EX']))
        print('e=%.1e  N_A(A)=%.4f  N_A(B, SU)=%.3e  ratio %.2e  P_B/P_theta(k*)=%.2e' % (e, case['N_A'], NB, NB / case['N_A'], PB[ik] / P_theta), flush=True)
    out['seconds'] = time.time() - t0
    json.dump(out, open(RES / 'sep_universe.json', 'w'), indent=1, default=float)
    print('seconds', out['seconds'])


if __name__ == '__main__':
    main()
