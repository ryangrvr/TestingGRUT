#!/usr/bin/env python3
"""inf_chain (Claude, own code): holder-induced inflaton -> Y pairs (Born) for any benchmark, from its own slice.

Same physics and formulas as inf_potential/code/yy_harmonics.py (couplings and thresholds verified by its skeptic [D]):
at the late vacuum the trilinear -(G phi_I Sigma* Y^2 + h.c.) couples phi_I to the real Y components with
|c_LL| = |c_HH| = g_I |cos alpha_I| |Sigma| and |c_LH| = 2 g_I |sin alpha_I| |Sigma|; a homogeneous oscillation
phi_I = sum_n phi_n cos(n w t) produces pairs at c^2 phi_n^2 beta/16 pi (identical) and c'^2 phi_n^2 beta/32 pi (L H).
Differences from the potential track: the background starts at the benchmark's slice (S1, a = 1 there, the convention of
this chain), the Gamma_I term acts from the slice as in the engine, and g_I, alpha_I are the chain's (Bench.G).
Each full period between same-side turning points is Fourier analysed (512 points, n = 1..12).
Output: dict with N_YL, N_YH (Y_L, Y_H quanta per H_i^3, a = 1 at the slice), tail estimate.  [C] rates; the
homogeneous pre-fragmentation background and the late-vacuum masses are the stated limits [E].
Usage: python3 -B yy_born.py KEY [tmax]  -> ../results/yy_born_<KEY>.json
"""
import json, math, sys
import numpy as np
import ch_common as ch
from benchmarks_reg import get_bench

ic = ch.ic
SIG = 40 / math.sqrt(2)
M_L = 22.24859546128699
M_H = 36.124783736376884


def beta2(m1, m2, W):
    x = 1 - (m1 + m2) ** 2 / W ** 2; y = 1 - (m1 - m2) ** 2 / W ** 2
    return math.sqrt(x * y) if x > 0 else 0.0


def count(b, tmax_Hi=12.0, nharm=12):
    m = b.m; bg = b.bg; hs = b.h_s
    Hi_MP = b.Hi / ch.MP
    pmin = b.pmin
    sol = ic.post_inflation(m, bg, b.gam * hs, tmax_Hi / hs, start_Nbefore=b.Nslice, rtol=1e-10, atol=1e-14, dense=True)
    t = np.linspace(0, tmax_Hi / hs, int(max(400001, 120 * b.mI * tmax_Hi)))
    lna, f, fd, rr = sol.sol(t)
    idx = np.where(np.sign(fd[1:]) != np.sign(fd[:-1]))[0]
    tps = t[idx] + fd[idx] / (fd[idx] - fd[idx + 1]) * (t[idx + 1] - t[idx])
    alpha = math.atan2(b.G.imag, b.G.real)
    cL = b.gI * abs(math.cos(alpha)) * SIG; cH = cL; cLH = 2 * b.gI * abs(math.sin(alpha)) * SIG
    NL = NH = 0.0; rows = []; last = None
    for j in range(0, len(tps) - 2, 2):
        t0, t1 = tps[j], tps[j + 2]
        T = (t1 - t0) * hs
        w = 2 * math.pi / T
        ts = np.linspace(t0, t1, 513)[:-1]
        Y = sol.sol(ts)
        ph = (Y[1] - pmin) / Hi_MP
        F = np.fft.rfft(ph) / len(ts)
        amp = 2 * np.abs(F[1:nharm + 1])
        a = math.exp(float(np.mean(Y[0])) + b.Nslice)            # a = 1 at the slice
        rL = rH = 0.0
        for n, An in enumerate(amp, start=1):
            W = n * w
            pLL = cL ** 2 * An ** 2 * beta2(M_L, M_L, W) / (16 * math.pi)
            pHH = cH ** 2 * An ** 2 * beta2(M_H, M_H, W) / (16 * math.pi)
            pLH = cLH ** 2 * An ** 2 * beta2(M_L, M_H, W) / (32 * math.pi)
            rL += 2 * pLL + pLH; rH += 2 * pHH + pLH
        NL += rL * a ** 3 * T; NH += rH * a ** 3 * T
        last = (rL + rH) * a ** 3
        if len(rows) < 8 or j % 40 == 0:
            rows.append(dict(t_Hi=float(t0 * hs), a=a, omega_Hi=w, harm_rel=[float(x) for x in amp[:6] / max(amp[0], 1e-300)],
                             phi1_Hi=float(amp[0]), rate_comoving=last))
    tail = last / b.gam if (b.gam > 0 and last is not None) else None
    return dict(g_I=b.gI, alpha_I=alpha, m_I_over_Hi=b.mI, N_YL=NL, N_YH=NH, tail_NY_estimate=tail, t_max_Hi=tmax_Hi, n_periods=len(tps) // 2,
                rows=rows, normalisation='a = 1 at the slice')


def main(key, tmax=12.0):
    b = get_bench(key)
    r = count(b, tmax)
    r['key'] = key; r['method'] = __doc__
    json.dump(r, open(ch.RES / ('yy_born_%s.json' % key), 'w'), indent=1, default=float)
    print('%s Born: N_YL=%.4g N_YH=%.4g tail=%s (a = 1 at the slice), m_I/H_i=%.2f g_I=%.3e' % (key, r['N_YL'], r['N_YH'], r['tail_NY_estimate'], b.mI, b.gI))


if __name__ == '__main__':
    main(sys.argv[1], float(sys.argv[2]) if len(sys.argv) > 2 else 12.0)
