#!/usr/bin/env python3
"""inf_potential (Claude, own code): Y-pair production by the oscillating inflaton through the holder trilinear itself.

-(g_I e^{i alpha_I} phi_I Sigma* Y^2 + h.c.) at the late vacuum (u = mu real, Sigma = h/sqrt2) couples phi_I to the real Y
components: phi_I a^2 (heavy, m_H = 36.12 H_i) and phi_I b^2 (light, m_L = 22.25 H_i) with |c| = g_I |cos alpha_I| |Sigma|,
and phi_I a b with |c'| = 2 g_I |sin alpha_I| |Sigma| [D].  For a homogeneous classical oscillation
phi_I(t) = sum_n phi_n cos(n w t + theta_n), the Born (perturbative) pair rate per volume is
    identical pair (b b or a a):  c^2 phi_n^2 beta / (16 pi),        a b:  c'^2 phi_n^2 beta / (32 pi),
with beta the two-body velocity factor at frequency n w (the n = 1 term is the usual perturbative decay rate times n_phi).
Each full period between turning points of the same side is Fourier-analysed (256 points, n = 1..10).  Narrow-resonance
checks: q = delta m^2 / (4 w^2) and q^2 w / H are printed (Born valid if q^2 w/H << 1 and q << 1).
Hidden heat [E]: each light Y deposits 0.6565 E_X (dq_heat ledger at K = 1.986); each heavy Y cascades to Y_L + A_L and
deposits 0.6565 + 1.79 (1e-6/e)^{1/2} E_X [Q dq_heat].  E_X = 188706 H_i^4 (M1r release, a = 1 at its slice) [Q].
Both holder calibrations: g_I fixed at S1 (1 e-fold before eps_H = 1) and at Send (eps_H = 1).
Control [D]: the M1r surrogate (quadratic, m_I = 3 H_i) gives zero (no harmonics, fundamental below every threshold);
[C] the LL channel is checked against linear light-Y modes on the same background (modes.py ycontrol).
Usage: python3 yy_harmonics.py T1 E1 HT6t   Output: ../results/yy_harmonics.json.   Labels: [C] Born rates on our
backgrounds; [E] the heat ledger and the homogeneous (pre-fragmentation) background.
"""
import json, math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import inf_common as ic
from exit_factor import model_from_bench

RES = Path(__file__).resolve().parent.parent / 'results'
ALPHA_I = math.atan2(10.0, -7.159456159)
SIG = 40 / math.sqrt(2)
M_L = 22.249
M_H = math.sqrt(900 + 2 * 7.159456159 * SIG)
DEP_L = 0.6565
DEP_AL = 1.7805
EX = 188706.29


def beta2(m1, m2, W):
    x = 1 - (m1 + m2) ** 2 / W ** 2; y = 1 - (m1 - m2) ** 2 / W ** 2
    return math.sqrt(x * y) if x > 0 else 0.0


def analyse(name, gam_Hi, tmax_Hi):
    m, rec = model_from_bench(name)
    bg = ic.background(m, 75.0)
    V0 = rec['V0_MP4']; Hi = rec['H_i_GeV']; hs = Hi / ic.MPL_GEV / math.sqrt(V0)
    Hi_MP = Hi / ic.MPL_GEV
    pmin = m.phi_min()
    sol = ic.post_inflation(m, bg, gam_Hi * hs, tmax_Hi / hs, rtol=1e-10, atol=1e-14, dense=True)
    t = np.linspace(0, tmax_Hi / hs, int(max(400001, 80 * rec['m_I_over_Hi'] * tmax_Hi)))
    lna, f, fd, rr = sol.sol(t)
    idx = np.where(np.sign(fd[1:]) != np.sign(fd[:-1]))[0]
    tps = t[idx] + fd[idx] / (fd[idx] - fd[idx + 1]) * (t[idx + 1] - t[idx])
    out = {}
    for lab, phiI_cal in (('S1', rec['phiI_S1_MP']), ('Send', rec['phiI_end_MP'])):
        gI = abs(10j - 7.159456159) * Hi_MP / abs(phiI_cal)
        cL = gI * abs(math.cos(ALPHA_I)) * SIG; cH = cL; cLH = 2 * gI * abs(math.sin(ALPHA_I)) * SIG
        NL = NH = NLL = 0.0; rows = []; last = None
        for j in range(0, len(tps) - 2, 2):
            t0, t1 = tps[j], tps[j + 2]
            T = (t1 - t0) * hs                                  # period in H_i^-1 (t_Hi = t_code * hs)
            w = 2 * math.pi / T
            ts = np.linspace(t0, t1, 257)[:-1]
            Y = sol.sol(ts)
            ph = (Y[1] - pmin) / Hi_MP                         # phi_I in H_i
            F = np.fft.rfft(ph) / len(ts)
            amp = 2 * np.abs(F[1:11])                          # phi_n, n = 1..10
            a = math.exp(float(np.mean(Y[0])))
            rate = 0.0; rL = rH = 0.0
            for n, An in enumerate(amp, start=1):
                W = n * w
                pLL = cL ** 2 * An ** 2 * beta2(M_L, M_L, W) / (16 * math.pi)
                pHH = cH ** 2 * An ** 2 * beta2(M_H, M_H, W) / (16 * math.pi)
                pLH = cLH ** 2 * An ** 2 * beta2(M_L, M_H, W) / (32 * math.pi)
                rL += 2 * pLL + pLH; rH += 2 * pHH + pLH
                NLL += 2 * pLL * a ** 3 * T
            NL += rL * a ** 3 * T; NH += rH * a ** 3 * T
            Hloc = math.sqrt((0.5 * np.mean(Y[2] ** 2) + np.mean([m.v(x) for x in Y[1]]) + np.mean(Y[3])) / 3) / hs
            q = 2 * cL * amp[0] / (4 * w * w)
            last = (rL + rH) * a ** 3
            if len(rows) < 12 or j % 20 == 0:
                rows.append(dict(t_Hi=float(t0 * hs), a=a, omega_Hi=w, H_Hi=Hloc, harm_rel=[float(x) for x in amp[:6] / amp[0]],
                                 phi1_Hi=float(amp[0]), q=q, q2w_over_H=q * q * w / Hloc, NY_rate_comoving=last))
        tail = last / (gam_Hi) if gam_Hi > 0 else None    # remaining production ~ e^{-Gamma t}
        R = NL * DEP_L + NH * (DEP_L + DEP_AL)
        out[lab] = dict(g_I=gI, N_YL=NL, N_YH=NH, N_YL_from_LL=NLL, R_YY_to_tmax=R, tail_NY_estimate=tail,
                        R_YY_with_tail=(R + (tail or 0) * (DEP_L + 0.5 * DEP_AL)), rows=rows)
    return dict(benchmark=name, gamma_over_Hi=gam_Hi, t_max_Hi=tmax_Hi, m_I_over_Hi=rec['m_I_over_Hi'], res=out)


def main(names):
    allr = dict(method=__doc__, runs=[])
    for nm in names:
        for gam, tm in ((0.5, 12.0), (0.0, 6.0)):
            r = analyse(nm, gam, tm)
            allr['runs'].append(r)
            for lab, v in r['res'].items():
                print('%-5s Gamma=%.1f g_I(%s)=%.3e: N_YL=%.3e N_YH=%.3e R_YY(t<%g)=%.3e R(+tail)=%.3e | first: w=%.1f q=%.3g q2w/H=%.3g harm=%s' % (
                    nm, gam, lab, v['g_I'], v['N_YL'], v['N_YH'], tm, v['R_YY_to_tmax'], v['R_YY_with_tail'],
                    v['rows'][0]['omega_Hi'], v['rows'][0]['q'], v['rows'][0]['q2w_over_H'], ['%.3f' % x for x in v['rows'][0]['harm_rel']]), flush=True)
            json.dump(allr, open(RES / 'yy_harmonics.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1:])
