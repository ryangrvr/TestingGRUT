#!/usr/bin/env python3
"""Skeptic of inf_higgs: the Higgs mass the S1 channel induces through the N-L-H loop, along the inflationary trajectory.

S1: -(y/2) phi_I Nbar N with m_N = M_N + y phi_I, seesaw y_nu L H~ N with y_nu^2 = mt M_N / v^2 (mt = 0.05 eV).
One-loop Coleman-Weinberg of the N-nu system [D]: mass matrix [[0, m_D], [m_D, M]], m_D^2 = y_nu^2 h^2/2,
M = |M_N + y phi_I|.  Expanding the heavy eigenvalue M + m_D^2/M in V = -(2/64 pi^2) sum m^4 (ln m^2/mu^2 - 3/2):
    V ⊃ -(1/8 pi^2) (ln(M^2/mu^2) - 1) M^2 m_D^2   =>   delta m_h^2(phi_I) = -(y_nu^2/8 pi^2) M(phi_I)^2 [ln(M^2/mu^2) - 1].
Expanding M^2 = M_N^2 + 2 y M_N phi_I + y^2 phi_I^2: the constant (EW-tuned), the R-odd A-term the track bounds
(A_ind ~ y y_nu^2 M_N / 4 pi^2 x log), and an R-EVEN piece lam_ind phi_I^2 |H|^2 ~ -(y^2 y_nu^2/8 pi^2) x log, larger than
the A piece by y phi_I / (2 M_N) (~ 65 at y_nat on HT6t's plateau) and of a sign fixed by the renormalisation scale mu,
not by the CP branch.  The track assumes it vanishes in every gravitational-only and plasma-only result.
  python3 -B sk_induced.py   -> results/sk_induced.json
"""
import sys
sys.dont_write_bytecode = True
import json
import math

import numpy as np
from scipy.optimize import minimize_scalar

import sk_common as C
import sk_infl as I

MNU, VEW = 0.05e-9, 174.1


def main():
    P = json.load(open(C.TRACKS / 'inf_higgs' / 'results' / 'plasma.json'))
    out = {}
    for name in ('HT6t', 'T1', 'E1'):
        pot, hist = I.history(name)
        Nb, H, R, eps, fI = hist.T
        sel = Nb <= 5.0
        Hg = pot.H_GeV
        mI = pot.m_I
        yn = P[name]['y_nat']
        o = dict(y_nat=yn)
        for MNf in (0.25, 0.0625):
            MN = MNf * mI                                    # H_i units
            ynu2 = MNU * MN * Hg / VEW ** 2
            c = ynu2 / (8 * math.pi ** 2)
            for y in (1e-5, 3e-5, 1e-4, yn):
                M = np.abs(MN + y * fI / pot.eta)            # H_i units, along N_b
                Ms = M[sel]; Hs = H[sel]
                f = lambda lmu: c * Ms ** 2 * (np.log(Ms ** 2) - 2 * lmu - 1) * (-1) / Hs ** 2
                # best renormalisation scale: minimise the largest |delta m^2|/H^2 over the last 5 e-folds
                r = minimize_scalar(lambda lmu: float(np.max(np.abs(f(lmu)))),
                                    bounds=(math.log(1e-3 * Ms.min()), math.log(1e3 * Ms.max())), method='bounded')
                best = float(r.fun)
                # least destabilising choice: largest -delta m^2 (most negative) over the window, minimised
                rS1 = int(np.argmin(np.abs(Nb - 1.0)))
                MS1 = M[rS1]
                dm_MN = -c * MS1 ** 2 * (math.log(MS1 ** 2 / MN ** 2) - 1)
                dm_MP = -c * MS1 ** 2 * (math.log(MS1 ** 2 / (1 / pot.eta) ** 2) - 1)
                o[f'MN={MNf} y={y:.3g}'] = dict(M_S1_over_Hi=float(MS1), M_end=float(M[0]), M_Nb5=float(Ms.max()),
                                               ynu2=ynu2, c_M2_S1=float(c * MS1 ** 2),
                                               best_mu_maxabs_dm2_over_H2_last5=best,
                                               best_mu_GeV=math.exp(r.x) * Hg,
                                               dm2_S1_mu_eq_MN=float(dm_MN), dm2_S1_mu_eq_MP=float(dm_MP),
                                               lam_ind_scale=float(y ** 2 * ynu2 / (8 * math.pi ** 2)),
                                               A_over_lam_piece_S1=float(2 * MN / max(abs(y * fI[rS1] / pot.eta), 1e-30)))
                print(name, MNf, '%.2e' % y, 'M(S1) %.1f H_i' % MS1, 'cM^2 %.3g' % (c * MS1 ** 2),
                      'best |dm2|/H2 over last 5 e-folds %.3g' % best, 'dm2(S1) mu=M_N %.3g mu=M_P %.3g' % (dm_MN, dm_MP),
                      flush=True)
        out[name] = o
    C.jdump(out, C.RES / 'sk_induced.json')


if __name__ == '__main__':
    main()
