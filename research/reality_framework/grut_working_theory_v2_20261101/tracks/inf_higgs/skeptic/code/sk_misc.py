#!/usr/bin/env python3
"""Skeptic of inf_higgs: small checks.
  (1) hard-gas composition: induced Higgs mass coefficient of a pure L + H gas (N -> L H) against the track's
      thermal-normalised c = 0.0228.  One-loop energy shift of a hard gas in a background h: dE = m_i^2(h)/2E per quantum.
      In the region of field h the 4 real Higgs-doublet quanta are 1 radial (m^2 = 3 lam h^2) + 3 longitudinal W+, W-, Z
      (equivalence theorem: m_W^2 = g^2 h^2/4, m_Z^2 = (g^2 + g'^2) h^2/4); leptons couple only through y_l (neglected).
      => m_ind^2 = (n_H / E) <c_i>,  n_H = n_SM / 2.   [A: physical-mode counting; gauge-independent at this order]
  (2) non-equilibrium Gaussian tail at the end of inflation vs M_t (mode-run variance, tree de Sitter barrier) -- the
      track's own 'conservative' alternative to Hawking-Moss, extended to the M_t range of its plasma claim.
  -> results/sk_misc.json
"""
import sys
sys.dont_write_bytecode = True
import math

import numpy as np

import sk_common as C
import sk_bench as B
import sk_infl as I


def composition():
    smi = C.sm_in(172.57)
    out = {}
    for mu in (1e11, 1e12, 1e13):
        lm = math.log(mu)
        g2 = float(smi.sp['G2'](lm)); gp2 = 0.6 * float(smi.sp['G1'](lm)); lam = float(smi.sp['lam'](lm))
        yt2 = float(smi.sp['T'](lm))
        kT = 3 * g2 / 16 + gp2 / 16 + yt2 / 4 + lam / 2
        cavg = (3 * lam + (3 * g2 + gp2) / 4) / 4            # per Higgs quantum
        c_LH = 0.5 * cavg                                    # n_H = n_SM/2
        # thermal normalisation of the track: c = kappa_T T^2 / (rho/E^2), E = 2.70 T, rho = pi^2 g* T^4/30
        c_track = 0.11 * 2.70 ** 2 / (math.pi ** 2 * 106.75 / 30)
        c_thermal_here = kT * 2.70 ** 2 / (math.pi ** 2 * 106.75 / 30)
        out[f'{mu:.0e}'] = dict(g2=g2, gp2=gp2, lam=lam, yt2=yt2, kappa_T=kT, c_avg_per_H=cavg, c_LH_gas=c_LH,
                                c_track=c_track, c_thermal_with_kappa_here=c_thermal_here, ratio_LH_to_track=c_LH / c_track)
        print('mu %.0e g2 %.3f gp2 %.3f lam %.4f yt2 %.3f kappa_T %.3f | c_LH %.4f track %.4f ratio %.2f' %
              (mu, g2, gp2, lam, yt2, kT, c_LH, c_track, c_LH / c_track))
    return out


def tail():
    pot, hist = I.history('HT6t')
    ev = B.load(B.key('HT6t', 'MD', 0))
    rows, spec, ks, wq = ev['rows'], ev['spec'], ev['ks'], ev['wq']
    te = ev['tau_end_inf']
    N0 = I.N0_MD('HT6t')
    out = {}
    j = int(np.argmin(np.abs(rows[:, 0] - te)))
    a, H, R = rows[j, 1], rows[j, 2], rows[j, 3]
    s2 = float(np.sum(wq[ks < a * H] * spec[ks < a * H, j]))
    Hg = pot.H_GeV
    for Mt in I.MTS:
        smi = C.sm_in(Mt)
        # tree dS barrier position at the end
        hs = np.exp(np.linspace(math.log(1e-3 * H * Hg), math.log(1e17), 3000))
        lmu = 0.5 * np.log(hs ** 2 + R * Hg ** 2)
        xi = 1 / 6 + (0.1 - 1 / 6) * np.exp(smi.sp['Ixi'](lmu) - smi.d['Ixi_Mt'])
        phi = np.exp(smi.sp['Gam'](lmu)) * hs
        V = 0.5 * xi * R * Hg ** 2 * phi ** 2 + 0.25 * smi.sp['lam'](lmu) * phi ** 4
        k = int(np.argmax(V))
        phib = phi[k] / Hg
        x = phib ** 2 / (2 * s2)
        l10 = (3 * (N0 - 0.0) + math.log1p(x) - x) / C.LN10
        out[str(Mt)] = dict(sigma2_kaH_end=s2, phi_b_end_over_H=phib, x=x, log10_patches_end=l10)
        print('Mt %.2f sigma2(k<aH) %.4f phi_b %.2f H_i  x %.0f  log10 patches %.0f' % (Mt, s2, phib, x, l10))
    return out


if __name__ == '__main__':
    o = dict(composition=composition(), gaussian_tail_end=tail())
    C.jdump(o, C.RES / 'sk_misc.json')
