#!/usr/bin/env python3
"""mf_higgs SKEPTIC: induced Higgs masses after the waterfall, recomputed on three histories with a different method
from the track's fluid (mh_common.fluid): the energy fractions are read from the histories themselves (mf_hybrid's
export and mf_chain's two rebuilt histories, whose stage 1 already carries the seesaw decays), and the quantum energies
are bounded, not evolved.  [E] throughout; every closure is a model choice.

Two medium masses for a region field phibar (forward scattering of on-shell quanta, m^2 = sum_i kappa_i n_i <1/E_i>,
with m_i^2(h) = m_i^2(0) + kappa_i h^2; the same formula gives the textbook thermal masses):
  (1) SM hard gas from N -> L H (the track's term): m^2 = c rho_SM / E_SM^2, c = 0.0228 (track), 0.026 (my composition
      count: [3 lam + (3 g^2 + g'^2)/4]/8 at mu ~ 1e14-1e15 GeV, cf. inf_higgs skeptic N11), E_SM = E_N/2;
  (2) NEW: the singlet-neutrino gas itself.  In a Higgs background the seesaw mixing raises the heavy mass,
      M_h^2 = M^2 + y_nu^2 h^2 (m_D = y_nu h/sqrt2), so kappa_N = y_nu^2 and m^2 = y_nu^2 rho_N / E_N^2 > 0.
      y_nu^2 = m_nu M / v^2 = 0.026 (N4) ... 0.086 (N2) at m-tilde = 0.05 eV.  This exists as soon as the N are produced
      (during the last chi infall), without any N decay, and does not depend on the SM gas composition.
E_N is bounded from above by its injection energy (3/4) mu5,max at the end of inflation, NOT redshifted (conservative:
m^2 ~ 1/E^2).  Output: results/sk_plasma.json
"""
import sys
sys.dont_write_bytecode = True
import math

import numpy as np

import sk_common as C

HP = {'EXP (mf_hybrid export)': C.TRACKS / 'mf_hybrid' / 'results' / 'bench_HFa.json',
      'MC_TR (mf_chain unblocked, 1.28e11)': C.TRACKS / 'mf_chain' / 'results' / 'histories' / 'bench_itHFa_track_it1.json',
      'MC_P0 (mf_chain Pauli, 2.40e11)': C.TRACKS / 'mf_chain' / 'results' / 'histories' / 'bench_itHFa_P0_it2.json'}
MNU, VEW = 0.05e-9, 174.1


def main():
    out = {}
    for lab, p in HP.items():
        d, h = C.load_hist(p)
        Hi = d['H_i_GeV']
        mP = C.MPL / Hi
        t, H, a = h['t'], h['H'], h['a']
        rho = 3 * mP ** 2 * H ** 2
        te = d['waterfall']['t_eps1']
        ie = int(np.argmin(np.abs(t - te)))
        f2 = 2000.0
        mu5 = math.sqrt(2 * h['f_chi'][ie] * rho[ie]) / (2 * f2)
        E_N = 0.75 * mu5
        Ms = (400.0, 200.0, 120.0, 250.0)
        y2 = {M: MNU * M * Hi / VEW ** 2 for M in Ms}
        fSM = np.nan_to_num(h['f_SM'])
        fN = np.maximum(np.nan_to_num(h['f_rad']) - fSM, 0.0)
        mN = min(y2.values()) * fN * rho / E_N ** 2
        mSM = 0.0228 * fSM * rho / (E_N / 2) ** 2
        # largest -V'' inside the abundance limit (phibar <= 3.4 H_i canonical) at M_t = 174.6, flat potential
        fl = C.Flat(174.6)
        eG = math.exp(fl.Gam(Hi))
        phis = np.linspace(0.05, 3.4, 60)
        V2max = max(-fl.V2_canon(ph * Hi / eG, Hi) for ph in phis)

        def first(x, thr):
            i = np.where((x >= thr) & (t > 0.1))[0]
            return float(t[i[0]]) if len(i) else None
        o = dict(H_i=Hi, t_eps1=te, mu5max_at_eps1=mu5, E_N_bound=E_N, y_nu2=y2,
                 note_EXP=('the export has f_SM = 0 before t = 0.6 (switch artefact, mf_hybrid skeptic); the SM column '
                           'is therefore only meaningful after 0.6 for EXP') if lab.startswith('EXP') else '',
                 minus_V2_max_within_abundance_limit_174p6=V2max)
        for nm, m2 in (('N_gas_ymin', mN), ('SM_hard_gas', mSM)):
            o[nm] = dict(t_m2_ge={f'{thr:g}': first(m2, thr) for thr in (1.0, 10.0, 100.0, 1e3, 1e4)},
                         m2_at={f'{x}': float(np.interp(x, t, m2)) for x in (0.18, 0.19, 0.195, 0.2, 0.25, 0.3, 0.6, 1.0, 3.0)},
                         m2_over_H2_at={f'{x}': float(np.interp(x, t, m2) / np.interp(x, t, H) ** 2) for x in (0.2, 1.0, 10.0)})
        o['f_rad_at'] = {f'{x}': float(np.interp(x, t, h['f_rad'])) for x in (0.18, 0.19, 0.2, 0.3, 1.0)}
        o['f_SM_at'] = {f'{x}': float(np.interp(x, t, fSM)) for x in (0.18, 0.19, 0.2, 0.3, 1.0)}
        out[lab] = o
        print(lab, 'E_N %.3g' % E_N, 'V2max %.3f' % V2max)
        for nm in ('N_gas_ymin', 'SM_hard_gas'):
            print('   ', nm, o[nm]['t_m2_ge'], {k: '%.3g' % v for k, v in o[nm]['m2_at'].items()})
    C.jdump(out, C.RES / 'sk_plasma.json')


if __name__ == '__main__':
    main()
