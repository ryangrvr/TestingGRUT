#!/usr/bin/env python3
"""mf_higgs SKEPTIC: post-inflation D with my own modes (sk_modes) on the bench export and mf_chain's rebuilt
histories; M_t,crit (gravitational only), criterion band, alpha_s and M_h variants, deadlines D(tau_th) = 0, and D with
the medium masses of sk_plasma (N gas, SM hard gas) scaled down by large factors.  [C on the argued criterion E]
Requires the runs saved by:  python3 -B sk_modes.py run EXP ; and the MC_* runs (see SKEPTIC.md).
Output: results/sk_postD.json
"""
import sys
sys.dont_write_bytecode = True
import json
import math

import numpy as np

import sk_common as C
import sk_modes as M

MNU, VEW = 0.05e-9, 174.1


def medium(o, hist_path, which):
    """medium mass^2 (H_i^2) at the saves of run o, from the fractions of the history at hist_path (sk_plasma model)"""
    d, h = C.load_hist(hist_path)
    Hi = d['H_i_GeV']
    mP = C.MPL / Hi
    t = h['t']
    rho = 3 * mP ** 2 * h['H'] ** 2
    ie = int(np.argmin(np.abs(t - d['waterfall']['t_eps1'])))
    E_N = 0.75 * math.sqrt(2 * h['f_chi'][ie] * rho[ie]) / 4000.0
    fSM = np.nan_to_num(h['f_SM'])
    fN = np.maximum(np.nan_to_num(h['f_rad']) - fSM, 0)
    if which == 'N':
        m2 = MNU * 120 * Hi / VEW ** 2 * fN * rho / E_N ** 2
    else:
        m2 = 0.0228 * fSM * rho / (E_N / 2) ** 2
    return np.interp(o['t'], t, m2)


def main():
    out = {}
    runs = {'EXP': 'EXP_base', 'MC_TR': 'MC_TR_Hi1.284e+11', 'MC_P0@1.3e11': 'MC_P0_Hi1.3e+11', 'MC_P0@2.4e11': 'MC_P0_Hi2.4e+11'}
    Ms = np.round(np.arange(172.70, 173.21, 0.02), 3)
    for lab, rn in runs.items():
        o = M.load(rn)
        r = dict(N_obs=o['N_obs'], K=o['K'], H_i=o['Hi'], xi=o['xi'])
        Ds = [M.danger(o, m)['D'] for m in Ms]
        r['Mt_crit'] = M.crit(Ms, Ds)
        r['D_at'] = {f'{m}': M.danger(o, m)['D'] for m in (172.28, 172.57, 172.86)}
        band = []
        for kw in (dict(cH=0.5), dict(cH=2.0), dict(kcf=0.5), dict(kcf=2.0)):
            band.append(M.crit(Ms, [M.danger(o, m, **kw)['D'] for m in Ms]))
        r['Mt_crit_band'] = [min(b for b in band if b), max(b for b in band if b)]
        r['Mt_crit_N60'] = M.crit(Ms, [M.danger(o, m, N=60.0)['D'] for m in Ms])
        Mw = np.round(np.arange(172.30, 173.61, 0.05), 3)
        for vl, kw in (('als=0.1171', dict(als=0.1171)), ('als=0.1189', dict(als=0.1189)), ('Mh=125.09', dict(Mh=125.09)),
                       ('Mh=125.31', dict(Mh=125.31))):
            r[f'Mt_crit {vl}'] = M.crit(Mw, [M.danger(o, m, **kw)['D'] for m in Mw])
        r['D(172.57, als=0.1171)'] = M.danger(o, 172.57, als=0.1171)['D']
        # deadlines: tau_th at which D(tau < tau_th) = 0
        dl = {}
        for Mt in (173.0, 173.15, 173.45, 173.75):
            taus = [0.3, 0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 30, 50, 80]
            Dv = [M.danger(o, Mt, tau_th=x)['D'] for x in taus]
            dd = None
            for i in range(len(taus) - 1):
                a, b = Dv[i], Dv[i + 1]
                if b >= 0 and (not np.isfinite(a) or a < 0):
                    dd = taus[i] + ((0 - a) * (taus[i + 1] - taus[i]) / (b - a) if np.isfinite(a) else 0.0)
                    break
            dl[f'{Mt}'] = dd
        r['deadlines'] = dl
        out[lab] = r
        print(lab, json.dumps({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k != 'D_at'}, default=str),
              {k: round(v, 1) for k, v in r['D_at'].items()}, flush=True)
    # medium masses in the criterion (EXP run; fractions from mf_chain's unblocked rebuild, which has the SM column
    # before t = 0.6; and the Pauli run with its own history)
    med = {}
    for lab, rn, hp in (('EXP run / MC_TR fractions', 'EXP_base', M.HPATH['MC_TR']),
                        ('MC_P0@2.4e11', 'MC_P0_Hi2.4e+11', M.HPATH['MC_P0'])):
        o = M.load(rn)
        for which in ('N', 'SM'):
            m2 = medium(o, hp, which)
            for fac in (1.0, 1e-2, 1e-4, 1e-5):
                med[f'{lab} | {which} x {fac:g}'] = {f'{Mt}': M.danger(o, Mt, mp2=fac * m2)['D'] for Mt in (172.86, 173.45, 174.05, 174.6)}
                print(lab, which, fac, {k: (round(v, 1) if np.isfinite(v) else 'no roll') for k, v in med[f'{lab} | {which} x {fac:g}'].items()}, flush=True)
    out['medium_in_criterion'] = med
    C.jdump(out, C.RES / 'sk_postD.json')


if __name__ == '__main__':
    main()
