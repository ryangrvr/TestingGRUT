#!/usr/bin/env python3
"""e2e_chain item 3 (Claude): CMB side at the self-consistent points.

At each (benchmark, H_i, K): N_* from mf_hybrid's Nstar_from_K (read-only) with the final history's own K; the tree
clock (c1, m2) re-calibrated by mf_hybrid's inflation.calibrate at the rescaled V0 (mf_hybrid derived(): lam v^4/4 =
3 M_P^2 H^2 with m_chi = 250, m_sig = 600 H_i fixed; kappa = 10; bias c = 1e5 H_i^3) to the benchmark's target n_s and
A_s at N_*, then the exact Mukhanov-Sasaki spectrum (inflation.ms, read-only) gives A_s, n_s, alpha_s, r; one retarget
pass as in mf_hybrid bench.run if MS misses the target (A_s by > 3e-3 or n_s by > 3e-4).  The e-fold capacity N_max
(largest N before sigma_c for which a retuned monotonic clock meets A_s and n_s, pivot within 0.2 M_P) is mf_chain
cmb_rescale.capacity's construction (the mf_hybrid skeptic's sk_structure.nefolds, read-only, its module mP patched).
n_s is put in by the target, not predicted.  Usage: python3 -B cmb.py NAME:H_GeV:K [...] -> results/cmb.json
"""
import json, math, sys
import numpy as np
import e2c_common as E

sys.path.insert(0, str(E.TRK / 'mf_hybrid' / 'code'))
sys.path.insert(0, str(E.TRK / 'mf_hybrid' / 'skeptic' / 'code'))
import mfh_common as MH                                   # noqa: E402 read-only
import inflation as I                                     # noqa: E402 read-only
import sk_common as SC                                    # noqa: E402 read-only
import sk_structure as SS                                 # noqa: E402 read-only

NS = {'HFa': MH.NS_PLANCK, 'HFb': MH.NS_ACT}


def capacity(H, kappa=10.0, ns=0.9649):
    mP = E.MP_GEV / H
    SC.mP = mP; SS.mP = mP
    P = SS.params(250.0, 600.0, kappa)
    best = (0.0, None)
    for dl in np.geomspace(1e-4, 0.2, 60):
        N = SS.nefolds(P, dl, ns)
        if N > best[0]:
            best = (N, dl)
    return dict(Nmax=best[0], dl_at_Nmax=best[1], X=P['g'] ** 2 * P['sc'])


def ms_at(name, H, Nst):
    ns_t = NS[name]
    p = MH.derived(H, 250.0, 600.0)
    mu2 = 10.0 * p['g'] ** 2 * (2 * p['sigma_c']) ** 2
    inf, o, _ = I.calibrate(p, ns_t, Nst, mu2=mu2, c_bias=1e5)
    msr, _ = I.ms(inf, Nst)
    first = dict(msr)
    if abs(msr['As'] / MH.AS - 1) > 3e-3 or abs(msr['ns'] - ns_t) > 3e-4:
        inf, o, _ = I.calibrate(p, ns_t + (ns_t - msr['ns']), Nst, mu2=mu2, c_bias=1e5, As_t=MH.AS ** 2 / msr['As'])
        msr, _ = I.ms(inf, Nst)
    return dict(ns_target=ns_t, slow_roll=dict(ns=o['ns'], As=o['As'], alpha_s=o.get('alpha_s'), r=o.get('r')),
                ms_first_pass=first, ms=msr, c1=o.get('c1'), m2=o.get('m2'), sigma_star_MP=o.get('sigma_star_MP'),
                Delta_sigma_MP=o.get('Delta_sigma_MP'))


if __name__ == '__main__':
    out = E.load(E.RES / 'cmb.json', {})
    for a in sys.argv[1:]:
        name, H, K = a.split(':'); H, K = float(H), float(K)
        Hs = H * (1 + 2.35e-6)
        Nst = MH.Nstar_from_K(Hs, H, K)
        r = 2 * Hs ** 2 / (math.pi ** 2 * E.A_S * E.MP_GEV ** 2)
        row = dict(name=name, H_i_GeV=H, K=K, Nstar=Nst, r_constH=r)
        try:
            row['capacity'] = capacity(H, 10.0, NS[name]); row['capacity_ok'] = row['capacity']['Nmax'] >= Nst
        except Exception as ex:
            row['capacity_error'] = repr(ex)
        try:
            row['MS'] = ms_at(name, H, Nst)
        except Exception as ex:
            row['MS_error'] = repr(ex)
        out[a] = row
        E.dump(E.RES / 'cmb.json', out)
        m = row.get('MS', {}).get('ms', {})
        print(a, 'N* %.3f r %.4g Nmax %s ok %s | MS ns %s As %s alpha_s %s r %s' % (Nst, r, row.get('capacity', {}).get('Nmax'), row.get('capacity_ok'),
              m.get('ns'), m.get('As'), m.get('alpha_s'), m.get('r')), row.get('MS_error', ''), row.get('capacity_error', ''), flush=True)
