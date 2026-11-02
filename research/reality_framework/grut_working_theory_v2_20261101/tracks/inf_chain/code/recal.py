#!/usr/bin/env python3
"""inf_chain (Claude, own driver): one-step consistency iteration -- re-calibrate a benchmark family at the H_i the
dark-matter calibration requires, using inf_potential's calibration routines (calibrate.py, imported read-only):
HT6t: (v, beta) for A_s = 2.1e-9 and n_s = 0.9649 with the Gamma_I = 0.5 history (H(S1) = H_i);  T1: alpha for A_s.
Output ../results/recal_<KEY>.json (shape, V0, N_*, n_s, r, alpha_s, H_*, m_I/H_i).  Usage: recal.py KEY family H_i [v0 beta0]
Labels: [C] given the model choices."""
import json, math, sys
import ch_common as ch
cal, ic = ch.cal, ch.ic


def main(key, fam, Hi, v0=None, b0=None):
    if fam == 'HT6t':
        m, r = cal.calibrate_hilltop_ns(Hi, ('gamma', 0.5), 6, x0=(v0, b0), verbose=True)
        params = dict(vev=m.vev, p=6, beta=m.beta)
    elif fam == 'T1':
        m, r = cal.calibrate_alpha('T', Hi, ('gamma', 0.5))
        params = dict(alpha=m.p['alpha'])
    else:
        raise KeyError(fam)
    mI = math.sqrt(m.m2_min() * r['V0']) * ic.MPL_GEV / Hi
    out = dict(key=key, family=fam, H_i_GeV=Hi, params=params, V0_MP4=r['V0'], Nstar=r['Nstar'], As=r['As'], ns=r['ns'], r=r['r'],
               alpha_s=r['alpha_s'], Hstar_GeV=r['Hstar_GeV'], m_I_over_Hi=mI, history='Gamma_I = 0.5 H_i', method=__doc__)
    json.dump(out, open(ch.RES / ('recal_%s.json' % key), 'w'), indent=1, default=float)
    print(json.dumps({k: v for k, v in out.items() if k != 'method'}, default=float))


if __name__ == '__main__':
    a = sys.argv[1:]
    main(a[0], a[1], float(a[2]), *(float(x) for x in a[3:5]))
