"""Item 3d of track higgs_reheat: the xi window edges (gravitational-only post-inflation criterion D_N = 0), and the
benchmark xi_EW = 0.1 at each (M_t, H_i).  Own code.  Writes ../results/edges.json (merged per call).
Run: python3 -B edges.py <Mt> <H_GeV>   (~3-5 min each).
xi here is xi(mu ~ H_i) (the value the post-inflation mode equations use; its running over H_i/10 ... 10 H_i is
< 1%).  Conversion to xi_EW = xi(M_t): (xi_EW - 1/6) = (xi_H - 1/6) exp(-int_{M_t}^{H_i} gamma_xi dln mu).
"""
import json
import os
import sys

import numpy as np

import modes as M

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')


def main(Mt, H_GeV):
    path = os.path.join(RES, 'edges.json')
    out = json.load(open(path)) if os.path.exists(path) else {}
    lam = M.Lam(Mt)
    sm = lam.sm
    toH = lambda xe: float(sm.xi(H_GeV, xe))
    c = sm.cpl(H_GeV); c0 = sm.cpl(Mt)
    fac = float(np.exp(c['Ixi'] - c0['Ixi']))
    toEW = lambda xh: float(1 / 6 + (xh - 1 / 6) / fac)
    cache = {}

    def D(xi, N):
        if xi not in cache:
            r, _ = M.run(xi, Mt=Mt, H_GeV=H_GeV, lam=lam, nsave=800, cHs=(0.5, 2.0))
            cache[xi] = r
        r = cache[xi]
        return r['danger']['best'][N]['D']

    def bisect(lo, hi, N, it=11):
        dlo, dhi = D(lo, N), D(hi, N)
        if np.sign(dlo) == np.sign(dhi):
            return None
        for _ in range(it):
            mid = np.sqrt(lo * hi)
            dm = D(mid, N)
            if np.sign(dm) == np.sign(dlo):
                lo, dlo = mid, dm
            else:
                hi, dhi = mid, dm
        return float(np.sqrt(lo * hi))
    res = dict(Mt=Mt, H_GeV=H_GeV, xi_H_over_xi_EW_running_factor=fac, h_max=lam.flat['h_max'])
    xb = toH(0.1)
    rb, _ = M.run(xb, Mt=Mt, H_GeV=H_GeV, lam=lam, nsave=800, cHs=(0.5, 2.0))
    res['benchmark'] = dict(xi_EW=0.1, xi_H=xb, D=rb['danger']['best'], D_cH=rb['danger_cH'],
                            max_2sigma_over_hmax=2 * rb['max_sigma2'] ** 0.5 / rb['h_max_over_H'])
    for N in (50, 60):
        lo = bisect(0.02, 0.3, N)
        hi = bisect(10.0, 60.0, N)
        res[f'N={N}'] = dict(xi_H_lower=lo, xi_H_upper=hi, xi_EW_lower=(toEW(lo) if lo else None),
                             xi_EW_upper=(toEW(hi) if hi else None))
        print(Mt, H_GeV, N, res[f'N={N}'], flush=True)
    # mean-field criterion 2 sigma_max = h_max (literature-style), from the cached runs
    xs = sorted(cache)
    mf = [(x, 2 * cache[x]['max_sigma2'] ** 0.5 / cache[x]['h_max_over_H']) for x in xs]
    up = [x for x, v in mf if v >= 1.0 and x > 1]
    res['mean_field_2sigma_eq_hmax_xi_H'] = min(up) if up else None
    res['scan'] = {f'{x:.5g}': dict(D60=cache[x]['danger']['best'][60]['D'],
                                    D60_cH05=cache[x]['danger_cH']['cH=0.5'][60]['D'],
                                    D60_cH2=cache[x]['danger_cH']['cH=2.0'][60]['D'],
                                    two_sigma_max_over_hmax=v) for x, v in mf}
    # cH sensitivity of the edges (N = 60) from the cached scan: sign changes
    for lab in ('cH=0.5', 'cH=2.0'):
        vals = [(x, cache[x]['danger_cH'][lab][60]['D']) for x in xs]
        lows = [x for x, d in vals if d < 0 and x < 1]
        highs = [x for x, d in vals if d >= 0 and x > 1]
        res[f'N=60 {lab} (from scan)'] = dict(smallest_safe_xi_H_below_1=(min(lows) if lows else None),
                                              smallest_unsafe_xi_H_above_1=(min(highs) if highs else None))
    out[f'Mt={Mt} H={H_GeV:.3g}'] = res
    json.dump(out, open(path, 'w'), indent=1, default=float)
    print(json.dumps({k: v for k, v in res.items() if k != 'scan'}, default=float, indent=1)[:3000])


if __name__ == '__main__':
    main(float(sys.argv[1]), float(sys.argv[2]))
