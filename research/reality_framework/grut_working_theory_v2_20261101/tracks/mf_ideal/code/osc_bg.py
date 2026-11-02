#!/usr/bin/env python3
"""mf_ideal sensitivity (Claude, own code): backgrounds whose matter-like interval is a quadratic oscillator of mass m_bg,
released at rest at the slice (so R(0) = 12 = de Sitter, then R = rho (1 + 3 cos 2 m_bg t) on average), decaying into
radiation at Gamma_bg tuned to K.  m_bg = 3, Gamma_bg = 0.5 is M1r's own background; m_bg -> infinity is the
averaged-fluid limit of the primary map (R steps from 12 to rho_m at the slice).

Background generation: rpmfi.c MODE 0 with the hidden fields at zero and G = 0 (an exact solution: the hidden sector
stays at zero), dumped every dt (MFI_DUMP), written as a MODE 5 table "t lna H R".  The holder is a spectator s
(MODE 5, table background, evolved s).  Every choice is a model choice.
Usage: python3 -B osc_bg.py MBG K [TON ...]     -> releases f_osc<MBG>_K<K>_m3_t<TON> (results/release/)
"""
import json, math, sys
import numpy as np
from scipy.optimize import brentq
import mfi_common as mc
from release_mfi import release

Z0 = [0j, 0j, 0j]


def osc_run(mbg, gam, end, nsamp, dump=False):
    return mc.run('osc_bg', 0, mbg, gam, Z0, 0j, end, nsamp, save=False, dump=dump)


def osc_K(mbg, gam):
    end = max(40.0 / gam, 20.0)
    o = osc_run(mbg, gam, end, 200)
    S = o['S']
    return float(S[-1, 1] ** 2 * S[-1, 2])


def gam_for(mbg, K):
    """smallest-width root of K(Gamma) = K (K(Gamma) has a minimum: a heavily damped oscillator released at rest rolls
    slowly); None if K is below the minimum."""
    G = np.geomspace(0.01, 5.0, 25)
    Ks = [osc_K(mbg, g) for g in G]
    for i in range(1, len(G)):
        if Ks[i] < K <= Ks[i - 1]:
            return float(math.exp(brentq(lambda lg: osc_K(mbg, math.exp(lg)) - K, math.log(G[i - 1]), math.log(G[i]), xtol=1e-10)))
    return None


def table(mbg, K):
    gam = gam_for(mbg, K)
    if gam is None:
        raise SystemExit('K = %g is below K_min of an m_bg = %g oscillator released at rest' % (K, mbg))
    dt = 0.001 if mbg > 12 else 0.002
    end = max(20.0 / gam, 30.0)
    o = osc_run(mbg, gam, end, int(round(end / dt)), dump=True)
    D = o['D']; S0 = o['S'][0]
    first = np.r_[0.0, 0.0, o['I0']['I0'], 0.0, np.zeros(16), S0[mc.COLS.index('H')], S0[mc.COLS.index('R')]]
    D = np.vstack([first, D])
    path = mc.SP / ('osc_m%g_K%g.txt' % (mbg, K))
    with open(path, 'w') as fh:
        fh.write('%d 4\n' % len(D))
        for r in D:
            fh.write('%.17g %.17g %.17g %.17g\n' % (r[0], r[1], r[20], r[21]))
    Kend = float(o['S'][-1, 1] ** 2 * o['S'][-1, 2])
    return path, gam, Kend, end


def main():
    mbg, K = float(sys.argv[1]), float(sys.argv[2])
    tons = [float(x) for x in sys.argv[3:]] or [0.1, 1.0]
    path, gam, Kend, end = table(mbg, K)
    print('osc background m_bg=%g K=%g: Gamma_bg=%.6g, K(table end %.0f)=%.6f, %s' % (mbg, K, gam, end, Kend, path), flush=True)
    meta = dict(m_bg=mbg, K=K, Gamma_bg=gam, K_table_end=Kend, table_end=end)
    for ton in tons:
        tag = 'f_osc%g_K%g_m3_t%g' % (mbg, K, ton)
        info = release(tag, K, 3.0, ton, bg='table:%s' % path, bcap=max(3.0, mbg))
        info['osc_background'] = meta
        json.dump(info, open(mc.RES / 'release' / ('release_%s.json' % tag), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
