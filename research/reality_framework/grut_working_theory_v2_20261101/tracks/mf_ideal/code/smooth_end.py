#!/usr/bin/env python3
"""mf_ideal sensitivity (Claude, own code): a SMOOTH end of inflation of duration t_end.

Background (model choice): at the slice (a = 1, H = H_i = 1, R = 12) the energy is vacuum-like, rho_v = 3; it converts at
rate 1/t_end into radiation (fraction 1 - f_d) and dust (f_d), the dust decays into radiation at Gamma_bg.  R = 4 rho_v +
rho_m (trace), H from 3H^2 = rho.  t_end -> 0 with f_d = 0 is the primary map's K = 1 history (R steps from 12 to 0).
The table "t lna H R" (dt = 0.0005 for t < 5, then 0.005) is fed to rpmfi.c MODE 5 with an evolved spectator.
Usage: python3 -B smooth_end.py TEND [TON ...]   (f_d = 0)   -> releases f_send<TEND>_m3_t<TON>
"""
import json, math, sys
import numpy as np
from scipy.integrate import solve_ivp
import mfi_common as mc
from release_mfi import release


def background(tend, fd=0.0, gam=1.0, tmax=60.0):
    gv = 1.0 / tend

    def rhs(t, y):
        la, rv, rm, rr = y
        H = math.sqrt((rv + rm + rr) / 3)
        return [H, -gv * rv, fd * gv * rv - (3 * H + gam) * rm, (1 - fd) * gv * rv + gam * rm - 4 * H * rr]
    t = np.r_[np.arange(0, 5, 0.0005), np.arange(5, tmax + 1e-9, 0.005)]
    s = solve_ivp(rhs, (0, tmax), [0.0, 3.0, 0.0, 0.0], t_eval=t, method='DOP853', rtol=1e-12, atol=1e-15)
    la, rv, rm, rr = s.y
    H = np.sqrt((rv + rm + rr) / 3)
    R = 4 * rv + rm
    return s.t, la, H, R


def table(tend):
    t, la, H, R = background(tend)
    path = mc.SP / ('send_%g.txt' % tend)
    with open(path, 'w') as fh:
        fh.write('%d 4\n' % len(t))
        for i in range(len(t)):
            fh.write('%.17g %.17g %.17g %.17g\n' % (t[i], la[i], H[i], R[i]))
    return path, float(np.exp(2 * la[-1]) * H[-1])


def main():
    tend = float(sys.argv[1])
    tons = [float(x) for x in sys.argv[2:]] or [0.1, 0.5, 1.0]
    path, K = table(tend)
    print('smooth end t_end=%g: K=%.5f' % (tend, K), flush=True)
    for ton in tons:
        tag = 'f_send%g_m3_t%g' % (tend, ton)
        info = release(tag, K, 3.0, ton, bg='table:%s' % path)
        info['smooth_end'] = dict(t_end=tend, f_d=0.0, K=K)
        json.dump(info, open(mc.RES / 'release' / ('release_%s.json' % tag), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
