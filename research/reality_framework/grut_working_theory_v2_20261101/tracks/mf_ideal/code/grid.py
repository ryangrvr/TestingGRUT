#!/usr/bin/env python3
"""mf_ideal (Claude, own code): the history grid (model choices) and a sequential runner of release_mfi.release.

Primary map ('P'): fluid background (abrupt end), K in {1.0, 1.3, 1.6, 2.0, 3.0}, m_s in {2, 3, 4}, t_on in {0.1, 1},
                   Gamma_s = 0 (no decay), s_* = 1e16 GeV, c0 = 1.
Small decay ('D'):  as P with m_s = 3 and Gamma_s = 0.05 H_i.
Fine scan ('F'):    K = 1, m_s in {2, 2.5, 3, 3.5, 4}, t_on in {0, 0.3, 0.5, 2, 3} (+ m_s = 2.5, 3.5 at t_on = 0.1, 1).
xi_s tracking ('X'): m_s^2(t) = m_s^2 (1 - R(t)/12) (G5 option b), all K and m_s.
Isocurvature ('I'): as P with m_s = 3, c0 = 1 +- 0.01 (fixed g_s; the preparation is re-solved at the shifted u).
Usage: python3 -B grid.py SET [SHARD NSHARD]   (SET = P, D, I or a comma list of tags); skips tags already in results.
"""
import sys
import mfi_common as mc
from release_mfi import release

KS = (1.0, 1.3, 1.6, 2.0, 3.0)
MS = (2.0, 3.0, 4.0)
TON = (0.1, 1.0)


def tag(K, ms, ton, gs=0.0, c0=1.0, extra=''):
    t = 'f_K%g_m%g_t%g' % (K, ms, ton)
    if gs:
        t += '_g%g' % gs
    if c0 != 1.0:
        t += '_c%.3f' % c0
    return t + extra


def grid(name):
    out = []
    if name == 'P':
        out = [dict(tag=tag(K, m, t), K=K, ms=m, ton=t) for K in KS for m in MS for t in TON]
    elif name == 'D':
        out = [dict(tag=tag(K, 3.0, t, gs=0.05), K=K, ms=3.0, ton=t, gs=0.05) for K in KS for t in TON]
    elif name == 'F':                                    # fine (m_s, t_on) scan at K = 1 (t_on = 0: instant turn-on)
        out = [dict(tag=tag(1.0, m, t), K=1.0, ms=m, ton=t) for m in (2.0, 2.5, 3.0, 3.5, 4.0) for t in (0.0, 0.3, 0.5, 2.0, 3.0)]
        out += [dict(tag=tag(1.0, m, t), K=1.0, ms=m, ton=t) for m in (2.5, 3.5) for t in TON]
    elif name == 'X':                                    # xi_s tracking (G5 b): m_s^2(t) = m_s^2 (1 - R(t)/12), no ramp
        out = [dict(tag=tag(K, m, 0.0, extra='_xi'), K=K, ms=m, ton=0.0, xis=-m * m / 12.0) for K in KS for m in MS]
    elif name == 'I':
        out = [dict(tag=tag(K, 3.0, t, c0=c0), K=K, ms=3.0, ton=t, c0=c0) for K in KS for t in TON for c0 in (0.99, 1.01)]
    return out


def main():
    name = sys.argv[1]
    shard = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    nsh = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    jobs = grid(name)[shard::nsh]
    for j in jobs:
        if (mc.RES / 'release' / ('release_%s.json' % j['tag'])).exists():
            print('skip', j['tag']); continue
        release(j['tag'], j['K'], j['ms'], j['ton'], gs=j.get('gs', 0.0), c0=j.get('c0', 1.0), xis=j.get('xis', 0.0))


if __name__ == '__main__':
    main()
