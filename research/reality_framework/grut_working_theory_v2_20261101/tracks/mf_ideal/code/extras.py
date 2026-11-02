#!/usr/bin/env python3
"""mf_ideal extras (Claude, own code): the short checks recorded in NOTES.md section 6, collected for reproduction.
  holder   : holder-amplitude scan c0 (|G s_*| relative to M1r's mu_eff = 10i calibration) at the best corners
  sstar    : s_* = 1e14, 1e18 GeV; mirror s_* < 0 (c0 = -1); fast decay Gamma_s = 0.5 (K = 1, m_s = 3, t_on = 0.5)
  iso05    : the c0 = 1 +- 0.01 pair at K = 1, m_s = 3, t_on = 0.5, and isocurvature.json update
  rtol C0  : rtol 1e-10 vs 1e-12 to tau = 3000 at K = 1, m_s = 3, t_on = 0.5, c0 = C0 -> results/rtol_check_c<C0>.json
  hartree  : direct relaxed-Hartree solves vs the S/Q closed form -> results/hartree_check.json
After 'holder', 'sstar', 'iso05' run vectors_mfi.py and chain_mfi.py on the new tags (printed).
Usage: python3 -B extras.py holder|sstar|iso05|hartree|rtol C0 [SHARD NSHARD]
"""
import json, sys
import mfi_common as mc
from release_mfi import release


def jobs_holder():
    out = [(1.0, 3.0, 0.5, c) for c in (0.9, 1.03, 1.05, 1.07, 1.1, 1.2)]
    out += [(K, m, t, c) for (K, m, t) in ((1.0, 4.0, 1.0), (1.0, 2.5, 0.5), (1.3, 3.0, 1.0)) for c in (1.03, 1.05, 1.07)]
    return [('f_K%g_m%g_t%g_c%.3f' % (K, m, t, c), dict(K=K, ms=m, ton=t, c0=c)) for K, m, t, c in out]


def jobs_sstar():
    return [('f_K1_m3_t0.5_s1e14', dict(K=1.0, ms=3.0, ton=0.5, sstar=1e14)), ('f_K1_m3_t0.5_s1e18', dict(K=1.0, ms=3.0, ton=0.5, sstar=1e18)),
            ('f_K1_m3_t0.5_mirror', dict(K=1.0, ms=3.0, ton=0.5, c0=-1.0)), ('f_K1_m3_t0.5_g0.5', dict(K=1.0, ms=3.0, ton=0.5, gs=0.5))]


def jobs_iso05():
    return [('f_K1_m3_t0.5_c%.3f' % c, dict(K=1.0, ms=3.0, ton=0.5, c0=c)) for c in (0.99, 1.01)]


def run_jobs(jobs, sh=0, n=1):
    for tag, kw in jobs[sh::n]:
        release(tag, kw['K'], kw['ms'], kw['ton'], gs=kw.get('gs', 0.0), sstar=kw.get('sstar', 1e16), c0=kw.get('c0', 1.0))
    print('tags:', ' '.join(t for t, _ in jobs))


def rtol_check(c0):
    out = {}
    u0 = mc.MU + (mc.U_PREP - mc.MU) * c0
    z = mc.prep(u0, 12.0)['z']; G = mc.U_PREP - mc.MU
    env, gam = mc.history_env(1.0, 3.0, 0.5, br=mc.br_of(1e16), c0=c0)
    for rt in (1e-10, 1e-12):
        o = mc.run('rt', 5, 3.0, gam, z, G, 3000.0, 600, env=env, binw=20.0, rtol=rt, atol=1e-13 * rt / 1e-10, save=False)
        out['%g' % rt] = {k: mc.bin_mean(o['B'], k, 2000, 3000) for k in ('NX', 'EX', 'NYL', 'NYH', 'Nsig')}
        out['%g' % rt]['END'] = o['END']
    r = {k: out['1e-12'][k] / out['1e-10'][k] - 1 for k in ('NX', 'EX', 'NYL', 'NYH', 'Nsig')}
    json.dump(dict(c0=c0, runs=out, rel=r), open(mc.RES / ('rtol_check_c%.3f.json' % c0), 'w'), indent=1, default=float)
    print(c0, r)


def hartree():
    import chain_mfi as cm
    out = {}
    for tag, Ws in (('f_K1_m3_t0.5', [16.0, 49.06, 68.6, 343.0, 1861.0]), ('f_K3_m4_t0.1', [30.0, 100.0, 1409.8, 3832.0, 6865.5])):
        o = json.load(open(mc.RES / 'release' / ('release_%s.json' % tag)))['result']
        out[tag] = cm.hartree_check(o['N_X'], o['E_X'], Ws)
    json.dump(out, open(mc.RES / 'hartree_check.json', 'w'), indent=1, default=float)


def iso_update():
    import chain_mfi as cm
    d = json.load(open(mc.RES / 'isocurvature.json'))
    d['K1_t0.5'] = cm.isocurvature(1.0, 0.5)
    json.dump(d, open(mc.RES / 'isocurvature.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    what = sys.argv[1]
    if what == 'rtol':
        rtol_check(float(sys.argv[2]))
    elif what == 'hartree':
        hartree()
    elif what == 'isoupdate':
        iso_update()
    else:
        sh = int(sys.argv[2]) if len(sys.argv) > 2 else 0
        n = int(sys.argv[3]) if len(sys.argv) > 3 else 1
        run_jobs({'holder': jobs_holder, 'sstar': jobs_sstar, 'iso05': jobs_iso05}[what](), sh, n)
