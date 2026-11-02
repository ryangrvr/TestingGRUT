#!/usr/bin/env python3
"""mf_chain: run mf_ideal's whole chain (bench_run.py, unchanged, read-only import) on one bench-format history, with the
outputs redirected to this track (mfc_common).  Optional isocurvature releases (c0 = 1 +- 0.01, preparation re-solved,
own vectors), exactly as mf_ideal's chain_mfi.isocurvature does for its grid but on the bench table.

Usage: python3 -B run_bench.py BENCH_JSON TAG [--R raw|trace] [--iso] [--tlong 1e4] [--gs G]
Labels: [C] given the model choices (history, ledger, closed form, p map are inherited [Q]).
"""
import argparse, json, math, sys
import numpy as np
import mfc_common as M
import bench_run as br                                   # mf_ideal, read-only


def iso_releases(tag, base_diag):
    """two extra releases at c0 = 0.99, 1.01 on the same MODE 5 table (s evolves with the bench's m_s^2(t))."""
    rel0 = json.load(open(M.ENG / 'release' / ('release_%s.json' % tag)))
    sp = rel0['spec']
    orig = M.rm.spec_env

    def spec_env(*a, **kw):
        env, gam = orig(*a, **kw)
        env['MFI_SM2TAB'] = '1'
        return env, gam
    M.rm.spec_env = spec_env
    out = {}
    try:
        for c0 in (0.99, 1.01):
            t2 = '%s_c%.2f' % (tag, c0)
            M.rm.release(t2, sp['K_target'], sp['m_s'], 0.0, gs=sp['Gamma_s'], sstar=sp['s_star_GeV'], c0=c0, bg=sp['bg'],
                         Tl=rel0['result']['T_long'], bcap=sp['step_cap_b'])
            M.vm.main(t2)
            o = M.cm.analyse(t2)
            M.dump(M.ENG / 'chain' / ('chain_%s.json' % t2), o)
            out[c0] = o
    finally:
        M.rm.spec_env = orig
    o1 = M.cm.analyse(tag)
    res = {}
    amp = math.sqrt(M.cm.BETA_ISO * M.cm.AS / (1 - M.cm.BETA_ISO))
    for lab in ('e1e-6', 'e3e-6'):
        v = {}
        for c0, o in out.items():
            Mm = o[lab]['M']
            v[c0] = dict(lnOmega=math.log(Mm * abs(o['N_X']) * o['K'] ** -1.5), S_per_Q=o[lab]['S_per_Q'], branch=o[lab]['branch'],
                         N_X=o['N_X'], M=Mm)
        D = (v[1.01]['lnOmega'] - v[0.99]['lnOmega']) / (math.log(1.01) - math.log(0.99))
        DN = (math.log(abs(v[1.01]['N_X'])) - math.log(abs(v[0.99]['N_X']))) / (math.log(1.01) - math.log(0.99))
        res[lab] = dict(D=D, D_charge_only=DN, pm=v, branch_c1=o1[lab]['branch'],
                        s_star_min_per_H_i=abs(D) / (2 * math.pi * amp))
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('bench'); ap.add_argument('tag')
    ap.add_argument('--R', default='raw'); ap.add_argument('--win', default='0.05')
    ap.add_argument('--iso', action='store_true'); ap.add_argument('--tlong', default='1e4'); ap.add_argument('--gs')
    a = ap.parse_args()
    argv = ['bench_run.py', a.bench, '--tag', a.tag, '--R', a.R, '--win', a.win, '--tlong', a.tlong]
    if a.gs is not None:
        argv += ['--gs', a.gs]
    sys.argv = argv
    br.main()
    if a.iso:
        iso = iso_releases(a.tag, None)
        M.dump(M.ENG / 'chain' / ('iso_%s.json' % a.tag), iso)
        for lab, r in iso.items():
            print('ISO %s %s: D = %.3f (charge only %.3f), s_* > %.3g H_i' % (a.tag, lab, r['D'], r['D_charge_only'], r['s_star_min_per_H_i']))


if __name__ == '__main__':
    main()
