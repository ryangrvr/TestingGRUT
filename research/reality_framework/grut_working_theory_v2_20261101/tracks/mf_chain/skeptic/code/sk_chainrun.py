#!/usr/bin/env python3
"""mf_chain SKEPTIC: full chain (mf_ideal's engine + vectors + ledger, read-only, outputs in the skeptic's scratch) on a
bench-format history, optionally with the isocurvature releases c0 = 0.99, 1.01.  Used for (i) the isocurvature at the
e = 3e-6 self-consistent points, which the track did not evaluate there, and (ii) the full chain on the skeptic's m_s
variants.  Usage: python3 -B sk_chainrun.py BENCH.json TAG [--iso] [--R raw|trace]
Labels: [C] given the model choices.
"""
import argparse, json, math, sys
import sk_common as S
import bench_run as br                                   # mf_ideal, read-only


def iso(tag):
    rel0 = json.load(open(S.ENG / 'release' / ('release_%s.json' % tag)))
    sp = rel0['spec']
    orig = S.rm.spec_env

    def spec_env(*a, **kw):
        env, gam = orig(*a, **kw)
        env['MFI_SM2TAB'] = '1'
        return env, gam
    S.rm.spec_env = spec_env
    out = {}
    try:
        for c0 in (0.99, 1.01):
            t2 = '%s_c%.2f' % (tag, c0)
            S.rm.release(t2, sp['K_target'], sp['m_s'], 0.0, gs=sp['Gamma_s'], sstar=sp['s_star_GeV'], c0=c0, bg=sp['bg'],
                         Tl=rel0['result']['T_long'], bcap=sp['step_cap_b'])
            S.vm.main(t2)
            out[c0] = S.cm.analyse(t2)
    finally:
        S.rm.spec_env = orig
    res = {}
    for lab in ('e1e-6', 'e3e-6'):
        lnO = {c0: math.log(o[lab]['M'] * abs(o['N_X']) * o['K'] ** -1.5) for c0, o in out.items()}
        D = (lnO[1.01] - lnO[0.99]) / (math.log(1.01) - math.log(0.99))
        DN = (math.log(abs(out[1.01]['N_X'])) - math.log(abs(out[0.99]['N_X']))) / (math.log(1.01) - math.log(0.99))
        res[lab] = dict(D=D, D_charge_only=DN, S_per_Q={str(c): o[lab]['S_per_Q'] for c, o in out.items()},
                        M={str(c): o[lab]['M'] for c, o in out.items()})
    return res


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('bench'); ap.add_argument('tag'); ap.add_argument('--iso', action='store_true')
    ap.add_argument('--R', default='raw'); ap.add_argument('--tlong', default='1e4')
    a = ap.parse_args()
    sys.argv = ['bench_run.py', a.bench, '--tag', a.tag, '--R', a.R, '--tlong', a.tlong]
    br.main()
    if a.iso:
        r = iso(a.tag)
        S.dump(S.ENG / 'chain' / ('iso_%s.json' % a.tag), r)
        for lab, x in r.items():
            print('ISO %s %s: D = %.3f (charge only %.3f)' % (a.tag, lab, x['D'], x['D_charge_only']), flush=True)
