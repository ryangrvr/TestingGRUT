#!/usr/bin/env python3
"""fx_holder track (Claude, own code) -- controls: rtol 1e-12 reruns of the (5, 1) Option 1A and trilinear histories, and
Option 1A in rptri.c against the unmodified m1_holder/code/rpholder.c binary (built in scratch). Output: ../results/controls.json.
Run from this directory after building rpholder_ref:  gcc -O3 -std=c11 -o $FX_SCRATCH/build/rpholder_ref ../../../../grut_formation_review_20261025/tracks/m1_holder/code/rpholder.c -lm
"""
import sys, json, subprocess; sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import numpy as np, fx_common as fc
out = {}
for lam in (1.0, 0.5):
    p = fc.params(lam); st = fc.stationary(fc.U_PREP, p); zh = st['z']
    G, I0, ph0 = fc.trilinear_G(zh, p, 5.0)
    for nm, lock, g in (('option1A', fc.LOCK_1A, 0j), ('trilinear', 0j, G)):
        a = fc.run(f'ctl_lam{lam:g}_b5_{nm}_r10', p, 5.0, 1.0, z0=zh, lock=lock, G=g, end=1000.0, nsamp=5000, save=False)
        b = fc.run(f'ctl_lam{lam:g}_b5_{nm}_r12', p, 5.0, 1.0, z0=zh, lock=lock, G=g, end=1000.0, nsamp=5000, rtol=1e-12, atol=1e-15, save=False)
        na, nb = fc.window_mean(a, 'NX', 400, 1000), fc.window_mean(b, 'NX', 400, 1000)
        out[f'lam{lam:g}_b5_g1_{nm}'] = dict(NX_rtol10=na, NX_rtol12=nb, rel=nb / na - 1, NYL10=fc.window_mean(a, 'NYL', 400, 1000), NYL12=fc.window_mean(b, 'NYL', 400, 1000), checks12=b['END'])
        print(nm, lam, na, nb, nb / na - 1, out[f'lam{lam:g}_b5_g1_{nm}']['NYL10'], out[f'lam{lam:g}_b5_g1_{nm}']['NYL12'], b['END'])
    # independent binary: m1_holder rpholder.c (unmodified) for Option 1A
    for b_, gam in ((5.0, 1.0), (3.0, 0.5)):
        gfmt = lambda v: '%.17g' % float(v)
        args = [str(fc.SP / 'build' / 'rpholder_ref'), '0', gfmt(b_), gfmt(gam), gfmt(p['xi']), gfmt(fc.ETA2), gfmt(p['cH']), gfmt(p['lx']), gfmt(p['ly']), gfmt(p['ls']),
                gfmt(p['kap']), gfmt(p['M']), gfmt(p['h']), gfmt(p['mu']), gfmt(fc.LOCK_1A.real), gfmt(fc.LOCK_1A.imag), gfmt(800.0), gfmt(1e-10), gfmt(1e-13), '4000', 'F'] + \
               [gfmt(v) for c in zh for v in (c.real, c.imag)] + ['0'] * 6
        o = subprocess.run(args, capture_output=True, text=True).stdout
        W80 = [l for l in o.split('\n') if l.startswith('W80')][0].split()
        ref = float(W80[2 + 4])
        mine = fc.run('ctl_1A_800', p, b_, gam, z0=zh, lock=fc.LOCK_1A, end=800.0, nsamp=4000, save=False)['W']['W80']['NX']
        out[f'lam{lam:g}_b{b_:g}_1A_vs_rpholder'] = dict(rpholder_W80=ref, rptri_W80=mine, rel=mine / ref - 1)
        print('rpholder check', lam, b_, ref, mine, mine / ref - 1)
json.dump(out, open(fc.RES / 'controls.json', 'w'), indent=1, default=float)
