"""qm_transfer numerical controls (own, Claude 2 Nov 2026) -> results/controls.json.  [C]
 c1 cold collisionless fraction (sigma = 1e-11, f = 0.03 / 0.3) must reproduce CDM (k = 1, 10, 68)
 c2 large-scale adiabaticity of the full kinetic model (P -> 1 at k = 0.3, 1): thermo and sudden, zeta = 1e-4
 c3 resolution: Nq/L = 24/10 (production) vs 36/14 and 16/8 at zeta = 1e-3 (k = 12, 48, 96)
 c4 approximation switches: K_rsa/K_end = 200/400 (production) vs 400/800, and rtol 1e-7 vs 1e-8 (zeta = 1e-3 and 3e-3)
usage: controls.py OUT.json [which]
"""
import json, os, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from qt_bg import Cosmo
from qt_hidden import Vapour, Frozen
from qt_boltz import Mode, CDM


def P(c, dm, k, **kw):
    kw = dict(dict(K_rsa=200., K_end=400., rtol=1e-7), **kw)
    return (Mode(c, dm, k, **kw).run()[1.0] / Mode(c, CDM(), k, **kw).run()[1.0])**2


def main():
    out_fn = sys.argv[1]; which = sys.argv[2].split(',') if len(sys.argv) > 2 else ['c1', 'c2', 'c3', 'c4']
    c = Cosmo()
    res = json.load(open(out_fn)) if os.path.exists(out_fn) else {}
    t0 = time.time()
    if 'c1' in which:
        r = {}
        for f in (0.03, 0.3):
            fr = Frozen(c, fraction=f, dist='MB', sigma=1e-11, Nq=20, Lv=10)
            r['f%g' % f] = {str(k): P(c, fr, k) for k in (1.0, 10.0, 68.0)}
        res['c1_cold_fraction'] = r; print('c1', r, '%.0fs' % (time.time() - t0), flush=True)
    if 'c2' in which:
        r = {}
        for mode in ('thermo', 'sudden'):
            v = Vapour(c, SQ=1.235, Q=1e12, s=1.0, mode=mode, Nq=24, Lv=10)
            r[mode] = {str(k): P(c, v, k) for k in (0.3, 1.0)}
        res['c2_adiabatic_lowk'] = r; print('c2', r, '%.0fs' % (time.time() - t0), flush=True)
    if 'c3' in which:
        r = {}
        for Nq, L in ((24, 10), (36, 14), (16, 8)):
            v = Vapour(c, SQ=1.235, Q=1e9, s=1.0, mode='thermo', Nq=Nq, Lv=L)
            r['Nq%d_L%d' % (Nq, L)] = {str(k): P(c, v, k) for k in (12.0, 48.0, 96.0)}
        res['c3_resolution_zeta1e-3'] = r; print('c3', r, '%.0fs' % (time.time() - t0), flush=True)
    if 'c4' in which:
        r = {}
        for z in (1e-3, 3e-3):
            v = Vapour(c, SQ=1.235, Q=z**-3, s=1.0, mode='thermo', Nq=24, Lv=10)
            r['zeta%g' % z] = {'K200_rtol1e-7': {str(k): P(c, v, k) for k in (8.0, 24.0, 68.0)},
                              'K400_rtol1e-7': {str(k): P(c, v, k, K_rsa=400., K_end=800.) for k in (8.0, 24.0, 68.0)},
                              'K200_rtol1e-8': {str(k): P(c, v, k, rtol=1e-8) for k in (8.0, 24.0, 68.0)}}
        res['c4_switches_tolerance'] = r; print('c4', r, '%.0fs' % (time.time() - t0), flush=True)
    json.dump(res, open(out_fn, 'w'), indent=1)


if __name__ == '__main__':
    main()
