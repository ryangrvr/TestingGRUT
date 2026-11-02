"""qm_transfer SKEPTIC: zeta_max on the physical M1r line S/Q = 1.235 (s = 1) from (i) the production grid
(K_rsa = 200, 'thermo'), (ii) converged K_rsa = 800 'thermo', (iii) K_rsa = 800 'sudden' (the right post-freeze-out
physics for Q < Q_c ~ 3e8-1e9, where the objects cannot act as a thermostat, sk_object_mc.py), and the Lyman-alpha
area criterion with k_up = 75 h/Mpc (track) or infinity.  Log-zeta interpolation between grid points.  [C on E proxies]
usage: sk_limits.py krsa_metrics.json OUT.json"""
import json, math, sys
d = json.load(open(sys.argv[1]))
import os
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', '..', 'code'))
import constraints as C
thr_inf = C.metrics(C.viel_T2(5.7, C.KG))['dA_inf']
def get(name, lab):
    for k, r in d['runs'].items():
        if k.rsplit('_K', 1)[0] == name and lab in r:
            return r[lab]['metrics']
    return None
def cross(pts, key, th, larger_bad):
    pts = [(z, m[key]) for z, m in pts if m is not None]
    for (z1, v1), (z2, v2) in zip(pts, pts[1:]):
        bad1 = v1 > th if larger_bad else v1 < th
        bad2 = v2 > th if larger_bad else v2 < th
        if not bad1 and bad2:
            f = (th - v1) / (v2 - v1)
            return z1 * (z2 / z1)**f
    return None
out = {}
zs = (1e-3, 3e-3, 1e-2)
th_lya = C.metrics(C.viel_T2(5.7, C.KG))['dA']
# satellites threshold: viel_6.5 N_gt_1e8
th_sat = C.metrics(C.viel_T2(6.5, C.KG))['N_gt_1e8']; th_lens = C.metrics(C.viel_T2(6.1, C.KG))['dn_1e7p5']
for tag, fam, lab in (('production_K200_thermo', 'P', 'prod_K200'), ('K800_thermo', 'P', 'K800'), ('K800_sudden', 'S', 'K800'),
                      ('production_K200_sudden', 'S', 'prod_K200')):
    pts = [(z, get('%s_SQ1.235_z%.0e' % (fam, z), lab)) for z in zs]
    out[tag] = dict(lya_5p7=cross(pts, 'dA', th_lya, True), lya_5p7_kup_inf=cross(pts, 'dA_inf', thr_inf, True))
    out[tag]['sat_6p5'] = cross(pts, 'N_gt_1e8', th_sat, False)
    out[tag]['lens_6p1'] = cross(pts, 'dn_1e7p5', th_lens, False)
    print(tag, {k: (None if v is None else float('%.3g' % v)) for k, v in out[tag].items()})
out['thresholds'] = dict(dA_5p7=th_lya, dA_inf_5p7=thr_inf, N8_6p5=th_sat, dn_6p1=th_lens)
json.dump(out, open(sys.argv[2], 'w'), indent=1)
