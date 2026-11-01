#!/usr/bin/env python3
"""fx_holder track (Claude, own code) -- scan of lambda_Y = lambda_Sigma across the perturbative band [0.24, 0.57]
(cf_scale) for the trilinear holder and the unheld original, (3, 0.5) history to H_i t = 1000: preparation (lightest
mass), late charge (400-1000 mean, 80-window spread), yield at the M1 H_i, recalibrated H_i, light-Y debris.
Output: ../results/lam_scan.json
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fx_common as fc


def main():
    rows = []
    for lam in (0.25, 0.3, 0.4, 0.5, 0.566, 0.75, 1.0):
        p = fc.params(lam)
        st = fc.stationary(fc.U_PREP, p)
        G, I0, ph0 = fc.trilinear_G(st['z'], p, 3.0)
        zo, _ = fc.prepare_held(p, 3.0)
        row = dict(lam=lam, prep_amps=st['amps'].tolist(), prep_delta=st['delta'], prep_lightest=float(min(v for v in st['eig'] if abs(v) > 1e-6)),
                   kappa_over_kappa_max=float(p['kap'] / (16 * p['lx'] ** 3 * p['ly'] / 27) ** 0.25),
                   lambda_Y_eff=float(p['ly'] - 27 * p['kap'] ** 4 / (16 * p['lx'] ** 3)))
        for nm, z0, g in (('original', zo, 0j), ('trilinear', st['z'], G)):
            rr = fc.run(f'scan_lam{lam:g}_{nm}', p, 3.0, 0.5, z0=z0, G=g, end=1000.0, nsamp=5000, save=False)
            ws = fc.window_stats(rr, 'NX', 400.0, 1000.0, 80.0)
            Y, ratio, Rr = fc.yield_over_target(rr, ws['mean'])
            led = fc.ledger(fc.window_mean(rr, 'NYH', 400, 1000), fc.window_mean(rr, 'NYL', 400, 1000), fc.window_mean(rr, 'Nsig', 400, 1000),
                            fc.window_mean(rr, 'EX', 400, 1000), Rr, p)
            row[nm] = dict(NX=ws['mean'], run80=[ws['run_min'], ws['run_max']], yield_over_YQ=ratio, H_new_GeV=fc.H_GEV * ratio ** (-2 / 3),
                           NYL=fc.window_mean(rr, 'NYL', 400, 1000), R_Y=led['R_Y'], R_central=led['R_central'], checks=rr['END'])
        row['trilinear_over_original'] = row['trilinear']['NX'] / row['original']['NX']
        rows.append(row)
        print('lam=%.3f k/kmax=%.3f lamYeff=%.4f lightest=%.3f amps=%s | orig Y/YQ=%.5f | tri Y/YQ=%.5f H=%.5e NYL=%.4f R_Y=%.4f | tri/orig=%.5f | checks %s'
              % (lam, row['kappa_over_kappa_max'], row['lambda_Y_eff'], row['prep_lightest'], np.round(row['prep_amps'], 4), row['original']['yield_over_YQ'],
                 row['trilinear']['yield_over_YQ'], row['trilinear']['H_new_GeV'], row['trilinear']['NYL'], row['trilinear']['R_Y'], row['trilinear_over_original'],
                 row['trilinear']['checks']))
    json.dump(dict(rows=rows, note='(3, 0.5) history, end 1000, late = 400-1000 mean; lambda = 0.25 sits at the kappa boundedness floor (0.241)'),
              open(fc.RES / 'lam_scan.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
