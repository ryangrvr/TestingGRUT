#!/usr/bin/env python3
"""mf_chain: isocurvature derivative D = d ln Omega_DM / d ln s_* on an existing bench run (release + vectors + chain
already done under TAG): two extra releases at c0 = 0.99, 1.01 (fixed g_s, preparation re-solved, same table, s evolves
with the bench's m_s^2(t)), own vectors, own chain; Omega_DM ~ M |N_X| K^{-3/2} at fixed H_i (mf_ideal's recipe).
Planck 2018 (arXiv:1807.06211): beta_iso < 0.038 (uncorrelated CDI, n = 1, k = 0.05/Mpc) -> s_* > |D| H_i/(2 pi sqrt(beta A_s/(1-beta))).
Usage: python3 -B iso.py TAG [TAG ...] -> results/engine/chain/iso_<TAG>.json"""
import json, math, sys
import mfc_common as M
from run_bench import iso_releases

if __name__ == '__main__':
    for tag in sys.argv[1:]:
        r = iso_releases(tag, None)
        rel = json.load(open(M.ENG / 'release' / ('release_%s.json' % tag)))
        sst = rel['spec']['s_star_GeV']
        bench = json.load(open(rel['spec']['bg'].split(':', 1)[1].replace('benchtab_', 'XX'))) if False else None
        M.dump(M.ENG / 'chain' / ('iso_%s.json' % tag), r)
        for lab, x in r.items():
            print('ISO %s %s: D = %.3f (charge only %.3f); branch %s; s_* > %.4g H_i' % (tag, lab, x['D'], x['D_charge_only'],
                                                                                       x['branch_c1'], x['s_star_min_per_H_i']), flush=True)
