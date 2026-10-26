#!/usr/bin/env python3
"""kin_rates: supplementary grid points for the 4->2 channels A and B around the opening of the on-shell s chain
(sqrt s = 4.25 .. 7.5 m_X), merged into the n->2 tables by build_interface.py.  usage: python3 -B run_supplement.py <out.json>"""
import json, sys, time
import numpy as np
import model as M
import kernels_n2 as K
PTS = [4.25, 4.35, 4.42, 4.47, 4.5, 4.55, 4.65, 4.75, 4.85, 5.0, 5.5, 6.0, 7.5]
out = {}
for ch in ('A_XXXbXb_XXb', 'B_XXXXb_XX', 'D_XXXb_Xs', 'E_XXbs_XXb'):
    rows = []
    rng = np.random.default_rng(abs(hash(('supp', ch))) % (2**32))
    t0 = time.time()
    for sq in PTS:
        if ch == 'D_XXXb_Xs' and sq < 3.5: continue
        for freg, fregX, N in ((0.02, 0.05, 60000), (0.04, 0.10, 30000)):
            r = K.kernel_n2(ch, sq * M.MX, N, rng, freg=freg, fregX=fregX)
            r['sqrt_s_over_mX'] = sq
            rows.append(r)
            print(ch, sq, freg, f"{r['G']:.3e} +- {r['err']:.1e} full {r['G_full']:.3e}", flush=True)
    out[ch] = dict(rows=rows, seconds=time.time() - t0)
json.dump(out, open(sys.argv[1], 'w'), indent=1)
