"""Locate the module's mid-H p_- offset: start at the module's start cell instead of v/H (verifier v_inflation)."""
import json, math
import fp1d_verify as F
out = {}
orig = F.V_GEV
for H, ycell in ((1000.0, 0.99), (1368.8745095370805, 0.73), (1670.0415360156744, 0.59), (2037.4685280397082, 0.49)):
    res = {}
    for lab, y0 in (('v_over_H', orig / H), ('module_start_cell', ycell)):
        F.V_GEV = y0 * H   # only moves the start node; potential uses running couplings (m2 from RG), not V_GEV
        r = F.p_minus_extrap(H, [50.0, 60.0], hgoal=0.002)
        res[lab] = {('N%d' % int(N)): r[N]['p'] for N in r}
    F.V_GEV = orig
    out['%.6g' % H] = res
    print(H, res, flush=True)
json.dump(out, open('results_start_offset.json', 'w'), indent=1)
