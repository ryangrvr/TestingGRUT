#!/usr/bin/env python3
"""inf_chain (Claude, own driver): check of dq_entropy's closed form S/Q(W) far above its fit range (W = 13-140) with the
copied relaxed-Hartree solver (condensed branch, continuation) and our no-condensate branch (hartree_nocond.py) above the
condensate's evaporation.  The relaxed state depends only on W [Q dq_entropy]; we use the HT6t release (N_X, E_X).
Output ../results/sq_check.json.  [C]"""
import json, sys, warnings
warnings.filterwarnings('ignore')
import hartree_nocond as hn, chain as cn
import ch_common as ch
rows = hn.scan(75913.2, 1.98078e6, [9.0, 30, 100, 300, 1000, 3000, 5000, 5400], [5400, 6000, 7562, 1e4, 2e4, 5e4, 1.3e5])
for r in rows:
    r['S_per_Q_closed'] = cn.SQ(r['W']); r['rel'] = r['S_per_Q'] / r['S_per_Q_closed'] - 1
json.dump(dict(method=__doc__, rows=rows, max_abs_rel_converged=max(abs(r['rel']) for r in rows if r['resid'] < 1e-10)),
          open(ch.RES / 'sq_check.json', 'w'), indent=1)
print('max |rel| (converged rows) = %.2e' % max(abs(r['rel']) for r in rows if r['resid'] < 1e-10))
