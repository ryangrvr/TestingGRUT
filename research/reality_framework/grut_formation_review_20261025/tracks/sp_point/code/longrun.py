"""sp_point: long integrations of the leading point (a >> 3000) to test that the charge fraction is final (own code).

Usage: python3 longrun.py ../results/longrun.json label method cpar a_end
"""
import json
import os
import sys
import time
import numpy as np
import spmodel as m

H = m.H_STAR
WINS = (1, 300, 600, 1000, 2000, 3000, 5000, 10000, 20000, 30000, 50000, 100000)


def main(out, label, method, cpar, aend):
    wins = tuple(w for w in WINS if w <= aend)
    z = m.initial_z(H, np.pi / 2)
    t0 = time.time()
    r = m.run(H, z, a_end=aend, method=method, cpar=cpar, windows=wins, ns=400)
    keep = ['epsX', 'NX', 'EX', 'JY', 'JYeig', 'NY', 'Neff', 'Neff_std', 'epsX_std', 'NX_std', 'YX_over_YQ',
            'eps_of_means', 'pairs_per_charge', 'EYexc']
    rec = dict(label=label, method=method, cpar=cpar, a_end=aend, runtime_s=time.time() - t0,
               controls=r['controls'], windows={k: {q: v[q] for q in keep} for k, v in r['windows'].items()},
               history=[{q: h[q] for q in ['a', 'NX', 'epsX', 'EX', 'JY', 'Neff', 'Gauss']} for h in r['history']])
    allr = json.load(open(out)) if os.path.exists(out) else {}
    allr[label] = rec
    json.dump(allr, open(out, 'w'), indent=1)
    for k, v in rec['windows'].items():
        print(label, k, 'eps=%.8f NX=%.6f EX=%.5f JY=%.6f JYeig=%.6f Neff=%.4f' % (v['epsX'], v['NX'], v['EX'], v['JY'], v['JYeig'], v['Neff']))
    print(label, rec['controls'], rec['runtime_s'])


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5]))
