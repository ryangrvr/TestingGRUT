#!/usr/bin/env python3
"""kin_rates: the 2->n stage of run_kernels.py on its own (same settings), usage: python3 -B run_2n_only.py <out_dir> [nproc]"""
import json, sys, time
from multiprocessing import Pool
import kernels_n2 as K
import run_kernels as RK

if __name__ == '__main__':
    nproc = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    twon = sorted(K.CHANNELS_2N, key=lambda c: -len(K.CHANNELS_2N[c]['legs']))
    t1 = time.time()
    with Pool(nproc) as pool:
        res_2n = dict(pool.imap_unordered(RK.run_2n, twon, chunksize=1))
    json.dump(dict(N6=RK.N6, N5=RK.N5, fregs=RK.FREGS, npts=RK.NPTS, channels=res_2n, seconds=time.time() - t1),
              open(f'{RK.OUT}/xsec_2n.json', 'w'), indent=1)
    print('2->n done', time.time() - t1, flush=True)
