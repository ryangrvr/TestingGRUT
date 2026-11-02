#!/usr/bin/env python3
"""dq_entropy item 1 (Claude, own code): one long run of the M1r release (lambda_Y = lambda_Sigma = 0.5, inflaton trilinear,
(m_I, Gamma_I) = (3, 0.5) H_i) with rpdq.c, to H_i t = END, with physical-time bins (20 wide) and conformal-time bins
(0.25 wide in eta = int dt/a).  Usage: long_run.py TAG CPORT_UNITS END [RTOL]
  CPORT_UNITS = c(M_Y) in units of 1/(4 pi^2): 0 = induced portals vanish at M_Y (benchmark); -13.31 = vanish at M_P.
Writes $DQ_SCRATCH/runs/TAG.npz (series) and $DQ_SCRATCH/runs/TAG_meta.json.
"""
import json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import dq_common as dq


def main():
    tag = sys.argv[1]; cu = float(sys.argv[2]); end = float(sys.argv[3])
    rtol = float(sys.argv[4]) if len(sys.argv) > 4 else 1e-10
    cport = cu * dq.CP_UNIT
    p, st, G, du = dq.prepare(0.5, 3.0, cport)
    t0 = time.time()
    nsamp = int(round(end / 5.0))
    rr = dq.run(tag, p, 3.0, 0.5, st['z'], G=G, cport=cport, end=end, rtol=rtol, atol=1e-13 * rtol / 1e-10, nsamp=nsamp, binw=20.0, cbin=0.25)
    meta = dict(tag=tag, cport_units=cu, cport=cport, end=end, rtol=rtol, prep_amps=st['amps'].tolist(), prep_delta=st['delta'],
                prep_eig=sorted(float(x) for x in st['eig']), prep_grad=st['grad'], du_prep=du, G=[G.real, G.imag], I0=rr['I0'],
                END=rr['END'], eta_end=rr['eta_end'], Yfinal=rr['Y'], seconds=time.time() - t0,
                last=dict(zip(dq.COLS, rr['S'][-1].tolist())))
    json.dump(meta, open(dq.SP / 'runs' / (tag + '_meta.json'), 'w'), indent=1, default=float)
    print(tag, 'done', meta['seconds'], meta['END'], 'eta_end', meta['eta_end'], 'a_end', meta['last']['a'])


if __name__ == '__main__':
    main()
