"""Compact late segment in the exact parent frequency/bin convention."""
from pathlib import Path
import argparse
import hashlib
import json
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent/'grut_m1_expansion_clock_20261028'
sys.path.insert(0,str(PARENT/'code'))
from source_ringdown import spectrum,compact_spectrum

def main():
    p=argparse.ArgumentParser();p.add_argument('--tag',choices=['h005','h0025'],required=True);a=p.parse_args()
    h=.05 if a.tag=='h005' else .025
    trace=ROOT/'.work'/f'late_{a.tag}_trace.csv';data=np.genfromtxt(trace,names=True,delimiter=',')
    assert abs(data['tau'][0])<1e-9 and abs(data['tau'][-1]-10000)<1e-9
    diagnostic,raw=spectrum(data,0,10000,h);compact=compact_spectrum(raw)
    dest=ROOT/f'late_{a.tag}_spectrum_10000_20000.npz';np.savez_compressed(dest,**compact)
    out={'source':a.tag,'relative_trace_window':[0,10000],'absolute_mX_t_window':[10000,20000],
         'trace_sha256':hashlib.sha256(trace.read_bytes()).hexdigest(),'radial_h':h,'spectrum':diagnostic,
         'scope':'Signed whole-window field flux and its parent-compatible compact bins. Not an independently certified particle source or all-future emission.'}
    (ROOT/f'late_{a.tag}_compact.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
