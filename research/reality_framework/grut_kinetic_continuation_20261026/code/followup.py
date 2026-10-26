"""Numerical refinements and labelled gas-only diagnostic source scan."""
from concurrent.futures import ProcessPoolExecutor,as_completed
import json
from bootstrap import ROOT
from moment_controls import Control,save

def task(s):
    label,R,form,rtol,method=s
    try:return label,Control(R,form,recycle=form!='gas',rtol=rtol,method=method).run()
    except Exception as e:return label,dict(error=repr(e),complete=False)

def main():
    specs=[('gas_band_coarse',41.707077,'gas',2e-7,'Radau'),('gas_band_BDF',41.707077,'gas',2e-8,'BDF'),
           ('instant_band_coarse',41.707077,'instant',2e-7,'Radau')]
    specs += [(f'gas_scan_R{r:g}',r,'gas',2e-8,'Radau') for r in [1.,5.,7.,10.,20.,40.5626228853,48.5,90.]]
    out={}
    with ProcessPoolExecutor(max_workers=3) as pool:
        fs=[pool.submit(task,s) for s in specs]
        for f in as_completed(fs):
            label,r=f.result();out[label]=r;save('FOLLOWUP_CONTROLS.json',dict(cases=out,physical_threshold=None))
            print('FOLLOWUP',label,r.get('final',r.get('error')),flush=True)

if __name__=='__main__':main()
