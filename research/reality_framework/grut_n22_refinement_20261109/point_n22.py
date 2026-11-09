"""Only N22, massless; unchanged B2 functional and physical solver."""
import argparse,json,time
import numpy as np
from shared import card,run_b2_point,write

def run(output,refined):
    cfg=card()
    if refined:cfg['eigensolver_seed_base']=cfg['refined_seed_base']
    if output.exists() or output.with_suffix('.npz').exists():raise FileExistsError('Preserve prior result')
    start=time.monotonic();result,checkpoint=run_b2_point(cfg,22,mass=0.,refined=refined)
    result['protocol']=cfg['id'];result['seed_base_used']=cfg['eigensolver_seed_base']
    result['empirical_low_eigenvalue_gap']=result['energies'][1]-result['energies'][0]
    result['gap_caveat']='The computed two-eigenvalue gap is diagnostic, not a certified spectral separation.'
    np.savez_compressed(output.with_suffix('.npz'),**checkpoint);write(output,result)
    print(json.dumps({'N':22,'refined':refined,'seconds':time.monotonic()-start}),flush=True)

if __name__=='__main__':
    from pathlib import Path
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--refined',action='store_true')
    a=p.parse_args();run(a.output,a.refined)
