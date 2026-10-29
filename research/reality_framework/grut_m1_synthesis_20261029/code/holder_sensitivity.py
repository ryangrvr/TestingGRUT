"""Homogeneous phase-response of the proposed curvature-holder trajectory.

Rotate Sigma at fixed field magnitudes and zero velocities; this is one chosen
direction, not a fully propagated gauge/metric/isocurvature calculation.
"""
from pathlib import Path
import json
import numpy as np
from holder_evolution import run
ROOT=Path(__file__).resolve().parents[1]

def main():
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--end',type=float,default=100);p.add_argument('--rtol',type=float,default=2e-9);a=p.parse_args()
 v=json.loads((ROOT/'phase_holder.json').read_text());state=next(v for v in v['states'] if v['xi']==4e8);base=np.array([complex(*v) for v in state['fields_over_H']]);rows=[];arrays={}
 for d in (0.,-.001,.001,-.0005,.0005):
  z=base.copy();z[2]*=np.exp(1j*d);name=f'delta_{d:+.4f}';print('start',name,flush=True)
  r,arr=run(name,z,4e8,a.end,a.rtol);r['delta_sigma']=d;rows.append(r);arrays[name]=arr;print('done',name,r['late_mean_NX'],flush=True)
 def slope(step):
  plus=next(r for r in rows if r['delta_sigma']==step)['late_mean_NX'];minus=next(r for r in rows if r['delta_sigma']==-step)['late_mean_NX'];return (np.log(abs(plus))-np.log(abs(minus)))/(2*step)
 out={'end_t_Hi':a.end,'rtol':a.rtol,'rows':rows,'dln_abs_NX_dtheta_sigma_step_0p001':float(slope(.001)),
      'dln_abs_NX_dtheta_sigma_step_0p0005':float(slope(.0005)),
      'scope':'Finite-difference homogeneous response at a fixed local minimum under one Sigma phase perturbation, in a spectator background. No full metric/gauge perturbation propagation, stochastic evolution, non-Gaussianity, or final charge-population mapping.'}
 np.savez_compressed(ROOT/'holder_sensitivity.npz',**arrays);(ROOT/'holder_sensitivity.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
if __name__=='__main__':main()
