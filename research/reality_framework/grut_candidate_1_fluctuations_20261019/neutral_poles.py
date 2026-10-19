#!/usr/bin/env python3
"""Strict scalar vacuum pole masses by full real-vertex tensor contraction."""
from pathlib import Path
import importlib.util
import hashlib
import json
import time
import numpy as np

HERE=Path(__file__).resolve().parent
PREV=HERE.parent/'grut_candidate_1_quantum_20261018'


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
    return mod


def main():
    start=time.perf_counter()
    protocol=(HERE/'neutral_pole_protocol.json').read_bytes()
    p=load('grut_pole_integrals',HERE/'pole.py')
    e=load('grut_potential_vertices',PREV/'effective_potential.py')
    q=np.array([1.,0,0,0]);H=e.tree_hessian(q);m2=np.diag(H)
    assert np.max(abs(H-np.diag(m2)))==0
    g=np.einsum('ijkb,b->ijk',e.K,q)
    value,_,tad,_=e.cw(q,True);dv=-np.linalg.solve(H,tad)
    prior=json.loads((PREV/'potential_results.json').read_text())
    xpole=json.loads((HERE/'pole_results.json').read_text())
    rows=[];quaderrors=[];cuts=[]
    for i,name in enumerate(['phi','x','y','D']):
        quartic=.5*sum(e.K[i,i,k,k]*p.A(m2[k],.1) for k in range(4))
        displacement=float(np.einsum('a,a->',g[i,i,:],dv))
        bubbles=[];bubble0=0.;bubbleshell=0.
        for k in range(4):
            for l in range(4):
                if g[i,k,l]==0:continue
                threshold=(np.sqrt(m2[k])+np.sqrt(m2[l]))**2
                assert m2[i]<threshold
                cuts.append(float(threshold-m2[i]))
                b0=p.B0(0,m2[k],m2[l],.1);bs=p.B0(m2[i],m2[k],m2[l],.1)
                quaderrors.extend([abs(b0-p.B0_gauss(0,m2[k],m2[l],.1,128)),abs(bs-p.B0_gauss(m2[i],m2[k],m2[l],.1,128))])
                v0=-.5*g[i,k,l]**2*b0/p.LOOP;vs=-.5*g[i,k,l]**2*bs/p.LOOP
                bubble0+=v0;bubbleshell+=vs
                bubbles.append({'internal_indices':[k,l],'vertex_over_v':float(g[i,k,l]),'pair_cut_over_v2':float(threshold),'B0_zero':float(b0),'B0_shell':float(bs)})
        curv=m2[i]+displacement+quartic+bubble0
        mass2=m2[i]+displacement+quartic+bubbleshell
        rows.append({'field':name,'tree_mass_squared_over_v2':float(m2[i]),
                     'vacuum_displacement_over_v2':displacement,'quartic_tadpole_over_v2':float(quartic),
                     'bubble_zero_over_v2':float(bubble0),'bubble_shell_over_v2':float(bubbleshell),
                     'strict_curvature_squared_over_v2':float(curv),'strict_pole_squared_over_v2':float(mass2),
                     'strict_pole_mass_over_v':float(np.sqrt(m2[i])+(mass2-m2[i])/(2*np.sqrt(m2[i]))),
                     'sqrt_truncated_pole_squared_over_v':float(np.sqrt(mass2)),
                     'pole_minus_curvature_squared_over_v2':float(mass2-curv),'cubic_channels_ordered':bubbles})
    checks={
        'all_four_curvatures_reproduce_potential':max(abs(r['strict_curvature_squared_over_v2']-prior['vacuum']['strict_one_loop_curvature_diagonal'][i]) for i,r in enumerate(rows))<1e-14,
        'x_y_degeneracy':abs(rows[1]['strict_pole_squared_over_v2']-rows[2]['strict_pole_squared_over_v2'])<1e-14,
        'independent_X_pole_formula_agrees':abs(rows[1]['strict_pole_squared_over_v2']-xpole['matching_result']['strict_pole_squared_over_v2'])<1e-14,
        'all_nonzero_cubic_channels_below_pair_cuts':min(cuts)>0,
        'adaptive_and_Gauss_bubbles_agree':max(quaderrors)<1e-11,
        'all_pole_squared_masses_positive':min(r['strict_pole_squared_over_v2'] for r in rows)>0,
    }
    checks={k:bool(v) for k,v in checks.items()}
    res={'protocol_sha256':hashlib.sha256(protocol).hexdigest(),'input_tensor_source_sha256':hashlib.sha256((PREV/'effective_potential.py').read_bytes()).hexdigest(),
         'results':rows,'maximum_quadrature_difference':float(max(quaderrors)),'checks':checks,'all_checks_pass':all(checks.values()),'elapsed_seconds':time.perf_counter()-start,
         'scope':'One-loop new-scalar vacuum pole thresholds, not a quantum finite-Q binding calculation.'}
    (HERE/'neutral_pole_results.json').write_text(json.dumps(res,indent=2)+'\n')
    (HERE/'neutral_pole_execution.jsonl').write_text(json.dumps({'command':'python neutral_poles.py','checks':checks,'elapsed_seconds':res['elapsed_seconds']})+'\n')
    print(json.dumps({'masses':[{k:r[k] for k in ['field','strict_pole_squared_over_v2','strict_pole_mass_over_v']} for r in rows],'checks':checks},indent=2))
    if not res['all_checks_pass']:raise SystemExit(1)


if __name__=='__main__':main()
