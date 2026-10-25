"""Run the full C solver and serialize explicitly named diagnostics."""
import json
import subprocess
import time
from pathlib import Path
import numpy as np
from model import parameters, minimum, initial_z, yield_factor, YQ
from build import BUILD

ROOT=Path(__file__).resolve().parents[1]
NAMES=['NX','NY','NS','Neff','EX','EY','ES','V_kappa','V_mu','epsX','NcircX',
       'EY_exc_quadratic','JY_approx','gauge_Higgs_norm_over_H','vector_mass_over_H',
       'conformal_Hamiltonian','physical_a4rho','Gauss','Hamiltonian_work_error','abs_aY']


def evolve(H,cH,delta,end=1000.,rtol=2e-10,atol=2e-13,ns=120,z0=None):
    point=minimum(H,cH);p=parameters(H)
    z=initial_z(point,delta) if z0 is None else np.asarray(z0,complex)
    args=[str(BUILD/'release')]+[repr(float(p[k])) for k in ['M','h','mu']]
    args += [repr(float(end)),repr(float(rtol)),repr(float(atol))]
    args += [repr(float(v)) for zz in z for v in [zz.real,zz.imag]]+[str(ns)]
    t=time.perf_counter()
    proc=subprocess.run(args,capture_output=True,text=True)
    if proc.returncode:raise RuntimeError(proc.stderr+'\n'+proc.stdout[-1000:])
    hist=[];windows={};final=None;control=None
    for line in proc.stdout.splitlines():
        fields=line.split();a=np.array(list(map(float,fields[1:])))
        if fields[0]=='S':hist.append(dict(a=a[0],**dict(zip(NAMES,a[1:]))))
        elif fields[0]=='W':
            lo,hi,span=a[:3];v=a[3:23]
            d=dict(zip(NAMES,v));d.update(a_lo=lo,a_hi=hi,span=span)
            for j,key in enumerate(['NX','Neff','epsX']):d[key+'_std']=float(np.sqrt(max(0,a[23+j]-d[key]**2)))
            factor=yield_factor(H)
            d.update(YX=d['NX']*factor,Yeff=d['Neff']*factor,
                     YX_over_YQ=d['NX']*factor/YQ,Yeff_over_YQ=d['Neff']*factor/YQ,
                     eps_from_mean_charge=d['NX']/d['NcircX'],
                     approximate_Y_decay_pairs_per_abs_X_charge=(1.5*d['JY_approx']/abs(d['NX']) if d['NX'] else None))
            windows[f'{lo:g}-{hi:g}']=d
        elif fields[0]=='F':final=a.tolist()
        elif fields[0]=='END':control=dict(steps=int(a[0]),rejects=int(a[1]),max_Gauss_over_initial_Ncirc=a[2],max_Hamiltonian_work_error_over_initial_energy=a[3])
    return dict(H_GeV=H,cH=cH,delta=delta,a_end=end,rtol=rtol,atol=atol,
                initial_point=point,initial_complex_fields_over_H=[[x.real,x.imag] for x in z],
                yield_factor=yield_factor(H),windows=windows,history=hist,final_canonical_state=final,
                controls=control,runtime_seconds=time.perf_counter()-t)


def write(path,obj):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');tmp.replace(path)


if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser();a.add_argument('--H',type=float,required=True);a.add_argument('--cH',type=float,default=1)
    a.add_argument('--delta',type=float,default=np.pi/2);a.add_argument('--end',type=float,default=1000)
    a.add_argument('--out',required=True);x=a.parse_args()
    r=evolve(x.H,x.cH,x.delta,end=x.end);write(x.out,r)
    w=r['windows'].get('600-1000',list(r['windows'].values())[-1])
    print(json.dumps(dict(H=x.H,cH=x.cH,time=r['runtime_seconds'],eps=w['epsX'],Y=w['YX_over_YQ'],controls=r['controls'])))
