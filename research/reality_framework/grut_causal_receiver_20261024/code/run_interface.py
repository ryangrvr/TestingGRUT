#!/usr/bin/env python3
import hashlib,json,subprocess,tempfile,time
from pathlib import Path
import numpy as np
from instrument import build
from run import ROOT,BASE,initial


def main():
    name='I025';c={**BASE,'dt':.015625};data=ROOT/'data';data.mkdir(exist_ok=True)
    source=build();prefix=data/name
    with tempfile.TemporaryDirectory(prefix='grut-cfs-interface-') as tmp:
        tmp=Path(tmp);binary=tmp/'radial';init=tmp/'initial.csv'
        subprocess.run(['g++','-O3','-std=c++17',str(source),'-o',str(binary)],check=True)
        preparation,_=initial(c,init)
        command=[str(binary),str(init),str(prefix)]+[str(c[k]) for k in ['h','dt','length','tmax','probe','absorber','sigma','alpha']]+['1',str(c['mass'])]
        print('START',name,c,flush=True);start=time.monotonic()
        with (data/f'{name}.log').open('w') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True)
    arrays={key:np.loadtxt(str(prefix)+f'_{key}.csv',delimiter=',',skiprows=1) for key in ['history','trace','final','memory']}
    hashes={}
    for key,array in arrays.items():
        hashes[key]=hashlib.sha256(array.tobytes()).hexdigest();np.savez_compressed(str(prefix)+f'_{key}.npz',**{key:array});Path(str(prefix)+f'_{key}.csv').unlink()
    history=arrays['history'];e0,q0=history[0,1:3];ei0,qi0=history[0,8:10]
    keep=np.unique(np.r_[np.arange(0,len(history),4),len(history)-1])
    errors=dict(energy=float(max(abs((history[:,8]+history[:,10])/ei0-1))),charge=float(max(abs((history[:,9]+history[:,11])/qi0-1))))
    baseline=json.loads((ROOT/'P025_time_RESULTS.json').read_text())
    identities={key:hashes[key]==baseline['raw_numeric_sha256'][key] for key in ['trace','final','memory']}
    late=history[:,0]>.8*c['tmax']
    out=dict(name=name,parameters=c,seconds=time.monotonic()-start,preparation=preparation,
        source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),base_source_sha256=hashlib.sha256((ROOT/'code/radial_cfs.cpp').read_bytes()).hexdigest(),
        raw_numeric_sha256=hashes,field_dynamics_bitwise_identical_to_P025_time=identities,
        exact_cut_radius=(round(c['probe']/c['h'])+.5)*c['h'],
        maximum_cut_budget_relative_error=errors,Ebar_initial=float(e0),Qbar_initial=float(q0),
        maximum_energy_plus_receiver_relative_error=float(max(abs((history[:,1]+history[:,3])/e0-1))),
        maximum_charge_plus_receiver_relative_error=float(max(abs((history[:,2]+history[:,4])/q0-1))),
        core60_late_mean_fraction=float(np.mean(history[late,5])/q0),core60_final_fraction=float(history[-1,5]/q0),
        maximum_layer_charged_mass2_shift=float(max(history[:,6])),
        final_net_Eflux=float(history[-1,10]),final_net_Qflux=float(history[-1,11]),
        history_columns=['tau','E','Q','E_receiver','Q_receiver','Q60','max_layer_mass2_shift','Fmax','E_inside_cut','Q_inside_cut','E_flux_cut','Q_flux_cut'],
        history=history[keep].tolist(),
        scope='Exact classical signed discrete energy/charge interface. This is not yet a positive asymptotic time/momentum particle spectrum.')
    (ROOT/f'{name}_RESULTS.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print('DONE',name,'cut errors',errors,'bitwise dynamics',identities,flush=True)


if __name__=='__main__':main()
