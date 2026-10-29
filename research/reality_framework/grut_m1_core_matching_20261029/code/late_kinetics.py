"""Same conditional expanding-core model with the extra late source included."""
from pathlib import Path
import argparse
import json
import sys
import time
import numpy as np
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parent/'grut_m1_expansion_clock_20261028'
sys.path.insert(0,str(PARENT/'code'))
import expansion
from hybrid import Hybrid

def append_late(e,source,tag):
    z=np.load(ROOT/f'late_{tag}_spectrum_10000_20000.npz');n=len(e)
    plus=source[:n].copy();minus=source[n:].copy();audit={}
    for species in ('plus','minus','neutral'):
        w,en=z[species+'_N_signed'],z[species+'_EN_signed'];positive=w>0
        energy=np.divide(en,w,out=np.zeros_like(w),where=positive)
        low=1 if species!='neutral' else np.sqrt(6)
        valid=positive&(energy>=low)&np.isfinite(energy);keep=valid&(energy>=e[0])&(energy<=e[-1])
        item={'input_signed_number':float(w.sum()),'negative_signed_number':float(w[w<0].sum()),
              'invalid_positive_number':float(w[positive&~valid].sum()),
              'outside_grid_positive_number':float(w[valid&~keep].sum()),
              'internal_negative_weight':float(z[species+'_N_negative_part'].sum())}
        if species=='neutral':
            missing=0.
            for ev,weight in zip(energy[keep],w[keep]):
                dd,f=expansion.daughters(e,ev);plus+=weight*f*dd;minus+=weight*f*dd;missing+=weight*(1-f)
            item['neutral_decay_unsupported_number']=float(missing)
        elif species=='plus':plus+=expansion.deposit(e,energy[keep],w[keep])
        else:minus+=expansion.deposit(e,energy[keep],w[keep])
        audit[species]=item
    return np.r_[plus,minus],audit

def run(tag,bins,amax):
    start=time.monotonic();e=expansion.grid(0,8.,'translated',bins)
    baseline,baseaudit=expansion.radial_source(e,tag);source,lateaudit=append_late(e,baseline,tag)
    eps=expansion.clock()['H0_mX2_over_ncore0'];model=Hybrid(e,'full',8,eps,.00025,True)
    q0,omega0,*_=model.terms(source);pulse=json.loads((PARENT/'sources'/f'ringdown_{tag}.json').read_text())['initial_pulse_energy_mX']
    source_energy=float(e@(source[:len(e)]+source[len(e):]))
    remainder=float(pulse-model.mass(q0)-source_energy)
    grid=np.unique(np.r_[0.,-np.log(1-np.array([.1,1,10,100,1000,10000,100000])*eps),np.log([1.01,1.1,2.,5.,amax])])
    grid=grid[(grid>=0)&(grid<=np.log(amax)+1e-14)]
    y0=np.r_[source,0.,0.,0.]
    sol=solve_ivp(model.rhs,(0,np.log(amax)),y0,method='BDF',jac=model.jac,rtol=2e-7,atol=1e-9,t_eval=grid)
    if not sol.success:raise RuntimeError(sol.message)
    rows=[model.observables_hybrid(x,y) for x,y in zip(sol.t,sol.y.T)]
    ledger=np.array([r['total_stationary_plus_gas_plus_work_mX']-r['interpolated_generator_energy_defect_mX'] for r in rows])
    final=rows[-1];full_percent=100*((float(model.mass(final['Qcore']))+final['gas_energy_mX']+remainder)/expansion.EREF-1)
    out={'source_tag':tag,'kinetic_nodes':len(e),'shift_bins':bins,'amax':amax,
         'initial_pulse_energy_mX':pulse,'initial_Q_core':float(q0),'initial_Omega':float(omega0),
         'initial_gas_energy_mX':source_energy,'held_unassigned_recomputed_mX':remainder,
         'initial_plus':float(source[:len(e)].sum()),'initial_minus':float(source[len(e):].sum()),
         'late_source_audit':lateaudit,'baseline_source_audit':baseaudit,
         'snapshots':rows,'full_modeled_final_excess_percent_with_remainder':full_percent,
         'checks':{'max_energy_plus_work_minus_generator_defect_drift_mX':float(np.max(abs(ledger-ledger[0]))),
                   'minimum_saved_population':float(sol.y[:2*len(e)].min()),'max_generator_defect':model.build_error},
         'solver':{'seconds':time.monotonic()-start,'nfev':sol.nfev,'njev':sol.njev,'nlu':sol.nlu},
         'scope':'Conditional parent homogeneous instantaneous-source model, now including the observed 10000–20000 radial segment and recomputing the held energy residual once. Not a self-consistent interacting-core ringdown, proof of mixing, all-future source, or cosmological 1% verdict.'}
    return out,sol,e

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',choices=['h005','h0025'],required=True)
    p.add_argument('--bins',type=int,default=64);p.add_argument('--amax',type=float,default=10.);p.add_argument('--tag',required=True);a=p.parse_args()
    out,sol,e=run(a.source,a.bins,a.amax)
    np.savez_compressed(ROOT/(a.tag+'.npz'),E=e,x=sol.t,Y=sol.y)
    (ROOT/(a.tag+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
