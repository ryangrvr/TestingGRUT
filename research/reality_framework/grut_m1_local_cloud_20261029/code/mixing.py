"""Collisionless first-arrival ledger at conditional inter-core distances.

This is a two-boundary single-source transport control. Reaching a sphere is
not proof of spatial homogenization and any daughter entering the inner
matching sphere is handed off, not assumed absorbed or scattered.
"""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.integrate import cumulative_trapezoid,simpson
from scipy.stats import qmc
from local_rates import pulse_series,MH,GAMMA
from daughter_cloud import reconstruct_times,BETA_STAR,ESTAR

ROOT=Path(__file__).resolve().parents[1]
D=(1.2141071422639381e-15)**(-1/3)


def primary_arrivals(z,tgrid,E,plus,minus,R,times):
    de=float(np.diff(z['energy_edges'])[0]);r0=float(z['probe_radius_mX_inverse'])
    p=np.sqrt(E*E-1);v=p/E;out={}
    for label,flux in [('plus',plus*de),('minus',minus)]:
        cum=cumulative_trapezoid(flux,x=tgrid,axis=0,initial=0)
        total=cum[-1];rows=[]
        for T in times:
            cut=T-(R-r0)/v
            # Each energy has its own retarded emission-time limit.
            counts=np.array([np.interp(cut[i],tgrid,cum[:,i],left=0,right=total[i]) for i in range(len(E))])
            rows.append({'arrived_number':float(counts.sum()),'arrived_energy_mX':float(counts@E),
                         'unarrived_number':float((total-counts).sum()),'unarrived_energy_mX':float((total-counts)@E)})
        out[label]={'initial_number':float(total.sum()),'initial_energy_mX':float(total@E),'snapshots':rows}
    return out


def samples(z,power,seed):
    h=z['neutral_positive'];weights=h.ravel();N=weights.sum();cum=np.cumsum(weights)/N
    U=qmc.Sobol(d=5,scramble=True,seed=seed).random_base2(power);n=len(U)
    idx=np.minimum(np.searchsorted(cum,U[:,0]),len(cum)-1);w=idx//h.shape[1];j=idx%h.shape[1]
    width=float(z['window_width_mX_t']);centers=z['window_center_mX_t'][w]
    temit=reconstruct_times(centers,width,U[:,1],float(z['window_center_mX_t'][-1]))
    edges=z['energy_edges'];lo=np.maximum(MH,edges[j]);hi=edges[j+1]
    if np.any(hi<=lo):raise ValueError('unphysical neutral sample')
    Eh=lo+(hi-lo)*U[:,2];gamma=Eh/MH;beta=np.sqrt(np.maximum(0,1-1/gamma**2))
    delay=-np.log(np.maximum(1-U[:,3],1e-300))*gamma/GAMMA
    r0=float(z['probe_radius_mX_inverse']);td=temit+delay;rd=r0+beta*delay;x=2*U[:,4]-1
    daughters=[]
    for xx in (x,-x):
        energy=gamma*ESTAR*(1+beta*BETA_STAR*xx);pr=gamma*ESTAR*(beta+BETA_STAR*xx)
        p=np.sqrt(np.maximum(energy*energy-1,0));v=p/energy;mu=np.divide(pr,p,out=np.zeros_like(p),where=p>0)
        daughters.append((energy,v,mu))
    return {'N':float(N),'n':n,'r0':r0,'temit':temit,'td':td,'rd':rd,'beta':beta,'Eh':Eh,'daughters':daughters}


def parent_flux(sample,R,times):
    Nh=sample['N'];n=sample['n'];r0=sample['r0'];td=sample['td'];rd=sample['rd'];Eh=sample['Eh'];temit=sample['temit'];beta=sample['beta']
    th=temit+np.divide(R-r0,beta,out=np.full_like(beta,np.inf),where=beta>0)
    daughter_times=[]
    for energy,v,mu in sample['daughters']:
        b=rd*np.sqrt(np.maximum(0,1-mu*mu));z0=rd*mu
        out=td+np.divide(np.sqrt(np.maximum(0,R*R-b*b))-z0,v,out=np.full_like(v,np.inf),where=v>0)
        canhit=(mu<0)&(b<r0)
        path=-np.sqrt(np.maximum(0,r0*r0-b*b))-z0
        hit=np.where(canhit,td+np.divide(path,v,out=np.full_like(v,np.inf),where=v>0),np.inf)
        daughter_times.append((energy,out,hit))
    rows=[]
    for T in times:
        h_out=(th<=T)&(th<td);h_inside=(td>T)&(th>T);decay_inside=(td<=T)&(td<th)
        energy_total=Eh*(h_out+h_inside);categories=[];charge_total=np.zeros(n)
        for sign,(energy,tout,thit) in zip((1,-1),daughter_times):
            escaped=decay_inside&(tout<=T)&(tout<thit)
            returned=decay_inside&(thit<=T)&(thit<tout)
            inside=decay_inside&~escaped&~returned
            energy_total+=energy*(escaped+returned+inside)
            charge_total+=sign*(escaped.astype(int)+returned.astype(int)+inside.astype(int))
            categories.append({'out_number':float(Nh*escaped.mean()),'out_energy_mX':float(Nh*np.mean(energy*escaped)),
                               'inner_hit_number':float(Nh*returned.mean()),'inner_hit_energy_mX':float(Nh*np.mean(energy*returned)),
                               'inside_number':float(Nh*inside.mean()),'inside_energy_mX':float(Nh*np.mean(energy*inside))})
        rows.append({'h_out_number':float(Nh*h_out.mean()),'h_out_energy_mX':float(Nh*np.mean(Eh*h_out)),
                     'h_inside_number':float(Nh*h_inside.mean()),'h_inside_energy_mX':float(Nh*np.mean(Eh*h_inside)),
                     'positive_daughters':categories[0],'negative_daughters':categories[1],
                     'per_parent_energy_closure_max_mX':float(np.max(abs(energy_total-Eh))),
                     'per_parent_charge_closure_max':int(np.max(abs(charge_total)))})
    return {'initial_h_number':Nh,'initial_h_energy_mX':float(Nh*Eh.mean()),'snapshots':rows}


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',default='source_t400');p.add_argument('--power',type=int,default=18)
    p.add_argument('--seed',type=int,default=19);p.add_argument('--dt',type=float,default=10.);p.add_argument('--tag',required=True);args=p.parse_args()
    z=np.load(ROOT/(args.source+'.npz'));tgrid,e,plus,minus=pulse_series(z,args.dt)
    T=np.array([10000.,20000.,50000.,100000.,200000.,500000.,1000000.]);s=samples(z,args.power,args.seed);rows=[]
    for fraction in (.5,1.):
        R=D*fraction
        rows.append({'radius_over_density_length':fraction,'radius_mX_inverse':R,
                     'primary':primary_arrivals(z,tgrid,e,plus,minus,R,T),
                     'neutral_parents_and_daughters':parent_flux(s,R,T)})
    out={'configuration':vars(args),'density_length_mX_inverse':D,'snapshot_times_mX_inverse':T.tolist(),
         'clock_optical_time_at_1e6':1.2141071422639381e-9,
         'clock_hubble_fraction_at_1e6':1.930099844452384e-15,'regions':rows,
         'scope':'First passage at two illustrative radii of a collisionless, single-source plume. No claim of actual spatial homogenization, core interaction, daughter interactions, overlapping plumes or formation history.'}
    (ROOT/(args.tag+'.json')).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
