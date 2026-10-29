"""Conditional ballistic distance from 2 MeV to equality; no perturbation solver."""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.integrate import simpson

ROOT=Path(__file__).resolve().parents[1]

def run(tag):
    d=json.loads((ROOT/(tag+'.json')).read_text())
    z=np.load(ROOT/(tag+'.npz')); e=z['E']; n=len(e)
    i=np.argmin(abs(z['x']-np.log(10)))
    if abs(z['x'][i]-np.log(10))>1e-10:raise ValueError('Need a=10 spectrum')
    raw=z['Y'][:n,i]+z['Y'][n:2*n,i]
    discarded=float(np.maximum(-raw,0).sum());weight=np.maximum(raw,0)
    # Standard cosmology used solely as a reference background.
    # Planck 2018 base LCDM H0=67.4, Omega_m=.315, z_eq=3402.
    H0=67.4; om=.315; ar=1/3403.;orr=om*ar
    Tcmb=2.7255*8.617333262145e-14  # GeV
    a0=Tcmb/.02*(3.909/10.75)**(1/3)
    ast=10*a0;p=np.sqrt(e*e-1)
    x=np.linspace(np.log(ast),np.log(ar),4001);a=np.exp(x)
    vel=p[:,None]*ast/np.sqrt(a[None,:]**2+(p[:,None]*ast)**2)
    integral=simpson(vel*a[None,:]/np.sqrt(orr+om*a[None,:]),x=x,axis=1)
    path_pc=299792.458/H0*1e6*integral
    quant=[];cum=np.cumsum(weight)/weight.sum()
    for q in (.5,.9,.99):quant.append(float(np.interp(q,cum,path_pc)))
    return {'tag':tag,'a_relative_start':10.,'a_absolute_start':ast,'a_equality':ar,
            'Nfree':float(weight.sum()),'discarded_negative_solver_population':discarded,
            'free_mass_fraction_of_reference':float(weight.sum()/187479.9981052459),
            'mean_path_pc':float(path_pc@weight/weight.sum()),'median_p90_p99_path_pc':quant,
            'floor_particles':float(weight[0]),'p_start_mean_mX':float(p@weight/weight.sum())}

def main():
    p=argparse.ArgumentParser();p.add_argument('tags',nargs='+');args=p.parse_args()
    result={'results':[run(t) for t in args.tags],
            'parameters':{'H0_km_s_Mpc':67.4,'Omega_m':.315,'z_eq':3402,
                          'T_CMB_K':2.7255,'gstar_s_start':10.75,'gstar_s_today':3.909},
            'scope':'Ballistic path only after the selected a=10 snapshot through equality. No subsequent interactions, core-emission injection, self-gravity, full transfer function, perturbations, or measured Lyman-alpha exclusion. Late radiation-density scaling is used all the way back to 2 MeV, so the e+e- entropy transition is approximate and early H is mismatched by about 10 percent.'}
    (ROOT/'FREEFLIGHT.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
