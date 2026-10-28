"""Fixed production and numerical-control histories. Cache never changes a run."""
import argparse
from common import *

CASES={
 'cold64':dict(N=64,L=48.,dt=.04,T=300.),
 'cold48':dict(N=48,L=48.,dt=.04,T=300.),
 'cold64_half':dict(N=64,L=48.,dt=.02,T=300.),
 'cold64_seed2':dict(N=64,L=48.,dt=.04,T=300.,number=2),
 'ellipse64':dict(N=64,L=48.,dt=.04,T=300.,ellipse=True),
 'object64':dict(N=64,L=64.,dt=.04,T=200.,kind='object'),
 'object64_perturbed':dict(N=64,L=64.,dt=.04,T=200.,kind='object',perturbation=.05),
 'object96_perturbed':dict(N=96,L=64.,dt=.04,T=200.,kind='object',perturbation=.05),
 'dilute64':dict(N=64,L=64.,dt=.04,T=600.,F=.5),
 'dilute48':dict(N=48,L=64.,dt=.04,T=600.,F=.5),
 'dilute64_seed2':dict(N=64,L=64.,dt=.04,T=600.,F=.5,number=2),
 'ellipse96':dict(N=96,L=48.,dt=.04,T=300.,ellipse=True,threads=4),
 'ellipse64_half':dict(N=64,L=48.,dt=.02,T=300.,ellipse=True),
}

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('cases',nargs='*');args=a.parse_args()
    for tag in args.cases or list(CASES):
        r=run(tag,**CASES[tag])
        print(json.dumps({'tag':tag,'seconds':r['seconds'],'energy_drift':r['energy_relative_drift'],
                          'charge_drift':r['charge_relative_drift'],'initial_E_over_Q':r['parameters']['initial_E_over_Q'],
                          'final':r['final']}),flush=True)
