"""Direct scattering holdouts at populated energy nodes and evolving charges."""
from pathlib import Path
import argparse
import json
import time
import numpy as np
from scipy.interpolate import PchipInterpolator
import expansion
import core_family
import core

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser();p.add_argument('--tau',type=float,required=True);p.add_argument('--count',type=int,default=60)
    p.add_argument('--tag',default='expanding96');args=p.parse_args();start=time.monotonic()
    d=json.loads((ROOT/(args.tag+'.json')).read_text())
    row=min(d['snapshots'],key=lambda r:abs(r['tau']-args.tau));q=row['Qcore']
    z=np.load(ROOT/(args.tag+'.npz'));idx=np.argmin(abs(z['x']-np.log(row['a_over_a0'])))
    e=z['E'];y=z['Y'][:,idx];n=len(e);sol=core_family.use_q(q);omega=float(sol.p[0])
    qnodes=np.arange(186000.,204001.,3000.);Qtables=[];ks=None
    for Q in qnodes:
        r=np.load(ROOT/f'core_Q{int(Q)}.npz');ks=r['k'];Qtables.append(r['W'])
    j=min(np.searchsorted(qnodes,q)-1,len(qnodes)-2);j=max(j,0)
    f=(q-qnodes[j])/(qnodes[j+1]-qnodes[j]);W=(1-f)*Qtables[j]+f*Qtables[j+1]
    interp=PchipInterpolator(ks,W,axis=0);old=np.load(ROOT.parent/'grut_m1_kinetics_20261028'/'CORE_KERNELS.npz')
    reference=PchipInterpolator(old['k'],old['W'],axis=0)
    def get(fun,k,k0):return fun(k0)*k/k0 if k<k0 else fun(k)
    candidates=[]
    for species in (0,1):
        for i,E in enumerate(e):
            N=float(y[species*n+i]);p_i=np.sqrt(E*E-1)
            ea=E if species==1 else E-2*omega
            if ea<=1 or N<=0:continue
            k=np.sqrt(ea*ea-1);eh=E+(omega if species==1 else -omega)
            if k>ks[-1]:continue
            w=get(interp,k,ks[0]);factor=N*np.pi/(p_i*E)
            channel=np.array([w[0,1], w[2,species] if eh>np.sqrt(6) else 0])*factor
            candidates.append({'species':species,'i':i,'E':float(E),'k':float(k),'eh':float(eh),'factor':float(factor),
                               'approx':channel,'rank':float(channel.sum())})
    candidates.sort(key=lambda a:a['rank'],reverse=True)
    selected=candidates[:args.count]
    # Always include populated bins close to the moving ell=4 resonance and
    # neutral threshold, even if interpolation itself underestimated them.
    used={(a['species'],a['i']) for a in selected}
    for r in candidates[args.count:]:
        if (0.31<r['k']<.41 or abs(r['eh']-np.sqrt(6))<.02) and (r['species'],r['i']) not in used:
            selected.append(r);used.add((r['species'],r['i']))
    sums={'direct':np.zeros((2,2)),'interpolated':np.zeros((2,2)),'frozen_reference':np.zeros((2,2))}
    details=[]
    for r in selected:
        w,u,rec,l=core.strengths(r['k'])
        oldw=get(reference,r['k'],old['k'][0])
        direct=np.array([w[0,1],w[2,r['species']] if r['eh']>np.sqrt(6) else 0])*r['factor']
        base=np.array([oldw[0,1],oldw[2,r['species']] if r['eh']>np.sqrt(6) else 0])*r['factor']
        for name,value in [('direct',direct),('interpolated',r['approx']),('frozen_reference',base)]:sums[name][r['species']]+=value
        details.append({'species':r['species'],'E':r['E'],'k':r['k'],'direct_rates':direct.tolist(),
                        'interpolated_rates':r['approx'].tolist(),'reference_rates':base.tolist(),
                        'unitarity_error':u,'reciprocity_error':rec,'ell_max':l})
    all_approx=sum(r['rank'] for r in candidates);covered=sum(r['rank'] for r in selected)
    out={'q':q,'omega':omega,'tau':row['tau'],'tag':args.tag,'count_selected':len(selected),
         'count_candidates':len(candidates),'coverage_by_interpolated_event_rate':covered/all_approx if all_approx else 1,
         'sum_rows':'species positive then negative; columns exchange then neutral',
         'rate_sums':{k:v.tolist() for k,v in sums.items()},'samples':details,'seconds':time.monotonic()-start,
         'scope':'Direct holdout only on selected populated nodes from a frozen-profile trajectory. Coverage is ranked with the interpolated model, so it is not a bound on unsampled resonances.'}
    name=f'KERNEL_HOLDOUT_{int(args.tau)}.json'
    (ROOT/name).write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(name,'Q',q,'selected',len(selected),'coverage',out['coverage_by_interpolated_event_rate'],
          'direct',sums['direct'],'interpolated',sums['interpolated'],'reference',sums['frozen_reference'],flush=True)

if __name__=='__main__':main()
