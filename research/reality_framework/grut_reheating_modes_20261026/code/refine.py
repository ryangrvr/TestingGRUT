"""Adaptive-follow-up band refinement; no background/phase parameter change."""
import json
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
from mode_common import ROOT,WORK,write
from modes import build,integrate,diagnostics

def main():
    source=json.loads((ROOT/'MODE_RESULTS.json').read_text());old=np.load(ROOT/'MODE_PROFILES.npz');times=old['times']
    fineks=np.geomspace(.8,6,65);exe=build();jobs=[]
    for r in source['rows']:
        for n in range(0,len(fineks),8):jobs.append((r,n,fineks[n:n+8]))
    cache={r['tag']:[] for r in source['rows']}
    with ThreadPoolExecutor(max_workers=3) as pool:
        tasks={pool.submit(integrate,exe,f"ref_{r['tag']}_{n}",r['beta'],r['Gamma_over_H'],kk,times):(r,n,kk) for r,n,kk in jobs}
        for future in as_completed(tasks):
            r,n,kk=tasks[future];bg,TT,control=future.result();dd,checks=diagnostics(r['beta'],bg,TT,kk)
            mean=abs(r['late_NX']);dd[:,:,1]/=mean;dd[:,:,3]/=mean**2
            cache[r['tag']].append((n,dd));print('REFINE',r['tag'],n,flush=True)
    rows=[];profiles=dict(times=times,ks=fineks)
    for r in source['rows']:
        tag=r['tag'];d=np.concatenate([q[1] for q in sorted(cache[tag],key=lambda q:q[0])],axis=1)
        profiles[tag]=d;intervals=[]
        for end in (3.,10.,40.,80.):
            mask=times<=end;dd=d[mask];peak=np.unravel_index(dd[:,:,0].argmax(),dd[:,:,0].shape)
            intervals.append(dict(end_tau=end,max_amplitude_response=float(dd[:,:,0].max()),max_charge_response=float(dd[:,:,1].max()),
                peak_k=float(fineks[peak[1]]),peak_tau=float(times[mask][peak[0]])))
        rows.append(dict(tag=tag,intervals=intervals))
    np.savez_compressed(ROOT/'REFINED_PROFILES.npz',**profiles)
    write('REFINED_RESULTS.json',dict(rows=rows,ks=fineks,reason='Resolved amplification band near k/H_i~2 on initial33-node main grid; fixed65-node log refinement from.8 to6, same backgrounds and seed definitions.'))
    print('REFINEMENT DONE',flush=True)

if __name__=='__main__':main()
