"""Overlap-defined persistence of candidate patterns; not particle worldlines."""
from analyze import *

def masks_for(path):
    N,L,t,a=snapshot(path);j=a[0]*a[4]-a[1]*a[3];dx=L/N
    smoothj=gaussian_filter(j.astype(float),1/dx,mode='wrap')
    smoothphi=gaussian_filter(np.asarray(a[2],float)**2,1/dx,mode='wrap')
    labs,_=periodic_labels((smoothj>4*j.mean())&(smoothphi<.25*C.NU**2))
    return labs

def main():
    runs=json.loads((ROOT/'RUNS.json').read_text())['cases'];analysis=json.loads((ROOT/'ANALYSIS.json').read_text())['cases'];out={}
    for tag in ['mixed_N96','onesign_N96','mixed_N128','mixed_dt_half']:
        if tag not in analysis:continue
        rows=analysis[tag]['rows'];paths=runs[tag]['snapshots'];tracks=[];active={};prev=None;previous_catalog={};nextid=0
        for path,row in zip(paths,rows):
            current={c['id']:c for c in row['catalogs'][1]['candidates']};labs=masks_for(path);matches={}
            if prev is not None and current and previous_catalog:
                ok=(prev>0)&(labs>0);pairs,counts=np.unique(np.stack([prev[ok],labs[ok]],axis=1),axis=0,return_counts=True)
                proposals=[]
                for (old,new),overlap in zip(pairs,counts):
                    if old not in previous_catalog or new not in current:continue
                    a=previous_catalog[old];b=current[new];dice=2*overlap/(a['cells']+b['cells']);ratio=b['Q_mask']/a['Q_mask']
                    if dice>.2 and .5<ratio<2:proposals.append((dice,int(old),int(new)))
                used=set()
                for dice,old,new in sorted(proposals,reverse=True):
                    if old not in active or old in used or new in matches:continue
                    matches[new]=active[old];used.add(old)
            newactive={}
            for ident,c in current.items():
                if ident in matches:tid=matches[ident]
                else:tid=nextid;nextid+=1;tracks.append(dict(id=tid,rows=[]))
                tracks[tid]['rows'].append(dict(t=row['t'],**c));newactive[ident]=tid
            active=newactive;previous_catalog=current;prev=labs
        persistent=[]
        for tr in tracks:
            rr=tr['rows'];duration=rr[-1]['t']-rr[0]['t']
            if duration>=59.9 and len(rr)>=4:
                persistent.append(dict(id=tr['id'],duration=duration,snapshots=len(rr),
                    Q_min=min(r['Q_mask'] for r in rr),Q_max=max(r['Q_mask'] for r in rr),
                    E_over_Q_min=min(r['E_mask_over_Qm'] for r in rr),E_over_Q_max=max(r['E_mask_over_Qm'] for r in rr),
                    omega_min=min(r['omega_from_current_over_m'] for r in rr),omega_max=max(r['omega_from_current_over_m'] for r in rr)))
        out[tag]=dict(tracks=tracks,persistent_summaries=persistent,
            definition='At least4 snapshots spanning>=59.9 m^-1, successive mask Dice>0.2 and charge ratio within(0.5,2). Greedy overlap matches patterns, not identities or stability eigenmodes.')
        print('TRACKED',tag,len(persistent),flush=True)
        save('PERSISTENCE.json',dict(cases=out,bound_object_certificate=False))

if __name__=='__main__':main()
