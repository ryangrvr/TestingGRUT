"""Audit the existing cold fields; do not infer masses or particle occupations."""
from cw_common import *

def resolve_path(tag,path):
    return P.WORK/tag/Path(path).name

def audit(path,catalog,include_spheres):
    N,L,t,arr=P.snapshot(path);q=np.asarray(arr[:3],float);p=np.asarray(arr[3:],float)
    f=ledger(q,p,L)
    # Reuse exactly the prior mask definition, with non-mutating analysis calls.
    computed,labs=A.catalog(q,p,L,4.)
    assert [x['id'] for x in computed['candidates']]==[x['id'] for x in catalog['candidates']]
    smooth=gaussian_filter(q[2],1/(L/N),mode='wrap')
    whole=region(q,p,f,np.ones((N,N,N),bool),L)
    records=[];compact_count=0
    for c in catalog['candidates']:
        entry=dict(candidate=c,mask=region(q,p,f,labs==c['id'],L))
        if include_spheres and c['compact'] and compact_count<3:
            compact_count+=1;rr=distances(N,L,c['centre']);spheres=[]
            for r in [2.,4.,6.,8.,10.]:
                spheres.append(dict(radius=r,**region(q,p,f,rr<r,L)))
            entry.update(spheres=spheres,neutral_sign_shell=sign_shell(smooth,rr,c['r90_mask_over_m']))
        records.append(entry)
    prior=float(A.local_energy(q,p,L).mean())
    return dict(t=t,N=N,L=L,candidates=records,whole_box=whole,
        global_energy_reconstruction_relative=abs(whole['E_total']/L**3/prior-1),
        global_charge_reconstruction_relative=abs(whole['Q']/L**3/f['j'].mean()-1),
        robust_neutral_sign_volume_fractions=dict(positive=float(np.mean(smooth>.5*C.NU)),negative=float(np.mean(smooth<-.5*C.NU)),weak=float(np.mean(abs(smooth)<=.5*C.NU))))

def main():
    runs=json.loads((PILOT/'RUNS.json').read_text())['cases']
    analysis=json.loads((PILOT/'ANALYSIS.json').read_text())['cases']
    final={};history={}
    for tag,run in runs.items():
        row=analysis[tag]['rows'][-1]
        final[tag]=audit(resolve_path(tag,run['snapshots'][-1]),row['catalogs'][1],True)
        save('FINAL_AUDIT.json',dict(cases=final,physical_bound_population=None,physical_f_sym=None))
        print('FINAL',tag,'candidates',len(final[tag]['candidates']),flush=True)
    for tag in ['mixed_N96','onesign_N96']:
        rows=[]
        for path,prior in zip(runs[tag]['snapshots'],analysis[tag]['rows']):
            if prior['t']<40:continue
            if prior['t']==final[tag]['t']:row=final[tag]
            else:row=audit(resolve_path(tag,path),prior['catalogs'][1],False)
            rows.append(row)
            print('HISTORY',tag,prior['t'],flush=True)
        history[tag]=rows
        save('HISTORY_AUDIT.json',dict(cases=history,physical_bound_population=None,physical_f_sym=None))
    allrows=list(final.values())+[r for rs in history.values() for r in rs]
    identity=max(max(c['mask']['kinetic_identity_relative'],c['mask']['rotation_fit_identity_relative']) for row in allrows for c in row['candidates'])
    energy=max(row['global_energy_reconstruction_relative'] for row in allrows)
    charge=max(row['global_charge_reconstruction_relative'] for row in allrows)
    v=read('VERIFICATION_RESULTS.json')
    v['groups']['production_reconstruction']=dict(passed=max(identity,energy,charge)<1e-11,
        kinetic_identity_relative=identity,energy_relative=energy,charge_relative=charge,
        unique_snapshots=len(final)+sum(len(x)-1 for x in history.values()))
    save('VERIFICATION_RESULTS.json',v)
    assert v['groups']['production_reconstruction']['passed']

if __name__=='__main__':main()
