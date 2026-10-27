"""Opposite-current candidate masks; no identification with antiparticles."""
from analyze import *

def negative_catalog(q,p,L,threshold=4.):
    N=q.shape[-1];dx=L/N;j=q[0]*p[1]-q[1]*p[0];reference=abs(float(j.mean()))
    sm=gaussian_filter(-j,1/dx,mode='wrap');ph=gaussian_filter(q[2]**2,1/dx,mode='wrap')
    mask=(sm>threshold*reference)&(ph<.25*C.NU**2);labs,n=periodic_labels(mask);rho=local_energy(q,p,L);rows=[]
    for ident in range(1,n+1):
        m=labs==ident;Q=float(-j[m].sum()*dx**3)
        if m.sum()<4 or Q<100:continue
        rows.append(dict(cells=int(m.sum()),opposite_Q_mask_magnitude=Q,E_mask_m=float(rho[m].sum()*dx**3)))
    rows.sort(key=lambda x:x['opposite_Q_mask_magnitude'],reverse=True)
    return dict(threshold=threshold,count=len(rows),candidates=rows,
        local_negative_current_magnitude_over_net= float(np.maximum(-j,0).sum()/abs(j.sum())))

def main():
    runs=json.loads((ROOT/'RUNS.json').read_text())['cases'];out={}
    for tag,r in runs.items():
        paths=r['snapshots'] if tag in ['mixed_N96','onesign_N96'] else r['snapshots'][-1:];rows=[]
        for path in paths:
            N,L,t,a=snapshot(path);q=np.asarray(a[:3]);p=np.asarray(a[3:]);rows.append(dict(t=t,catalogs=[negative_catalog(q,p,L,k) for k in [2.,4.,8.]]))
        out[tag]=rows
    save('OPPOSITE_CURRENT.json',dict(cases=out,physical_antiparticle_interpretation=False,
        scope='Same phi-depletion/current thresholds with current sign reversed and normalized to |mean net current|. Final states for all7 runs and histories for the two main96³ controls. Counts are opposite-current candidates, not anti-Q-balls or an antiparticle abundance.'))
    print([(tag,[(c['threshold'],c['count'],[x['opposite_Q_mask_magnitude'] for x in c['candidates'][:3]]) for c in rows[-1]['catalogs']]) for tag,rows in out.items()],flush=True)

if __name__=='__main__':main()
