"""Periodic catalogs and positive Hamiltonian diagnostics; no physical relic count."""
from common import *
from scipy.ndimage import gaussian_filter,label,generate_binary_structure
from scipy.optimize import linear_sum_assignment
from analysis_io import read_analysis,write_analysis

def positive_energy(q,p,L):
    """Same global FD4 Hamiltonian as the engine, locally nonnegative gradient.

    Per axis: .5/dx²*[.5*(forward²+backward²)+second_difference²/12].
    Its integral equals the engine's signed-link Hamiltonian exactly on the torus.
    Regional assignments differ by boundary terms; shell diagnostics remain necessary.
    """
    dx=L/q.shape[-1];ge=np.zeros(q.shape[-3:])
    for f in q:
        for ax in range(3):
            forward=np.roll(f,-1,ax)-f;backward=f-np.roll(f,1,ax)
            ge+=.5/dx**2*(.5*(forward**2+backward**2)+(forward-backward)**2/12)
    return .5*np.sum(p*p,axis=0)+potential(q)+ge

def labels(mask):
    lab,n=label(mask,structure=generate_binary_structure(3,1))
    if not n:return lab,0
    parent=np.arange(n+1)
    def find(a):
        while parent[a]!=a:parent[a]=parent[parent[a]];a=parent[a]
        return a
    for ax in range(3):
        a=np.take(lab,0,axis=ax);b=np.take(lab,-1,axis=ax);m=(a>0)&(b>0)
        for x,y in np.unique(np.stack([a[m],b[m]],axis=1),axis=0):
            x=find(x);y=find(y)
            if x!=y:parent[max(x,y)]=min(x,y)
    roots=np.array([find(i) for i in range(n+1)]);out=np.zeros(n+1,dtype=int)
    for i,r in enumerate(np.unique(roots[1:]),1):out[roots==r]=i
    return out[lab],int(out.max())

def sphere_diagnostics(q,p,L,centre,radii=(4.,6.,8.,10.,12.)):
    N=q.shape[-1];dx=L/N;dv=dx**3;axis=np.arange(N)*dx
    delta=[((axis-centre[i]+L/2)%L)-L/2 for i in range(3)]
    rr=np.sqrt(delta[0][:,None,None]**2+delta[1][None,:,None]**2+delta[2][None,None,:]**2)
    j=q[0]*p[1]-q[1]*p[0];rho=positive_energy(q,p,L);r2=q[0]**2+q[1]**2
    rows=[]
    for r in radii:
        m=rr<=r;sh=(rr>r)&(rr<=r+2)
        Q=float(j[m].sum()*dv);E=float(rho[m].sum()*dv)
        rows.append({'radius':r,'Q':Q,'E':E,'E_over_Q':E/Q if Q>0 else None,
                     'shell_X_rms_over_v':float(np.sqrt(r2[sh].mean())/V),
                     'shell_phi_deviation_rms_over_v':float(np.sqrt(((q[2][sh]-V)**2).mean())/V),
                     'shell_mean_charge_over_global':float(j[sh].mean()/j.mean()),
                     'shell_energy_per_volume':float(rho[sh].mean())})
    return rows

def catalog(q,p,L,threshold=4.,smooth=1.):
    N=q.shape[-1];dx=L/N;dv=dx**3;j=q[0]*p[1]-q[1]*p[0]
    sj=gaussian_filter(j,smooth/dx,mode='wrap');sp=gaussian_filter(q[2]**2,smooth/dx,mode='wrap')
    mask=(sj>threshold*j.mean())&(sp<.25*V**2)
    labs,n=labels(mask);rho=positive_energy(q,p,L);rows=[]
    for i in range(1,n+1):
        m=labs==i
        if m.sum()<4:continue
        Q=float(j[m].sum()*dv)
        if Q<100:continue
        where=np.array(np.where(m)).T;w=np.maximum(j[m],0);centre=[];strength=[]
        for ax in range(3):
            z=np.sum(w*np.exp(2j*np.pi*where[:,ax]/N))/w.sum()
            centre.append(float(np.angle(z)%(2*np.pi)*L/(2*np.pi)));strength.append(float(abs(z)))
        delta=((where*dx-centre+L/2)%L)-L/2;rr=np.linalg.norm(delta,axis=1);order=np.argsort(rr)
        r90=float(rr[order][np.searchsorted(np.cumsum(w[order]),.9*w.sum())])
        omega=float(j[m].sum()/(q[0][m]**2+q[1][m]**2).sum())
        rows.append({'Q_mask':Q,'E_mask':float(rho[m].sum()*dv),
                     'E_mask_over_Q':float(rho[m].sum()*dv/Q),'centre':centre,
                     'r90_mask':r90,'centre_strength':min(strength),'omega_current':omega,
                     'cells':int(m.sum()),'compact':bool(r90<L/4 and min(strength)>.5)})
    rows.sort(key=lambda x:x['Q_mask'],reverse=True)
    return {'threshold':threshold,'count':len(rows),'mask_charge_fraction':float(j[mask].sum()/j.sum()),
            'mask_volume_fraction':float(mask.mean()),'objects':rows},labs

def analyze(path,object_run=False,spheres=False):
    N,L,t,a=snapshot(path);q=np.asarray(a[:3]);p=np.asarray(a[3:]);j=q[0]*p[1]-q[1]*p[0]
    rho=positive_energy(q,p,L);dv=(L/N)**3
    cats=[]
    for threshold in (2.,4.,8.):
        c,_=catalog(q,p,L,threshold);cats.append(c)
    r={'N':N,'L':L,'t':t,'E':float(rho.sum()*dv),'Q':float(j.sum()*dv),
       'charge_contrast':float(j.std()/j.mean()),'negative_local_current_over_net':float(np.maximum(-j,0).sum()/j.sum()),
       'negative_current_is_particle_count':False,'phi_negative_volume_fraction':float(np.mean(q[2]<0)),
       'phi_mean_over_v':float(q[2].mean()/V),'phi_rms_over_v':float(np.sqrt(np.mean(q[2]**2))/V),
       'catalogs':cats,'positive_vs_signed_energy_relative':float(abs(rho.sum()/local_energy(q,p,L).sum()-1))}
    # Fraction of resolved kinetic+gradient energy near the lattice edge.
    f=np.fft.fftfreq(N)*2*np.pi;fr=np.fft.rfftfreq(N)*2*np.pi;dx=L/N
    one=lambda z:4*np.sin(z/2)**2*(1+np.sin(z/2)**2/3)/dx**2
    k2=one(f)[:,None,None]+one(f)[None,:,None]+one(fr)[None,None,:]
    qk=rfftn(q,axes=(-3,-2,-1),norm='forward');pk=rfftn(p,axes=(-3,-2,-1),norm='forward')
    w=np.full(N//2+1,2.);w[[0,-1]]=1
    ek=.5*np.sum(abs(pk)**2+k2*abs(qk)**2,axis=0)*w[None,None,:]
    high=np.maximum(np.maximum(abs(f)[:,None,None],abs(f)[None,:,None]),abs(fr)[None,None,:])>.75*np.pi
    r['axial_UV75_fraction_kinetic_gradient']=float(ek[high].sum()/ek.sum())
    if object_run:
        r['central_spheres']=sphere_diagnostics(q,p,L,[L/2]*3,(8.,12.,16.,20.,24.))
        r['max_X_over_v']=float(np.sqrt((q[0]**2+q[1]**2).max())/V)
        r['min_phi_over_v']=float(q[2].min()/V)
    if spheres:
        for c in cats:
            for obj in c['objects'][:20]:obj['spheres']=sphere_diagnostics(q,p,L,obj['centre'])
    return r

def collect():
    # Merge per-run records after execution; concurrent histories write only their own raw files.
    from run import CASES
    runs={}
    for tag in CASES:
        p=WORK/tag/'result.json'
        if p.exists():runs[tag]=json.loads(p.read_text())
    save('RUNS.json',{'cases':runs})
    out={}
    if (ROOT/'ANALYSIS.json').exists() or (ROOT/'ANALYSIS.json.gz').exists():out=read_analysis()['cases']
    for tag,run in runs.items():
        if tag in out:continue
        rows=[]
        for p in run['snapshots']:
            t=snapshot(p)[2]
            rows.append(analyze(p,object_run=run['args']['kind']=='object',spheres=t>=run['args']['T']-60))
        out[tag]={'rows':rows,'physical_eta':None,'physical_pair_ratio':None,'physical_f_sym':None}
        write_analysis({'cases':out,'scope':__doc__})
        print(json.dumps({'tag':tag,'last':{k:rows[-1][k] for k in ('t','charge_contrast','phi_negative_volume_fraction','axial_UV75_fraction_kinetic_gradient')},
                          'catalog4':{k:rows[-1]['catalogs'][1][k] for k in ('count','mask_charge_fraction')}}),flush=True)

if __name__=='__main__':collect()
