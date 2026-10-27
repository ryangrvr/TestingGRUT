"""Conservative field diagnostics; no particle or bound-core interpretation by fiat."""
from scipy.ndimage import gaussian_filter,label,generate_binary_structure
from scipy.fft import rfftn
import hashlib
from pilot import *

def periodic_labels(mask):
    lab,count=label(mask,structure=generate_binary_structure(3,1))
    if not count:return lab,0
    parent=np.arange(count+1)
    def find(x):
        while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
        return x
    for ax in range(3):
        a=np.take(lab,0,axis=ax);b=np.take(lab,-1,axis=ax);ok=(a>0)&(b>0)
        pairs=np.unique(np.stack([a[ok],b[ok]],axis=1),axis=0)
        for x,y in pairs:
            x=find(x);y=find(y)
            if x!=y:parent[max(x,y)]=min(x,y)
    roots=np.array([find(i) for i in range(count+1)]);unique=np.unique(roots[1:]);mapping=np.zeros(count+1,dtype=int)
    for i,r in enumerate(unique,1):mapping[roots==r]=i
    return mapping[lab],len(unique)

def local_energy(q,p,L):
    N=q.shape[-1];dx=L/N;r2=q[0]**2+q[1]**2
    pot=C.LP/4*(q[2]**2-C.NU**2)**2+C.G**2/2*q[2]**2*r2+C.LX/4*r2*r2
    grad=np.zeros_like(r2)
    for f in q:
        for ax in range(3):
            grad+=.5/dx**2*((4/3)*(np.roll(f,-1,ax)-f)**2-(np.roll(f,-2,ax)-f)**2/12)
    return .5*np.sum(p*p,axis=0)+pot+grad

def catalog(q,p,L,threshold=4.,smooth=1.):
    N=q.shape[-1];dx=L/N;dv=dx**3;j=q[0]*p[1]-q[1]*p[0];mean=float(j.mean())
    sj=gaussian_filter(j,smooth/dx,mode='wrap');sp=gaussian_filter(q[2]**2,smooth/dx,mode='wrap')
    mask=(sj>threshold*mean)&(sp<.25*C.NU**2);labs,n=periodic_labels(mask)
    rho=local_energy(q,p,L);rows=[];centres=[]
    ids=np.arange(N);coord=(ids+.0)*dx
    for index in range(1,n+1):
        m=labs==index
        if m.sum()<4:continue
        charge=float(j[m].sum()*dv)
        if charge<100:continue
        weights=np.maximum(j[m],0);where=np.array(np.where(m)).T
        centre=[];coherence=[]
        for ax in range(3):
            val=np.sum(weights*np.exp(2j*np.pi*where[:,ax]/N))/weights.sum();coherence.append(abs(val));centre.append(float(np.angle(val)%(2*np.pi)/(2*np.pi)*L))
        distances=np.minimum(abs(where*dx-centre),L-abs(where*dx-centre));rr=np.linalg.norm(distances,axis=1);order=np.argsort(rr)
        r90=float(rr[order][np.searchsorted(np.cumsum(weights[order]),.9*weights.sum())])
        r2=q[0][m]**2+q[1][m]**2;omega=float(j[m].sum()/max(r2.sum(),1e-300))
        rot=.5*j[m]**2/np.maximum(r2,1e-30);kx=.5*(p[0][m]**2+p[1][m]**2)
        energy_mask=float(rho[m].sum()*dv)
        rows.append(dict(id=index,cells=int(m.sum()),Q_mask=charge,E_mask_m=energy_mask,E_mask_over_Qm=energy_mask/charge,
            centre=centre,circular_centre_strength=min(coherence),r90_mask_over_m=r90,
            omega_from_current_over_m=omega,radial_X_kinetic_fraction=float(max(0.,(kx.sum()-rot.sum())/max(kx.sum(),1e-300))),
            compact=bool(r90<L/4 and min(coherence)>.5),
            scope='Superlevel-mask candidate. Its charge omits the tail; energy includes background/crossing edges. Not a bound-object mass.'))
    rows.sort(key=lambda r:r['Q_mask'],reverse=True)
    return dict(threshold=threshold,smoothing_length_over_m=smooth,candidate_count=len(rows),
        Q_fraction_in_masks=float(j[mask].sum()/j.sum()),mask_volume_fraction=float(mask.mean()),candidates=rows),labs

def spectra(q,p,L):
    N=q.shape[-1];vol=L**3;k2=symbols(N,L);om=np.sqrt(1+k2)
    qh=rfftn(q,axes=(-3,-2,-1),norm='forward');ph=rfftn(p,axes=(-3,-2,-1),norm='forward')
    weights=np.full(N//2+1,2.);weights[[0,-1]]=1.;weights=weights[None,None,:]
    kinetic=.5*np.sum(abs(ph)**2,axis=0);grad=.5*k2*np.sum(abs(qh)**2,axis=0)
    nref=(np.sum(abs(ph[:2])**2,axis=0)+om**2*np.sum(abs(qh[:2])**2,axis=0))/(2*om)
    qref=(np.conj(qh[0])*ph[1]-np.conj(qh[1])*ph[0]).real
    nplus=.5*np.sum(weights*(nref+qref));nminus=.5*np.sum(weights*(nref-qref));Q=float(np.mean(q[0]*p[1]-q[1]*p[0]))
    eref=float(np.sum(weights*om*nref));r2=q[0]**2+q[1]**2
    residual=float(np.mean(.5*abs(C.G**2*q[2]**2-1)*r2+C.LX/4*r2*r2)/eref)
    freq=np.fft.fftfreq(N,d=L/N)*2*np.pi;last=np.fft.rfftfreq(N,d=L/N)*2*np.pi
    physicalk=np.sqrt(freq[:,None,None]**2+freq[None,:,None]**2+last[None,None,:]**2)
    uv=physicalk>.5*np.pi/(L/N);dynamical=float(np.sum(weights*(kinetic+grad)))
    j=q[0]*p[1]-q[1]*p[0];jh=rfftn(j-j.mean(),norm='forward');power=weights*abs(jh)**2/Q**2
    edges=np.r_[0.,np.geomspace(2*np.pi/L*.999,physicalk.max()*1.001,33)];bins=np.clip(np.digitize(physicalk.ravel(),edges)-1,0,len(edges)-2)
    chargepower=np.bincount(bins,weights=power.ravel(),minlength=len(edges)-1)
    dynpower=np.bincount(bins,weights=(weights*(kinetic+grad)).ravel(),minlength=len(edges)-1)
    maxaxis=np.maximum(np.maximum(abs(freq)[:,None,None],abs(freq)[None,:,None]),abs(last)[None,None,:])
    edgeshell={str(f):float(np.sum((weights*(kinetic+grad))[maxaxis>f*np.pi/(L/N)])/dynamical) for f in [.5,.75,.9]}
    return dict(k_edges=edges,charge_contrast_variance_per_bin=chargepower,
        kinetic_plus_gradient_per_bin=dynpower,
        charge_Parseval_error=float(abs(power.sum()-np.var(j)/Q**2)),
        UV_fraction_of_kinetic_plus_gradient=float(np.sum((weights*(kinetic+grad))[uv])/dynamical),
        UV_definition='Euclidean continuum |k| > half the axial Nyquist momentum; potential energy is not assigned to k.',
        boundary_shell_energy_fractions=edgeshell,
        boundary_shell_definition='max(|k_x|,|k_y|,|k_z|) exceeds the given fraction of axial Nyquist; uses kinetic+gradient only.',
        vacuum_reference=dict(nplus_per_m3=float(nplus),nminus_per_m3=float(nminus),charge_identity_relative=float(abs((nplus-nminus)/Q-1)),
            energy_X_ref_per_m4=eref,absolute_interaction_residual_over_ref=residual,
            physical_particle_interpretation=False,
            scope='Vacuum-mass classical reference projection only. Interaction residual and coherent/cores/domain decomposition prevent automatic particle interpretation.'))

def analyze_snapshot(path,details=False):
    N,L,t,raw=snapshot(path);q=np.array(raw[:3],dtype=float);p=np.array(raw[3:],dtype=float);j=q[0]*p[1]-q[1]*p[0]
    cats=[];mainlabs=None
    for threshold in [2.,4.,8.]:
        cat,labs=catalog(q,p,L,threshold);cats.append(cat)
        if threshold==4:mainlabs=labs
    spec=spectra(q,p,L);rho=local_energy(q,p,L)
    row=dict(t=t,N=N,L=L,nQ=float(j.mean()),rho=float(rho.mean()),charge_contrast_RMS=float(j.std()/j.mean()),
        negative_current_absolute_fraction=float(np.maximum(-j,0).sum()/np.maximum(j,0).sum()),
        negative_current_scope='Local classical current reversal, not an antiparticle abundance.',
        phi_mean_over_v=float(q[2].mean()/C.NU),phi_rms_over_v=float(np.sqrt(np.mean(q[2]**2))/C.NU),
        phi_negative_fraction=float(np.mean(q[2]<0)),catalogs=cats,spectra=spec)
    if details:
        return row,q,p,j,mainlabs
    return row

def main():
    runs=json.loads((ROOT/'RUNS.json').read_text())['cases'];out={}
    codehash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    if (ROOT/'ANALYSIS.json').exists():
        cached=json.loads((ROOT/'ANALYSIS.json').read_text())
        if cached.get('analysis_code_sha256')==codehash:out=cached['cases']
    for tag,r in runs.items():
        if 'error' in r or tag in out:continue
        rows=[]
        for path in r['snapshots']:
            rows.append(analyze_snapshot(path))
        out[tag]=dict(rows=rows,scope='Cold classical periodic pilot; candidate masks and reference spectra, no physical f_sym.')
        save('ANALYSIS.json',dict(cases=out,analysis_code_sha256=codehash,physical_f_sym=None,physical_formation_kernel=None))
        print('ANALYZED',tag,[(x['t'],x['catalogs'][1]['candidate_count'],x['catalogs'][1]['Q_fraction_in_masks'],x['spectra']['UV_fraction_of_kinetic_plus_gradient']) for x in rows],flush=True)

if __name__=='__main__':main()
