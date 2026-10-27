"""Exact active-sector ledgers and explicitly non-particle reference spectra."""
from en_common import *

def spectral(q,p,L):
    N=q.shape[-1];volume=L**3;k2=symbol(N,L);k=np.sqrt(k2)
    qh=rfftn(q,axes=(-3,-2,-1),norm='forward');ph=rfftn(p,axes=(-3,-2,-1),norm='forward')
    weights=np.full(N//2+1,2.);weights[[0,-1]]=1.;weights=weights[None,None,:]
    kinetic=.5*np.sum(abs(ph)**2,axis=0);grad=.5*k2*np.sum(abs(qh)**2,axis=0)
    dynamical=weights*(kinetic+grad);physical=np.fft.fftfreq(N)*2*np.pi*N/L;last=np.fft.rfftfreq(N)*2*np.pi*N/L
    maxaxis=np.maximum(np.maximum(abs(physical)[:,None,None],abs(physical)[None,:,None]),abs(last)[None,None,:])
    shell=float(dynamical[maxaxis>.75*np.pi*N/L].sum()/dynamical.sum())
    kphys=np.sqrt(physical[:,None,None]**2+physical[None,:,None]**2+last[None,None,:]**2)
    bins=np.linspace(0,kphys.max()*1.00001,51);ids=np.clip(np.digitize(kphys.ravel(),bins)-1,0,len(bins)-2)
    power=np.bincount(ids,weights=dynamical.ravel(),minlength=len(bins)-1)
    r2=q[0]**2+q[1]**2;hf=2*LX*np.mean(r2)+G2*np.mean(q[2]**2)
    refs=[];Qm=(np.conj(qh[0])*ph[1]-np.conj(qh[1])*ph[0]).real;nonzero=k2>0
    Qtotal=float(np.mean(q[0]*p[1]-q[1]*p[0]));Qmean=float(q[0].mean()*p[1].mean()-q[1].mean()*p[0].mean())
    for mass2,name in [(0.,'massless'),(hf,'isotropic_Hartree_reference')]:
        omega=np.sqrt(k2+mass2);safe=np.where(nonzero,omega,1.)
        nr=(np.sum(abs(ph[:2])**2,axis=0)+omega**2*np.sum(abs(qh[:2])**2,axis=0))/(2*safe)
        npart=float(np.sum((weights*(nr+Qm))[nonzero])/2);nanti=float(np.sum((weights*(nr-Qm))[nonzero])/2)
        identity=abs((npart-nanti)-(Qtotal-Qmean))/max(abs(Qtotal),1.)
        refs.append(dict(name=name,mass2=mass2,raw_nonzero_nplus=npart,raw_nonzero_nminus=nanti,
            raw_nminus_over_net_current=nanti/Qtotal,charge_identity_relative=identity,
            physical_particle_interpretation=False,
            scope='Classical finite-band field reference. No zero-point subtraction, quasiparticle certification or full quantum/kinetic population.'))
    potential=LX/4*r2*r2+G2/2*q[2]**2*r2+LP/4*q[2]**4
    E=energy(q,p,L);recon=float(abs((dynamical.sum()+potential.mean())/E-1))
    return dict(boundary_shell_075=shell,kinetic_gradient_k_edges=bins,kinetic_gradient_spectrum=power,
        spectra_definition='Forward FFT and exact FD4 symbol. Potential energy not assigned to individual k. Boundary shell uses max physical momentum component>.75 axial Nyquist.',
        reference_projections=refs,energy_Parseval_reconstruction=recon,
        scalar_Hartree_mass2=G2*np.mean(r2)+3*LP*np.mean(q[2]**2),
        light_Hartree_mass2=hf,physical_particle_distribution=None)

def trace_quantities(trace):
    rho=trace['rho'];meanr2=trace['q0']**2+trace['q1']**2
    mix0=G2/2*trace['qphi']**2*meanr2
    EX0=.5*(trace['p0']**2+trace['p1']**2)+LX/4*meanr2**2+mix0
    EP0=.5*trace['pphi']**2+LP/4*trace['qphi']**4
    Ephi=trace['Kphi']+trace['Gphi']+trace['Uphi'];EX=rho-Ephi
    meanQ=trace['q0']*trace['p1']-trace['q1']*trace['p0']
    Emin=.75*LX**(1/3)*abs(trace['nQ'])**(4/3)
    contrast=np.sqrt(np.maximum(trace['charge_RMS']**2-trace['nQ']**2,0))/abs(trace['nQ'])
    fluct=np.maximum(trace['X2']-meanr2,0)/np.maximum(trace['X2'],1e-300)
    integral=np.r_[0.,np.cumsum(.5*(trace['work_phi_to_X'][1:]+trace['work_phi_to_X'][:-1])*np.diff(trace['eta']))]
    work_error=float(np.max(abs((EX-EX[0])-integral))/rho[0])
    def first(v,level):
        ix=np.flatnonzero(v>=level);return float(trace['tau'][ix[0]]) if len(ix) else None
    return dict(tau=trace['tau'],eta=trace['eta'],rho=rho,nQ=trace['nQ'],
        Ephi_fraction=Ephi/rho,phi_plus_interaction_energy_fraction=(Ephi+trace['Umix'])/rho,
        spatial_mean_energy_fraction=(EX0+EP0)/rho,
        energy_outside_spatial_mean_fraction=1-(EX0+EP0)/rho,
        spatial_mean_charge_fraction=meanQ/trace['nQ'],
        X_field_fluctuation_variance_fraction=fluct,charge_contrast_RMS=contrast,
        minimum_energy_at_total_charge_fraction=Emin/rho,
        scalar_exchange_work_identity_relative=work_error,
        first_X_variance_half=first(fluct,.5),first_charge_contrast_one=first(contrast,1.),
        first_phi_fraction_001=first(Ephi/rho,.01),first_phi_fraction_01=first(Ephi/rho,.1),
        first_phi_fraction_05=first(Ephi/rho,.5))

def main():
    runs=read('RUNS.json')['cases'];out={};curves={}
    for tag,r in runs.items():
        if 'error' in r:raise RuntimeError(r)
        trace=np.genfromtxt(WORK/tag/'state.tsv',names=True,delimiter='\t');v=trace_quantities(trace)
        N,L,s,raw=snapshot(r['snapshots'][-1]);q=np.array(raw[:3]);p=np.array(raw[3:]);sp=spectral(q,p,L)
        A=r['parameters']['clock'];a=A['a0']+A['C']*trace['eta'];mL2=30.**2-np.sqrt(2)*M.MU*40.
        frozen_est=3*M.KAP**2*trace['max_X2']/(4*LX*mL2*a*a)
        row=dict(tag=tag,N=N,L=L,initial_charge_offset=r['parameters']['sample_charge_offset'],
            final={key:float(value[-1]) for key,value in v.items() if isinstance(value,np.ndarray)},
            thresholds={key:value for key,value in v.items() if key.startswith('first_')},
            scalar_exchange_work_identity_relative=v['scalar_exchange_work_identity_relative'],
            frozen_Sigma_quadratic_Y_force_estimate_max=float(frozen_est.max()),
            frozen_Sigma_quadratic_Y_force_estimate_final=float(frozen_est[-1]),
            approximation_scope='Frozen-Sigma static quadratic-Y estimate only; not a nonlinear full-model error bound.',
            spectra=sp,conformal_light_crossings=float(s/L),
            physical_f_sym=None,physical_nbar_per_net_charge=None)
        out[tag]=row
        stride=1 if tag=='N48' else (2 if tag=='halfstep48' else 10)
        for key,value in v.items():
            if isinstance(value,np.ndarray):curves[tag+'_'+key]=value[::stride]
        save('FIELD_RESULTS.json',dict(cases=out,physical_f_sym=None,physical_bound_population=None,
            physical_particle_distribution=None,full_model_nonlinear=False))
        print('ANALYSIS',tag,'phi',row['final']['Ephi_fraction'],'meanQ',row['final']['spatial_mean_charge_fraction'],
              'refanti',[x['raw_nminus_over_net_current'] for x in sp['reference_projections']],'UV',sp['boundary_shell_075'],flush=True)
    np.savez_compressed(ROOT/'FIELD_CURVES.npz',**curves)
    # Exactly labelled, block-averaged field views; not used for analysis.
    r=runs['N48'];N,L,s,raw=snapshot(r['snapshots'][-1]);q=np.array(raw[:3]);p=np.array(raw[3:])
    j=q[0]*p[1]-q[1]*p[0];factor=N//24
    def block(x):return x.reshape(24,factor,24,factor,24,factor).mean(axis=(1,3,5))
    np.savez_compressed(ROOT/'VIEW_FIELDS.npz',charge=block(j),phi=block(q[2]),X2=block(q[0]**2+q[1]**2),L=L,eta=s)

if __name__=='__main__':main()
