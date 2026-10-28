"""Classical free-wave occupation outside a supplied isolated core.

This is a test of the object perturbation, not a primordial/cosmological population.
Interaction and taper-gradient sensitivity are always reported. Periodic recurrences
remain part of the experiment; an exterior count is not an infinite-time escaped yield.
"""
from common import *
from analyze import positive_energy

def measure(path,cut,width):
    N,L,t,a=snapshot(path);q=np.asarray(a[:3]);p=np.asarray(a[3:]);dx=L/N
    axis=(np.arange(N)-N//2)*dx
    rr=np.sqrt(axis[:,None,None]**2+axis[None,:,None]**2+axis[None,None,:]**2)
    u=np.clip((rr-cut)/width,0,1);W=np.sin(np.pi*u/2)**2
    qw=q[:2]*W;pw=p[:2]*W
    qk=rfftn(qw,axes=(-3,-2,-1),norm='forward');pk=rfftn(pw,axes=(-3,-2,-1),norm='forward')
    ff=np.fft.fftfreq(N)*2*np.pi;fr=np.fft.rfftfreq(N)*2*np.pi
    one=lambda z:4*np.sin(z/2)**2*(1+np.sin(z/2)**2/3)/dx**2
    k2=one(ff)[:,None,None]+one(ff)[None,:,None]+one(fr)[None,None,:];om=np.sqrt(1+k2)
    weight=np.full(N//2+1,2.);weight[[0,-1]]=1;weight=weight[None,None,:]
    actions=(np.sum(abs(pk)**2,axis=0)+om**2*np.sum(abs(qk)**2,axis=0))/(2*om)
    current=(np.conj(qk[0])*pk[1]-np.conj(qk[1])*pk[0]).real
    totalN=float(np.sum(weight*actions)*L**3);totalQ=float(np.sum(weight*current)*L**3)
    energy=float(np.sum(weight*om*actions)*L**3)
    r2=q[0]**2+q[1]**2
    interaction=np.abs(.5*(G2*q[2]**2-1)*r2+LX/4*r2*r2)
    # Compare the tapered free energy with the weighted untapered free energy.
    ge=np.zeros_like(r2)
    for z in q[:2]:
        for ax in range(3):
            forward=np.roll(z,-1,ax)-z;backward=z-np.roll(z,1,ax)
            ge+=.5/dx**2*(.5*(forward**2+backward**2)+(forward-backward)**2/12)
    untapered=.5*np.sum(p[:2]**2,axis=0)+.5*r2+ge
    weighted_energy=float(np.sum(W*W*untapered)*dx**3)
    realQ=float(np.sum(W*W*(q[0]*p[1]-q[1]*p[0]))*dx**3)
    return {'t':t,'cut_radius':cut,'taper_width':width,'Nplus_reference':(totalN+totalQ)/2,
            'Nminus_reference':(totalN-totalQ)/2,'Q_tapered':totalQ,'E_free_tapered':energy,
            'Nminus_over_original_Q':(totalN-totalQ)/400000,
            'absolute_interaction_over_free_energy':float(np.sum(W*W*interaction)*dx**3/energy),
            'taper_gradient_energy_fraction':(energy-weighted_energy)/energy,
            'Fourier_real_space_charge_relative_error':abs(totalQ-realQ)/max(1,abs(realQ))}

def main():
    out={}
    for tag in ('object64','object64_perturbed','object96_perturbed'):
        path=WORK/tag/'result.json'
        if not path.exists():continue
        r=json.loads(path.read_text());rows=[]
        for p in r['snapshots']:
            if snapshot(p)[2]<100:continue
            for cut,width in ((12.,4.),(16.,4.),(16.,8.),(20.,4.)):
                rows.append(measure(p,cut,width))
        out[tag]={'rows':rows}
    save('EXTERIOR_RESULTS.json',{'scope':__doc__,'cases':out,'physical_cosmic_p':None})
    for tag,v in out.items():
        x=v['rows'][-4:]
        print(json.dumps({'tag':tag,'last':x}),flush=True)

if __name__=='__main__':main()
