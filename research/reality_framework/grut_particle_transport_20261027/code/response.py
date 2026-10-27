"""Matrix canonical response on actual nonlinear fields, with phase/probe controls."""
from pt_common import *
from pt_production import mapping,selected
from concurrent.futures import ProcessPoolExecutor,as_completed

def run(s0,seed):
    N,L=48,12.;tag=f's{s0:g}_seed{seed}';z=np.load(WORK/'base48'/f'state_s{s0:g}.npy')
    lib=library();set_box(L);rng=np.random.default_rng(seed);source=rng.normal(size=(N,N,N))
    tangent=np.zeros((3,6,N,N,N))
    for b in range(3):tangent[b,3+b]=source
    where,half,flip,kc,kl=mapping(N,L)
    src=fftn(source,norm='forward').ravel()[where];sid=np.floor(kc/(2*np.pi/L)+.5).astype(int)
    shells=np.unique(sid);shells=shells[shells>0];indices=[np.flatnonzero(sid==n) for n in shells]
    denominators=[np.sum(abs(src[ix])**2) for ix in indices]
    times=np.arange(151)*.04;R=np.zeros((len(times),len(shells),3,3),complex);previous=0.
    for j,t in enumerate(times):
        if j:step(lib,z,t-previous,.004,tangent=tangent,threads=2)
        previous=t
        for b in range(3):
            fh=rfftn(tangent[b,:3],axes=(-3,-2,-1),norm='forward');qh=selected(fh,half,flip)
            for a in range(3):
                for k,ix in enumerate(indices):R[j,k,a,b]=np.sum(qh[a,ix]*src[ix].conj())/denominators[k]
    normal=.5*(R[:,:,0,0]+R[:,:,1,1]+1j*(R[:,:,1,0]-R[:,:,0,1]))
    anomalous=.5*(R[:,:,0,0]-R[:,:,1,1]+1j*(R[:,:,1,0]+R[:,:,0,1]))
    w=np.linspace(-20,20,2001);dw=w[1]-w[0];taper=np.cos(np.pi*times/(2*times[-1]))**2
    trap=np.full(len(times),.04);trap[[0,-1]]*=.5
    transform=np.exp(-1j*w[:,None]*times[None,:]);Gr=transform@(normal*(trap*taper)[:,None]);rho=-2*Gr.imag
    rows=[]
    for k,shell in enumerate(shells):
        rp=rho[:,k];positive=w>0;negative=w<0
        ip=np.flatnonzero(positive)[np.argmax(rp[positive])];im=np.flatnonzero(negative)[np.argmax(-rp[negative])]
        zp=float(np.sum((w*rp)[positive])*dw/(2*np.pi));zm=float(np.sum((w*rp)[negative])*dw/(2*np.pi))
        rows.append(dict(shell=int(shell),k_center=float(shell*2*np.pi/L),modes=len(indices[k]),
            initial_slope=float(normal[1,k].real/.04),Z_positive=zp,Z_negative=zm,Z_sum=zp+zm,
            response_peak_positive=float(w[ip]),response_peak_negative=float(w[im]),
            anomalous_over_normal_RMS=float(np.linalg.norm(anomalous[:,k])/np.linalg.norm(normal[:,k]))))
    np.savez_compressed(ROOT/f'RESPONSE_{tag}.npz',s0=s0,delay=times,shell=shells,normal=normal,anomalous=anomalous,
        phi_to_X=R[:,:,:2,2],X_to_phi=R[:,:,2,:2],omega=w,rho=rho)
    result=dict(s0=s0,probe_seed=seed,T=6.,ds=.04,maximum_step=.004,rows=rows,
        interpretation='Classical canonical tangent on a nonstationary background. A stochastic spatial source estimates shell traces. The f-sum is a canonical numerical check, not a proof of quasiparticles. Spectral resolution and finite-probe/phase errors are retained; no FDT or on-shell occupation is assumed.',physical_nbar=None)
    save(f'RESPONSE_{tag}.json',result);print('RESPONSE',tag,'low-k',[(r['k_center'],r['response_peak_positive'],r['anomalous_over_normal_RMS']) for r in rows[:4]],flush=True)
    return result

def main():
    cases=[(20.,7),(20.44,7),(20.,11)]
    with ProcessPoolExecutor(max_workers=int(os.environ.get('GRUT_RESPONSE_WORKERS','1'))) as pool:
        fs=[pool.submit(run,*c) for c in cases]
        for f in as_completed(fs):f.result()

if __name__=='__main__':main()
