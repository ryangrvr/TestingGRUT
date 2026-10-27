"""Streaming nonlinear histories; unique module name avoids legacy collisions."""
from pt_common import *
from prepare import initialize
from concurrent.futures import ProcessPoolExecutor,as_completed

CASES=[
 dict(tag='base48',N=48,L=12.,seed=1,width=None,hmax=.002),
 dict(tag='base48_half',N=48,L=12.,seed=1,width=None,hmax=.001),
 dict(tag='IR12_s1_w1',N=32,L=12.,seed=1,width=1,hmax=.004),
 dict(tag='IR12_s2_w1',N=32,L=12.,seed=2,width=1,hmax=.004),
 dict(tag='IR12_s1_w3',N=32,L=12.,seed=1,width=3,hmax=.004),
 dict(tag='IR12_fine_s1_w1',N=48,L=12.,seed=1,width=1,hmax=.004),
 dict(tag='IR24_s1_w1',N=64,L=24.,seed=1,width=1,hmax=.004),
 dict(tag='IR24_s2_w1',N=64,L=24.,seed=2,width=1,hmax=.004)]

def mapping(N,L,kmax=8.):
    kc2,kl=geometry(N,L);where=np.flatnonzero(kc2.ravel()<=kmax*kmax+1e-12)
    x,y,z=np.unravel_index(where,(N,N,N));flip=z>N//2
    xx=np.where(flip,(-x)%N,x);yy=np.where(flip,(-y)%N,y);zz=np.where(flip,(-z)%N,z)
    half=np.ravel_multi_index((xx,yy,zz),(N,N,N//2+1))
    return where,half,flip,np.sqrt(kc2.ravel()[where]),kl.ravel()[where]

def selected(fh,half,flip):
    a=fh.reshape(fh.shape[0],-1)[:,half].copy();a[:,flip]=a[:,flip].conj();return a

def observe(z,L,kh2,weights,bins,nbins):
    f=rfftn(z,axes=(-3,-2,-1),norm='forward');qh,ph=f[:3],f[3:]
    means=z.mean(axis=(1,2,3));vec=means[[0,1,3,4,2,5]]
    r2=z[0]**2+z[1]**2;phi=z[2];j=z[0]*z[4]-z[1]*z[3]
    U=np.array([LX/4*np.mean(r2*r2),LP/4*np.mean(phi**4),G2/2*np.mean(phi*phi*r2)])
    K=.5*np.sum(weights*np.sum(abs(ph)**2,axis=0));G=.5*np.sum(weights*kh2*np.sum(abs(qh)**2,axis=0));rho=K+G+U.sum()
    Q=j.mean();Q0=means[0]*means[4]-means[1]*means[3]
    Kphi=.5*np.sum(weights*abs(ph[2])**2);Gphi=.5*np.sum(weights*kh2*abs(qh[2])**2)
    Em=.5*np.sum(means[3:]**2)+LX/4*(means[0]**2+means[1]**2)**2+LP/4*means[2]**4+G2/2*means[2]**2*(means[0]**2+means[1]**2)
    Xq=abs(qh[0])**2+abs(qh[1])**2;Xp=abs(ph[0])**2+abs(ph[1])**2
    Qh=(qh[0].conj()*ph[1]-qh[1].conj()*ph[0]).real
    hf=2*LX*np.mean(r2)+G2*np.mean(phi*phi);W=weights.copy();W[0,0,0]=0.
    refs=[];tails=[];kcont2,_=geometry(z.shape[-1],L,half=True)
    for m2 in [0.,7.,12.,20.,26.,hf]:
        om=np.sqrt(np.where(kh2>0,kh2+m2,1.));nm=.5*((Xp/om+om*Xq)/2-Qh)
        refs.append(float(np.sum(W*nm)/Q));tails.append(float(np.sum(W*nm*(kcont2>64))/Q))
    row=np.r_[rho,Q,Kphi+Gphi+U[1],U[2],1-Em/rho,Q0/Q,hf,
        np.sqrt(max(0,np.mean(j*j)-Q*Q))/abs(Q),1-(means[0]**2+means[1]**2)/np.mean(r2),means,refs,tails]
    v=f[[0,1,3,4,2,5]];cov=np.zeros((nbins,6,6));covk=cov.copy()
    for a in range(6):
      for b in range(a,6):
        val=W*np.real(v[a].conj()*v[b]);cv=np.bincount(bins,weights=val.ravel(),minlength=nbins)
        ck=np.bincount(bins,weights=(val*kh2).ravel(),minlength=nbins)
        cov[:,a,b]=cov[:,b,a]=cv;covk[:,a,b]=covk[:,b,a]=ck
    return f,row,cov,covk,np.outer(vec,vec),U

NAMES=['rho','nQ','Ephi','Emix','outside_mean','Q0_over_Q','Hartree_m2','charge_contrast','X_variance_fraction',
 'q0','q1','qphi','p0','p1','pphi',*[f'ref_{x}' for x in ['0','7','12','20','26','H']],*[f'tail_{x}' for x in ['0','7','12','20','26','H']]]

def run(case):
    start=time.monotonic();tag=case['tag'];N,L=case['N'],case['L'];directory=WORK/tag;directory.mkdir(exist_ok=True)
    z,params=initialize(N,L,case['seed'],case['width']);z=np.ascontiguousarray(z);rho0,Q0=invariants(z,L)
    params.update(case);threads=int(os.environ.get('GRUT_PARTICLE_THREADS','3'));params['threads']=threads
    set_box(L);lib=library();np.save(directory/'initial.npy',z)
    ds=.04;samples=np.arange(851)*ds
    where,half,flip,kc,kl=mapping(N,L)
    psi=np.lib.format.open_memmap(directory/'psi.npy',mode='w+',dtype=np.complex64,shape=(len(samples),len(where)))
    kc2,kh2=geometry(N,L,half=True);nbins=96;edges=np.linspace(0,24.,nbins+1)
    bins=np.clip(np.digitize(np.sqrt(kc2).ravel(),edges)-1,0,nbins-1)
    weights=np.full(N//2+1,2.);weights[[0,-1]]=1.;weights=np.broadcast_to(weights[None,None,:],kh2.shape).copy()
    # Covariances retain the zero mode separately as an averaged outer product.
    averages={};windows=[(14.,24.),(24.,34.)]
    for lo,hi in windows:
      for rule in ['conformal','cosmic']:
        averages[f'{lo:g}_{hi:g}_{rule}']=dict(cov=np.zeros((nbins,6,6)),covk=np.zeros((nbins,6,6)),mean_outer=np.zeros((6,6)),
          potential=np.zeros(3),weight=0.,rho=0.,Q=0.)
    traces=[];snap_s={14.,20.,20.24,20.44,20.68,24.,34.};previous=0.
    for index,s in enumerate(samples):
        if index:step(lib,z,s-previous,case['hmax'],threads=threads)
        previous=s
        f,row,cv,ck,mo,U=observe(z,L,kh2,weights,bins,nbins);v=selected(f[:2],half,flip)
        psi[index]=v[0]+1j*v[1];traces.append(row)
        if any(abs(s-x)<1e-9 for x in snap_s):np.save(directory/f'state_s{s:g}.npy',z)
        for lo,hi in windows:
            if s<lo-1e-10 or s>hi+1e-10:continue
            trap=.5 if min(abs(s-lo),abs(s-hi))<1e-10 else 1.
            for rule in ['conformal','cosmic']:
                w=trap*(1. if rule=='conformal' else 8.869218738806+1.9857962750634526*s)
                a=averages[f'{lo:g}_{hi:g}_{rule}'];a['cov']+=w*cv;a['covk']+=w*ck;a['mean_outer']+=w*mo
                a['potential']+=w*U;a['weight']+=w;a['rho']+=w*row[0];a['Q']+=w*row[1]
        if index in [0,250,500,750,850]:print('HISTORY',tag,'s',s,'phi+',float((row[2]+row[3])/row[0]),'seconds',round(time.monotonic()-start,1),flush=True)
    psi.flush();trace=np.array(traces);np.savez_compressed(directory/'trace.npz',s=samples,tau=clock(samples),values=trace,names=NAMES)
    np.savez_compressed(directory/'modes.npz',where=where,kcont=kc,kfd2=kl,L=L,N=N,s=samples,nQ=Q0)
    public={};recon=[]
    for key,a in averages.items():
        w=a.pop('weight');a={k:v/w for k,v in a.items()};C=a['cov'];K=a['covk'];M=a['mean_outer'];U=a['potential']
        Q=C[:,0,3].sum()-C[:,1,2].sum()+M[0,3]-M[1,2]
        energy=.5*(C[:,2,2]+C[:,3,3]+C[:,5,5]).sum()+.5*(K[:,0,0]+K[:,1,1]+K[:,4,4]).sum()+.5*(M[2,2]+M[3,3]+M[5,5])+U.sum()
        recon.append(dict(window=key,charge_relative=abs(Q/a['Q']-1),energy_relative=abs(energy/a['rho']-1)))
        for name,value in a.items():public[key+'_'+name]=value
    public['k_edges']=edges;np.savez_compressed(ROOT/f'AVERAGE_{tag}.npz',**public)
    result=dict(case=case,initial=params,seconds=time.monotonic()-start,output_ds=ds,end_s=34.,end_tau=float(clock(34.)),
        energy_max_relative_drift=float(np.max(abs(trace[:,0]/rho0-1))),charge_max_relative_drift=float(np.max(abs(trace[:,1]/Q0-1))),
        export_checks=recon,statistical_band_maximum_k=8.,temporal_field_storage='complex64 selected full-space psi Fourier coefficients; no vacuum subtraction',
        names=NAMES,initial_file_sha256=sha(directory/'initial.npy'),final_file_sha256=sha(directory/'state_s34.npy'),
        physical_nbar=None,physical_f_sym=None)
    save(f'RUN_{tag}.json',result);return result

def main():
    wanted=set(sys.argv[1:]);cases=[c for c in CASES if not wanted or c['tag'] in wanted]
    with ProcessPoolExecutor(max_workers=int(os.environ.get('GRUT_PARTICLE_WORKERS','2'))) as pool:
        futures={pool.submit(run,c):c['tag'] for c in cases}
        for f in as_completed(futures):
            r=f.result();print('DONE',r['case']['tag'],r['seconds'],r['energy_max_relative_drift'],r['charge_max_relative_drift'],flush=True)

if __name__=='__main__':main()
