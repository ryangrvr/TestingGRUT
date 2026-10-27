"""Explicit pure Gaussian completions of the previously omitted IR band."""
from pt_common import *
from concurrent.futures import ThreadPoolExecutor,as_completed
import neutral as NE

def keys_for(L,cutoff=8.):
    n=int(np.ceil(cutoff*L/(2*np.pi)));a=np.arange(-n,n+1)
    s=(a[:,None,None]**2+a[None,:,None]**2+a[None,None,:]**2).ravel()
    keys,count=np.unique(s[(s>0)&((2*np.pi/L)**2*s<=cutoff**2+1e-12)],return_counts=True)
    return keys,count

def factor(TT,PT,ks,lam,U,phi_mass,width):
    ww2=ks[:,None]**2+lam
    # No ground state exists for the negative-curvature oscillators.
    # Their widths are supplied positive initial-state data, not frequencies.
    w=np.sqrt(np.where(ww2>0,ww2,ks[:,None]**2+width*width))
    S=np.zeros((len(ks),12,12));S[:,:6,:6]=U[None,:,:]/np.sqrt(2*w[:,None,:])
    S[:,6:,6:]=U[None,:,:]*np.sqrt(w[:,None,:]/2)
    evolved=TT@S;wp=np.sqrt(ks*ks+phi_mass);sp=np.zeros((len(ks),2,2))
    sp[:,0,0]=1/np.sqrt(2*wp);sp[:,1,1]=np.sqrt(wp/2);ep=PT@sp
    out=np.zeros((len(ks),6,14));out[:,[0,1,3,4],:12]=evolved[:,[0,1,6,7],:]
    out[:,2,12:]=ep[:,0,:];out[:,5,12:]=ep[:,1,:]
    return out

def prepare():
    M=E.M;exe=M.build();ne=exe.parent/'neutral_pt'
    subprocess.run(['gcc','-O3','-std=c11','-DNEUTRAL_SECTOR','-I',str(exe.parent),str(E.T.MODE/'code/modes.c'),'-lm','-o',str(ne)],check=True)
    old=np.load(PARENT/'MATCH_L12.npz');lam=old['initial_mass_eigenvalues'];U=old['initial_mass_eigenvectors']
    v0,_=M.prepared(3.);phi_mass=NE.mass(v0,3.);rows=[]
    J=np.block([[np.zeros((6,6)),np.eye(6)],[-np.eye(6),np.zeros((6,6))]])
    for L in [12.,24.]:
        target=ROOT/f'PREPARATION_L{L:g}.npz';keys,deg=keys_for(L);ks=2*np.pi/L*np.sqrt(keys)
        F1=np.zeros((len(keys),6,14));F3=F1.copy();known=np.zeros(len(keys),bool)
        if L==12.:
            look={int(k):i for i,k in enumerate(old['keys'])}
            for j,key in enumerate(keys):
                if key not in look:continue
                ix=look[int(key)];f=np.zeros((6,14));a=old['evolved_scalar'][ix];p=old['evolved_phi'][ix]
                f[[0,1,3,4],:12]=a[[0,1,6,7]];f[2,12:]=p[0];f[5,12:]=p[1]
                F1[j]=F3[j]=f;known[j]=True
        todo=np.flatnonzero(~known);checks=[]
        def batch(start):
            ix=todo[start:start+96];kk=ks[ix]
            bg,T,c=M.integrate(exe,f'pt_L{L:g}_{start}',3.,.5,kk,np.array([0.,20.]),rtol=2e-12,e=0.)
            pb,P,pc=NE.run(ne,dict(tag=f'ptL{L:g}_{start}',beta=3.,Gamma_over_H=.5),kk,np.array([0.,20.]),rtol=2e-12,suffix='_pt')
            assert np.max(abs(bg-pb)/(1+abs(bg)))<1e-6
            return ix,T[-1],P[-1],c
        with ThreadPoolExecutor(max_workers=2) as pool:
            fs=[pool.submit(batch,start) for start in range(0,len(todo),96)]
            for future in as_completed(fs):
                ix,T,P,c=future.result();kk=ks[ix]
                F1[ix]=factor(T,P,kk,lam,U,phi_mass,1.);F3[ix]=factor(T,P,kk,lam,U,phi_mass,3.)
                sym=max(np.max(abs(t.T@J@t-J)) for t in T);det=max(abs(np.linalg.det(P)-1))
                checks.append(dict(kmin=kk.min(),kmax=kk.max(),shells=len(ix),symplectic_absolute=sym,phi_determinant_error=det,control=c))
                print('IR_MATCH',L,len(ix),'sym',sym,flush=True)
        np.savez_compressed(target,keys=keys,ks=ks,degeneracy=deg,factor_width1=F1,factor_width3=F3,
            light_mean_q=old['light_mean_q'],light_mean_p=old['light_mean_p'],a_match=old['a_match'],C_match=old['C_match'],
            initial_mass_eigenvalues=lam,initial_mass_eigenvectors=U)
        expected=[]
        for width,F in [(1.,F1),(3.,F3)]:
            cov=np.einsum('k,kia,kja->ij',deg/L**3,F,F)
            expected.append(dict(width=width,light_covariance_at_match=cov,
                charged_q_RMS_over_mean=float(np.sqrt(cov[0,0]+cov[1,1])/np.linalg.norm(old['light_mean_q'][:2])),
                charged_p_RMS_over_mean=float(np.sqrt(cov[3,3]+cov[4,4])/np.linalg.norm(old['light_mean_p'][:2]))))
        rows.append(dict(L=L,shells=len(keys),new_shells=len(todo),minimum_k=ks.min(),cutoff=8.,
            unstable_initial_oscillators=int(np.sum(ks[:,None]**2+lam<=0)),checks=checks,expected=expected))
    result=dict(rows=rows,old_initial_mass_eigenvalues=lam,old_ground_state_obstruction_k=np.sqrt(-lam[lam<0]),
        initial_rule='Stable oscillators retain omega=sqrt(k^2+eigenvalue). Negative-curvature oscillators use supplied positive widths sqrt(k^2+sigma^2), sigma=1 or3. Full scalar plus phi state is a pure minimum-uncertainty Gaussian in the e0 system, not an instantaneous ground state or selected cosmic state.',
        later_state='The light marginal need not be pure. Its classical realization has no vacuum subtraction or selected quantum completion.')
    save('PREPARATION_RESULTS.json',result)

def initialize(N,L,seed,width=None):
    if width is None:
        from initial import initialize as previous
        q,p,params=previous(N,L,seed);return np.r_[q,p],params
    d=np.load(ROOT/f'PREPARATION_L{L:g}.npz');lookup={int(x):i for i,x in enumerate(d['keys'])}
    F=d[f'factor_width{int(width)}'];coeff=np.zeros((6,N,N,N),complex);nmax=int(np.ceil(8*L/(2*np.pi)));count=0
    for i in range(-nmax,nmax+1):
      for j in range(-nmax,nmax+1):
       for k in range(-nmax,nmax+1):
        n=(i,j,k);first=next((a for a in n if a),0)
        if first<=0:continue
        key=i*i+j*j+k*k
        if key not in lookup:continue
        assert max(abs(i),abs(j),abs(k))<N/2
        rng=np.random.default_rng(np.random.SeedSequence([seed,int(L*100),i+128,j+128,k+128]))
        z=(rng.normal(size=14)+1j*rng.normal(size=14))/np.sqrt(2)
        light=F[lookup[key]]@z/np.sqrt(L**3)
        coeff[:,i%N,j%N,k%N]=light;coeff[:,(-i)%N,(-j)%N,(-k)%N]=light.conjugate();count+=2
    coeff[:,0,0,0]=np.r_[d['light_mean_q'],d['light_mean_p']]
    z=ifftn(coeff,axes=(-3,-2,-1),norm='forward');assert np.max(abs(z.imag))<1e-10;z=np.ascontiguousarray(z.real)
    rho,Q=invariants(z,L);mean=np.r_[d['light_mean_q'],d['light_mean_p']]
    return z,dict(N=N,L=L,seed=seed,IR_width=width,nonzero_modes=count,nQ=Q,rho=rho,
        q_RMS_over_mean=np.sqrt(np.mean((z[0]-mean[0])**2+(z[1]-mean[1])**2))/np.linalg.norm(mean[:2]),
        p_RMS_over_mean=np.sqrt(np.mean((z[3]-mean[3])**2+(z[4]-mean[4])**2))/np.linalg.norm(mean[3:5]),
        state_rule='Explicit positive widths for unstable initial IR modes; no ground-state claim or charge/energy retuning.')

if __name__=='__main__':prepare()
