"""Propagate a specified finite-band Gaussian through all scalar linear modes."""
from en_common import *
from concurrent.futures import ThreadPoolExecutor,as_completed
import neutral as N

def shell_keys(L,cutoff=10.):
    n=int(np.ceil(cutoff*L/(2*np.pi)));a=np.arange(-n,n+1)
    s=(a[:,None,None]**2+a[None,:,None]**2+a[None,None,:]**2).ravel();k=2*np.pi/L*np.sqrt(s)
    return np.unique(s[(k>=1.8)&(k<=cutoff)])

def main():
    exe=M.build();neutral=exe.parent/'neutral_match'
    subprocess.run(['gcc','-O3','-std=c11','-DNEUTRAL_SECTOR','-I',str(exe.parent),str(T.MODE/'code/modes.c'),'-lm','-o',str(neutral)],check=True)
    for L,cut in [(12.,10.),(16.,8.)]:
        keys=shell_keys(L,cut);ks=2*np.pi/L*np.sqrt(keys);times=np.array([0.,20.]);parts={}
        def batch(start):
            kk=ks[start:start+20]
            bg,TT,control=M.integrate(exe,f'en_match_L{L:g}_{start}',3.,.5,kk,times,e=0.,rtol=2e-12)
            pb,PT,pc=N.run(neutral,dict(tag=f'en_L{L:g}_{start}',beta=3.,Gamma_over_H=.5),kk,times,rtol=2e-12,suffix='_match')
            assert np.max(abs(bg-pb)/(1+abs(bg)))<1e-6
            return start,bg,TT[-1],PT[-1],control,pc
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures=[pool.submit(batch,start) for start in range(0,len(ks),20)]
            for f in as_completed(futures):
                start,bg,TT,PT,control,pc=f.result();parts[start]=(TT,PT,control,pc)
        TT=np.concatenate([parts[i][0] for i in sorted(parts)]);PT=np.concatenate([parts[i][1] for i in sorted(parts)])
        v0,_=M.prepared(3.);M0,_,_,_,_=M.matrices(v0,3.);lam,U=np.linalg.eigh(M0)
        assert np.min(ks[:,None]**2+lam)>0
        w=np.sqrt(ks[:,None]**2+lam);S=np.zeros((len(ks),12,12))
        S[:,:6,:6]=U[None,:,:]/np.sqrt(2*w[:,None,:]);S[:,6:,6:]=U[None,:,:]*np.sqrt(w[:,None,:]/2)
        evolved=TT@S
        wphi=np.sqrt(ks**2+N.mass(v0,3.));SP=np.zeros((len(ks),2,2));SP[:,0,0]=1/np.sqrt(2*wphi);SP[:,1,1]=np.sqrt(wphi/2)
        ep=PT@SP
        _,_,_,chi,chip=M.matrices(bg[-1],3.);a=np.exp(bg[-1,15]);C=a*a*bg[-1,16]
        np.savez_compressed(ROOT/f'MATCH_L{L:g}.npz',keys=keys,ks=ks,evolved_scalar=evolved,evolved_phi=ep,
            background=bg,initial_mass_eigenvalues=lam,initial_mass_eigenvectors=U,
            light_mean_q=np.r_[chi[:2],0.],light_mean_p=np.r_[chip[:2],0.],a_match=a,C_match=C,tau_match=20.)
        J=np.block([[np.zeros((6,6)),np.eye(6)],[-np.eye(6),np.zeros((6,6))]])
        sym=max(np.max(abs(t.T@J@t-J)) for t in TT);det=max(abs(np.linalg.det(PT)-1))
        record=dict(L=L,seed_upper_limit=cut,mode_shells=len(ks),tau_match=20.,a_match=a,C_match=C,
            scalar_transfer_symplectic=float(sym),phi_determinant=float(det),
            initial_minimum_frequency2=float(np.min(ks[:,None]**2+lam)),
            batch_controls=[dict(start=i,scalar=parts[i][2],neutral=parts[i][3]) for i in sorted(parts)],
            scope='Full X/Y/Sigma and phi linear propagation in the e0 scalar limit to20; only the light marginal is nonlinearly continued. Finite seed band, conditional instantaneous conformal Gaussian, not a selected full quantum vacuum.')
        current=read('MATCH_RESULTS.json') if (ROOT/'MATCH_RESULTS.json').exists() else dict(rows=[])
        current['rows']=[r for r in current['rows'] if r['L']!=L]+[record];save('MATCH_RESULTS.json',current)
        print('MATCH',L,len(ks),'a',a,'C',C,'sym',sym,'phi',det,flush=True)

if __name__=='__main__':main()
