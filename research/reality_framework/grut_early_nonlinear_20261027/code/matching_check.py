"""Measure the scalar-limit matching error before nonlinear breakup."""
from en_common import *
import neutral as N

def main():
    exe=M.build();ne=exe.parent/'neutral_en_check'
    subprocess.run(['gcc','-O3','-std=c11','-DNEUTRAL_SECTOR','-I',str(exe.parent),str(T.MODE/'code/modes.c'),'-lm','-o',str(ne)],check=True)
    times=np.r_[0.,np.arange(20.,50.00001,.1)];ks=np.array([2.13,7.1])
    bg,tt,ctl=M.integrate(exe,'en_reduced_check',3.,.5,ks,times,e=0.,rtol=2e-12)
    pb,pt,pc=N.run(ne,dict(tag='en_reduced_check',beta=3.,Gamma_over_H=.5),ks,times,rtol=2e-12,suffix='_matched')
    _,_,_,chi,chip=M.matrices(bg[1],3.);a=np.exp(bg[1,15]);C=a*a*bg[1,16]
    seed=T.phase_seed(3.);full=np.einsum('tkij,j->tki',tt,seed)
    state=np.r_[chi[:2],chip[:2]];starts=[]
    for j in range(2):starts.extend(np.r_[full[1,j,:2],full[1,j,6:8],pt[1,j].ravel()])
    init=np.r_[state,starts]
    def rhs(s,y):
        q=y[:2];p=y[2:4];r2=q@q;out=np.empty_like(y);out[:4]=np.r_[p,-LX*r2*q]
        for j,k in enumerate(ks):
            off=4+8*j;vec=y[off:off+4];phi=y[off+4:off+8].reshape(2,2)
            mass=(k*k+LX*r2)*np.eye(2)+2*LX*np.outer(q,q)
            out[off:off+4]=np.r_[vec[2:],-mass@vec[:2]]
            out[off+4:off+8]=np.vstack([phi[1],-(k*k+G2*r2)*phi[0]]).ravel()
        return out
    tau=times[1:];elapsed=2*(tau-20)/(a+np.sqrt(a*a+2*C*(tau-20)))
    sol=solve_ivp(rhs,(0,elapsed[-1]),init,t_eval=elapsed,method='DOP853',rtol=2e-12,atol=2e-13);assert sol.success
    actual=sol.y.T;fq=[]
    for v in bg[1:]:fq.append(M.matrices(v,3.)[3][:2])
    fq=np.array(fq);orbit=apo(chi[:2],chip[:2]);berr=float(np.max(np.linalg.norm(actual[:,:2]-fq,axis=1))/orbit)
    rows=[]
    for j,k in enumerate(ks):
        off=4+8*j;target=full[1:,j,:2];approx=actual[:,off:off+2]
        err=float(np.max(np.linalg.norm(approx-target,axis=1))/np.max(np.linalg.norm(target,axis=1)))
        phi=actual[:,off+4:off+8].reshape(-1,2,2);pherr=float(np.max(abs(phi-pt[1:,j]))/np.max(abs(pt[1:,j])))
        rows.append(dict(k=k,phase_tangent_envelope_relative_error=err,neutral_transfer_relative_error=pherr))
    _,_,pot,F,_,_,R,_=M.geometric(bg[1],3.,M.XI,M.ETA2)
    x=bg[1,2]+1j*bg[1,3]
    fullforce=M.realvec(F)[:2]*a**3+(M.XI-1/6)*a*a*R*chi[:2]
    fbare=LX*(chi[:2]@chi[:2])*chi[:2]
    ferror=float(np.linalg.norm(fullforce-fbare)/np.linalg.norm(fbare))
    save('MATCHING_CHECKS.json',dict(tau_match=20.,through_tau=50.,background_orbit_scaled_error=berr,
        initial_charged_force_relative_mismatch=ferror,momenta=rows,
        full_curvature_R_at_match=float(R),fixed_radiation_C=C,
        scope='Measured homogeneous and linear approximation errors. Not a nonlinear heavy/gauge error bound, physical matching certificate or quantum renormalization.'))
    print((ROOT/'MATCHING_CHECKS.json').read_text(),flush=True)

if __name__=='__main__':main()
