"""Necessary numerical checks on the new potential, state and reused engine."""
from en_common import *
from run import run
from initial import initialize

def final_arrays(row):
    _,_,_,raw=snapshot(row['snapshots'][-1]);return np.array(raw[:3]),np.array(raw[3:])

def main():
    groups={};rng=np.random.default_rng(271026);N=6;L=5.
    q=rng.normal(size=(3,N,N,N));p=rng.normal(size=q.shape)
    force=lap(q,L)-gradient(q);ward=float(abs(np.sum(q[0]*force[1]-q[1]*force[0]))/np.sum(q*q))
    def rhs(t,y):
        qq=y[:q.size].reshape(q.shape);pp=y[q.size:].reshape(q.shape)
        return np.r_[pp.ravel(),(lap(qq,L)-gradient(qq)).ravel()]
    sol=solve_ivp(rhs,(0,.03),np.r_[q.ravel(),p.ravel()],method='DOP853',rtol=2e-12,atol=2e-14)
    row=run('verify_random',N=N,L=L,end_tau=.03,sample=.01,hmax=.0005,threads=1,arrays=(q,p),clock=dict(tau0=0.,a0=1.,C=0.),snap=1.)
    qc,pc=final_arrays(row);actual=np.r_[qc.ravel(),pc.ravel()]
    err=float(np.max(abs(actual-sol.y[:,-1])))
    groups['genuine_3D_C_vs_FFT_DOP853_and_Ward']=dict(passed=ward<1e-12 and err<1e-8,trajectory_absolute=err,global_U1_Ward=ward)

    base=np.array([10.,3.,1.]);vel=np.array([-.1,3.,.2]);Tend=2.88
    def homogeneous(t,y):return np.r_[y[3:],-gradient(y[:3])]
    exact=solve_ivp(homogeneous,(0,Tend),np.r_[base,vel],method='DOP853',rtol=2e-13,atol=2e-14).y[:,-1]
    q=np.broadcast_to(base[:,None,None,None],(3,4,4,4)).copy();p=np.broadcast_to(vel[:,None,None,None],q.shape).copy();errors=[]
    for h in [.08,.04,.02]:
        r=run(f'verify_h{h:g}',N=4,L=4.,end_tau=Tend,sample=.24,hmax=h,threads=1,arrays=(q,p),clock=dict(tau0=0.,a0=1.,C=0.),snap=10.)
        a,b=final_arrays(r);errors.append(float(np.max(abs(np.r_[a[:,0,0,0],b[:,0,0,0]]-exact))))
    groups['homogeneous_fourth_order_refinement']=dict(passed=errors[-1]<1e-6 and errors[0]/errors[1]>10 and errors[1]/errors[2]>10,steps=[.08,.04,.02],state_errors=errors)

    q1,p1,r1=initialize(32,12.);q2,p2,r2=initialize(64,12.)
    shared=float(max(np.max(abs(q1-q2[:,::2,::2,::2])),np.max(abs(p1-p2[:,::2,::2,::2]))))
    groups['same_continuous_Gaussian_realization_across_grids']=dict(passed=shared<1e-10,shared_grid_absolute_error=shared,
        scope='Identical finite-band Fourier functions; numerical gradient energy varies with the lattice operator.')
    mc=read('MATCH_RESULTS.json')['rows'];sym=max(r['scalar_transfer_symplectic'] for r in mc);det=max(r['phi_determinant'] for r in mc)
    groups['full_scalar_and_phi_linear_preparation']=dict(passed=sym<1e-7 and det<1e-7,scalar_symplectic=sym,phi_determinant=det,
        all_four_fields_in_preparation=True,full_four_field_nonlinear_evolution=False)
    save('ENGINE_CHECKS.json',dict(groups=groups));print((ROOT/'ENGINE_CHECKS.json').read_text(),flush=True)
    assert all(g['passed'] for g in groups.values())

if __name__=='__main__':main()
