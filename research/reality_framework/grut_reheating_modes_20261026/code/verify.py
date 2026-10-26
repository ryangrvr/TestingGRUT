"""Independent finite differences, Hamiltonian algebra and Radau controls."""
import numpy as np
from scipy.integrate import solve_ivp
from mode_common import ROOT,WORK,write
from modes import (matrices,generator,prepared,potential_force,geometric,rhs,
                   realvec,TCH,J,XI,ETA2,build,integrate,diagnostics)

def algebra():
    rng=np.random.default_rng(261026);hess=ward=legendre=canon=0.
    for _ in range(40):
        v,_=prepared(3.)
        v[2:8]+=rng.normal(size=6)*3;v[8:14]=rng.normal(size=6);v[15]=rng.uniform(0,2)
        v[16]=rng.uniform(.1,1);M,g,gp,chi,chip=matrices(v,3.)
        z,vel,V,F,f,fd,R,hd=geometric(v,3.,XI,ETA2);a=np.exp(v[15]);r=realvec(z)
        def force(q):
            zz=(q[::2]+1j*q[1::2])/np.sqrt(2);_,ff=potential_force(zz)
            out=realvec(ff);out[:2]+=XI*R*q[:2]
            return out
        fdmat=np.empty((6,6))
        for j in range(6):
            step=np.zeros(6);step[j]=1e-3
            fdmat[:,j]=(-force(r+2*step)+8*force(r+step)-8*force(r-step)+force(r-2*step))/(12e-3)
        expected=a*a*(fdmat-R/6*np.eye(6));hess=max(hess,float(np.max(abs(M-expected))/max(1,np.max(abs(M)))))
        gf=a**3*TCH@realvec(F)-a*a*R/6*g;gf[:2]+=a*a*XI*R*g[:2]
        ward=max(ward,float(np.linalg.norm(M@g-gf)/max(1,np.linalg.norm(gf))))
        k=.003;e=1e-7;G=generator(v,3.,k,e)
        canon=max(canon,float(np.max(abs(G.T@J+J@G))))
        q=rng.normal(size=6);p=rng.normal(size=6)
        D=k*k+e*e*(g@g);bar=g@p-gp@q;u=p+e*e/k**2*g*bar
        Om=M+k*k*np.eye(6);current=g@u-gp@q
        L=.5*(u@u-q@Om@q)-.5*e*e/D*current*current
        H=.5*(p@p+q@Om@q)+.5*e*e/k**2*bar*bar
        legendre=max(legendre,float(abs(p@u-L-H)/max(1,abs(H))))
    return dict(passed=hess<1e-7 and ward<1e-9 and legendre<1e-11 and canon<1e-11,
                samples=40,Hessian_relative=hess,Ward_relative=ward,Legendre_relative=legendre,canonical_generator_absolute=canon)

def independent(exe):
    rng=np.random.default_rng(141026);rows=[]
    for b,gam,k in [(3.,.5,2.3112763174514255),(3.,.5,1e-4),(10.,.1,2.3112763174514255)]:
        v,_=prepared(b);seed=rng.normal(size=(12,2))
        def fun(t,y):
            out=np.empty_like(y);out[:19]=rhs(t,y[:19],b,gam,XI,ETA2)
            out[19:]=(generator(y[:19],b,k)@y[19:].reshape(12,2)).ravel()
            return out
        sol=solve_ivp(fun,(0,3),np.r_[v,seed.ravel()],method='Radau',rtol=2e-9,atol=2e-12)
        assert sol.success,sol.message
        bg,tt,ctl=integrate(exe,f'radau_b{b:g}_k{k:g}',b,gam,[k],[0,3],rtol=2e-12)
        target=tt[-1,0]@seed;actual=sol.y[19:,-1].reshape(12,2)
        err=float(np.max(abs(actual-target)/(1+abs(target))))
        rows.append(dict(beta=b,Gamma=gam,k=k,scaled_state_error=err,nfev=sol.nfev))
    return dict(passed=max(r['scaled_state_error'] for r in rows)<1e-6,rows=rows)

def main():
    exe=build();groups={'Hessian_Ward_and_A0_Hamiltonian':algebra()}
    groups['independent_Radau_modes']=independent(exe)
    data=np.load(ROOT/'MODE_PROFILES.npz');times=data['times'];ks=data['ks']
    tags=[t for t in data.files if t.startswith('b') and not t.endswith(('_bg','_transfer'))]
    sym=[]
    import json
    rows=json.loads((ROOT/'MODE_RESULTS.json').read_text())['rows']
    groups['canonical_transfer_symplecticity']=dict(passed=max(r['controls']['symplectic_absolute_max'] for r in rows)<1e-6,
        maximum_absolute=max(r['controls']['symplectic_absolute_max'] for r in rows),
        maximum_scaled=max(r['controls']['symplectic_scaled_max'] for r in rows))
    # Tight primary method control at the resolved late amplification band.
    kk=[1e-5,.01,1.783559058474993,2.3112763174514255,40.]
    bg,t1,_=integrate(exe,'tight_control',3.,.5,kk,times,rtol=2e-12)
    bg,t0,_=integrate(exe,'normal_control',3.,.5,kk,times,rtol=2e-10)
    d1,_=diagnostics(3.,bg,t1,kk);d0,_=diagnostics(3.,bg,t0,kk)
    err=float(np.max(abs(d1[:,:,:2]-d0[:,:,:2])/(1+abs(d1[:,:,:2]))))
    bg,te,_=integrate(exe,'ungauged_control',3.,.5,kk,times,e=0.)
    de,_=diagnostics(3.,bg,te,kk)
    gauge_main=float(np.max(abs(de[:,1:,:2]-d1[:,1:,:2])/(1+abs(d1[:,1:,:2]))))
    # The k=0, e=0 global gauge tangent must track (g,g') exactly.
    bg,tg,_=integrate(exe,'global_Ward_control',3.,.5,[0.],[0,1,3,10],e=0.,rtol=2e-12)
    _,g,gp,_,_=matrices(bg[0],3.);vec=np.r_[g,gp];ward_error=0.
    for i,v in enumerate(bg):
        _,g,gp,_,_=matrices(v,3.);expected=np.r_[g,gp]
        ward_error=max(ward_error,float(np.linalg.norm(tg[i,0]@vec-expected)/np.linalg.norm(expected)))
    groups['tolerance_gauge_and_global_tangent_controls']=dict(passed=err<1e-6 and gauge_main<1e-3 and ward_error<1e-8,
        tighter_response_relative_error=err,gauge_effect_in_control_main_band=gauge_main,global_gauge_tangent_relative_error=ward_error,
        infrared_gauge_effect=float(np.max(abs(de[:,0,:2]-d1[:,0,:2])/(1+abs(d1[:,0,:2])))))
    write('VERIFICATION_RESULTS.json',dict(groups=groups,scope='Spectator six-real-scalar/longitudinal linear evolution; no metric/inflaton perturbations, nonlinear run or physical cosmic seed covariance.'))
    print(json.dumps(json.loads((ROOT/'VERIFICATION_RESULTS.json').read_text()),indent=2),flush=True)
    if not all(g['passed'] for g in groups.values()):raise SystemExit(1)

if __name__=='__main__':main()
