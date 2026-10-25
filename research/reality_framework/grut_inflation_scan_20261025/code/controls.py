"""Mathematical controls and independently written physical-time ODE checks."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from model import minimum,parameters,potential,force,initial_z
from run import evolve,write,ROOT


def independent_physical_time(H,cH,delta,a_end,z0=None,rtol=2e-12):
    p=parameters(H);point=minimum(H,cH)
    z=initial_z(point,delta) if z0 is None else z0
    def rhs(t,w):
        a=np.sqrt(1+2*t);positions=w[:3];vel=w[3:]
        return np.r_[vel,-3*vel/a**2-force(positions,p)]
    sol=solve_ivp(rhs,[0,(a_end*a_end-1)/2],np.r_[z,np.zeros(3,complex)],
                  method='DOP853',rtol=rtol,atol=2e-14)
    assert sol.success
    pos,vel=sol.y[:3,-1],sol.y[3:,-1]
    chi=a_end*pos;pi=pos+a_end*a_end*vel
    return np.ravel([[z.real,z.imag] for z in chi]+[[z.real,z.imag] for z in pi]),sol.nfev


def main():
    rng=np.random.default_rng(20261025);fd=[];ward=[]
    for _ in range(48):
        H=10**rng.uniform(11,13.67);p=parameters(H)
        z=rng.normal(size=3)+1j*rng.normal(size=3);z[0]*=20
        f=force(z,p)
        grad=[]
        for j in range(3):
            for direction in [1,1j]:
                dz=np.zeros(3,complex);dz[j]=direction*1e-5
                grad.append((potential(z+dz,p)-potential(z-dz,p))/2e-5)
        exact=np.ravel([[2*x.real,2*x.imag] for x in f])
        fd.append(float(np.linalg.norm(exact-grad)/max(np.linalg.norm(exact),1)))
        ward.append(float(abs(np.sum(np.array([1,3,6])*np.imag(np.conj(z)*f)))/max(np.linalg.norm(z)*np.linalg.norm(f),1)))
    minima=[]
    for H,cH in [(3e11,1.),(4.7e13,1.),(3e11/np.sqrt(10),10.)]:
        point=minimum(H,cH);p=parameters(H)
        ref=np.array(point['amplitudes_over_H']);scale=max(1.,abs(point['V_inflation_over_H4']))
        best=[]
        for j in range(8):
            x0=ref*np.exp(rng.normal(0,.6,3))
            def v(a):return potential(np.array([a[0],-a[1],-a[2]],complex),p,cH)/scale
            def g(a):
                f=force(np.array([a[0],-a[1],-a[2]],complex),p,cH)
                return np.array([2*f[0].real,-2*f[1].real,-2*f[2].real])/scale
            r=minimize(v,x0,jac=g,bounds=[(1e-9,None)]*3,method='L-BFGS-B',
                       options=dict(ftol=1e-15,gtol=1e-12,maxiter=5000))
            best.append(float(abs(r.fun*scale-point['V_inflation_over_H4'])/scale))
        minima.append(dict(H=H,cH=cH,max_8_start_minimum_energy_relative_shift=max(best)))
    short=[]
    for H,cH,delta,a_end in [(3e11,1.,np.pi/2,20.),(4.7e13,1.,.3,30.),(1e12,10.,.2,15.)]:
        r=evolve(H,cH,delta,end=a_end,rtol=2e-12,atol=2e-14)
        other,nfev=independent_physical_time(H,cH,delta,a_end)
        state=np.array(r['final_canonical_state'])
        err=float(np.linalg.norm(state-other)/np.linalg.norm(other))
        short.append(dict(H=H,cH=cH,delta=delta,a_end=a_end,physical_time_vs_C_state_relative_error=err,
                          independent_nfev=nfev,C_controls=r['controls']))
    H=4.7e13;pt=minimum(H,1.);z=initial_z(pt,0);z[0]*=np.exp(-1j*np.pi/24)
    old=evolve(H,1.,0,end=np.sqrt(201),rtol=2e-12,atol=2e-14,z0=z)
    q=old['history'][-1];oldcharge=[q['NX'],q['NY'],q['NS']]
    reference=np.array([-16.2613874661,6.39590444861,-.4877209799485])
    golden=max(abs(np.array(oldcharge)-reference))
    aligned=evolve(H,1.,0,end=300,rtol=2e-12,atol=2e-14)
    plus=evolve(3e11,1.,.3,end=100,rtol=2e-12,atol=2e-14)
    minus=evolve(3e11,1.,-.3,end=100,rtol=2e-12,atol=2e-14)
    cp=max(abs(plus['history'][-1][k]+minus['history'][-1][k]) for k in ['NX','NY','NS'])
    seed=evolve(H,1.,1e-30,end=200,rtol=2e-12,atol=1e-38,ns=1000)
    hist=seed['history']
    env1=max(abs(x['NX']) for x in hist if 30<x['a']<40)
    env2=max(abs(x['NX']) for x in hist if 140<x['a']<160)
    out=dict(author='Codex, same-session numerical controls; not outside replication.',
             finite_difference_force_max_relative_error=max(fd),gauge_Ward_max_relative_error=max(ward),
             multiple_start_minima=minima,physical_time_cross_checks=short,
             old_short_release=dict(charge=oldcharge,reference=reference.tolist(),max_abs_difference=golden),
             CP_conjugate_charge_sum_max=cp,
             exact_aligned_max_abs_charge=max(abs(x[k]) for x in aligned['history'] for k in ['NX','NY','NS']),
             linear_seed=dict(delta=1e-30,atol=1e-38,window_30_40_NX_envelope=env1,
                              window_140_160_NX_envelope=env2,envelope_growth_ratio=env2/env1,
                              diagnostic='Charge envelope of a homogeneous infinitesimal Sigma-phase seed, not a finite-k spectrum.'))
    write(ROOT/'CONTROL_RESULTS.json',out)
    write(ROOT/'CONTROL_HISTORIES.json',dict(linear_seed=seed,aligned=aligned,old_short_release=old))
    assert max(fd)<1e-7 and max(ward)<1e-13
    assert max(x['physical_time_vs_C_state_relative_error'] for x in short)<1e-7
    assert golden<1e-7 and cp<1e-9 and out['exact_aligned_max_abs_charge']==0
    assert env2/env1>1e8
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
