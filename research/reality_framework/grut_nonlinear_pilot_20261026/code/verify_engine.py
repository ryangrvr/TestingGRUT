"""Check the implemented 3D Hamiltonian before nonlinear interpretation."""
import numpy as np
from scipy.integrate import solve_ivp
from pilot import *
sys.path.insert(0,str(ROOT.parent/'grut_formation_onset_20261026'/'code'))
from elliptic import orbit

def raw_lap(f,L):
    d=L/f.shape[-1];out=-7.5*f
    for ax in range(-3,0):out+=(4/3)*(np.roll(f,1,ax)+np.roll(f,-1,ax))-(np.roll(f,2,ax)+np.roll(f,-2,ax))/12
    return out/d**2

def projection(data,cos):return 2*np.mean(np.asarray(data)*cos,axis=(-3,-2,-1))

def main():
    groups={};rng=np.random.default_rng(3801);q=rng.normal(size=(3,12,12,12));L=18.
    a=raw_lap(q,L);b=lap(q,L);err=float(np.max(abs(a-b))/(1+np.max(abs(a))))
    ward=abs(np.mean(q[0]*a[1]-q[1]*a[0]))
    groups['spatial_operator_and_U1']=dict(passed=err<1e-12 and ward<1e-12,FFT_vs_stencil=err,charge_Ward=ward)
    # Exercise all three C++ neighbour directions on a genuinely random3D
    # state, using an independent FFT evaluation of the same discrete H.
    N=8;ql=rng.normal(size=(3,N,N,N))*.1;ql[2]+=C.NU;pl=rng.normal(size=ql.shape)*.05
    def lattice_rhs(t,y):
        qq=y[:3*N**3].reshape(3,N,N,N);pp=y[3*N**3:].reshape(qq.shape)
        return np.r_[pp.ravel(),(lap(qq,12.)-C.gradient(qq)).ravel()]
    iv=np.r_[ql.ravel(),pl.ravel()]
    ref=solve_ivp(lattice_rhs,(0,.02),iv,method='DOP853',rtol=2e-13,atol=2e-15).y[:,-1]
    cr=run('verify_random3D',N=N,L=12.,T=.02,dt=.002,snap=.02,sample=.02,bits=64,threads=2,initial_arrays=(ql,pl))
    _,_,_,actual=snapshot(cr['snapshots'][-1]);cerr=float(np.max(abs(np.asarray(actual).ravel()-ref)))
    groups['spatial_operator_and_U1']['C_vs_FFT_time_reference_absolute']=cerr
    groups['spatial_operator_and_U1']['passed']=bool(groups['spatial_operator_and_U1']['passed'] and cerr<1e-10)
    init=[]
    for N in [24,32,48]:
        x,p,d=initialize(N,L=48.,seed=1,neutral='mixed');init.append([abs(d['energy']/(C.EX+C.V0)-1),abs(d['nQ']/C.Q0-1)])
    groups['matched_initial_state']=dict(passed=max(max(z) for z in init)<1e-12,errors=init)

    bg,t0,P,angle,diag=orbit();v=diag['initial'];initial=np.array([v[0],v[1],C.NU*1e-8,v[2],v[3],0.])
    def force(t,y):return np.r_[y[3:],-C.gradient(y[:3])]
    ref=solve_ivp(force,(0,20),initial,method='DOP853',rtol=2e-13,atol=2e-15).y[:,-1]
    hom=[]
    for dt in [.08,.04,.02]:
        N=8;q=np.broadcast_to(initial[:3,None,None,None],(3,N,N,N)).copy();p=np.broadcast_to(initial[3:,None,None,None],q.shape).copy()
        r=run('verify_hom_'+str(dt),N=N,L=16.,T=20.,dt=dt,snap=20.,sample=1.,bits=64,threads=2,initial_arrays=(q,p))
        _,_,_,a=snapshot(r['snapshots'][-1]);computed=np.mean(a,axis=(1,2,3));e=np.linalg.norm(computed-ref)/np.linalg.norm(ref)
        hom.append(dict(dt=dt,state_relative=float(e),energy=r['energy_relative_drift'],charge=r['charge_relative_drift']))
    groups['homogeneous_time_convergence']=dict(passed=hom[-1]['state_relative']<1e-6 and hom[-1]['state_relative']<hom[-2]['state_relative']/8,
        rows=hom)

    N=16;L=2*np.pi/.341484601885;k2=float(symbols(N,L)[1,0,0]);cos=np.broadcast_to(np.cos(np.arange(N)*2*np.pi/N)[:,None,None],(N,N,N));keff=np.sqrt(k2)
    linear=[]
    for kind in ['X','phi']:
        amp=1e-6 if kind=='X' else 1e-9
        q0=np.array([v[0],v[1],0.]);p0=np.array([v[2],v[3],0.]);q=np.broadcast_to(q0[:,None,None,None],(3,N,N,N)).copy();p=np.broadcast_to(p0[:,None,None,None],q.shape).copy()
        j=0 if kind=='X' else 2;q[j]+=amp*cos
        def f(t,y):
            qq=bg.sol(t)[:2];H=C.hessian([*qq,0.])+k2*np.eye(3);return np.r_[y[3:],-H@y[:3]]
        perturb=np.zeros(6);perturb[j]=amp
        refsol=solve_ivp(f,(0,10),perturb,method='DOP853',rtol=2e-13,atol=1e-23,dense_output=True)
        for dt in [.04,.02]:
            r=run('verify_linear_'+kind+'_'+str(dt),N=N,L=L,T=10.,dt=dt,snap=1.,sample=1.,bits=64,threads=2,initial_arrays=(q,p))
            times=[];computed=[]
            for path in r['snapshots']:
                _,_,t,a=snapshot(path);times.append(t);computed.append(projection(a,cos))
            pred=refsol.sol(times).T;comp=np.array(computed);err=float(np.max(abs(comp-pred))/np.max(abs(pred)))
            linear.append(dict(kind=kind,dt=dt,relative=err,k_effective=keff,energy=r['energy_relative_drift'],charge=r['charge_relative_drift']))
    groups['charged_and_neutral_linear_modes']=dict(passed=all(r['relative']<2e-5 for r in linear if r['dt']==.02),rows=linear,
        scope='Plane modes in a genuine 3D grid, compared with continuum time ODE at the exact discrete momentum.')

    save('ENGINE_CHECKS.json',dict(groups=groups,physical_formation_derived=False,physical_f_sym=None))
    print(groups,flush=True)
    if not all(r['passed'] for r in groups.values()):raise SystemExit(1)

if __name__=='__main__':main()
