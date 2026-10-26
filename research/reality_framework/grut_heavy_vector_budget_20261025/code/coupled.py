"""Initial common quantum mode through six-scalar Coulomb-gauge evolution.

Only one initially independent mode is propagated. Late longitudinal
matching is checked at four times; no full initial covariance is claimed.
"""
import json,math,subprocess
from multiprocessing import Pool
import numpy as np
from scipy.integrate import solve_ivp
from background import ROOT,WORK,MH,HH,MU,KAP,prepared,write
from build import BUILD,build
from vector_modes import load_profile


def inflation(k,e,rtol=2e-11):
    _,prep=prepared();mi=e*prep['gauge_Higgs_norm_over_H'];a0=k/2000.
    am2=(a0*mi)**2;z2=am2/(k*k+am2);z=math.sqrt(z2)
    omega=math.sqrt(k*k+am2);u=np.exp(2000j)/math.sqrt(2*omega)
    initial=np.array([u/z,z*(-1j*omega-a0*k*k/(k*k+am2))*u])
    def rhs(l,state):
        a=math.exp(l);am2=(a*mi)**2;inertia=am2/(k*k+am2)
        return [state[1]/(a*inertia),-am2*state[0]/a]
    result=solve_ivp(rhs,[math.log(a0),0],initial,method='DOP853',rtol=rtol,
                     atol=[abs(initial[0])*1e-14,abs(initial[1])*1e-14])
    assert result.success
    B,Pi=result.y[:,-1];zz=mi*mi/(k*k+mi*mi)
    return B,Pi/(zz*B)


def late(k,e,a0,B,Bprime):
    aa,ff,interpolator,tail=load_profile()
    def amass(a):
        if a<=aa[-1]:return a*a*e*e*float(interpolator(math.log(a)))
        return (6*e*HH*a)**2+e*e*tail
    z2=amass(a0)/(k*k+amass(a0));Pi=z2*Bprime
    def rhs(l,state):
        a=math.exp(l);mass=amass(a);inertia=mass/(k*k+mass)
        return [a*state[1]/inertia,-a*mass*state[0]]
    end=20/math.sqrt(6*e*HH)
    initial=np.array([B,Pi]);state=initial
    boundaries=[math.log(a0)]+([math.log(aa[-1])] if a0<aa[-1] else [])+[math.log(end)]
    for left,right in zip(boundaries[:-1],boundaries[1:]):
        if right-left<1e-14:continue
        result=solve_ivp(rhs,[left,right],state,method='DOP853',rtol=2e-10,
                         atol=[abs(initial[0])*1e-13,max(abs(initial[1])*1e-13,1e-30)])
        assert result.success;state=result.y[:,-1]
    Bf,Pf=state;am2=amass(end);zz=am2/(k*k+am2)
    return float((abs(Pf)**2/zz+am2*abs(Bf)**2)/(2*math.sqrt(k*k+am2))-.5)


def build_coupled():
    build();exe=BUILD/'coupled'
    subprocess.run(['gcc','-O3','-std=c11','-Wall','-Wextra','-I',str(BUILD),
                    str(ROOT/'code/coupled.c'),'-lm','-o',str(exe)],check=True)
    return exe


def run_batch(exe,xs,tag,rtol=2e-11,matches=(10.,20.,30.,40.)):
    e=1e-7;z,_=prepared();ks=[x*math.sqrt(6*e*HH) for x in xs]
    inf=[inflation(k,e) for k in ks]
    args=[MH,HH,MU,KAP,40.,rtol,rtol*.001]+[v for c in z for v in [c.real,c.imag]]+[200,e]
    for k,(B,ratio) in zip(ks,inf):args += [k,ratio.real,ratio.imag]
    run=subprocess.run([str(exe)]+[format(v,'.17g') for v in args],capture_output=True,text=True,check=True)
    (WORK/(tag+'.txt')).write_text(run.stdout)
    projections=[];finish=None
    for line in run.stdout.splitlines():
        parts=line.split()
        if parts[0]=='P':
            values=[float(v) for v in parts[1:]];a=values[0]
            if not any(abs(a-t)<1e-10 for t in matches):continue
            for j,(k,x,(B0,ratio)) in enumerate(zip(ks,xs,inf)):
                theta=complex(*values[1+4*j:3+4*j]);tp=complex(*values[3+4*j:5+4*j])
                n=late(k,e,a,B0*theta,B0*tp)
                projections.append(dict(a_match=a,k_over_kstar=x,k_over_H=k,
                                        theta=[theta.real,theta.imag],theta_prime=[tp.real,tp.imag],
                                        late_common_mode_occupation=n,
                                        dN_dlnk=k**3*n/(2*math.pi**2)))
        elif parts[0]=='END':
            finish=dict(steps=int(parts[1]),rejects=int(parts[2]),max_gauss=float(parts[3]),max_work=float(parts[4]))
    return dict(tag=tag,xs=xs,rtol=rtol,finish=finish,projections=projections)


def band_job(args):return run_batch(*args)


def main():
    exe=build_coupled()
    control=run_batch(exe,[.3,1.,3.],'common_control')
    tight=run_batch(exe,[.3,1.,3.],'common_tight',2e-12)
    write(ROOT/'COUPLED_CONTROL.json',dict(rows=[control,tight],scope=__doc__))
    xx=np.geomspace(.2,5.,49);jobs=[]
    for start in range(0,len(xx),3):
        batch=xx[start:start+3].tolist()
        while len(batch)<3:batch.append(batch[-1])
        jobs.append((exe,batch,'common_band_%02d'%start,2e-11,(40.,)))
    with Pool(4) as pool:
        runs=[]
        for i,result in enumerate(pool.imap(band_job,jobs),1):
            runs.append(result);print('COUPLED_BAND_PROGRESS',i,len(jobs),flush=True)
    unique={}
    for run in runs:
        for row in run['projections']:unique[row['k_over_kstar']]=row
    points=[unique[x] for x in sorted(unique)]
    integral=float(np.trapezoid([p['dN_dlnk'] for p in points],np.log([p['k_over_H'] for p in points])))
    result=dict(nodes=len(points),N_A_common_band=integral,match_a=40.,points=points,runs=runs,
                scope='One initially independent common BD mode propagated with the full six-real-scalar Coulomb-gauge Hamiltonian. Late Proca matching; other initial physical modes not calculated.')
    write(ROOT/'COUPLED_BAND.json',result)
    print('COUPLED_BAND',integral,len(points),flush=True)


if __name__=='__main__':main()
