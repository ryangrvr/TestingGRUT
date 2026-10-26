"""Reconstruct and evolve the supplied heavy spot with the author's own engine."""
from pathlib import Path
import json, math, os, subprocess, hashlib
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.optimize import root
from build import build, WORK

ROOT=Path(__file__).resolve().parents[1]
LX,LY,LS,CH=.001,1.,1.,10.
MH,HH,RR,FRAC=30.,40.,.45,.5
MU=RR*MH**2/(math.sqrt(2)*HH)
KAP=MH*LX*math.sqrt(FRAC*4*(1-RR**2)/(12*CH))
HBASE=434867798206.6179
MRED,YQ=2.435e18,4.8294009076988175e-9


def native(x):
    if isinstance(x,np.generic):return x.item()
    if isinstance(x,np.ndarray):return x.tolist()
    raise TypeError(type(x).__name__)


def write(path,data):path.write_text(json.dumps(data,indent=2,allow_nan=False,default=native)+'\n')


def force(z,cH=0.):
    x,y,s=z
    return np.array([2*LX*abs(x)**2*x+3*KAP*y*np.conj(x)**2-cH*x,
                     MH*MH*y+2*LY*abs(y)**2*y+KAP*x**3+2*MU*s*np.conj(y),
                     2*LS*(abs(s)**2-HH*HH/2)*s+MU*y*y])


def prepared():
    direction=-1j
    def equations(v):
        x,yr,yi,s=v;ff=force([x,yr+1j*yi,direction*s],CH)
        return np.array([ff[0].real,ff[1].real,ff[1].imag,(ff[2]*np.conj(direction)).real])
    v=root(equations,[76.4,-2.11,-.94,28.293],tol=1e-12).x
    z=np.array([complex(v[0]),v[1]+1j*v[2],direction*v[3]])
    residual=float(np.linalg.norm(equations(v)))
    assert residual<1e-8
    eigcols=[]
    for j in range(4):
        d=np.zeros(4);d[j]=1e-5
        eigcols.append((equations(v+d)-equations(v-d))/(2e-5))
    hessian=np.array(eigcols).T
    eig=np.linalg.eigvalsh((hessian+hessian.T)/2)
    assert min(eig)>0
    return z,dict(values_H_units=[[a.real,a.imag] for a in z],conditional_gradient_norm=residual,
                  conditional_hessian_eigenvalues=eig.tolist(),
                  gauge_Higgs_norm_over_H=float(np.sqrt(2*np.dot([1,9,36],abs(z)**2))),
                  nonzero_X_torque=float(-6*np.imag(KAP*np.conj(z[1])*z[0]**3)),
                  physical_velocities='zero; conformal chi_prime=chi at a=1',
                  scope='Full conditional minimum at supplied X/Sigma phases; no phase-holder action derived.')


def execute(exe,z,tag,end,tol,samples):
    vals=[MH,HH,MU,KAP,end,tol,tol*1e-3]+[v for q in z for v in [q.real,q.imag]]+[samples]
    run=subprocess.run([str(exe)]+[format(v,'.17g') for v in vals],capture_output=True,text=True,check=True)
    raw=WORK/(tag+'.txt');raw.write_text(run.stdout)
    history=[];windows=[];finish=None;state=None
    for line in run.stdout.splitlines():
        parts=line.split()
        if parts[0]=='S':history.append([float(x) for x in parts[1:]])
        elif parts[0]=='W':
            vv=[float(x) for x in parts[1:]]
            if vv[1]-vv[0]>1e-10:
                windows.append(dict(start=vv[0],end=vv[1],span=vv[2],observables=vv[3:30]))
        elif parts[0]=='END':finish=dict(steps=int(parts[1]),rejects=int(parts[2]),max_gauss=float(parts[3]),max_work=float(parts[4]))
        elif parts[0]=='F':state=[float(v) for v in parts[1:]]
    assert finish is not None
    aa=np.array(history)
    np.savez_compressed(WORK/(tag+'.npz'),samples=aa)
    return dict(tag=tag,a_end=end,rtol=tol,atol=tol*1e-3,sample_count=samples,
                windows=windows,finish=finish,final_state=state,
                raw_sha256=hashlib.sha256(raw.read_bytes()).hexdigest()),aa


def entropy(H):
    T=(90*MRED*MRED*H*H/(math.pi**2*106.75))**.25
    return 2*math.pi**2/45*106.75*T**3


def main():
    WORK.mkdir(exist_ok=True,parents=True)
    z,prep=prepared();exe=build()
    specs=[('central',40.,2e-12,20000),('tight',40.,2e-13,5000),('long',60.,2e-12,10000),
           ('window100',100.,2e-12,5000)]
    with ThreadPoolExecutor(3) as pool:
        futures=[pool.submit(execute,exe,z,*s) for s in specs]
        runs=[f.result() for f in futures]
    for run,history in runs:
        if run['tag']=='window100':
            parts=[w for w in run['windows'] if w['start']>=30 and w['end']<=100]
            span=sum(w['span'] for w in parts)
            run['combined_window_30_100']=[sum(w['observables'][j]*w['span'] for w in parts)/span for j in range(27)]
    central,profile=runs[0]
    final_window=next(w for w in central['windows'] if w['start']==30 and w['end']==40)
    ob=final_window['observables'];nx,ex,eps,nyh,nyl,ns=ob[0],ob[4],ob[9],ob[20],ob[21],ob[22]
    raw_yield=abs(nx)*HBASE**3/entropy(HBASE)
    scale=(YQ/raw_yield)**(2/3);H=HBASE*scale
    np.savez_compressed(ROOT/'HIGGS_PROFILE.npz',a=profile[:,0],F_over_H=profile[:,14])
    physical=dict(H_GeV=H,M_Y_GeV=MH*H,h_GeV=HH*H,mu_GeV=MU*H,kappa=KAP,
                  lambda_X=LX,lambda_Y=LY,lambda_Sigma=LS,cH=CH,e=1e-7,
                  base_H_GeV=HBASE,uniform_scale=scale,base_yield_over_YQ=raw_yield/YQ)
    result=dict(parameters=physical,dimensionless=dict(M=MH,h=HH,mu=MU,kappa=KAP,r=RR),
                preparation=prep,late=dict(N_X=nx,E_X=ex,epsilon=eps,N_YH=nyh,N_YL=nyl,N_Sigma=ns),
                runs=[r for r,a in runs],
                source='Original Codex b448f25 release.c with runtime kappa and new windows only; conditional preparation solved afresh.',
                profile_columns='a,F_over_H (physical Higgs norm); sample table at array index14 includes leading a column.')
    write(ROOT/'BACKGROUND_RESULTS.json',result)
    print(json.dumps(dict(late=result['late'],parameters=physical,prep_residual=prep['conditional_gradient_norm'],
                          conservation=[r['finish'] for r,a in runs]),default=native),flush=True)


if __name__=='__main__':main()
