"""Fixed M1 parameters, continuum seeds, Hamiltonian and executable interface."""
from pathlib import Path
import os, json, hashlib, struct, subprocess, time
import numpy as np
from scipy.fft import rfftn, irfftn
from scipy.optimize import brentq
from scipy.interpolate import CubicHermiteSpline

ROOT=Path(__file__).resolve().parents[1]
RF=ROOT.parent
WORK=Path(os.environ.get('GRUT_M1_DYNAMICS_WORK','/workspace/scratch/0c6d6cf183c6/m1_dynamics_runs'))
WORK.mkdir(parents=True,exist_ok=True)
G=.1; G2=G*G; LX=.001; LP=.03; EPS=.002; V=10.
PARENT=RF/'grut_matched_candidate_20261027'

def save(name,data):
    def default(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    (ROOT/name).write_text(json.dumps(data,indent=2,default=default,allow_nan=False)+'\n')

def build():
    src=ROOT/'code/field3d.cpp'; key=hashlib.sha256(src.read_bytes()).hexdigest()[:16]
    exe=WORK/('field3d_'+key)
    if not exe.exists():subprocess.run(['g++','-std=c++17','-O3','-ffp-contract=off','-fopenmp',str(src),'-o',str(exe)],check=True)
    return exe

def potential(q):
    r2=q[0]**2+q[1]**2;s=q[2]
    return LP/4*(s*s-V*V)**2-EPS*V*(s-V)**3+G2/2*s*s*r2+LX/4*r2*r2

def potential_gradient(q):
    a,b,s=q;r2=a*a+b*b;f=G2*s*s+LX*r2
    return np.array([f*a,f*b,(LP*(s*s-V*V)+G2*r2)*s-3*EPS*V*(s-V)**2])

def lap(q,L):
    out=np.zeros_like(q);dx=L/q.shape[-1]
    for ax in (-3,-2,-1):
        out+=(-np.roll(q,2,ax)+16*np.roll(q,1,ax)-30*q+16*np.roll(q,-1,ax)-np.roll(q,-2,ax))/(12*dx*dx)
    return out

def local_energy(q,p,L):
    dx=L/q.shape[-1];ge=np.zeros(q.shape[-3:])
    for a in q:
        for ax in (0,1,2):
            ge+=.5/dx**2*((4/3)*(np.roll(a,-1,ax)-a)**2-(np.roll(a,-2,ax)-a)**2/12)
    return .5*np.sum(p*p,axis=0)+potential(q)+ge

def totals(q,p,L):
    dv=(L/q.shape[-1])**3
    return {'E':float(local_energy(q,p,L).sum()*dv),'Q':float((q[0]*p[1]-q[1]*p[0]).sum()*dv)}

def seed(N,L,number,lo=.18,hi=.55,count=32):
    rng=np.random.default_rng(number);maxn=int(np.ceil(hi*L/(2*np.pi)));modes=[];seen=set()
    while len(modes)<count:
        n=rng.integers(-maxn,maxn+1,3);k=2*np.pi*np.linalg.norm(n)/L
        if not lo<=k<=hi:continue
        if next((v for v in n if v),0)<0:n=-n
        key=tuple(n)
        if key in seen:continue
        seen.add(key);modes.append((n,float(rng.uniform(0,2*np.pi))))
    x=np.arange(N)*2*np.pi/N;f=np.zeros((N,N,N))
    for n,phase in modes:f+=np.cos(n[0]*x[:,None,None]+n[1]*x[None,:,None]+n[2]*x[None,None,:]+phase)
    f*=np.sqrt(2/count)
    return f,[{'n':n.tolist(),'phase':phase} for n,phase in modes]

def uniform(N,L,number=1,F=1.,amplitude=.01,ellipse=False):
    S=brentq(lambda s:LP/G2*(s*s-1)*s-3*EPS/G2*(s-1)**2+F*F*s,0,1)
    om=np.sqrt(S*S+LX/G2*F*F)
    a,am=seed(N,L,number);b,bm=seed(N,L,number+1000)
    c,cm=seed(N,L,number+2000,.1,.65)
    radius=V*F*(1+amplitude*a);theta=amplitude*b
    q=np.array([radius*np.cos(theta),radius*np.sin(theta),V*(S+amplitude*c)])
    p=np.array([-om*q[1],om*q[0],np.zeros_like(a)])
    ratio=(9+np.sqrt(77))/2 if ellipse else 1.
    if ellipse:
        # Canonical fixed-charge stretch. No physical particle interpretation.
        q[0]*=np.sqrt(ratio);q[1]/=np.sqrt(ratio)
        p[0]*=np.sqrt(ratio);p[1]/=np.sqrt(ratio)
    target=om*(V*F)**2*L**3
    p[:2]*=target/totals(q,p,L)['Q']
    params={'kind':'elliptic' if ellipse else 'cold','supplied_F':F,'S':S,'Omega':om,
            'seed':number,'relative_seed_rms':amplitude,'ellipse_axis_ratio':ratio,
            'harmonic_shape_action_reference':1.75 if ellipse else 0.,
            'modes':{'amplitude':am,'phase':bm,'neutral':cm}}
    return q,p,params

def object_state(N,L,perturbation=0.):
    a=np.loadtxt(PARENT/'PROFILE_M1.csv',delimiter=',',skiprows=1)
    F=CubicHermiteSpline(a[:,0],a[:,1],a[:,2]);S=CubicHermiteSpline(a[:,0],a[:,3],a[:,4])
    om=json.loads((PARENT/'DEFORMATION_RESULTS.json').read_text())['primary']['Omega']
    x=(np.arange(N)-N//2)*L/N
    xx,yy,zz=np.meshgrid(x,x,x,indexing='ij');r=np.sqrt(xx*xx+yy*yy+zz*zz)
    q=np.array([V*F(r),np.zeros_like(r),V*S(r)])
    p=np.array([np.zeros_like(r),om*q[0],np.zeros_like(r)])
    raw=totals(q,p,L)
    if perturbation:
        # Smooth l=0/l=2 deformation localized on the core. Same total Q.
        shape=np.exp(-(r/12)**2)*(1+(xx*xx-yy*yy)/(r*r+4))
        q[0]*=1+perturbation*shape
        q[2]+=perturbation*V*(1-S(r))*shape
    p[1]=om*q[0]
    p[:2]*=200000/totals(q,p,L)['Q']
    return q,p,{'kind':'object','perturbation_amplitude':perturbation,'Omega_reference':om,
                'reference_Q':200000,'unadjusted_interpolation_totals':raw,
                'periodic_tail_max_abs_X':float(np.max(abs(q[0][[0,-1],:,:])))}

def write_initial(path,q,p,L):
    with open(path,'wb') as f:
        f.write(b'GRUTM1D\n');f.write(struct.pack('<id5d',q.shape[-1],L,V,G2,LX,LP,EPS))
        np.asarray(q,dtype='<f8').tofile(f);np.asarray(p,dtype='<f8').tofile(f)

def snapshot(path):
    with open(path,'rb') as f:N,L,t,bits=struct.unpack('<iddi',f.read(24))
    arr=np.memmap(path,dtype='<f8' if bits==64 else '<f4',mode='r',offset=24,shape=(6,N,N,N))
    return N,L,t,arr

def run(tag,N,L,dt,T,snap=20.,sample=.5,threads=2,number=1,ellipse=False,kind='cold',perturbation=0.,F=1.,amplitude=.01):
    directory=WORK/tag;directory.mkdir(exist_ok=True);prefix=directory/'state';initial=directory/'initial.bin'
    args=locals().copy();args={k:v for k,v in args.items() if k not in ('directory','prefix','initial')}
    resultpath=directory/'result.json'
    codehash=hashlib.sha256((ROOT/'code/field3d.cpp').read_bytes()+(ROOT/'code/common.py').read_bytes()).hexdigest()
    if resultpath.exists():
        old=json.loads(resultpath.read_text())
        if old['code_sha256']==codehash and old['args']==args:return old
        raise RuntimeError('Refuse stale run collision: '+tag)
    q,p,params=object_state(N,L,perturbation) if kind=='object' else uniform(N,L,number,F,amplitude,ellipse)
    params.update(totals(q,p,L));params['initial_E_over_Q']=params['E']/params['Q']
    write_initial(initial,q,p,L)
    params['initial_sha256']=hashlib.sha256(initial.read_bytes()).hexdigest()
    del q,p
    (directory/'parameters.json').write_text(json.dumps(params,indent=2)+'\n')
    start=time.time()
    cmd=[str(build()),str(initial),str(prefix),str(dt),str(T),str(sample),str(snap),'4',str(threads),'64']
    with open(directory/'run.log','w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    data=np.genfromtxt(str(prefix)+'.tsv',names=True,delimiter='\t')
    result={'args':args,'code_sha256':codehash,'parameters':params,'seconds':time.time()-start,
            'energy_relative_drift':float(np.max(abs(data['rho']/data['rho'][0]-1))),
            'charge_relative_drift':float(np.max(abs(data['nQ']/data['nQ'][0]-1))),
            'final':{n:float(data[n][-1]) for n in data.dtype.names},
            'snapshots':[str(f) for f in sorted(directory.glob('state_s*.bin'))]}
    resultpath.write_text(json.dumps(result,indent=2)+'\n')
    return result

