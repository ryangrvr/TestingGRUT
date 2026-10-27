"""Reproducible Fourier seed functions, matched data and C++ orchestration."""
from pathlib import Path
import os,sys,json,struct,subprocess,hashlib,time
import numpy as np
from scipy.fft import rfftn,irfftn
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[2]
sys.path.insert(0,str(ROOT.parent/'grut_formation_onset_20261026'/'code'))
import common as C
WORK=Path(os.environ.get('GRUT_3D_WORKDIR','/tmp/grut_nonlinear_pilot_20261026'));WORK.mkdir(exist_ok=True,parents=True)

def save(name,data):
    def cv(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,Path):return str(v)
        raise TypeError(type(v).__name__)
    (ROOT/name).write_text(json.dumps(data,indent=2,default=cv,allow_nan=False)+'\n')

def build():
    src=ROOT/'code/field3d.cpp';key=hashlib.sha256(src.read_bytes()).hexdigest()[:16];exe=WORK/('field3d_'+key)
    if not exe.exists():subprocess.run(['g++','-std=c++17','-O3','-ffp-contract=off','-fopenmp',str(src),'-o',str(exe)],check=True)
    return exe

def symbols(N,L):
    dx=L/N;f=np.fft.fftfreq(N)*2*np.pi;r=np.fft.rfftfreq(N)*2*np.pi
    def one(u):return 4*np.sin(u/2)**2*(1+np.sin(u/2)**2/3)/dx**2
    return one(f)[:,None,None]+one(f)[None,:,None]+one(r)[None,None,:]

def lap(q,L):
    k2=symbols(q.shape[-1],L)
    return irfftn(-k2*rfftn(q,axes=(-3,-2,-1)),s=q.shape[-3:],axes=(-3,-2,-1))

def energy(q,p,L):
    r2=q[0]**2+q[1]**2
    pot=C.LP/4*(q[2]**2-C.NU**2)**2+C.G**2/2*q[2]**2*r2+C.LX/4*r2*r2
    grad=-.5*np.sum(q*lap(q,L),axis=0)
    return float(np.mean(.5*np.sum(p*p,axis=0)+pot+grad))

def field_modes(N,L,seed,kind,count=32):
    rng=np.random.default_rng(seed);lo,hi=(.18,.55) if kind!='phi' else (.10,.65)
    maxn=int(np.ceil(hi*L/(2*np.pi)));modes=[];seen=set()
    while len(modes)<count:
        n=rng.integers(-maxn,maxn+1,3);k=2*np.pi*np.linalg.norm(n)/L
        if not lo<=k<=hi:continue
        first=next((v for v in n if v),0)
        if first<0:n=-n
        key=tuple(n)
        if key in seen:continue
        seen.add(key);modes.append((n,float(rng.uniform(0,2*np.pi))))
    # Continuous Fourier functions are identical on all resolved grids.
    axis=np.arange(N)*2*np.pi/N;out=np.zeros((N,N,N))
    for n,phase in modes:
        out+=np.cos(n[0]*axis[:,None,None]+n[1]*axis[None,:,None]+n[2]*axis[None,None,:]+phase)
    out*=np.sqrt(2/count)
    return out,[dict(n=n.tolist(),phase=phase) for n,phase in modes]

def write_initial(path,q,p,L):
    N=q.shape[-1]
    with open(path,'wb') as f:
        f.write(b'GRUT3D1\n');f.write(struct.pack('<id4d',N,L,C.NU,C.G**2,C.LX,C.LP));np.asarray(q,dtype='<f8').tofile(f);np.asarray(p,dtype='<f8').tofile(f)

def initialize(N=64,L=48.,seed=1,neutral='mixed',amplitude=1e-6):
    amp,am=field_modes(N,L,seed,'amp');phase,pm=field_modes(N,L,seed+1000,'phase');ph,sm=field_modes(N,L,seed+2000,'phi')
    radius=C.RC*(1+amplitude*amp);theta=amplitude*phase
    q=np.array([radius*np.cos(theta),radius*np.sin(theta),C.NU*amplitude*(ph if neutral=='mixed' else 1+.1*ph)])
    p=np.array([-C.Q0/radius*np.sin(theta),C.Q0/radius*np.cos(theta),np.zeros_like(radius)])
    E_before=energy(q,p,L);pr=-np.sqrt(2*(C.EX+C.V0-E_before))
    p[0]+=pr*np.cos(theta);p[1]+=pr*np.sin(theta)
    initial=dict(N=N,L=L,seed=seed,neutral=neutral,relative_seed=amplitude,energy=energy(q,p,L),
        nQ=float(np.mean(q[0]*p[1]-q[1]*p[0])),radial_momentum=pr,
        modes=dict(amplitude=am,phase=pm,neutral=sm),
        scope='Cold microcanonical classical seed control; no vector heat, cosmic covariance or quantum vacuum samples.')
    return q,p,initial

def snapshot(path):
    with open(path,'rb') as f:N,L,t,bits=struct.unpack('<iddi',f.read(24))
    arr=np.memmap(path,dtype='<f8' if bits==64 else '<f4',mode='r',offset=24,shape=(6,N,N,N))
    return N,L,t,arr

def run(tag,N=64,L=48.,seed=1,neutral='mixed',dt=.06,T=200.,sample=1.,snap=40.,order=4,threads=4,bits=32,amplitude=1e-6,initial_arrays=None):
    directory=WORK/tag;directory.mkdir(exist_ok=True);prefix=directory/'state';init=directory/'initial.bin'
    if initial_arrays is None:q,p,params=initialize(N,L,seed,neutral,amplitude)
    else:q,p=initial_arrays;params=dict(N=N,L=L,verification_initial_state=True)
    write_initial(init,q,p,L);params.update(dt_requested=dt,T=T,sample=sample,snapshot_interval=snap,order=order,threads=threads,snapshot_bits=bits)
    (directory/'parameters.json').write_text(json.dumps(params,indent=2)+'\n');del q,p
    cmd=[str(build()),str(init),str(prefix),str(dt),str(T),str(sample),str(snap),str(order),str(threads),str(bits)]
    start=time.time()
    with open(directory/'run.log','w') as log:subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)
    data=np.genfromtxt(str(prefix)+'.tsv',names=True,delimiter='\t')
    result=dict(tag=tag,parameters=params,seconds=time.time()-start,
        charge_relative_drift=float(max(abs(data['nQ']-data['nQ'][0]))/max(abs(data['nQ'][0]),1.)),
        energy_relative_drift=float(max(abs(data['rho']/data['rho'][0]-1))),
        final={name:float(data[name][-1]) for name in data.dtype.names},
        snapshots=[str(p) for p in sorted(directory.glob('state_s*.bin'))])
    (directory/'run_result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result

if __name__=='__main__':
    r=run('benchmark_N48',N=48,L=48.,T=10.,snap=10.,threads=4);save('BENCHMARK.json',r);print(r,flush=True)
