"""Read-only interfaces to the reviewed scalar pilot; unique module name."""
from pathlib import Path
import os,sys,json,hashlib,subprocess,time,ctypes,struct
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np
from scipy.fft import fftn,ifftn,rfftn,irfftn,fft
from scipy.integrate import solve_ivp,quad
ROOT=Path(__file__).resolve().parents[1]; REPO=ROOT.parents[2]
PARENT=ROOT.parent/'grut_early_nonlinear_20261027'
WORK=Path(os.environ.get('GRUT_PARTICLE_WORKDIR',str(REPO.parent/'grut_particle_transport_work')))
WORK.mkdir(parents=True,exist_ok=True)
os.environ.setdefault('GRUT_EARLY_WORKDIR',str(WORK/'parent'))
os.environ.setdefault('GRUT_MODE_WORKDIR',str(WORK/'modes'))
os.environ.setdefault('GRUT_REHEATING_WORKDIR',str(WORK/'reheating'))
sys.path.insert(0,str(PARENT/'code'))
import en_common as E
LX,G2,LP=E.LX,E.G2,E.LP

def save(name,obj):
    def convert(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    (ROOT/name).write_text(json.dumps(obj,indent=2,default=convert,allow_nan=False)+'\n')

def read(name):return json.loads((ROOT/name).read_text())
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def geometry(N,L,half=False):
    m=np.fft.fftfreq(N)*N;z=np.fft.rfftfreq(N)*N if half else m
    k2=(2*np.pi/L)**2*(m[:,None,None]**2+m[None,:,None]**2+z[None,None,:]**2)
    x=2*np.pi*m/N;y=2*np.pi*z/N
    f=lambda u:4*np.sin(u/2)**2*(1+np.sin(u/2)**2/3)/(L/N)**2
    kl=f(x)[:,None,None]+f(x)[None,:,None]+f(y)[None,None,:]
    return k2,kl

def library():
    src=ROOT/'code/fieldstep.cpp'; out=WORK/('step_'+sha(src)[:16]+'.so')
    if not out.exists():
        temporary=out.with_name(out.name+f'.{os.getpid()}.tmp')
        subprocess.run(['g++','-std=c++17','-O3','-ffp-contract=off','-fopenmp','-fPIC','-shared',str(src),'-o',str(temporary)],check=True)
        os.replace(temporary,out)
    lib=ctypes.CDLL(str(out));dp=ctypes.POINTER(ctypes.c_double)
    lib.evolve.argtypes=[dp,dp,ctypes.c_int,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_int]
    lib.evolve.restype=ctypes.c_int
    return lib

def step(lib,state,duration,hmax=.002,tangent=None,threads=2):
    assert state.dtype==np.float64 and state.flags.c_contiguous
    nsteps=int(np.ceil(duration/hmax-1e-12));h=duration/nsteps if nsteps else 0.
    ptr=ctypes.POINTER(ctypes.c_double)
    dp=tangent.ctypes.data_as(ptr) if tangent is not None else ptr()
    nt=0 if tangent is None else len(tangent)
    assert tangent is None or (tangent.dtype==np.float64 and tangent.flags.c_contiguous)
    rc=lib.evolve(state.ctypes.data_as(ptr),dp,nt,state.shape[-1],float(step.L),h,nsteps,threads)
    if rc:raise RuntimeError(('evolve',rc))

def set_box(L):step.L=float(L)

def clock(s):
    a0=8.869218738806; C=1.9857962750634526
    return 20+a0*np.asarray(s)+.5*C*np.asarray(s)**2

def invariants(z,L):
    return E.energy(z[:3],z[3:],L),float(np.mean(z[0]*z[4]-z[1]*z[3]))
