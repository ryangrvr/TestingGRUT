"""Pinned early histories and a regulated, explicit reduced scalar model."""
from pathlib import Path
import os,sys,json,hashlib,struct,subprocess,time
import numpy as np
from scipy.integrate import solve_ivp
from scipy.fft import rfftn,irfftn
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[2]
WORK=Path(os.environ.get('GRUT_EARLY_WORKDIR','/tmp/grut_early_nonlinear_20261027'));WORK.mkdir(parents=True,exist_ok=True)
TIMING=ROOT.parent/'grut_source_timing_20261026'
sys.path.insert(0,str(TIMING/'code'))
import timing_common as T
M=T.M
LX=.001;G2=.01;LP=.03

def save(name,obj):
    def cv(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        if isinstance(v,Path):return str(v)
        raise TypeError(type(v).__name__)
    (ROOT/name).write_text(json.dumps(obj,indent=2,default=cv,allow_nan=False)+'\n')

def read(name):return json.loads((ROOT/name).read_text())

def symbol(N,L):
    a=2*np.pi*np.fft.fftfreq(N);b=2*np.pi*np.fft.rfftfreq(N);dx=L/N
    def f(x):return 4*np.sin(x/2)**2*(1+np.sin(x/2)**2/3)/dx**2
    return f(a)[:,None,None]+f(a)[None,:,None]+f(b)[None,None,:]

def lap(q,L):return irfftn(-symbol(q.shape[-1],L)*rfftn(q,axes=(-3,-2,-1)),s=q.shape[-3:],axes=(-3,-2,-1))

def gradient(q):
    r=q[0]**2+q[1]**2;s=q[2]
    return np.array([(LX*r+G2*s*s)*q[0],(LX*r+G2*s*s)*q[1],(LP*s*s+G2*r)*s])

def energy(q,p,L):
    r=q[0]**2+q[1]**2
    return float(np.mean(.5*np.sum(p*p-q*lap(q,L),axis=0)+LX/4*r*r+G2/2*q[2]**2*r+LP/4*q[2]**4))

def apo(q,p):
    E=.5*np.sum(p*p,axis=-1)+LX/4*np.sum(q*q,axis=-1)**2
    Q=q[...,0]*p[...,1]-q[...,1]*p[...,0]
    arg=np.clip(-3*Q**2/(4*E)*np.sqrt(3*LX/(4*E)),-1,1)
    y=2*np.sqrt(4*E/(3*LX))*np.cos(np.arccos(arg)/3)
    return np.sqrt(y)

def build():
    src=ROOT/'code/early3d.cpp';key=hashlib.sha256(src.read_bytes()).hexdigest()[:16];exe=WORK/('early3d_'+key)
    if not exe.exists():subprocess.run(['g++','-std=c++17','-O3','-ffp-contract=off','-fopenmp',str(src),'-o',str(exe)],check=True)
    return exe

def write_initial(path,q,p,L):
    with open(path,'wb') as f:
        f.write(b'GRUT3D1\n');f.write(struct.pack('<id4d',q.shape[-1],L,0.,G2,LX,LP))
        np.asarray(q,dtype='<f8').tofile(f);np.asarray(p,dtype='<f8').tofile(f)

def snapshot(path):
    with open(path,'rb') as f:N,L,s,bits=struct.unpack('<iddi',f.read(24))
    a=np.memmap(path,dtype='<f8',offset=24,mode='r',shape=(6,N,N,N));return N,L,s,a
