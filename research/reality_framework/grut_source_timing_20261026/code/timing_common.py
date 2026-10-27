"""Pinned quiet histories, linear equations and late-width source clocks."""
from pathlib import Path
import os,sys,json,hashlib,subprocess,tempfile
import numpy as np
from scipy.integrate import solve_ivp,quad
from scipy.special import gammainc,gammaincc,gamma
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[2]
WORK=Path(os.environ.get('GRUT_TIMING_WORKDIR','/tmp/grut_source_timing_20261026'));WORK.mkdir(parents=True,exist_ok=True)
MODE=ROOT.parent/'grut_reheating_modes_20261026'
sys.path.insert(0,str(MODE/'code'))
import modes as M
sys.path.insert(0,str(ROOT.parent/'grut_kinetic_continuation_20261026'/'code'))
import ks_common as K
KS=np.array([.02,.05,.1,.2,.4,.8,1.2,1.6,2.,2.13,2.4,2.8,3.2,4.,5.6,5.8,6.,8.])

def save(name,value):
    def cv(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,Path):return str(x)
        raise TypeError(type(x).__name__)
    (ROOT/name).write_text(json.dumps(value,indent=2,default=cv,allow_nan=False)+'\n')

def read(name):return json.loads((ROOT/name).read_text())

def phase_seed(beta,step=1e-5):
    plus,_=M.prepared(beta,delta=step);minus,_=M.prepared(beta,delta=-step)
    q=(plus[2:8]-minus[2:8])/(2*step)*np.sqrt(2)
    return np.r_[q,q]

def build():
    from scipy.integrate._ivp import dop853_coefficients as dc
    src=ROOT/'code/response.c';key=hashlib.sha256(src.read_bytes()).hexdigest()[:12]
    directory=WORK/('build_'+key);directory.mkdir(exist_ok=True)
    def fmt(name,arr):
        a=np.asarray(arr)
        if a.ndim==1:return 'static const double '+name+'['+str(len(a))+']={'+','.join(format(x,'.17g') for x in a)+'};\n'
        return 'static const double '+name+'['+str(a.shape[0])+']['+str(a.shape[1])+']={'+','.join('{'+','.join(format(x,'.17g') for x in row)+'}' for row in a)+'};\n'
    header=''.join(fmt(n,a) for n,a in [('AA',dc.A[:12,:12]),('CC',dc.C[:12]),('BB',dc.B),('EE3',dc.E3),('EE5',dc.E5)])
    (directory/'dop853.h').write_text(header)
    exe=directory/'response'
    if not exe.exists():subprocess.run(['gcc','-O3','-std=c11','-Wall','-Wextra','-I',str(directory),str(src),'-lm','-o',str(exe)],check=True)
    return exe

def integrate(exe,tag,beta,decay,ks,times,e=1e-7,rtol=2e-10):
    import time
    v,res=M.prepared(beta);seed=phase_seed(beta)
    params=np.array([beta,decay,e,M.XI,M.ETA2,M.MU,M.KAP,times[-1],rtol],dtype='<f8')
    inp=WORK/(tag+'.in');out=WORK/(tag+'.bin')
    with inp.open('wb') as f:
        np.array([len(ks),len(times)],dtype='<i4').tofile(f);params.tofile(f);np.asarray(v,dtype='<f8').tofile(f)
        np.asarray(ks,dtype='<f8').tofile(f);np.asarray(times,dtype='<f8').tofile(f)
        np.tile(seed,(len(ks),1)).astype('<f8').tofile(f)
    start=time.monotonic();proc=subprocess.run([str(exe),str(inp),str(out)],capture_output=True,text=True,check=True)
    data=np.fromfile(out,dtype='<f8').reshape(len(times),20+12*len(ks))
    assert np.allclose(data[:,0],times,rtol=0,atol=1e-10)
    return data[:,1:20],data[:,20:].reshape(len(times),len(ks),12),dict(
        check=proc.stdout.strip(),seconds=time.monotonic()-start,raw_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),preparation_residual=res)

def responses(beta,bg,vectors):
    a=np.exp(bg[:,15]);z=bg[:,2:8:2]+1j*bg[:,3:8:2];vel=bg[:,8:14:2]+1j*bg[:,9:14:2]
    chi=np.sqrt(2)*a[:,None]*z;chip=np.sqrt(2)*a[:,None]**2*(vel+bg[:,16,None]*z)
    q=np.empty((len(bg),6));p=q.copy();q[:,::2]=chi.real;q[:,1::2]=chi.imag;p[:,::2]=chip.real;p[:,1::2]=chip.imag
    radius=np.sum(q[:,:2]**2,axis=1)
    amplitude=np.einsum('tki,ti->tk',vectors[:,:,:2],q[:,:2])/radius[:,None]
    dq=vectors[:,:,0]*p[:,None,1]-vectors[:,:,1]*p[:,None,0]+q[:,None,0]*vectors[:,:,7]-q[:,None,1]*vectors[:,:,6]
    nx=q[:,0]*p[:,1]-q[:,1]*p[:,0]
    # Exactly the old per-radian charge normalization, with a mean through64--80.
    return amplitude,dq,nx

def crossings(times,a,values,ks,level):
    maximum=np.max(abs(values),axis=1);where=np.flatnonzero(maximum>=level)
    if not len(where):return None
    i=int(where[0]);j=int(np.argmax(abs(values[i])))
    return dict(level=level,k=float(ks[j]),tau_bracket=[float(times[max(0,i-1)]),float(times[i])],
                a_bracket=[float(a[max(0,i-1)]),float(a[i])],sampled_value=float(maximum[i]))
