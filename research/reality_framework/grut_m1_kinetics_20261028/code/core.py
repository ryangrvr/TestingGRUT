"""Full reciprocal angular strengths; same FD operator, batched incident columns."""
from setup import *
import subprocess,ctypes,time
from scipy.linalg import solve_banded
from scipy.interpolate import PchipInterpolator

def library():
    src=ROOT/'code/exterior.cpp';lib=WORK/('exterior_'+sha(src)[:16]+'.so')
    if not lib.exists():subprocess.run(['g++','-O3','-std=c++17','-shared','-fPIC',str(src),'-o',str(lib)],check=True)
    obj=ctypes.CDLL(str(lib));ptr=ctypes.POINTER(ctypes.c_double)
    obj.ratios.argtypes=[ptr,ctypes.c_int,ctypes.c_double,ctypes.c_double,ctypes.c_double,ptr]
    return obj

LIB=library()
def matrix(k,ell,h=.05,R=60.,Rfar=500.):
    n=round(R/h);r=np.arange(1,n+1)*h;f,s=profile(r)
    ea=np.sqrt(1+k*k);nu=ea+OMEGA;energies=np.array([ea+2*OMEGA,ea,nu])
    k2=energies**2-np.array([1.,1.,6.]);raw=np.empty(6)
    ptr=ctypes.POINTER(ctypes.c_double)
    LIB.ratios(k2.ctypes.data_as(ptr),ell,h,R,Rfar,raw.ctypes.data_as(ptr))
    q=raw[::2]+1j*raw[1::2];opened=k2>0;flux=np.where(opened,q.imag/h,0)
    diag=np.c_[s*s+2*A*f*f-energies[0]**2,s*s+2*A*f*f-ea*ea,B*(3*s*s-1)-6*C*(s-1)+f*f-nu*nu]
    diag+=ell*(ell+1)/r[:,None]**2+2/h**2
    ab=np.zeros((7,3*n),complex);ab[3]=diag.ravel();ab[0,3:]=-1/h**2;ab[6,:-3]=-1/h**2
    for left,right,v in [(0,1,A*f*f),(0,2,np.sqrt(2)*f*s),(1,2,np.sqrt(2)*f*s)]:
        ii=np.arange(n)*3+left;jj=np.arange(n)*3+right;ab[3+left-right,jj]=v;ab[3+right-left,ii]=v
    ab[3,-3:]-=q/h**2
    channels=np.flatnonzero(opened & (flux>1e-250));rhs=np.zeros((3*n,len(channels)),complex)
    for j,c in enumerate(channels):rhs[-3+c,j]=(np.conj(q[c])-q[c])/h**2/np.sqrt(flux[c])
    solution=solve_banded((3,3),ab,rhs,check_finite=False)[-3:]
    for j,c in enumerate(channels):solution[c,j]-=1/np.sqrt(flux[c])
    smat=solution[channels]*np.sqrt(flux[channels,None]);p=np.zeros((3,3))
    p[np.ix_(channels,channels)]=abs(smat)**2
    return p,float(np.max(abs(smat.conj().T@smat-np.eye(len(channels))))),float(np.max(abs(smat-smat.T)))

def strengths(k,h=.05,l_extra=0):
    # The old packet explicitly checked ell32 through k=3. For larger k,
    # increase the angular cutoff with the size of the actual core.
    lmax=max(32,int(np.ceil(12*k)))+l_extra
    w=np.zeros((3,3));unit=rec=0.
    for ell in range(lmax+1):
        p,u,r=matrix(float(k),ell,h=h);w+=(2*ell+1)*p;unit=max(unit,u);rec=max(rec,r)
    return w,unit,rec,lmax

def main():
    # Resolve the previously identified ell4 resonance and its shoulders.
    ks=np.unique(np.r_[np.geomspace(.002,.1,18),np.linspace(.1,1.2,111),
                       np.linspace(.326,.355,60),np.linspace(1.21,3.3,106),
                       np.linspace(3.35,12.0,72),[.4,.7,1.,1.25,1.5,2.,3.]])
    cache_stamp=hashlib.sha256((sha(M1/'PROFILE_M1.csv')+sha(ROOT/'code/core.py')+sha(ROOT/'code/exterior.cpp')+repr(OMEGA)).encode()).hexdigest()[:16]
    rows=[];checks=[];t=time.monotonic()
    for i,k in enumerate(ks):
        path=WORK/('core_'+cache_stamp+'_'+format(float(k),'.12g')+'.npz')
        if path.exists():
            with np.load(path) as d: w=d['w'];u=float(d['u']);rec=float(d['rec']);lm=int(d['lm'])
        else:
            w,u,rec,lm=strengths(k);np.savez(path,w=w,u=u,rec=rec,lm=lm)
        rows.append(w);checks.append((u,rec,lm))
        if (i+1)%40==0: print('core',i+1,'/',len(ks),'seconds',round(time.monotonic()-t,1),flush=True)
    rows=np.array(rows);np.savez_compressed(ROOT/'CORE_KERNELS.npz',k=ks,W=rows,checks=np.array(checks))
    from scattering import scatter
    old=[]
    for k in (.4,1.,1.5,3.):
        p,_,_=matrix(k,0);a=scatter(k,return_amplitudes=True)
        old.append({'k':k,'max_probability_difference':float(np.max(abs(p[:,1]-[a['P_X'],a['P_antiparticle'],a['P_phi']])) )})
    refined=[]
    for k in (1.,1.5,4.,8.):
        base,u,r,l=strengths(k);fine,*_=strengths(k,h=.025);tail,*_=strengths(k,l_extra=24)
        off=~np.eye(3,dtype=bool)
        refined.append({'k':k,'lmax':l,'base_off_diagonal':base[off],
                        'fine_off_diagonal':fine[off],'tail_off_diagonal':tail[off]})
    write_json('CORE_CHECKS.json',{'same_operator_comparison':old,'refinement':refined,
        'max_unitarity_error':max(x[0] for x in checks),'max_reciprocity_error':max(x[1] for x in checks),
        'sampled_k_count':len(ks),'canonical_neutral_coupling':'sqrt(2) F S'})

if __name__=='__main__':main()
