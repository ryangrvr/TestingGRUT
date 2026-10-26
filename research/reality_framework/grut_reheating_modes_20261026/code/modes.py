"""Finite-time scalar/longitudinal tangent map; no periodic Floquet assumption."""
import sys,math,json,subprocess,hashlib,tempfile
from concurrent.futures import ThreadPoolExecutor,as_completed
import numpy as np
from scipy.integrate._ivp import dop853_coefficients as dc
from mode_common import ROOT,WORK,PARENT,write
sys.path.insert(0,str(PARENT/'code'))
from reheating import prepared,potential_force,geometric,rhs,XI,ETA2,MU,KAP,H0,MRED,INP

J=np.block([[np.zeros((6,6)),np.eye(6)],[-np.eye(6),np.zeros((6,6))]])
TCH=np.zeros((6,6))
for i,q in enumerate([1,3,6]):TCH[2*i,2*i+1]=-q;TCH[2*i+1,2*i]=q

def realvec(z):
    out=np.empty(6);out[::2]=np.sqrt(2)*z.real;out[1::2]=np.sqrt(2)*z.imag
    return out

def matrices(v,beta,xi=XI,eta2=ETA2):
    z,vel,V,F,f,fd,R,hd=geometric(v,beta,xi,eta2);x,y,s=z;a=np.exp(v[15])
    A=np.zeros((3,3),complex);B=A.copy()
    A[0,0]=.004*abs(x)**2+xi*R;B[0,0]=.002*x*x+6*KAP*y*x.conjugate()
    A[0,1]=3*KAP*x.conjugate()**2;A[1,0]=A[0,1].conjugate()
    A[1,1]=4*abs(y)**2+900;B[1,1]=2*y*y+2*MU*s
    A[1,2]=2*MU*y.conjugate();A[2,1]=A[1,2].conjugate()
    A[2,2]=4*abs(s)**2-1600;B[2,2]=2*s*s
    M=np.empty((6,6));M[::2,::2]=(A+B).real;M[::2,1::2]=-(A-B).imag
    M[1::2,::2]=(A+B).imag;M[1::2,1::2]=(A-B).real
    M=a*a*(M-R/6*np.eye(6))
    chi=a*realvec(z);chip=a*a*realvec(vel+v[16]*z)
    return M,TCH@chi,TCH@chip,chi,chip

def generator(v,beta,k,e=1e-7):
    M,g,gp,_,_=matrices(v,beta);c=e*e/k**2 if e else 0.;a=np.exp(v[15])
    return np.block([[-c*np.outer(g,gp),np.eye(6)+c*np.outer(g,g)],
                     [-k*k*np.eye(6)-M-c*np.outer(gp,gp),c*np.outer(gp,g)]])/a

def build():
    out=__import__('pathlib').Path(tempfile.mkdtemp(prefix='build_',dir=WORK))
    def fmt(name,arr):
        a=np.asarray(arr)
        if a.ndim==1:return 'static const double '+name+'['+str(len(a))+']={'+','.join(format(x,'.17g') for x in a)+'};\n'
        return 'static const double '+name+'['+str(a.shape[0])+']['+str(a.shape[1])+']={'+','.join('{'+','.join(format(x,'.17g') for x in row)+'}' for row in a)+'};\n'
    (out/'dop853.h').write_text(''.join(fmt(n,a) for n,a in [('AA',dc.A[:12,:12]),('CC',dc.C[:12]),('BB',dc.B),('EE3',dc.E3),('EE5',dc.E5)]))
    exe=out/'modes';subprocess.run(['gcc','-O3','-std=c11','-Wall','-Wextra','-I',str(out),str(ROOT/'code/modes.c'),'-lm','-o',str(exe)],check=True)
    return exe

def integrate(exe,tag,beta,decay,ks,times,rtol=2e-10,e=1e-7):
    v,res=prepared(beta);pars=np.array([beta,decay,e,XI,ETA2,MU,KAP,times[-1],rtol],dtype='<f8')
    inp=WORK/(tag+'.in');out=WORK/(tag+'.bin')
    with inp.open('wb') as f:
        np.array([len(ks),len(times)],dtype='<i4').tofile(f);pars.tofile(f);np.asarray(v,dtype='<f8').tofile(f)
        np.asarray(ks,dtype='<f8').tofile(f);np.asarray(times,dtype='<f8').tofile(f)
    proc=subprocess.run([str(exe),str(inp),str(out)],capture_output=True,text=True,check=True)
    q=np.fromfile(out,dtype='<f8').reshape(len(times),20+144*len(ks))
    assert np.allclose(q[:,0],times,atol=1e-12,rtol=0)
    return q[:,1:20],q[:,20:].reshape(len(times),len(ks),12,12),dict(raw_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),check=proc.stdout.strip(),preparation_residual=res)

def diagnostics(beta,bg,trans,ks):
    nt,nk=trans.shape[:2];out=np.zeros((nt,nk,8));M0,g0,gp0,chi0,_=matrices(bg[0],beta)
    vals,ev=np.linalg.eigh(M0)
    scales=np.repeat([np.linalg.norm(chi0[2*i:2*i+2]) for i in range(3)],2)
    D=np.diag(scales);sym_abs=sym_scaled=0.
    for j,k in enumerate(ks):
        omega=np.sqrt(k*k+np.maximum(vals,0)+1)
        W=(ev*omega)@ev.T;Wi=(ev/omega)@ev.T
        cov=np.block([[Wi,np.zeros((6,6))],[np.zeros((6,6)),W]])/2
        Ainv=np.eye(6)-1e-14*np.outer(g0,g0)/(k*k+1e-14*(g0@g0))
        init=np.block([[D,np.zeros((6,6))],[D,Ainv@D]])
        for n,v in enumerate(bg):
            M,g,gp,chi,chip=matrices(v,beta);TT=trans[n,j]
            amp=np.zeros(12);amp[:2]=chi[:2]/(chi[:2]@chi[:2])
            charge=np.zeros(12);charge[:2]=[-chip[1],chip[0]];charge[6:8]=[chi[1],-chi[0]]
            nx=-2*np.imag((v[2]+1j*v[3]).conjugate()*(v[8]+1j*v[9]))*np.exp(3*v[15])
            ar=amp@TT;nr=charge@TT
            pref=k**3/(2*np.pi**2)
            # Charge variance is in comoving N_X units; normalization is
            # applied after the background's late mean has been determined.
            out[n,j]=[np.linalg.norm(ar@init),np.linalg.norm(nr@init),
                      pref*(ar@cov@ar),pref*(nr@cov@nr),np.linalg.eigvalsh(M)[0],nx,
                      np.exp(v[15]),v[16]]
            if n in (0,nt//2,nt-1):
                err=float(np.max(abs(TT.T@J@TT-J)));sym_abs=max(sym_abs,err)
                sym_scaled=max(sym_scaled,err/max(1.,float(np.linalg.norm(TT,ord=2)**2)))
    return out,dict(symplectic_absolute_max=sym_abs,symplectic_scaled_max=sym_scaled)

def main():
    exe=build();ks=np.r_[1e-5,3e-5,1e-4,3e-4,.001,.003,np.geomspace(.01,40,33)]
    times=np.unique(np.r_[np.linspace(0,3,61),np.geomspace(3,80,81),[10.,40.]])
    points=[(3.,.5),(5.,1.),(10.,.1),(10.,3.),(12.,.1),(2.,.1)]
    jobs=[]
    for b,gam in points:
        tag=f'b{b:g}_g{gam:g}'
        for start in range(0,len(ks),8):jobs.append((tag,b,gam,start,ks[start:start+8]))
    assembled={tag:[] for tag,_,_,_,_ in jobs}
    with ThreadPoolExecutor(max_workers=3) as pool:
        tasks={pool.submit(integrate,exe,f'{tag}_k{start}',b,gam,kk,times):(tag,b,gam,start,kk) for tag,b,gam,start,kk in jobs}
        for future in as_completed(tasks):
            tag,b,gam,start,kk=tasks[future];bg,TT,control=future.result()
            np.savez_compressed(WORK/f'{tag}_k{start}.npz',bg=bg,trans=TT,ks=kk,times=times)
            assembled[tag].append((start,bg,TT,control));print('BATCH',tag,start,control['check'],flush=True)
    summaries=[];curves={}
    for b,gam in points:
        tag=f'b{b:g}_g{gam:g}';parts=sorted(assembled[tag],key=lambda x:x[0]);bg=parts[0][1]
        TT=np.concatenate([r[2] for r in parts],axis=1);diag,check=diagnostics(b,bg,TT,ks)
        nx=diag[:,0,5];late=times>=64;mean=float(np.trapezoid(nx[late],times[late])/(times[late][-1]-times[late][0]))
        diag[:,:,1]/=abs(mean);diag[:,:,3]/=mean*mean
        rho_rad=bg[-1,14]*MRED*MRED*H0*H0;temp=(30*rho_rad/(np.pi**2*106.75))**.25
        sent=2*np.pi**2/45*106.75*temp**3;yieldratio=abs(mean)*H0**3/(np.exp(3*bg[-1,15])*sent*INP['benchmark']['YQ'])
        intervals=[]
        for end in (3.,10.,40.,80.):
            mask=times<=end;main=ks>=.01;sub=diag[mask][:,main]
            peak=np.unravel_index(np.argmax(sub[:,:,0]),sub[:,:,0].shape)
            vk=np.trapezoid(sub[:,:,2],np.log(ks[main]),axis=1)
            vn=np.trapezoid(sub[:,:,3],np.log(ks[main]),axis=1)
            intervals.append(dict(end_tau=end,max_amplitude_response=float(sub[:,:,0].max()),
                max_charge_response=float(sub[:,:,1].max()),peak_k=float(ks[main][peak[1]]),peak_tau=float(times[mask][peak[0]]),
                conditional_Gaussian_rms_fraction_amplitude=float(np.sqrt(vk.max())),
                conditional_Gaussian_rms_fraction_charge=float(np.sqrt(vn.max()))))
        summaries.append(dict(tag=tag,beta=b,Gamma_over_H=gam,yield_over_YQ=yieldratio,late_NX=mean,intervals=intervals,
                              controls=check,batches=[r[3] for r in parts]))
        curves[tag]=diag;curves[tag+'_bg']=bg
        # Preserve representative exact transfer matrices for review/control.
        sel=np.unique(np.searchsorted(times,[0,1,3,10,40,80]));ik=[0,2,6,14,22,30,38]
        curves[tag+'_transfer']=TT[np.ix_(sel,ik,np.arange(12),np.arange(12))]
        print('RESULT',tag,yieldratio,intervals,check,flush=True)
    curves['ks']=ks;curves['times']=times
    np.savez_compressed(ROOT/'MODE_PROFILES.npz',**curves)
    write('MODE_RESULTS.json',dict(rows=summaries,ks=ks,times=times,
        profile_columns=['amplitude_relative_response','charge_relative_response','Gaussian_amplitude_power','Gaussian_charge_power','lowest_conformal_mass_eigenvalue','NX','a','H'],
        scope=__doc__,physical_cosmological_initial_covariance_derived=False))

if __name__=='__main__':main()
