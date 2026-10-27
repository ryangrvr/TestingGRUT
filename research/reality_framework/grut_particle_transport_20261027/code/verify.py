"""Tests of the new streaming map and its exact discrete tangents."""
from pt_common import *

def main():
    N=8;L=5.;set_box(L);lib=library();rng=np.random.default_rng(712)
    z=rng.normal(size=(6,N,N,N));initial=z.copy();T=.08
    step(lib,z,T,.002,threads=2)
    def rhs(t,v):
        a=v.reshape(6,N,N,N);return np.r_[a[3:],E.lap(a[:3],L)-E.gradient(a[:3])].ravel()
    sol=solve_ivp(rhs,(0,T),initial.ravel(),method='DOP853',rtol=2e-12,atol=2e-13)
    separate=float(np.max(abs(z.ravel()-sol.y[:,-1])))
    directory=WORK/'short_reviewed_engine';directory.mkdir(exist_ok=True)
    E.write_initial(directory/'initial.bin',initial[:3],initial[3:],L)
    subprocess.run([str(E.build()),str(directory/'initial.bin'),str(directory/'state'),'.002',str(T),'.02',str(T),'4','2','64','0','1','0'],check=True,capture_output=True)
    _,_,_,old=E.snapshot(sorted(directory.glob('state_s*.bin'))[-1]);reviewed=float(np.max(abs(z-old)))
    d=rng.normal(size=(2,6,N,N,N));d0=d.copy();a=initial.copy()
    step(lib,a,T,.002,tangent=d,threads=2)
    def omega(x,y):return float(np.sum(x[:3]*y[3:]-x[3:]*y[:3]))
    symplectic=abs(omega(d[0],d[1])-omega(d0[0],d0[1]))/max(abs(omega(d0[0],d0[1])),1.)
    fd=[]
    for eps in [1e-3,5e-4]:
        p=initial+eps*d0[0];m=initial-eps*d0[0]
        step(lib,p,T,.002,threads=2);step(lib,m,T,.002,threads=2)
        fd.append(float(np.max(abs((p-m)/(2*eps)-d[0]))))
    mean=np.load(PARENT/'MATCH_L12.npz');y=np.r_[mean['light_mean_q'],mean['light_mean_p']]
    def hom(t,x):return np.r_[x[3:],-E.gradient(x[:3,None,None,None]).ravel()]
    exact=solve_ivp(hom,(0,1.),y,rtol=2e-13,atol=1e-12,method='DOP853').y[:,-1]
    errors=[]
    for h in [.02,.01,.005]:
        a=np.broadcast_to(y[:,None,None,None],(6,N,N,N)).copy();step(lib,a,1.,h,threads=2)
        errors.append(float(np.linalg.norm(a[:,0,0,0]-exact)))
    ratios=np.array(errors[:-1])/errors[1:]
    passed=separate<2e-9 and reviewed<2e-11 and max(fd)<2e-7 and fd[1]<.4*fd[0] and symplectic<1e-12 and np.all((ratios>14)&(ratios<19))
    result=dict(passed=bool(passed),FFT_DOP853_max_absolute=separate,reviewed_engine_max_absolute=reviewed,
        tangent_finite_difference_errors=fd,tangent_symplectic_relative=symplectic,homogeneous_errors=errors,homogeneous_refinement_ratios=ratios,
        scope='Same-author independent formulas and reviewed-engine comparison; not external replication.')
    save('ENGINE_CHECKS.json',result);print(json.dumps(result,default=lambda x:x.tolist(),indent=2));assert passed

if __name__=='__main__':main()
