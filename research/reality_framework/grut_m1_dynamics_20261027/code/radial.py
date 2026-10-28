"""Open-boundary nonlinear survival, with an explicit absorbing receiver."""
from common import *

CASES={
 'radial_main':dict(h=.2,dt=.025,L=200.,ab=100.,amp=.05),
 'radial_unperturbed':dict(h=.2,dt=.025,L=200.,ab=100.,amp=0.),
 'radial_fine':dict(h=.1,dt=.0125,L=200.,ab=100.,amp=.05),
 'radial_half':dict(h=.2,dt=.0125,L=200.,ab=100.,amp=.05),
 'radial_boundary':dict(h=.2,dt=.025,L=300.,ab=200.,amp=.05),
}

def run_radial(tag,params):
    src=ROOT/'code/radial.cpp';key=hashlib.sha256(src.read_bytes()).hexdigest()[:16];exe=WORK/('radial_'+key)
    if not exe.exists():subprocess.run(['g++','-std=c++17','-O3','-ffp-contract=off',str(src),'-o',str(exe)],check=True)
    codehash=hashlib.sha256(src.read_bytes()+Path(__file__).read_bytes()).hexdigest()
    d=WORK/tag;d.mkdir(exist_ok=True);result=d/'result.json'
    if result.exists():
        out=json.loads(result.read_text())
        if out['code_sha256']==codehash:return out
        raise RuntimeError('Stale radial result')
    h,dt,L,ab,amp=[params[k] for k in ('h','dt','L','ab','amp')];N=round(L/h);r=(np.arange(N)+.5)*h
    a=np.loadtxt(PARENT/'PROFILE_M1.csv',delimiter=',',skiprows=1)
    fs=CubicHermiteSpline(a[:,0],a[:,1],a[:,2]);ss=CubicHermiteSpline(a[:,0],a[:,3],a[:,4])
    f=fs(np.minimum(r,160));s=ss(np.minimum(r,160));f[r>160]=0;s[r>160]=1
    shape=np.exp(-(r/12)**2);f*=1+amp*shape;s+=amp*(1-s)*shape
    om=json.loads((PARENT/'DEFORMATION_RESULTS.json').read_text())['primary']['Omega']
    vol=4*np.pi/3*((r+h/2)**3-(r-h/2)**3)
    py=-om*f;py*=2000/np.sum(-f*py*vol)
    initial=d/'initial.csv';prefix=d/'state'
    np.savetxt(initial,np.c_[r,f,np.zeros_like(r),s,np.zeros_like(r),py,np.zeros_like(r)],delimiter=',',
               header='rho,x,y,s,px,py,ps',comments='')
    begin=time.time()
    cmd=[str(exe),str(initial),str(prefix),str(h),str(dt),str(L),'1000','60',str(ab),'.2']
    with open(d/'run.log','w') as log:subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT)
    hist=np.genfromtxt(str(prefix)+'_history.csv',names=True,delimiter=',')
    E=hist['E']+hist['E_abs'];Q=hist['Q']+hist['Q_abs']
    out={'tag':tag,'params':params,'code_sha256':codehash,'seconds':time.time()-begin,
         'receiver_energy_balance_relative_error':float(np.max(abs(E/E[0]-1))),
         'receiver_charge_balance_relative_error':float(np.max(abs(Q/Q[0]-1))),
         'initial_E_over_Q':float(hist['E'][0]/hist['Q'][0]),
         'final_E_over_Q':float(hist['E'][-1]/hist['Q'][-1]),
         'final_core20_E_over_Q':float(hist['E20'][-1]/hist['Q20'][-1]),
         'final_core20_Q':float(hist['Q20'][-1]/G2),
         'absorbed_Q_over_original':float(hist['Q_abs'][-1]/hist['Q'][0]),
         'absorbed_energy_over_initial':float(hist['E_abs'][-1]/hist['E'][0])}
    result.write_text(json.dumps(out,indent=2)+'\n');return out

def main():
    out={};arrays={}
    for tag,params in CASES.items():
        r=run_radial(tag,params);out[tag]=r;print(json.dumps(r),flush=True)
        arrays[tag+'_history']=np.genfromtxt(WORK/tag/'state_history.csv',names=True,delimiter=',')
        arrays[tag+'_final']=np.genfromtxt(WORK/tag/'state_final.csv',names=True,delimiter=',')
    save('RADIAL_RESULTS.json',{'scope':__doc__,'cases':out,'physical_cosmic_p':None})
    np.savez_compressed(ROOT/'RADIAL_DATA.npz',**arrays)

if __name__=='__main__':main()
