"""Regenerate prior pulse traces from published inputs; verify archived histories."""
from common import *

def main():
    src=DYN/'code/radial.cpp';key=hashlib.sha256(src.read_bytes()).hexdigest()[:16];exe=WORK/('radial_'+key)
    if not exe.exists():subprocess.run(['g++','-std=c++17','-O3','-ffp-contract=off',str(src),'-o',str(exe)],check=True)
    primary=json.loads((DYN/'RADIAL_RESULTS.json').read_text())['cases']
    extra=json.loads((DYN/'EXCITATION_RESULTS.json').read_text())['cases'];cases={**primary,**extra}
    archives={}
    for name in ('RADIAL_DATA.npz','EXCITATION_DATA.npz'):
        with np.load(DYN/name) as d:
            for k in d.files:archives[k]=d[k]
    result=[]
    for tag in ('radial_main','radial_fine','pulse10','pulse20','pulse40_fine'):
        par=cases[tag]['params'];h,dt,L,ab,amp=[par[k] for k in ('h','dt','L','ab','amp')]
        n=round(L/h);r=(np.arange(n)+.5)*h;f,s=profile(r);shape=np.exp(-(r/12)**2)
        f*=1+amp*shape;s+=amp*(1-s)*shape
        vol=4*np.pi/3*((r+h/2)**3-(r-h/2)**3)
        py=-OMEGA*f;py*=2000/np.sum(-f*py*vol)
        d=WORK/'source'/tag;d.mkdir(parents=True,exist_ok=True)
        initial=d/'initial.csv';np.savetxt(initial,np.c_[r,f,np.zeros_like(r),s,np.zeros_like(r),py,np.zeros_like(r)],
                delimiter=',',header='rho,x,y,s,px,py,ps',comments='')
        subprocess.run([str(exe),str(initial),str(d/'state'),str(h),str(dt),str(L),'1000','60',str(ab),'.2'],check=True,stdout=subprocess.DEVNULL)
        hist=np.genfromtxt(d/'state_history.csv',names=True,delimiter=',');ref=archives[tag+'_history']
        err=max(float(np.max(abs(hist[n]-ref[n])))/max(1,float(np.max(abs(ref[n])))) for n in hist.dtype.names)
        result.append({'tag':tag,'archived_history_max_scaled_difference':err,
                       'regenerated_trace_sha256':hashlib.sha256((d/'state_trace.csv').read_bytes()).hexdigest()})
        if err>1e-10:raise RuntimeError('Published history not reproduced: '+tag)
    save('SOURCE_REPRODUCTION.json',{'rows':result,'previous_scratch_trace_truncated':'pulse20',
                                   'original_inputs_not_modified':True});print(json.dumps(result,indent=2))

if __name__=='__main__':main()
