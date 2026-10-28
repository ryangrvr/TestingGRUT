"""Archive compact exact traces and visualization projections from raw histories."""
from common import *
import platform

def main():
    from run import CASES
    runs=json.loads((ROOT/'RUNS.json').read_text())['cases']
    if set(runs)!=set(CASES):raise RuntimeError('Production histories incomplete')
    traces={};projections={};manifest={}
    for tag,r in runs.items():
        raw=WORK/tag/'state.tsv';data=np.genfromtxt(raw,names=True,delimiter='\t')
        traces[tag]=data
        paths=[WORK/tag/'initial.bin',raw]+[Path(x) for x in r['snapshots']]
        manifest[tag]=[{'name':p.name,'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in paths]
        if tag not in ('cold64','ellipse64','dilute64','object96_perturbed'):continue
        N,L,t,a=snapshot(r['snapshots'][-1]);q=np.asarray(a[:3]);p=np.asarray(a[3:]);j=q[0]*p[1]-q[1]*p[0]
        projections[tag+'_charge_max_projection']=j.max(axis=2)
        projections[tag+'_phi_min_projection_over_v']=q[2].min(axis=2)/V
        projections[tag+'_geometry']=np.array([N,L,t])
        if tag=='object96_perturbed':
            axis=(np.arange(N)-N//2)*L/N
            rr=np.sqrt(axis[:,None,None]**2+axis[None,:,None]**2+axis[None,None,:]**2)
            edges=np.linspace(0,24,49);bins=np.digitize(rr.ravel(),edges)-1
            f=np.sqrt(q[0]**2+q[1]**2)/V;s=q[2]/V
            ok=(bins>=0)&(bins<len(edges)-1);count=np.bincount(bins[ok],minlength=len(edges)-1)
            fbar=np.bincount(bins[ok],weights=f.ravel()[ok],minlength=len(edges)-1)/np.maximum(count,1)
            sbar=np.bincount(bins[ok],weights=s.ravel()[ok],minlength=len(edges)-1)/np.maximum(count,1)
            np.savetxt(ROOT/'OBJECT_FINAL_RADIAL.csv',np.c_[.5*(edges[:-1]+edges[1:]),fbar,sbar,count],
                       delimiter=',',header='radius_over_mX,F_amplitude_average,S_average,cells',comments='')
    np.savez_compressed(ROOT/'TRACE_DATA.npz',**traces)
    np.savez_compressed(ROOT/'FINAL_PROJECTIONS.npz',**projections)
    save('RAW_MANIFEST.json',{'status':'reproducible intermediate histories in scratch, not durable full-volume deliverables',
                             'files':manifest,'full_data_saved_in_git':False,
                             'compact_data':'TRACE_DATA.npz, FINAL_PROJECTIONS.npz and measured catalogs',
                             'reproduce':'code/run.py regenerates each listed initial state and history'})
    save('ENVIRONMENT.json',{'python':platform.python_version(),'platform':platform.platform(),
         'numpy':np.__version__,'compiler':subprocess.check_output(['g++','--version'],text=True).splitlines()[0],
         'production_runs':len(runs),'AI_subagents':0,'independent_review':False})
    print(json.dumps({'archived_runs':len(runs),'trace_bytes':(ROOT/'TRACE_DATA.npz').stat().st_size,
                      'projection_bytes':(ROOT/'FINAL_PROJECTIONS.npz').stat().st_size}))

if __name__=='__main__':main()
