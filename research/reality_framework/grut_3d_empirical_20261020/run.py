"""Build, freeze, and run the Cartesian campaign; sequential to limit memory.

python run.py preflight | freeze | run [case IDs ...]
No fits to bulk frequencies and no tuning against the empirical target.
"""
from pathlib import Path
import hashlib, json, os, subprocess, sys, time
import numpy as np

HERE=Path(__file__).resolve().parent
RAW=HERE/'raw'; RAW.mkdir(exist_ok=True)
BUILD=Path(os.environ.get('GRUT_BUILD_DIR','/workspace/scratch/85cb4751a782/grut3d_build'))
BUILD.mkdir(exist_ok=True)
EXE=BUILD/'cartesian3d'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def build():
    cmd=['g++','-O3','-march=native','-fopenmp','-ffp-contract=off','-std=c++17','-Wall','-Wextra',str(HERE/'cartesian.cpp'),'-o',str(EXE)]
    subprocess.run(cmd,check=True)
    return cmd
def command(c,out):
    opts={'--case':c['case'],'--h':c['h'],'--L':c['L'],'--Lz':c.get('Lz',c['L']),'--T':c.get('T',160),'--dt':c.get('dt',.05),'--sample':c.get('sample',5),'--eps':c.get('eps',0),'--angle':c.get('angle',0),'--sym':c.get('sym','none'),'--threads':c.get('threads',8),'--profile':HERE/'profile.dat','--modes':HERE/'modes.dat','--out':out}
    opts.update(c.get('extra',{}))
    return [str(EXE)]+[str(z) for item in opts.items() for z in item]
def execute(c):
    out=RAW/(c['id']+'.csv');cmd=command(c,out)
    print('START',c['id'],flush=True);begin=time.time()
    p=subprocess.run(cmd,text=True,capture_output=True)
    rec={'case':c,'command':cmd,'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr,'wall_seconds':time.time()-begin,'source_sha256':sha(HERE/'cartesian.cpp')}
    (RAW/(c['id']+'.execution.json')).write_text(json.dumps(rec,indent=2,default=str)+'\n')
    print('END',c['id'],p.returncode,p.stdout.strip(),flush=True)
    if p.returncode: raise RuntimeError(p.stderr)
    return out
def preflight():
    compile_command=build()
    base={'case':'single','h':.75,'L':12,'T':3,'sample':3,'dt':.05,'threads':2}
    u=execute(dict(base,id='control_u1',sym='u1'))
    x=np.genfromtxt(u,delimiter=',',names=True)
    results={'compile_command':compile_command,'global_phase_relative_residual':float(max(x['res0_norm']/x['gen0_norm']))}
    f=execute(dict(base,id='control_fd0',sym='fd',extra={'--state-out':BUILD/'fd0.bin'}))
    z=np.fromfile(BUILD/'fd0.bin',dtype=np.float64).reshape(2,6,-1);errors=[]
    for d in [.002,.001]:
        for s in [-1,1]: execute(dict(base,id=f'control_fd_{d}_{s}',sym='fd',extra={'--fd-offset':s*d,'--state-out':BUILD/f'fd{s}.bin'}))
        difference=(np.fromfile(BUILD/'fd1.bin',dtype=np.float64).reshape(2,6,-1)[0]-np.fromfile(BUILD/'fd-1.bin',dtype=np.float64).reshape(2,6,-1)[0])/(2*d)
        errors.append(float(np.linalg.norm(difference-z[1])/np.linalg.norm(z[1])))
    results.update(tangent_fd_errors=errors,tangent_fd_error_reduction=errors[0]/errors[1],source_sha256=sha(HERE/'cartesian.cpp'))
    # Rectangular indexing agrees with cube at the center for a causal short run.
    cube=execute(dict(base,id='control_cube',L=18))
    rect=execute(dict(base,id='control_rect',L=18,Lz=24))
    a=np.genfromtxt(cube,delimiter=',',names=True);b=np.genfromtxt(rect,delimiter=',',names=True)
    results['rectangular_core_Q8_relative_difference']=float(abs(a[-1]['c0_Q8']-b[-1]['c0_Q8'])/a[-1]['c0_Q8'])
    results['checks']={'global_phase':results['global_phase_relative_residual']<1e-10,'tangent_linearization':errors[1]<1e-5 and 3.5<errors[0]/errors[1]<4.5,'rectangular_indexing':results['rectangular_core_Q8_relative_difference']<1e-7}
    (HERE/'PREFLIGHT.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2),flush=True)
    if not all(results['checks'].values()):raise RuntimeError('preflight failed')
def freeze():
    if (HERE/'protocol.json').exists():raise RuntimeError('freeze already exists; use additive amendments')
    rng=np.random.default_rng(20261020)
    modes=[]
    for i in range(-6,7):
      for j in range(-6,7):
       for k in range(-6,7):
        t=(i,j,k)
        if t>(0,0,0) and .2<=np.linalg.norm(t)*np.pi/32<=.6:modes.append(t)
    modes=np.array(modes)[rng.choice(len(modes),24,replace=False)]
    amps=rng.normal(size=(24,3));amps/=np.sqrt(np.sum(amps**2,axis=0)/2)
    np.savetxt(HERE/'modes.dat',np.column_stack([modes,amps]),fmt='%.17g')
    cases=[]
    for tag,h in [('coarse',.75),('medium',.5),('fine',.375)]:
      for eps,name in [(0,'reference'),(.05,'generic')]:
        cases.append(dict(id=f'single_{name}_{tag}',case='single',h=h,L=18,T=160,eps=eps,sym='space' if eps else 'none'))
      cases.append(dict(id=f'pair_reference_{tag}',case='pair',h=h,L=24,Lz=48,T=160,sym='all'))
      cases.append(dict(id=f'pair_generic_{tag}',case='pair',h=h,L=24,Lz=48,T=160,eps=.05))
    cases.extend([
      dict(id='single_negative_fine',case='single',h=.375,L=18,T=160,eps=-.05),
      dict(id='single_generic_box',case='single',h=.5,L=27,T=160,eps=.05,sym='space'),
      dict(id='single_generic_dt',case='single',h=.5,L=18,T=160,dt=.025,eps=.05,sym='space'),
      dict(id='pair_reference_box',case='pair',h=.5,L=30,Lz=60,T=160,sym='all'),
      dict(id='pair_reference_dt',case='pair',h=.5,L=24,Lz=48,T=160,dt=.025,sym='all'),
      dict(id='pair_rotated_medium',case='pair',h=.5,L=24,Lz=48,T=160,angle=.3),
      dict(id='pair_rotated_fine',case='pair',h=.375,L=24,Lz=48,T=160,angle=.3),
      dict(id='pair_negative_fine',case='pair',h=.375,L=24,Lz=48,T=160,eps=-.05)])
    for tag,h in [('coarse',1),('medium',2/3),('fine',.5)]:
      cases.append(dict(id=f'production_{tag}',case='condensate',h=h,L=32,T=240,eps=.005,sample=10))
    spec={'local_freeze_unix':time.time(),'public_preregistration':False,'base_commit':'9ae004a06868674487bca679697bb073a10b7c75','scope':'Leading classical D=0, flat local microphysics; full Cartesian 3D. Finite perturbation family and finite observation interval, not all-perturbation or eternal stability.','pilot_disclosure':'Development pilots T=3 or 6 only, U(1), tangent finite differences, rectangular indexing and throughput; no T>6 object outcome inspected before this freeze.','cases':cases,'checks_before_science':json.loads((HERE/'PREFLIGHT.json').read_text()),'criteria':{'single_retention':'Every sampled tau=120..160: dominant R12 ball Q >= .90 initial total Q; local E/abs(Q)<1; one dominant 2%-of-peak connected component at endpoint. Report R8 and 5/10% sensitivity.','pair_retention':'Every sampled tau=120..160: two disjoint R12 balls contain >= .90 initial total Q in sum and each >= .15 initial total Q. Report each laboratory E/Q; this is sufficient but not necessary for binding.','spatial_convergence':'Report actual three-grid differences, inferred order where meaningful, finest change <= 2% of initial total Q for late core charges and <=5% for phase-space tangent norm/initial norm. Require decreasing differences; no pass for nonmonotonic unresolved sequences.','symmetry':'R36 canonical translation and rotation tangent residual relative to instantaneous generator: report every grid/time, seek decreasing O(h^2) residual. Never demand a symmetry tangent have constant positive norm. Full event rotation .3 rad must approach same core charges on refinement.','boundary_time':'Larger box and halved timestep changes < 1% initial Q in late core charge and <5% for normalized relative-shift tangent norm; otherwise continuum conclusion limited.','conservation':'max relative total Q error <1e-9 and Hamiltonian drift <1e-3. Boundary is conservative fixed vacuum, no sponge.','production':'Periodic mixed branch F=.5,S=sqrt(.75),Omega=sqrt(.775), zero-mean fully 3D perturbations, Qbar/volume exactly .25 sqrt(.775); count endpoint components at thresholds .02,.05,.10. Fragmentation does not certify a frozen cosmological mass function.'},'empirical_point':{'a':1,'b':.1,'gX':.1,'mX_GeV':100,'v_GeV':1000,'primordial_charge_yield':'YQ=YB inferred from D+3nu BBN, not fitted to Planck dark density; an additional initial-state postulate.','entropy':'No later SM entropy injection; hidden U1 exact in adopted EFT.','observable':'Lower bound on present hidden matter omega=Omega h^2, compared with Planck 2018 TT,TE,EE+lowE+lensing base-LambdaCDM omega_c.','blindness':'Published target already known; held out from parameter calibration, not an unseen/blind or prospective prediction.'},'files_sha256':{p.name:sha(p) for p in [HERE/'cartesian.cpp',HERE/'run.py',HERE/'profile.dat',HERE/'modes.dat',HERE/'PREFLIGHT.json']}}
    (HERE/'protocol.json').write_text(json.dumps(spec,indent=2)+'\n');print('Frozen',len(cases),'cases',flush=True)
if __name__=='__main__':
    os.chdir(HERE)
    if sys.argv[1]=='preflight':preflight()
    elif sys.argv[1]=='freeze':freeze()
    elif sys.argv[1]=='run':
      spec=json.loads((HERE/'protocol.json').read_text())
      for name,digest in spec['files_sha256'].items():
        if sha(HERE/name)!=digest:raise RuntimeError('frozen source changed '+name)
      build()
      requested=set(sys.argv[2:])
      for c in spec['cases']:
        if requested and c['id'] not in requested:continue
        if (RAW/(c['id']+'.execution.json')).exists():raise RuntimeError('do not overwrite '+c['id'])
        execute(c)
    else:raise SystemExit('preflight | freeze | run [ids]')
