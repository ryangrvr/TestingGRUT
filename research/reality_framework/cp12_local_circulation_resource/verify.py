"""Fresh reproducibility gate; does not overwrite committed result files."""
from pathlib import Path
import argparse,hashlib,json,math,os,platform,subprocess,sys,tempfile,time
import numpy as np
import scipy

P=Path(__file__).resolve().parent
DEPENDENCIES={
 'cp12_ground_connectivity/frozen_model.py':'e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a',
 'cp12_ground_connectivity/exact_destination.py':'1a1ccb4f06602179f14c2834dea3566fada022ca6a23f58370abe0d4f1e5b35b',
 'cp12_holonomy_locality_gate/check.py':'77e89f7240638107f04d013d42cb61dd50fcbde62c637260a1e1122f16282263',
}

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def compare(a,b,path='root'):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys()==b.keys(),path
        for k in a:
            if k!='seconds':compare(a[k],b[k],path+'.'+k)
    elif isinstance(a,list):
        assert isinstance(b,list) and len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):compare(x,y,path+f'[{i}]')
    elif isinstance(a,float):
        assert isinstance(b,(int,float)) and math.isfinite(a) and math.isfinite(b),path
        assert abs(a-b)<=5e-10*max(1.,abs(a),abs(b)),(path,a,b)
    else:assert a==b,(path,a,b)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true');args=ap.parse_args()
    for name,h in DEPENDENCIES.items():assert digest(P.parent/name)==h,name
    manifest=P/'MANIFEST.json'
    if manifest.exists():
        for name,h in json.loads(manifest.read_text())['files_sha256'].items():assert digest(P/name)==h,name
    elif not args.record:raise RuntimeError('Missing MANIFEST.json')
    rows=[];start=time.monotonic()
    jobs=[(P.parent/'cp12_holonomy_locality_gate/check.py','CP12_HOLONOMY_LOCALITY_GRAPH_CHECK_PASS',None),
          (P/'check.py','CP12_LOCAL_CIRCULATION_GRAPH_CONTROLS_PASS','graph_results.json'),
          (P/'wilson_check.py','CP12_HIGHER_DENSITY_WILSON_MOMENTS_PASS','wilson_results.json'),
          (P/'quantum_check.py','CP12_LOCAL_CIRCULATION_QUANTUM_CONTROLS_PASS','quantum_results.json')]
    with tempfile.TemporaryDirectory(prefix='cp12-local-circulation-') as tmp:
        env=os.environ.copy();env.update(CP12_OUTPUT_DIR=tmp,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        for script,marker,result in jobs:
            t=time.monotonic();r=subprocess.run([sys.executable,'-I',str(script)],env=env,text=True,capture_output=True)
            if r.returncode or marker not in r.stdout:
                print(r.stdout);print(r.stderr,file=sys.stderr);raise AssertionError(script)
            if result:compare(json.loads((P/result).read_text()),json.loads((Path(tmp)/result).read_text()))
            row=dict(script=str(script.relative_to(P.parent)),returncode=r.returncode,marker=marker,
                     fresh_execution=True,saved_result_comparison=bool(result),
                     stdout_sha256=hashlib.sha256(r.stdout.encode()).hexdigest(),seconds=time.monotonic()-t)
            rows.append(row);print(json.dumps(row),flush=True)
    out=dict(status='CP12_LOCAL_CIRCULATION_RESOURCE_VERIFIED',fresh_checkers=len(rows),checks=rows,
             parent_commit='ae91c9beb0d8dda0f7692ea8bab0c5f808355530',
             dependency_sha256=DEPENDENCIES,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,
             floating_record_tolerance='5e-10 * max(1,abs(a),abs(b)); elapsed time excluded',
             seconds=time.monotonic()-start,
             scope='Fresh numerical/algebraic checks. Analytical infinite-rotor and large-N proof requires mathematical audit.')
    if args.record:(P/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
