"""Fresh small controls; validate saved expensive-response bounds separately."""
from pathlib import Path
import argparse,hashlib,json,os,platform,subprocess,sys,time
import numpy as np
import scipy

P=Path(__file__).resolve().parent

def manifest():
    paths=[P/f for f in ['control.py','probe.py','check.py','verify.py','summarize.py']]
    paths+=[P.parent/'cp12_unmasked_formation'/f for f in ['run.py','check.py','verify.py']]
    paths+=[P.parent/'cp12_ground_connectivity'/f for f in ['frozen_model.py','invariant_evolution.py','exact_destination.py']]
    return {str(f.relative_to(P.parent)):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths}

def main(record):
    start=time.monotonic();m=manifest()
    if record:(P/'MANIFEST.json').write_text(json.dumps(m,indent=2)+'\n')
    assert json.loads((P/'MANIFEST.json').read_text())==m
    fresh=[]
    for script in [P.parent/'cp12_unmasked_formation'/'verify.py',P/'check.py']:
        r=subprocess.run([sys.executable,str(script)],capture_output=True,text=True,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
        assert r.returncode==0,(r.stdout[-1500:],r.stderr[-1500:])
        fresh.append(dict(script=str(script.relative_to(P.parent)),returncode=0,stdout_sha256=hashlib.sha256(r.stdout.encode()).hexdigest()))
        print(str(script.relative_to(P.parent)),'PASS',flush=True)
    records=[]
    for n,b in [(5,10),(6,6)]:
        a=json.loads((P/f'N{n}_B{b}_S64.json').read_text())
        assert a['initial_energy']==n and a['initial_variance']==4*n
        assert a['wait']==.1 and a['lag']==.1 and len(a['responses'])==3
        rs=a['responses'];assert abs(rs[0]['projected_response']+rs[1]['projected_response']-rs[2]['projected_response'])<1e-12
        for row in rs:
            ex=row['forward_x_error_upper'];ey=row['forward_y_error_upper']
            bound=2*(ey*row['dual_y_error_upper']+ex*row['dual_x_error_upper']+ex*ey)
            assert row['detector_norm']==1.
            assert abs(bound-row['observable_absolute_error_upper'])<1e-14
            assert bound<row['global_absolute_error_upper']
            assert row['sign_resolved']==(abs(row['projected_response'])>bound)
        assert a['numerical_norm_error']<1e-10 and a['numerical_baseline_energy_error']<1e-9
        records.append(dict(N=n,resolved_channels=[row['channel'] for row in rs if row['sign_resolved']]))
    checks=json.loads((P/'CHECKS.json').read_text())
    assert len(checks['dense'])==24 and len(checks['physical'])==3
    assert checks['hostile_control']['applicability_rejection_required']
    out=dict(status='CP12_OBSERVABLE_ERROR_CHECKS_PASS',fresh=fresh,saved_expensive_runs=records,expensive_runs_rerun_by_this_verifier=False,
             bounds_evaluated_in_ordinary_floating_point=True,
             runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
             formation_gate='Not established; diameter-two controls do not give extended volume growth, a causal cone, or autonomous state selection.',seconds=time.monotonic()-start)
    (P/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--record-manifest',action='store_true');a=ap.parse_args();main(a.record_manifest)
