"""Fresh controls plus explicit, limited validation of saved pilot records."""
from pathlib import Path
import argparse,hashlib,json,os,platform,subprocess,sys,time
import numpy as np
import scipy

P=Path(__file__).resolve().parent
DEP=P.parent/'cp12_ground_connectivity'
SOURCES=['run.py','check.py','verify.py']
DEPENDENCIES=['frozen_model.py','invariant_evolution.py','exact_destination.py']

def manifest():
    return dict(sources={f:hashlib.sha256((P/f).read_bytes()).hexdigest() for f in SOURCES},dependencies={f:hashlib.sha256((DEP/f).read_bytes()).hexdigest() for f in DEPENDENCIES})

def main(record):
    start=time.monotonic();current=manifest()
    if record:(P/'MANIFEST.json').write_text(json.dumps(current,indent=2)+'\n')
    assert json.loads((P/'MANIFEST.json').read_text())==current
    fresh=[]
    for script in [P.parent/'cp12_topology_sensitive_preparation'/'verify.py',P/'check.py']:
        result=subprocess.run([sys.executable,str(script)],capture_output=True,text=True,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'})
        assert result.returncode==0,(str(script),result.stdout[-2000:],result.stderr[-2000:])
        fresh.append(dict(script=str(script.relative_to(P.parent)),returncode=0,stdout_sha256=hashlib.sha256(result.stdout.encode()).hexdigest()))
        print(script.name,'PASS',flush=True)
    load=lambda f:json.loads((P/f).read_text())
    a=load('probe_N5_Q10_B10.json');b=load('probe_N6_Q12_B6.json');c=load('N7_Q14_B3.json');d=load('N6_Q12_B6.json')
    assert abs(a['response_by_final_distance']['1'])>a['each_response_absolute_error_upper']
    assert abs(b['response_by_final_distance']['2'])<b['each_response_absolute_error_upper']
    assert c['rows'][-1]['state_error_upper']>1
    assert d['rows'][-1]['T']==.2 and d['rows'][-1]['state_error_upper']<.001
    assert max(row['energy_error'] for row in c['rows'])<1e-10
    assert max(row['norm_error'] for row in c['rows'])<1e-10
    assert abs(d['rows'][-1]['observables']['ball_2']-6)<1e-10
    out=dict(status='CP12_UNMASKED_BENCHMARK_CHECKS_PASS_FORMATION_GATE_NOT_PASSED',fresh=fresh,
             runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
             saved_pilot_checks=dict(five_event_control_sign_resolved=True,six_event_distance_two_error_contains_zero=True,seven_event_time_two_error_inconclusive=True,six_event_intrinsic_ball_saturation=True),
             expensive_pilots_rerun_by_this_verifier=False,
             prototype_records_without_later_joint_fields=['probe_N5_Q10_B8.json','probe_N6_Q12_B6.json','N7_Q14_B3.json'],
             gate='Size and regulator controls unresolved; no autonomous preparation or spacetime formation derived.',seconds=time.monotonic()-start)
    (P/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--record-manifest',action='store_true');a=ap.parse_args();main(a.record_manifest)
