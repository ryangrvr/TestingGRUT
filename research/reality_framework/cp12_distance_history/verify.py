"""Fresh controls and explicit validation of separately saved large runs."""
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse,hashlib,json,os,platform,subprocess,sys,time
import numpy as np
import scipy

P=Path(__file__).resolve().parent

def main(record):
    start=time.monotonic()
    sources=list(sorted(P.glob('*.py')))
    sources += [P.parent/'cp12_ground_connectivity'/f for f in ['frozen_model.py','invariant_evolution.py','exact_destination.py']]
    source_hashes={str(f.relative_to(P.parent)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sources}
    files=sorted(P.glob('N*_Q*_B*_W*_T*.json'));assert len(files)==6
    record_hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in files}
    m=dict(source_hashes=source_hashes,evolution_record_hashes=record_hashes)
    if record:(P/'MANIFEST.json').write_text(json.dumps(m,indent=2)+'\n')
    assert json.loads((P/'MANIFEST.json').read_text())==m
    fresh=[]
    with TemporaryDirectory() as cache:
        env={**os.environ,'CP12_FORMATION_CACHE':cache,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
        for script in [P.parent/'cp12_response_error_gate'/'verify.py',P/'check.py',P/'edge_speed.py',P/'attribution.py',P/'pulse.py']:
            r=subprocess.run([sys.executable,str(script)],capture_output=True,text=True,env=env)
            assert r.returncode==0,(str(script),r.stdout[-1800:],r.stderr[-1800:])
            fresh.append(dict(script=str(script.relative_to(P.parent)),returncode=0,stdout_sha256=hashlib.sha256(r.stdout.encode()).hexdigest()))
            print(str(script.relative_to(P.parent)),'PASS',flush=True)
    summary=[]
    for f in files:
        a=json.loads(f.read_text());n=a['N'];q=a['Q'];cs=a['classes'];C=np.array(a['joint_responses']);Q=np.array(a['geometry_transitions'])
        pr=np.array([r['probability'] for r in a['past']]);assert abs(pr.sum()-1)<1e-10
        assert q==n-1 and a['norm_error']<1e-10
        assert np.max(abs(Q.sum(axis=1)-pr))<1e-10
        assert np.max(abs(C.sum(axis=1)-np.array([r['value'] for r in a['past']])))<1e-12
        assert np.max(abs(C.sum(axis=0)-np.array([r['value'] for r in a['post_recorded'][:-1]])))<1e-12
        assert abs(C.sum()-a['post_recorded'][-1]['value'])<1e-12
        sourceQ=np.array(a['source_filtered_geometry_transitions'])
        assert np.max(abs(sourceQ.sum(axis=1)-pr*np.array([r['source_density_at_wait'] for r in a['past']])))<1e-10
        drift=np.array([r['value'] for r in a['past']])-np.array([r['value'] for r in a['post_recorded'][:-1]])
        assert np.all(abs(drift)<=np.array(a['past_post_difference_upper_from_transition_weights'])+1e-11)
        assert abs(-np.sum(pr*np.log2(pr))-a['recording_entropy_bits'])<1e-10
        assert a.get('coherent_branch_sum_error',0)<1e-10
        for row in a['past']+a['post_recorded']+a['post_unrecorded']:
            e=2*(row['forward_y']*row['dual_y']+row['forward_x']*row['dual_x']+row['forward_x']*row['forward_y'])
            assert abs(e-row['error'])<1e-14
            assert abs(row['interval'][0]-(row['value']-e))<1e-14
            assert abs(row['interval'][1]-(row['value']+e))<1e-14
        far=next(r for r in a['past'] if r['disconnected'])
        assert abs(far['value'])>far['error']
        M=n*(n-1)//2;record_work_error=2*.1*q*(M-q)/M*(2*a['wait_error']+a['wait_error']**2)
        summary.append(dict(file=f.name,N=n,lag=a['lag'],disconnected_response=far['value'],analytic_truncation_bound=far['error'],conditional_interval=far['conditional_interval'],record_work_error_upper=record_work_error))
    s=json.loads((P/'STATE_COUPLING.json').read_text())
    for model in s['models']:
        for row in model['records']:
            assert row['coherent_cross_norm_square']<=row['global_operator_norm_square']+1e-10
            assert row['coherent_cross_norm_square_interval'][0]<=row['coherent_cross_norm_square']<=row['coherent_cross_norm_square_interval'][1]
    speed=json.loads((P/'EDGE_SPEED.json').read_text());assert len(speed['incidence_controls'])==9 and len(speed['history_speed_checks'])==6
    out=dict(status='CP12_DISTANCE_HISTORY_CHECKS_PASS',fresh=fresh,finite_size_history_summary=summary,
             large_evolutions_rerun_by_this_verifier=False,state_coupling_evolutions_rerun_by_this_verifier=False,
             floating_point_qualification='Analytically justified error bounds numerically evaluated; no interval-arithmetic certificate of roundoff.',
             locality_ruling='Resolved rewiring-mediated response between initially disconnected components in the tested auxiliary preparation. No growing intrinsic causal cone or state-selected spacetime established.',
             pool_norm_scope='Exact full-Hamiltonian collective capacity; not a proof every state realizes it or a uniform finite-time large-N extrapolation.',
             runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),seconds=time.monotonic()-start)
    (P/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--record-manifest',action='store_true');a=ap.parse_args();main(a.record_manifest)
