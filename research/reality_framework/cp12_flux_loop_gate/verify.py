"""Fresh independent small checks plus explicit larger-result validation."""
from pathlib import Path
from tempfile import TemporaryDirectory
from fractions import Fraction as F
import argparse,hashlib,json,os,platform,subprocess,sys,time
import numpy as np
import scipy
from loops import P,REL,UP
from onset_uniform import formulas


def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()


def load(name):return json.loads((P/name).read_text())


def main(record=False):
    start=time.monotonic()
    dependencies=[REL/'analyze.py',P.parent/'cp12_unmasked_formation'/'run.py',
        UP/'dynamics.py',P.parent/'cp12_distance_history'/'evolution.py']
    dependencies += [P.parent/'cp12_ground_connectivity'/name for name in
                     ['frozen_model.py','invariant_evolution.py','exact_destination.py']]
    dependencies += [P.parent/'cp12_protected_loop_onset'/name for name in ['REPORT.md','check.py']]
    inputs=[REL/'MANIFEST.json',UP/'MANIFEST.json',REL/'MEAN_ENERGY_CONTROL.json',UP/'N4_Q4_B256_T20.json']
    inputs += list(sorted(REL.glob('RELATIONAL_N*_B*.json')))
    files=sorted(f for f in P.iterdir() if f.is_file() and f.name not in ['MANIFEST.json','VERIFICATION.json'])
    manifest=dict(packet_hashes={f.name:sha(f) for f in files},
        dependency_hashes={str(f.relative_to(P.parent)):sha(f) for f in dependencies},
        input_hashes={str(f.relative_to(P.parent)):sha(f) for f in inputs})
    if not record:assert load('MANIFEST.json')==manifest
    previous=json.loads((REL/'MANIFEST.json').read_text())
    for name,h in previous['packet_hashes'].items():assert sha(REL/name)==h,name
    for name,h in previous['dependency_hashes'].items():assert sha(P.parent/name)==h,name
    fresh=[]
    with TemporaryDirectory() as root:
        root=Path(root);cache=root/'cache';cache.mkdir();snapshots=root/'snapshots';snapshots.mkdir()
        env={**os.environ,'CP12_FORMATION_CACHE':str(cache),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
        # Freshly build the small matrix, independent of the private cache.
        code='''from pathlib import Path
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply
from loops import model
sub,H,B,states,sizes,psi=model(4,4,256)
A=H-2*eye(H.shape[0],format='csr')
for t in [0.,1.,3.,5.,10.]:
 v=expm_multiply(-1j*t*A,psi,traceA=-1j*t*float(A.diagonal().sum())) if t else psi.copy()
 np.savez_compressed(Path(SNAPSHOT_ROOT)/f'N4_Q4_B256_T{t:g}.npz',wave=v)
'''.replace('SNAPSHOT_ROOT',repr(str(snapshots)))
        result=subprocess.run([sys.executable,'-c',code],cwd=P,env=env,capture_output=True,text=True)
        assert result.returncode==0,(result.stdout,result.stderr)
        fresh.append(dict(script='fresh N4 matrix and snapshot reconstruction',returncode=0))
        jobs=[('onset_exact','run()'),('controls','run()'),('references','run()'),('early_time','run()'),
              ('populations',f'run(4,4,256,[0.,1.,3.,5.,10.],{str(snapshots)!r})'),
              ('persistence',f'run({str(snapshots)!r},quick=True)'),
              ('measurement_controls',f'run({str(snapshots)!r})'),
              ('transport',f'run({str(snapshots)!r})')]
        for module,call in jobs:
            dest=root/module;dest.mkdir()
            code=f'from pathlib import Path;import {module} as m;m.P=Path({str(dest)!r});m.{call}'
            r=subprocess.run([sys.executable,'-c',code],cwd=P,env=env,capture_output=True,text=True)
            assert r.returncode==0,(module,r.stdout[-1000:],r.stderr[-3000:])
            if module=='onset_exact':
                d=json.loads((dest/'ONSET_EXACT.json').read_text())
                for key in ['protected_triangle_probability_leading_coefficient','occupied_flux_triangle_probability_leading_coefficient','hole_bearing_circulation_probability_T14_coefficient']:
                    assert d[key]==load('ONSET_EXACT.json')[key]
            if module=='references':
                d=json.loads((dest/'ENERGY_REFERENCES.json').read_text())
                for a,b in zip(d['rows'],load('ENERGY_REFERENCES.json')['rows']):
                    assert max(abs(a['values'][key]-b['values'][key]) for key in a['values'])<1e-10
            if module=='populations':
                d=json.loads((dest/'CHANNELS_N4_Q4_B256.json').read_text())
                for a,b in zip(d['rows'],load('CHANNELS_N4_Q4_B256.json')['rows']):
                    assert max(abs(a['values'][key]-b['values'][key]) for key in a['values'])<1e-10
            if module=='persistence':
                d=json.loads((dest/'PERSISTENCE_QUICK.json').read_text())['cases'][0]
                saved=next(c for c in load('PERSISTENCE.json')['cases'] if c['preparation']==d['preparation'])
                for a,b in zip(d['channels'],saved['channels']):
                    for row in a['rows']:
                        ref=min(b['rows'],key=lambda x:abs(x['lag']-row['lag']))
                        assert abs(row['destinations']['protected']['conditional_probability']-ref['destinations']['protected']['conditional_probability'])<1e-10
            if module=='transport':
                d=json.loads((dest/'CONDITIONAL_DENSITY.json').read_text())
                for a,b in zip(d['rows'],load('CONDITIONAL_DENSITY.json')['rows']):
                    assert abs(a['density_susceptibility']-b['density_susceptibility'])<1e-10
                    assert abs(a['density_susceptibility_conditioned_on_final_circulation']-b['density_susceptibility_conditioned_on_final_circulation'])<1e-10
            fresh.append(dict(script=module+'.py',returncode=0,stdout_sha256=hashlib.sha256(r.stdout.encode()).hexdigest()))
            print(module,'PASS',flush=True)
    # Validate all saved positive-score intervals and channel partitions.
    for file in sorted(P.glob('CHANNELS_N*.json')):
        d=json.loads(file.read_text());assert d['maximum_previous_record_difference']<1e-9
        for row in d['rows']:
            vals=row['values'];e=row['state_error_upper']
            assert abs(vals['protected_cycle_any']-vals['hole_only_cycle_event']-vals['flux_only_cycle_event']-vals['both_cycle_channels'])<1e-10
            assert vals['directed_flux_cycle_any']<=vals['electric_circulation_cycle_any']+1e-12<=vals['protected_cycle_any']+1e-10
            for key,p in vals.items():
                expected=[max(0.,np.sqrt(max(0.,p))-e)**2,min(1.,(np.sqrt(max(0.,p))+e)**2)]
                assert max(abs(a-b) for a,b in zip(expected,row['positive_score_intervals'][key]))<1e-14
            if d['Q']==d['N']:
                assert vals['electric_circulation_network_spanning']<=vals['full_graph_spanning_ring']+1e-12
    onset=load('ONSET_EXACT.json')
    assert F(onset['occupied_flux_triangle_probability_leading_coefficient'])==F(2,225)
    assert F(onset['vacuum_flux_triangle_geometric_loss_probability_T6_coefficient'])==F(1,3)
    for row in load('UNIFORM_ONSET.json')['rows']:
        expected=formulas(row['N'],row['Q'])
        assert row['exact_formulas']=={key:str(v) for key,v in expected.items()}
        assert max(abs(row['numerical_matrix_coefficients'][key]-float(v)) for key,v in expected.items())<1e-9
    for filename in ['PERSISTENCE.json','CIRCULATION_SURVIVAL.json']:
        for c in load(filename)['cases']:
            for channel in c['channels']:
                p0=channel['initial_probability']
                for row in channel['rows']:
                    ds=row['destinations']
                    assert abs(ds['protected']['sequential_joint_probability']-ds['hole']['sequential_joint_probability']-ds['occupied_flux']['sequential_joint_probability'])<1e-11
                    for data in ds.values():
                        assert -1e-10<=data['conditional_probability']<=1+1e-10
                        lo,hi=data['conditional_interval'];assert lo-1e-10<=data['conditional_probability']<=hi+1e-10
                    if 'electric_circulation' in ds:
                        assert ds['electric_circulation']['conditional_probability']<=ds['protected']['conditional_probability']+1e-10
                    change=row['noninterventional_projector_change']
                    assert -1e-10<=change<=1+1e-10
    response=load('CONDITIONAL_DENSITY.json')
    for c in response['rows']:
        lo,hi=c['density_susceptibility_interval'];assert lo*c['density_susceptibility']>0 and hi*c['density_susceptibility']>0
        a,b=c['density_susceptibility_conditioned_on_final_circulation_interval'];assert a*b>0
        assert abs(c['density_susceptibility']-c['commutator_vector_susceptibility'])<1e-10
        assert c['vacuum_flux_density_phase_response_exact']==0 and c['no_hopping_after_same_preparation_response_exact']==0
        assert min(p['same_triangle_protection_interval'][0] for p in c['paths'])>.993
        assert min(p['same_triangle_electric_circulation_interval'][0] for p in c['paths'])>.66
    dense=load('DENSE_CIRCULATION_BOUND.json');p=dense['projected_probability'];e=dense['forward_error_upper'];d=dense['backward_error_upper']
    assert abs(dense['positive_dual_interval'][0]-max(0.,p-e*d)**2/p)<1e-14
    assert dense['combined_interval'][0]==0
    for reference in load('ENERGY_REFERENCES.json')['rows']:assert abs(reference['full_mean_energy']-2)<1e-10
    audit=load('SOURCE_AUDIT.json');assert audit['onset_checker_returncode']==0 and 'PASS' in audit['onset_checker_stdout']
    out=dict(status='CP12_FLUX_LOOP_GATE_VERIFIED',fresh=fresh,
        saved_N5_N6_and_dense_dual_hash_validated=True,all_mixed_state_persistence_curves_silently_rerun=False,
        main_same_branch_density_interval=response['rows'][0]['density_susceptibility_interval'],
        main_final_circulation_conditioned_density_interval=response['rows'][0]['density_susceptibility_conditioned_on_final_circulation_interval'],
        ruling='Finite circulation, protection transfer and conditional compatibility with matter-density response established in the stated tests. Autonomous persistent three-dimensional geometry and causal propagation are not established.',
        qualification='Exact rational onset controls are unregulated. Dynamical full-law cutoff bounds are analytically justified and numerically evaluated in ordinary floating point, without a formal roundoff interval certificate.',
        runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        seconds=time.monotonic()-start)
    if record:(P/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    assert load('MANIFEST.json')==manifest
    (P/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True);return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--record-manifest',action='store_true')
    main(ap.parse_args().record_manifest)
