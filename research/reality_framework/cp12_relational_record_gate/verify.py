"""Reproducible small checks; honest distinction from saved larger results."""
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse,hashlib,json,os,platform,subprocess,sys,time
import numpy as np
import scipy
from analyze import catalog,UP

P=Path(__file__).resolve().parent


def sha(f):return hashlib.sha256(f.read_bytes()).hexdigest()


def load(name):return json.loads((P/name).read_text())


def main(record=False):
    start=time.monotonic()
    dependencies=[P.parent/'cp12_ground_connectivity'/name for name in
                  ['frozen_model.py','invariant_evolution.py','exact_destination.py']]
    dependencies += [P.parent/'cp12_unmasked_formation'/'run.py',
                     P.parent/'cp12_distance_history'/'evolution.py',UP/'dynamics.py',
                     P.parent/'cp12_exchangeability_selection'/'REPORT.md']
    inputs=[UP/'MANIFEST.json']
    for f in sorted(P.glob('RELATIONAL_*.json')):
        d=json.loads(f.read_text());src=UP/d['input_trajectory']
        assert sha(src)==d['input_trajectory_sha256'];inputs.append(src)
    files=sorted(f for f in P.iterdir() if f.is_file() and f.name not in ['MANIFEST.json','VERIFICATION.json'])
    manifest=dict(packet_hashes={f.name:sha(f) for f in files},
        dependency_hashes={str(f.relative_to(P.parent)):sha(f) for f in dependencies},
        source_trajectory_hashes={f.name:sha(f) for f in sorted(set(inputs))})
    if not record:assert load('MANIFEST.json')==manifest
    # Independently respect the previous packet's dependency attestations.
    upstream=json.loads((UP/'MANIFEST.json').read_text())
    for rel,h in upstream['source_hashes'].items():
        assert sha(P.parent/rel)==h,rel
    fresh=[]
    with TemporaryDirectory() as temp:
        root=Path(temp);cache=root/'cache';cache.mkdir()
        env={**os.environ,'CP12_FORMATION_CACHE':str(cache),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
        for module,call in [('checks','run()'),('analyze','analyze(4,4,256,[0.,1.,3.,5.,10.])'),('energy_control','run()')]:
            dest=root/module;dest.mkdir()
            code=f'from pathlib import Path; import {module} as m; m.P=Path({str(dest)!r}); m.{call}'
            result=subprocess.run([sys.executable,'-c',code],cwd=P,env=env,capture_output=True,text=True)
            assert result.returncode==0,(module,result.stdout[-1500:],result.stderr[-2500:])
            if module=='analyze':
                now=json.loads((dest/'RELATIONAL_N4_Q4_B256.json').read_text())
                saved=load('RELATIONAL_N4_Q4_B256.json')
                for a,b in zip(now['rows'],saved['rows']):
                    assert max(abs(a['means'][k]-b['means'][k]) for k in a['means'])<1e-10
                    for label in a['controls']:
                        assert abs(a['controls'][label]['core_shape_TV']-b['controls'][label]['core_shape_TV'])<1e-10
            if module=='energy_control':
                now=json.loads((dest/'MEAN_ENERGY_CONTROL.json').read_text());saved=load('MEAN_ENERGY_CONTROL.json')
                assert max(abs(now['means'][k]-saved['means'][k]) for k in now['means'])<1e-10
                assert abs(now['maximum_entropy_reference']['means']['has_cycle']-saved['maximum_entropy_reference']['means']['has_cycle'])<1e-10
            fresh.append(dict(script=module+'.py',returncode=0,stdout_sha256=hashlib.sha256(result.stdout.encode()).hexdigest()))
            print(module,'PASS',flush=True)
    summaries=[]
    for f in sorted(P.glob('RELATIONAL_*.json')):
        d=json.loads(f.read_text());n=d['N'];q=d['Q'];C=catalog(n)
        assert d['core_shape_canonical_masks']==C['keys']
        src=json.loads((UP/d['input_trajectory']).read_text())
        assert d['maximum_replay_count_difference']<1e-9
        assert d['maximum_replay_norm_drift']<1e-10 and d['maximum_replay_energy_drift']<1e-8
        for row in d['rows']:
            ref=min(src['rows'],key=lambda a:abs(a['T']-row['T']))
            assert abs(ref['T']-row['T'])<1e-9
            assert row['state_error_upper']==ref['state_error_upper']
            e=row['state_error_upper'];assert abs(row['bounded_operator_error_upper']-(2*e+e*e))<1e-13
            hist=np.array(row['core_shape_distribution'])
            assert abs(hist.sum()-1)<1e-10 and np.min(hist)>-1e-12
            expect=hist@C['scores']
            assert max(abs(row['means'][k]-x) for k,x in zip(d['metrics'],expect))<1e-11
            count=np.array(row['protected_edge_count_distribution'])
            meanp=float(count@np.arange(len(count)))
            assert abs(meanp-q*(1-row['removable_fraction']))<1e-9
            assert abs(row['means']['ball_1_fraction']-(1+2*meanp/n)/n)<1e-11
            for label,control in row['controls'].items():
                h=np.array(control['core_shape_distribution'])
                assert abs(h.sum()-1)<1e-10 and np.min(h)>-1e-12
                assert abs(control['core_shape_TV']-np.sum(abs(hist-h))/2)<1e-12
                expected=h@C['scores']
                for key,x in zip(d['metrics'],expected):
                    assert abs(control['means'][key]-x)<1e-11
                    assert abs(control['contrasts'][key]-(row['means'][key]-x))<1e-11
                assert abs(control['contrasts']['ball_1_fraction'])<1e-11
            assert abs(sum(row['shortcut']['contrasts_per_Q']))<1e-10
            assert abs(sum(row['shortcut']['neutral_edges_per_Q'])-row['removable_fraction'])<1e-10
            a=.1*(C['M']-q)/C['M']
            assert row['matched_resource_diagonal_mixture_energy']>=q/2+a*q-1e-10
        if d['budget'] in [256,128] and n==q:
            final=d['rows'][-1];tv=final['controls']['gauss_parent_resources']['core_shape_TV']
            assert tv>final['bounded_operator_error_upper']
            summaries.append(dict(N=n,Q=q,budget=d['budget'],T=final['T'],
                gauss_parent_core_TV=tv,core_TV_lower=max(0.,tv-final['bounded_operator_error_upper']),
                core_connected_probability=final['means']['connected']))
    low=load('RELATIONAL_N4_Q4_B128.json');high=load('RELATIONAL_N4_Q4_B256.json')
    cutoffdiff=max(abs(a['controls']['gauss_parent_resources']['core_shape_TV']-b['controls']['gauss_parent_resources']['core_shape_TV']) for a,b in zip(low['rows'],high['rows']))
    assert cutoffdiff<1e-9
    dense=load('REFINED_DENSE_GEOMETRY.json');source=load('RELATIONAL_N6_Q12_B10.json')['rows'][-1]
    for row in dense['rows']:
        err=2*dense['forward_error_upper']*row['backward_error_upper']+dense['forward_error_upper']**2
        assert abs(err-row['expectation_error_upper'])<1e-13
    signs=np.array(dense['TV_witness_class_signs'])
    a=np.array(source['core_shape_distribution']);b=np.array(source['controls']['gauss_parent_resources']['core_shape_distribution'])
    witness=float(signs@(a-b)/2)
    assert abs(witness-dense['projected_TV'])<1e-12
    assert dense['rows'][-1]['interval'][0]>.0144
    assert dense['rows'][1]['interval'][1]<.163
    mechanism=load('MECHANISM.json')
    for d in mechanism['models']:
        assert d['rewiring_core_shape_commutator_residual']<1e-12
        M=d['N']*(d['N']-1)/2;gap=.1*d['Q']*(M-d['Q'])/M
        assert abs(d['configuration_diagonal_energy_gap']-gap)<1e-12
    shortcut=load('SHORTCUT_CONTROLS.json')
    for m in shortcut['models']:
        for row in m['rows']:
            for control in row['controls'].values():
                assert abs(sum(control['contrasts_per_Q']))<1e-10
    for m in shortcut['models']:
        if m['N']==m['Q'] and m['N'] in [5,6]:
            last=m['rows'][-1]
            assert last['controls']['gauss_parent_resources']['contrasts_per_Q'][2]>last['bounded_score_or_contrast_error_upper']
    out=dict(status='CP12_RELATIONAL_RECORD_GATE_VERIFIED',fresh=fresh,
        saved_larger_results_hash_validated=True,larger_replay_or_dense_dual_silently_rerun=False,
        finite_relational_summaries=summaries,N4_core_TV_cutoff_comparison_max_difference=cutoffdiff,
        dense_core_TV_witness_lower=dense['rows'][-1]['interval'][0],
        ruling='Finite relational shape correlations differ from the specified controls. Extended persistent local geometry, a causal cone and spacetime formation are not established. This is not a theory-wide no-go result.',
        qualification='Analytic full-law cutoff and operator bounds evaluated in ordinary floating point. No formal interval-arithmetic certificate for roundoff.',
        runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        seconds=time.monotonic()-start)
    if record:(P/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    assert load('MANIFEST.json')==manifest
    (P/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True);return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--record-manifest',action='store_true')
    main(ap.parse_args().record_manifest)
