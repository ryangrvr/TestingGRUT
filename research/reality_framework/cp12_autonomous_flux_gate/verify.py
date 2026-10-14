"""Fresh small controls and explicit validation of separately saved evolutions."""
from pathlib import Path
from tempfile import TemporaryDirectory
import argparse,hashlib,json,math,os,platform,subprocess,sys,time
import numpy as np
import scipy
from energy_audit import MAIN,audit
from deletion_audit import audit as deletion_audit

P=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main(record=False):
    started=time.monotonic()
    sources=list(sorted(P.glob('*.py')))
    sources += [P.parent/'cp12_ground_connectivity'/f for f in ['frozen_model.py','invariant_evolution.py','exact_destination.py']]
    sources += [P.parent/'cp12_distance_history'/f for f in ['evolution.py','check.py','basis.py']]
    sources += [P.parent/'cp12_unmasked_formation'/'run.py',
                P.parent/'cp12_flux_protection_tradeoff'/'REPORT.md',
                P.parent/'cp12_flux_protection_tradeoff'/'check.py']
    records=sorted(P.glob('N*_Q*_B*_T*.json'))+sorted(P.glob('RECORDS_N*.json'))+[P/'MOBILITY.json',P/'ENERGY_AUDIT.json',P/'DELETION_AUDIT.json',P/'REFINED_COUNTS.json']
    manifest=dict(source_hashes={str(f.relative_to(P.parent)):sha(f) for f in sources},
                  evolution_and_resource_record_hashes={f.name:sha(f) for f in records})
    if not record:assert json.loads((P/'MANIFEST.json').read_text())==manifest
    fresh=[]
    with TemporaryDirectory() as root:
        root=Path(root);cache=root/'cache';cache.mkdir()
        env={**os.environ,'CP12_FORMATION_CACHE':str(cache),'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'1'}
        flux=P.parent/'cp12_flux_protection_tradeoff'/'check.py'
        result=subprocess.run([sys.executable,str(flux)],capture_output=True,text=True,env=env)
        assert result.returncode==0 and 'CP12_FLUX_PROTECTION_EXACT_LOCAL_IDENTITY_PASS' in result.stdout,(result.stdout,result.stderr)
        extension=dict(commit='4e2b7a876a3a45ed554408ed803b572dd9cab580',
                       checker_sha256=sha(flux),stdout=result.stdout,returncode=0)
        (P/'INDEPENDENT_EXTENSION_CHECK.json').write_text(json.dumps(extension,indent=2)+'\n')
        fresh.append(dict(script='cp12_flux_protection_tradeoff/check.py',returncode=0,stdout_sha256=hashlib.sha256(result.stdout.encode()).hexdigest()))
        print('Independent extension checker PASS',flush=True)
        jobs=[('exact_controls','main()'),('leakage_controls','main()'),('mobility','main()'),
              ('records','run(4,4,256,[0.,1.,3.,5.,10.])'),
              ('dynamics','run(4,4,256,20.,.1,8)')]
        for module,call in jobs:
            dest=root/module;dest.mkdir()
            code=f'from pathlib import Path; import {module} as m; m.P=Path({str(dest)!r}); m.{call}'
            result=subprocess.run([sys.executable,'-c',code],cwd=P,capture_output=True,text=True,env=env)
            assert result.returncode==0,(module,result.stdout[-1500:],result.stderr[-2500:])
            fresh.append(dict(script=module+'.py',returncode=0,stdout_sha256=hashlib.sha256(result.stdout.encode()).hexdigest()))
            if module=='mobility':
                actual=json.loads((dest/'MOBILITY.json').read_text());saved=json.loads((P/'MOBILITY.json').read_text())
                for a,b in zip(actual['rows'],saved['rows']):
                    assert abs(a['susceptibility']-b['susceptibility'])<1e-10
            if module=='records':
                actual=json.loads((dest/'RECORDS_N4_Q4_B256.json').read_text())
                saved=json.loads((P/'RECORDS_N4_Q4_B256.json').read_text())
                for a,b in zip(actual['rows'],saved['rows']):assert abs(a['graph_entropy_bits']-b['graph_entropy_bits'])<1e-9
            if module=='dynamics':
                actual=json.loads((dest/'N4_Q4_B256_T20.json').read_text())
                saved=json.loads((P/'N4_Q4_B256_T20.json').read_text())
                for a,b in zip(actual['rows'],saved['rows']):
                    for key in ['removable_fraction','flux_protected_fraction','graph_TV_from_initial']:
                        assert abs(a[key]-b[key])<1e-10
            print(module,'PASS',flush=True)
    summaries=[]
    for f in sorted(P.glob('N*_Q*_B*_T*.json')):
        d=json.loads(f.read_text());q=d['Q'];rs=d['rows'];errors=[]
        for r in rs:
            assert abs(r['removable_fraction']+r['flux_protected_fraction']+r['hole_protected_fraction']-1)<1e-8
            e=r['state_error_upper'];bound=min(1.,2*e+e*e)
            assert abs(bound-r['bounded_fraction_error_upper'])<1e-12
            assert r['norm_error']<1e-8 and r['energy_error_projected']<1e-7
            assert r['electric_square_projected']+1e-9>=q*(2*r['flux_protected_fraction']+r['hole_protected_fraction'])
            errors.append(e)
        assert min(np.diff(errors),default=0)>=-1e-12
        assert d['rewiring_removability_commutator_norm']<1e-9
        if f.name in MAIN[:4]:
            h=d['time_step'];late=[r for r in rs if r['T']>=3]
            low=max(0.,min(r['removable_fraction']-r['bounded_fraction_error_upper'] for r in rs)-h)
            flux_low=max(0.,min(r['flux_protected_fraction']-r['bounded_fraction_error_upper'] for r in late)-h)
            assert abs(low-d['continuous_time_bounds']['removable_fraction_lower_over_entire_window'])<1e-12
            assert abs(flux_low-d['continuous_time_bounds']['flux_protected_fraction_lower_over_late_window'])<1e-12
            assert low>.1 and flux_low>0
            final=next(r for r in rs if abs(r['T']-10)<1e-9)
            assert final['vacuum_flux_probability']>final['bounded_fraction_error_upper']
            summaries.append(dict(file=f.name,removable_fraction_continuous_lower=low,
                                  persistent_flux_fraction_continuous_lower=flux_low,
                                  final_state_error=errors[-1]))
    for f in P.glob('RECORDS_N*.json'):
        d=json.loads(f.read_text());G=d['graph_count']
        for r in d['rows']:
            assert -1e-9<=r['graph_entropy_bits']<=r['graph_classical_entropy_bits']+1e-8<=math.log2(G)+1e-7
            assert 1/G-1e-9<=r['graph_purity']<=1+1e-8
            assert r['classical_graph_uncertainty_given_full_configuration_bits']<=r['classical_conditional_entropy_support_bound_bits']+1e-8
            assert r['zero_removability_graph_offdiagonal_max']<1e-12
            assert 0<=r['graph_coherence_l1_normalized']<=1+1e-8
    mobility=json.loads((P/'MOBILITY.json').read_text())
    source_free=json.loads((P/'N4_Q4_B256_T20.json').read_text())
    for r in mobility['rows']:
        assert abs(r['susceptibility'])>r['susceptibility_error_upper']
        assert abs(r['no_hopping_after_same_preparation_response'])<1e-11
        ref=next(x for x in source_free['rows'] if abs(x['T']-r['wait'])<1e-9)
        assert abs(ref['removable_fraction']-r['before_source_removable_fraction'])<1e-10
    energy=json.loads((P/'ENERGY_AUDIT.json').read_text())
    for m in energy['models']:
        expected=audit(json.loads((P/m['file']).read_text()))
        assert {k:v for k,v in m.items() if k!='file'}==expected
        for r in m['rows']:
            assert r['electric_energy']+1e-9>=r['electric_tradeoff_rhs_projected']
            assert r['electric_energy_interval'][0]>=0
    deletion=json.loads((P/'DELETION_AUDIT.json').read_text())
    for m in deletion['models']:
        expected=deletion_audit(json.loads((P/m['file']).read_text()))
        assert {k:v for k,v in m.items() if k!='file'}==expected
    refined=json.loads((P/'REFINED_COUNTS.json').read_text())
    for m in refined['models']:
        original=json.loads((P/f"N{m['N']}_Q{m['Q']}_B{m['budget']}_T{m['T']:g}.json").read_text())
        r=next(r for r in original['rows'] if abs(r['T']-m['T'])<1e-9)
        for obs in m['rows']:
            assert abs(obs['value']-r[obs['detector']])<1e-10
            assert abs(obs['expectation_error_upper']-(2*m['forward_error']*obs['backward_error']+m['forward_error']**2))<1e-12
            if obs['detector']=='flux_protected_fraction' and m['budget']==6:assert obs['interval'][0]>0
    dense=json.loads((P/'N6_Q12_B10_T1.json').read_text())
    assert dense['continuous_time_bounds']['removable_fraction_lower_over_entire_window']>.1
    convergence=[]
    for low,high in [('N4_Q4_B32_T20.json','N4_Q4_B256_T20.json'),
                     ('N4_Q4_B128_T20.json','N4_Q4_B256_T20.json'),
                     ('N5_Q6_B12_T4.json','N5_Q6_B128_T10.json'),
                     ('N5_Q6_B32_T10.json','N5_Q6_B128_T10.json'),
                     ('N6_Q6_B48_T4.json','N6_Q6_B128_T10.json'),
                     ('N6_Q12_B6_T0.4.json','N6_Q12_B10_T1.json')]:
        a=json.loads((P/low).read_text());b=json.loads((P/high).read_text())
        common={round(r['T'],7):r for r in b['rows']};worst=0.
        for r in a['rows']:
            s=common[round(r['T'],7)]
            bound=r['bounded_fraction_error_upper']+s['bounded_fraction_error_upper']
            for key in ['removable_fraction','flux_protected_fraction','hole_protected_fraction']:
                difference=abs(r[key]-s[key]);assert difference<=bound+1e-10
                worst=max(worst,difference)
        convergence.append(dict(low=low,high=high,maximum_normalized_count_difference=worst))
    if record:(P/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    assert json.loads((P/'MANIFEST.json').read_text())==manifest
    out=dict(status='CP12_AUTONOMOUS_FLUX_GATE_VERIFIED',fresh=fresh,
             continuous_time_summaries=summaries,
             regulator_convergence=convergence,
             larger_N5_N6_evolutions_rerun_by_this_verifier=False,
             N5_record_unfolding_rerun_by_this_verifier=False,
             ruling='Autonomous internal graph records and sustained partial protection established in tested finite preparations; low mean removability, spatial selection and emergent spacetime not established.',
             floating_point_qualification='Analytic truncation and continuity bounds numerically evaluated; no interval-arithmetic certificate of roundoff.',
             runtime=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
             seconds=time.monotonic()-started)
    (P/'VERIFICATION.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--record-manifest',action='store_true')
    main(ap.parse_args().record_manifest)
