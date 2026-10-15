"""Fresh small/medium reproduction plus explicit hash validation of saved work."""
from pathlib import Path
import argparse,hashlib,json,os,platform,tempfile,time
import numpy as np
import scipy
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply
import audit
import short_lag
import independent_checks

P=Path(__file__).resolve().parent


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def dependencies():
    old=json.loads((P.parent/'cp12_flux_loop_gate/MANIFEST.json').read_text())
    paths=set(old['dependency_hashes'])
    paths.update(['cp12_flux_loop_gate/'+f for f in ['loops.py','references.py','persistence.py']])
    for name,expected in old['dependency_hashes'].items():
        assert digest(P.parent/name)==expected,('parent dependency changed',name)
    return {name:digest(P.parent/name) for name in sorted(paths)}


def inputs():
    paths=['cp12_flux_loop_gate/MANIFEST.json','cp12_flux_loop_gate/REPORT.md',
           'cp12_flux_loop_gate/CONDITIONAL_DENSITY.json',
           'cp12_relational_record_gate/RELATIONAL_N4_Q4_B256.json',
           'cp12_relational_record_gate/RELATIONAL_N5_Q5_B128.json',
           'cp12_relational_record_gate/RELATIONAL_N6_Q6_B128.json',
           'cp12_relational_record_gate/RELATIONAL_N6_Q12_B10.json']
    return {name:digest(P.parent/name) for name in paths}


def manifest():
    files=[p for p in P.iterdir() if p.is_file() and p.name not in ['MANIFEST.json','VERIFICATION.json']]
    return dict(parent_commit='7d414961a56a49df1265b9501cd305331354ffc6',
        packet_hashes={p.name:digest(p) for p in sorted(files)},
        dependency_hashes=dependencies(),input_hashes=inputs())


def check_manifest(m):
    for name,expected in m['packet_hashes'].items():assert digest(P/name)==expected,('packet hash',name)
    for section in ['dependency_hashes','input_hashes']:
        for name,expected in m[section].items():assert digest(P.parent/name)==expected,(section,name)


def difference(a,b):
    """Compare scientific values, not replay timestamps or zip-byte provenance."""
    worst=0.
    for field in ['no_initial_filter']:
        for phase in ['zero','plus','minus']:
            for k,v in a[field]['paths'][phase]['values'].items():
                worst=max(worst,abs(v-b[field]['paths'][phase]['values'][k]))
        for k,v in a[field]['susceptibilities'].items():
            worst=max(worst,abs(v['value']-b[field]['susceptibilities'][k]['value']))
        worst=max(worst,abs(a[field]['density_response_inside_final_circulation']['value']-
                            b[field]['density_response_inside_final_circulation']['value']))
    for x,y in zip(a['record_comparisons'],b['record_comparisons']):
        assert x['record']==y['record']
        worst=max(worst,abs(x['initial_probability']-y['initial_probability']),
                  abs(x['nonselective_recording_work']-y['nonselective_recording_work']))
        for k,v in x['response_interference'].items():
            worst=max(worst,abs(v['value']-y['response_interference'][k]['value']))
        for branch in ['selected_branch','complement_branch','nonselective_record']:
            for phase in ['zero','plus','minus']:
                for k,v in x[branch]['paths'][phase]['values'].items():
                    worst=max(worst,abs(v-y[branch]['paths'][phase]['values'][k]))
    return worst


def run(record=False):
    start=time.monotonic();mp=P/'MANIFEST.json'
    if mp.exists() and not record:check_manifest(json.loads(mp.read_text()))
    worst=0.;fresh=[]
    with tempfile.TemporaryDirectory(prefix='cp12-unfiltered-verify-') as temp:
        tmp=Path(temp);waves=tmp/'waves';waves.mkdir();out=tmp/'outputs';out.mkdir()
        oldcache=os.environ.get('CP12_FORMATION_CACHE');os.environ['CP12_FORMATION_CACHE']=str(tmp/'fresh-models')
        oldp=audit.P;audit.P=out;independent_checks.P=out
        try:
            for n,budget,waits in [(4,256,[5.,10.]),(5,128,[10.])]:
                sub,H,B,states,sizes,psi=audit.model(n,n,budget)
                A=H-n/2*eye(H.shape[0],format='csr')
                for wait in waits:
                    v=expm_multiply(-1j*wait*A,psi,traceA=-1j*wait*float(A.diagonal().sum()))
                    np.savez(waves/f'N{n}_Q{n}_B{budget}_T{wait:g}.npz',wave=v)
                    result=audit.main_case(n,n,budget,wait,.5,waves)
                    old=json.loads((P/f'UNFILTERED_N{n}_Q{n}_B{budget}_W{wait:g}_L0.5.json').read_text())
                    error=difference(result,old);worst=max(worst,error);assert error<1e-9
                    assert result['no_initial_filter']['susceptibilities']['density']['interval'][0]>0
                    assert result['no_initial_filter']['density_response_inside_final_circulation']['interval'][0]>0
                    fresh.append(dict(N=n,Q=n,budget=budget,wait=wait,lag=.5,
                                      response_dimension=result['dimension'],maximum_saved_difference=error))
            physical=short_lag.independent_operator_control()
            coeff=short_lag.main_coefficients(waves)
            saved=json.loads((P/'SHORT_LAG_CONTROLS.json').read_text())['main_state_coefficients']
            coefficient_difference=max(abs(x['projected_circulation_weighted_density_T2_coefficient']-
                y['projected_circulation_weighted_density_T2_coefficient']) for x,y in zip(coeff,saved))
            assert coefficient_difference<1e-10
            assert all(x['full_law_circulation_weighted_density_T2_coefficient_interval'][0]>0 for x in coeff)
            dilation=independent_checks.run()
        finally:
            audit.P=oldp
            if oldcache is None:os.environ.pop('CP12_FORMATION_CACHE',None)
            else:os.environ['CP12_FORMATION_CACHE']=oldcache
    capacity=json.loads((P/'CYCLE_CAPACITY.json').read_text())
    assert capacity['status']=='CP12_CIRCULATION_CAPACITY_PASS'
    for case in capacity['cases']:
        assert case['controls']['minimum_inequality_margin']>=0
        for d in case['probabilities'].values():
            lo,hi=d['full_law_interval'];assert 0<=lo<=d['projected']+1e-12 and d['projected']<=hi+1e-12<=1+1e-12
    refs=json.loads((P/'ENERGY_REFERENCE_RESPONSE.json').read_text())
    assert len(refs['rows'])==2
    for row in refs['rows']:
        assert abs(row['full_mean_energy']-2)<1e-10
        assert row['response']['density_response_inside_final_circulation']['interval'][0]>0
    if record:mp.write_text(json.dumps(manifest(),indent=2)+'\n')
    check_manifest(json.loads(mp.read_text()))
    result=dict(status='CP12_UNFILTERED_RECORD_GATE_PASS',
        fresh_model_cache=True,fresh_snapshot_reconstruction=True,
        freshly_replayed_cases=fresh,maximum_main_saved_difference=worst,
        independent_physical_control=physical,
        fresh_weighted_coefficient_maximum_difference=coefficient_difference,
        dense_unitary_dilations=dilation,
        hash_validated_without_default_rerun=['Both N4 mixed-reference finite-lag evolutions',
            'All four larger capacity scans, including N6_Q6 and N6_Q12',
            'Saved N4 short-lag display samples'],
        source_and_input_hashes_checked=True,
        Python=platform.python_version(),NumPy=np.__version__,SciPy=scipy.__version__,
        seconds=time.monotonic()-start,
        scope='Fresh main-state and independent checks plus explicitly separate saved-file validation. Analytic full-law electric bounds evaluated in ordinary floating point; no formal roundoff interval certificate. No test of a macroscopic geometry or a causal cone.')
    (P/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],maximum_difference=worst,seconds=result['seconds'])),flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--record-manifest',action='store_true')
    run(ap.parse_args().record_manifest)
