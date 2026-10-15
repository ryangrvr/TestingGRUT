"""Check saved evidence, independent identities and pinned source contents.

This fast verifier does not claim to rerun the large dynamics. The explicit
reproduction commands and their actual execution scope are recorded separately.
"""
from pathlib import Path
import hashlib,json

P=Path(__file__).resolve().parent


def read(name):return json.loads((P/name).read_text())


def main():
    checked=[]
    for name in ['CHECKS.json','N6_SYMMETRY_CHECK.json','CURRENT_ACTION_CHECK.json',
                 'LOCAL_ONSET.json','ONSET_N5.json','ONSET_N6.json',
                 'RAW_ONSET_N5_PYTHON.json','RAW_ONSET_N6.json']:
        assert read(name)['status']=='PASS',name;checked.append(name)
    local=read('LOCAL_ONSET.json')
    for n in [5,6]:
        onset=read(f'ONSET_N{n}.json')
        raw=read('RAW_ONSET_N5_PYTHON.json' if n==5 else 'RAW_ONSET_N6.json')
        row=next(x for x in local['cases'] if x['N']==n)
        assert onset['leading_probability_coefficient']==raw['leading_probability_coefficient']==row['coefficient']
        assert onset['integer_squared_norm']==raw['integer_squared_norm']==row['integer_squared_norm']
    n5=read('N5_Q10_B50_A0.25_V0.25.json');n6=read('N6_Q12_B12_A0.25_V0.25.json')
    def row(r,t):return next(x for x in r['rows'] if x['time']==t)
    assert row(n5,5.)['probability_interval'][0]>.39
    assert row(n6,1.)['probability_interval'][0]>.0004
    assert row(n6,2.)['probability_interval'][0]==0
    for r in [n5,n6,read('N5_Q10_B30_A0.25_V0.25.json'),read('N6_Q12_B8_A0.25_V0.25.json')]:
        assert r['status']=='PASS' and r['full_law_omitted_states']
        for point in r['rows']:
            assert point['norm_error']<1e-9
            if 'current' in point:
                c=point['current'];assert abs(c['total']-c['initial_source']-c['coherent_return'])<1e-12
        h=r['two_time_histories']
        assert abs(h['unmeasured_later_probability']-h['initial_flag_recorded_later_probability']-h['interference_in_later_probability'])<1e-9
        assert 0<=h['conditional_interval'][0]<=h['conditional_interval'][1]<=1
        assert r['finite_difference_checks'][1]['error']<=.3*r['finite_difference_checks'][0]['error']+1e-12
    low=read('N6_Q12_B8_A0.25_V0.25.json')
    assert row(low,.5)['probability_interval'][0]==0
    # Both cutoffs must have overlapping controlled probability intervals.
    for high,low in [(n5,read('N5_Q10_B30_A0.25_V0.25.json')),(n6,low)]:
        for t in sorted(set(x['time'] for x in high['rows'])&set(x['time'] for x in low['rows'])):
            a,b=row(high,t)['probability_interval'],row(low,t)['probability_interval']
            assert max(a[0],b[0])<=min(a[1],b[1]),(high['N'],t,a,b)
    response=read('RESPONSE_N5_B50.json')
    assert response['status']=='PASS' and response['charged_initial_event_response']['zero_excluded']
    assert response['unfiltered_response']['zero_excluded']
    assert abs(response['charged_event_probability']-row(n5,2.)['event_probability'])<1e-10
    assert max(abs(x['difference_from_commutator']) for x in response['finite_pulse_controls'])<1e-9
    audit=read('SOURCE_AUDIT.json')
    for relative,digest in audit.items():
        assert hashlib.sha256((P.parent/relative).read_bytes()).hexdigest()==digest,relative
    if (P/'MANIFEST.json').exists():
        for f in read('MANIFEST.json')['files']:
            path=P/f['path'];assert path.stat().st_size==f['bytes'] and hashlib.sha256(path.read_bytes()).hexdigest()==f['sha256'],str(path)
    print(json.dumps(dict(status='PASS',independent_evidence_files=len(checked),
        full_law_positive_N5_N6=True,cutoff_intervals_overlap=True,
        source_current_identity=True,history_interference_identity=True,
        conditional_response_resolved=True,independent_integer_onsets_match=True,
        pinned_sources_match=True,scope='Saved-evidence consistency and source-content validation; not a fresh rerun of large simulations or external mathematical certification.'),indent=2))


if __name__=='__main__':main()
