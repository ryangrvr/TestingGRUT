#!/usr/bin/env python3
"""Verify execution identity, preserved outcomes and packet integrity."""
import hashlib
import json
import re
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    raw=json.loads((HERE/'results_initial.json').read_text())
    corrected=json.loads((HERE/'results_with_control.json').read_text())
    attacks=json.loads((HERE/'challenge_results.json').read_text())
    algebra=json.loads((HERE/'algebra_results.json').read_text())
    writer=json.loads((HERE/'writer_results_refined.json').read_text())
    initial=json.loads((HERE/'writer_results_initial.json').read_text())
    assert raw['protocol_sha256']==digest(HERE/'protocol.json')
    assert raw['source_sha256']==digest(HERE/'calculate_original.py')
    assert corrected['source_sha256']==digest(HERE/'calculate.py')
    assert corrected['original_source_sha256']==digest(HERE/'calculate_original.py')
    assert corrected['protocol_sha256']==digest(HERE/'protocol.json')
    assert corrected['control_protocol_sha256']==digest(HERE/'control_protocol.json')
    review=HERE.parent/'grut_memory_candidate_20261019_review'
    assert corrected['review_source_sha256']==digest(review/'uniform_field_control.py')
    assert attacks['initial_results_sha256']==digest(HERE/'results_initial.json')
    assert attacks['challenge_protocol_sha256']==digest(HERE/'challenge_protocol.json')
    assert attacks['challenge_source_sha256']==digest(HERE/'challenge.py')
    assert algebra['source_sha256']==digest(HERE/'algebra.py')
    assert writer['protocol_sha256']==digest(HERE/'writer_protocol.json')
    assert writer['original_source_sha256']==digest(HERE/'autonomous_writer.py')
    assert initial['source_sha256']==digest(HERE/'autonomous_writer.py')
    assert writer['source_sha256']==digest(HERE/'writer_replay.py')
    assert writer['amendment_sha256']==digest(HERE/'execution_amendment.json')
    assert writer['original_results_sha256']==digest(HERE/'writer_results_initial.json')
    assert len(raw['runs'])==7 and len(writer['runs'])==3
    assert len(corrected['runs'])==10
    assert all(corrected['control_method_checks'].values())
    assert not corrected['mechanism_claims']['facilitation_improves_minimum_retention_all_sizes']
    assert corrected['mechanism_claims']['matched_uniform_control_retains_better_on_registered_samples']
    replay_error=0.
    for name,old in raw['runs'].items():
        new=corrected['runs'][name]
        for experiment in ['retention','readout']:
            error=float(np.max(abs(np.array(old[experiment]['means'])-new[experiment]['means'])))
            replay_error=max(replay_error,error)
    assert replay_error<1e-10
    assert raw['all_registered_acceptance_pass'] and all(raw['checks'].values())
    assert all(algebra['checks'].values())
    assert sum(attacks['checks'].values())==3 and len(attacks['checks'])==9
    assert not attacks['checks']['autonomous_write_and_hold_L4']
    assert not attacks['checks']['local_independent_bit_retains_to_t10_L5']
    assert all(writer['checks'].values())
    for name,run in corrected['runs'].items():
        p=np.array(run['retention']['success'])
        means=np.array(run['retention']['means'])
        assert np.max(abs(p-(.5+abs(means[0]-means[1])/4)))<1e-13
    pattern=r'^(.*?)\s+L=(\d+): min retention t<=10 = ([0-9.]+); at t=40 = ([0-9.]+); min over \[0,40\] = ([0-9.]+); best readout = ([0-9.]+) at t=([0-9.]+)$'
    review_rows=re.findall(pattern,(review/'control_results.txt').read_text(),re.M)
    assert len(review_rows)==6
    for label,L,short,last,minimum,readout,when in review_rows:
        name=('candidate_L' if label.startswith('M1') else 'uniform_floor_L')+L
        r=corrected['runs'][name]; ts=np.array(r['retention']['times']); ps=np.array(r['retention']['success'])
        values=[ps[ts<=10].min(),ps[-1],ps.min(),r['readout']['best_success']]
        assert max(abs(a-float(b)) for a,b in zip(values,[short,last,minimum,readout]))<5.01e-7
        assert abs(r['readout']['best_registered_time']-float(when))<1e-12
    # Check the final plot/data summary against actual refined results.
    summary=json.loads((HERE/'summary.json').read_text())
    assert abs(summary['writer_minimum_memory_t10']-writer['runs'][-1]['minimum_all_cell_success_through_t10'])<1e-14
    assert summary['matched_comparison']==corrected['matched_comparison']
    assert summary['mechanism_claims']==corrected['mechanism_claims']
    for p in HERE.iterdir():
        if p.is_file() and p.suffix=='.py':compile(p.read_text(),str(p),'exec')
        if p.is_file() and p.suffix=='.json':json.loads(p.read_text())
        if p.is_file() and p.suffix=='.md':
            for dest in re.findall(r'\[[^\]]*\]\(([^)]+)\)',p.read_text()):
                if not dest.startswith(('http://','https://','#')):
                    assert (p.parent/dest.split('#')[0]).exists(),(p.name,dest)
    paths=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.json')
    manifest=dict(parent_commit='f28b49d977f4a3e2bbe8625b768ff9a1ae16e44f',
        amended_packet_commit='1997d205830875e73437c2bcd50ca621bf4bd30c',
        review_commit='cc249673086b984f090d170f394d4d635e9c451d',
        files=[dict(path=p.name,bytes=p.stat().st_size,sha256=digest(p)) for p in paths],
        registered_core_conditions=raw['checks'],retained_adversarial_conditions=attacks['checks'],
        writer_conditions=writer['checks'],
        control_method_checks=corrected['control_method_checks'],mechanism_claims=corrected['mechanism_claims'],
        original_replay_max_conditional_mean_error=replay_error,
        classification='Finite Ising repetition register with matter feedback and programmed writer; facilitation benefit rejected by matched control. Universal GRUT absent.',
        scope='Syntax, byte identity, links, numerical arithmetic and recorded outcomes; not an external validation or all-times scattering certificate.')
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(files=len(paths),core_pass=sum(raw['checks'].values()),
          attack_conditions_pass=sum(attacks['checks'].values()),attack_conditions_fail=6,
          writer_pass=sum(writer['checks'].values()),matched_control_sizes=3,
          original_replay_error=replay_error,facilitation_benefit_supported=False,
          classification=manifest['classification'])))


if __name__=='__main__':main()
