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
    attacks=json.loads((HERE/'challenge_results.json').read_text())
    algebra=json.loads((HERE/'algebra_results.json').read_text())
    writer=json.loads((HERE/'writer_results_refined.json').read_text())
    initial=json.loads((HERE/'writer_results_initial.json').read_text())
    assert raw['protocol_sha256']==digest(HERE/'protocol.json')
    assert raw['source_sha256']==digest(HERE/'calculate.py')
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
    assert raw['all_registered_acceptance_pass'] and all(raw['checks'].values())
    assert all(algebra['checks'].values())
    assert sum(attacks['checks'].values())==3 and len(attacks['checks'])==9
    assert not attacks['checks']['autonomous_write_and_hold_L4']
    assert not attacks['checks']['local_independent_bit_retains_to_t10_L5']
    assert all(writer['checks'].values())
    for name,run in raw['runs'].items():
        p=np.array(run['retention']['success'])
        means=np.array(run['retention']['means'])
        assert np.max(abs(p-(.5+abs(means[0]-means[1])/4)))<1e-13
    # Check the final plot/data summary against actual refined results.
    summary=json.loads((HERE/'summary.json').read_text())
    assert abs(summary['writer_minimum_memory_t10']-writer['runs'][-1]['minimum_all_cell_success_through_t10'])<1e-14
    for p in HERE.iterdir():
        if p.is_file() and p.suffix=='.py':compile(p.read_text(),str(p),'exec')
        if p.is_file() and p.suffix=='.json':json.loads(p.read_text())
        if p.is_file() and p.suffix=='.md':
            for dest in re.findall(r'\[[^\]]*\]\(([^)]+)\)',p.read_text()):
                if not dest.startswith(('http://','https://','#')):
                    assert (p.parent/dest.split('#')[0]).exists(),(p.name,dest)
    paths=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.json')
    manifest=dict(parent_commit='f28b49d977f4a3e2bbe8625b768ff9a1ae16e44f',
        files=[dict(path=p.name,bytes=p.stat().st_size,sha256=digest(p)) for p in paths],
        registered_core_conditions=raw['checks'],retained_adversarial_conditions=attacks['checks'],
        writer_conditions=writer['checks'],
        classification='Narrower working model; universal gravity, infinite-time memory and empirical validation absent.',
        scope='Syntax, byte identity, links, numerical arithmetic and recorded outcomes; not an external validation or all-times scattering certificate.')
    (HERE/'MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(files=len(paths),core_pass=sum(raw['checks'].values()),
          attack_conditions_pass=sum(attacks['checks'].values()),attack_conditions_fail=6,
          writer_pass=sum(writer['checks'].values()),classification=manifest['classification'])))


if __name__=='__main__':main()
