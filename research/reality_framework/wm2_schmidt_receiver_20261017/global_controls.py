"""G1-only tests, executed after its separate charter and original workers."""
from pathlib import Path
from fractions import Fraction as Q
from math import sqrt
import json
import numpy as np
import controls
from bounds import S5
from allocation import global_keep
from mps_global import Receiver

BASE=Path(__file__).resolve().parent


def source_check():
    source=(BASE/'mps.py').read_text()
    expected=source.replace('import gates as G\n','import gates as G\nfrom allocation import global_keep\n')
    expected=expected.replace('rank_cap=16,total_rank_cap=96','rank_cap=96,total_rank_cap=96')
    old="keep={q:min(self.rank_cap,len(d['s'])) for q,d in decomposed.items()}\n        if sum(keep.values())>self.total_rank_cap:raise Limit('total Schmidt rank ceiling')"
    new="keep=global_keep({q:d['s'] for q,d in decomposed.items()},self.total_rank_cap)"
    assert old in expected
    expected=expected.replace(old,new)
    assert expected==(BASE/'mps_global.py').read_text()
    worker=(BASE/'worker.py').read_text().replace('from mps import Receiver','from mps_global import Receiver').replace("method='S8'","method='S8-G1'")
    assert worker==(BASE/'worker_global.py').read_text()
    return dict(status='PASS',only_rank_allocation_differs=True,arithmetic_and_integration_unchanged=True)


def allocation_checks():
    witnesses=[]
    for counts in [[2,2],[2,3]]:
        words=[u for u in S5.type_words(counts) if u];n=len(words);q=sum(counts)-1
        spectrum=np.full(n,1/sqrt(n));keep=global_keep({q:spectrum},96)
        assert keep[q]==n
        weight=Q(max(0,n-16),n)
        witnesses.append(dict(counts=counts,nonempty_histories=n,single_charge=q,
            conditional_history_Schmidt_rank=1,joint_Schmidt_rank=n,
            original_16_discarded_weight=str(weight),
            original_16_discarded_norm_display=sqrt(float(weight)),
            G1_kept=keep[q],G1_discarded_weight='0',actual_prepared_state_claim=False))
    rng=np.random.default_rng(89108);cases=[]
    for charges,cap in [(3,17),(6,96),(5,31)]:
        spectra={q:np.sort(rng.random(45))[::-1] for q in range(charges)}
        keep=global_keep(spectra,cap)
        actual=sum(float((s[keep[q]:]**2).sum()) for q,s in spectra.items())
        ordered=sorted([float(x) for s in spectra.values() for x in s],reverse=True)
        optimal=sum(x*x for x in ordered[cap:])
        assert sum(keep.values())==cap and abs(actual-optimal)<2e-13
        cases.append(dict(charges=charges,total_cap=cap,retained=sum(keep.values()),
            actual_tail_squared=actual,global_sorted_tail_squared=optimal))
    return dict(witnesses=witnesses,global_optimal_proposals=cases,
        scope='rank-capacity controls, not formation observations')


def run():
    result=dict(source=source_check(),allocation=allocation_checks())
    controls.Receiver=Receiver
    result['physical_and_arithmetic']=controls.run()
    result['status']='PASS';return result


if __name__=='__main__':
    result=run();(BASE/'G1_CONTROLS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
