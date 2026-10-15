"""Independent legacy-integer controls on large-basis outgoing destinations."""
from pathlib import Path
import json
import numpy as np
from formation import model
from invariant_evolution import Invariant
from exact_destination import transitions


def main():
    sub,H,B,states,sizes,psi=model(6,12,12)
    reference=Invariant.__new__(Invariant)
    reference.__dict__.update({k:v for k,v in sub.__dict__.items() if k not in ['engine_arrays','engine_ids']})
    reference.cache={}
    indices=np.random.default_rng(20261015).choice(len(states),256,replace=False)
    count=zeros=0
    for i in indices:
        for target,c in transitions(6,states[i],600):
            fast=sub.canonical(target);slow=reference.canonical(target)
            assert fast==slow,(target,fast,slow)
            count+=1;zeros+=fast[1]==0
    result=dict(status='PASS',N=6,budget=12,sampled_parent_states=256,
        integer_canonical_comparisons=count,cancelled_orbits_checked=zeros,
        scope='Independent legacy Python canonicalization on random large-basis columns and every unregulated destination of those columns.')
    (Path(__file__).resolve().parent/'N6_SYMMETRY_CHECK.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)


if __name__=='__main__':main()
