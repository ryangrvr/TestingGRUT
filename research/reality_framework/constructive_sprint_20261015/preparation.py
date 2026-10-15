"""Constructive, topology-neutral input circuit accounting and small checks.

This is an external preparation protocol, not autonomous preparation by CP12.
Only graph occupancy is simulated; fermion sea and zero electric fields are
simple additional supplied basis resources.
"""
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json

P=Path(__file__).resolve().parent


def probabilities(M,Q):
    states={(0,Q):F(1)}
    for j in range(M):
        out={}
        for (mask,remaining),p in states.items():
            left=M-j
            if remaining<left:out[mask,remaining]=p*F(left-remaining,left)
            if remaining:out[mask|(1<<j),remaining-1]=p*F(remaining,left)
        states=out
    assert all(r==0 for _,r in states)
    assert len(states)==comb(M,Q)
    assert set(states.values())=={F(1,comb(M,Q))}
    assert sum(states.values())==1
    return dict(M=M,Q=Q,final_graphs=len(states),every_exact_probability=str(F(1,comb(M,Q))),
        conditional_rotation_stages=M,maximum_distinct_rotation_angles=M*(Q+1),
        graph_phase_mask=False,external_controller=True)


def run():
    out=dict(status='TOPOLOGY_NEUTRAL_PREPARATION_CHECK_PASS',
        rows=[probabilities(M,Q) for M,Q in [(6,4),(10,10),(15,12)]],
        scope='Exact branch-probability check of a coherent sequential controlled-rotation circuit: at edge j rotate sqrt(r/(M-j)) into occupation and decrement a reversible remaining-count register. Positive amplitudes give the uniform Q-edge Dicke ray. This does not derive a physical controller, quantum gates, hbar, initialization energy or preparation by the frozen Hamiltonian.')
    (P/'PREPARATION.json').write_text(json.dumps(out,indent=2)+'\n');return out


if __name__=='__main__':print(json.dumps(run(),indent=2))
