"""Coherent physical-basis check of the fixed neutral-exterior identity.

Supplied complex vectors validate an analytical identity. They are not
substituted for the neutral-origin formation calculations.
"""
from collections import defaultdict
from itertools import combinations
from math import lcm
from pathlib import Path
import json
import numpy as np
from boundary import gauss,transitions,hopping_edge,flow,outside_weight
from pauli_algebra import family
from exact_destination import hop

P=Path(__file__).resolve().parent


def state(n,occupied):
    s,C,_=family(n);g,bits,fields=s;E=list(fields);es=list(combinations(range(n),2))
    for i in range(n):
        if occupied>>i&1:continue
        j=(i+1)%n;e=es.index(tuple(sorted((i,j))));side=0 if i<j else 1
        bits,_=hop(bits,n+e,i);E[2*e+side]-=1
    s=g,bits,tuple(E);gauss(n,s);assert flow(n,s)[0]==C
    return s,C


def norm2(v):return float(sum(abs(a)**2 for a in v.values()))


def main():
    rng=np.random.default_rng(7684601);rows=[]
    for n in range(5,9):
        amplitudes=rng.normal(size=1<<n)+1j*rng.normal(size=1<<n)
        amplitudes/=np.linalg.norm(amplitudes)
        exterior=defaultdict(complex);moving=defaultdict(complex);expected=0.;checked=0
        sc=lcm(8,10*n*(n-1)//2)
        for mask,a in enumerate(amplitudes):
            s,C=state(n,mask);info=flow(n,s);expected+=2*(n-mask.bit_count())*abs(a)**2
            for t,c in transitions(n,s,sc):
                e=hopping_edge(n,s,t)
                if e is None:continue
                gauss(n,t);checked+=1
                if not(C>>e&1):exterior[t]+=a*c/sc
                moving[t]+=a*c/sc*outside_weight(n,s,t,e,info)
        rest={t:moving[t]-exterior.get(t,0.) for t in moving}
        cross=sum(a.conjugate()*rest.get(t,0.) for t,a in exterior.items())
        measured=norm2(exterior);total=norm2(moving)
        assert abs(measured-expected)<1e-11 and abs(cross)<1e-11 and total>=measured-1e-11
        rows.append(dict(N=n,Q=2*n,coherent_source_states=len(amplitudes),raw_hops_checked=checked,
            neutral_exterior_action_square=measured,endpoint_vacancy_formula=expected,
            moving_outside_action_square=total,neutral_exterior_other_outside_cross_abs=abs(cross)))
    result=dict(status='PASS',rows=rows,
        identity='For a fixed recorded core C with every non-core edge fermion filled and its fields zero, ||sum_(e outside C) K_e psi||^2 = <sum_(e outside C) g_e(2-n_i-n_j)>. Exterior-hole tags make the outputs orthogonal, for arbitrary coherent physical amplitudes.',
        moving_core_lower_bound='If C is also the complete initial through-cycle core, its moving-outside action is at least the neutral-exterior value: outside and core-hop outputs have different exterior-hole patterns.',
        consequence='Within this explicitly stated neutral-exterior sector, destructive interference cannot cancel the first exterior hopping action. Zero action requires full Pauli filling at every active exterior endpoint. This does not exclude later interference or dressed transport in the unrestricted state.',
        scope='Mathematical identity plus supplied complex-vector validation, not a neutral-origin formation result, finite-time probability bound, spectral-gap result or large-size locality theorem.')
    (P/'NEUTRAL_EXTERIOR.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
