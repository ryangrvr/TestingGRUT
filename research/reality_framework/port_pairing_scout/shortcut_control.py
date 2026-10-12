#!/usr/bin/env python3
"""All-pairs reconnection produces a two-step geometry+fermion shortcut.

For a connected N-cycle choose any two distinct bonds separated by many
ring steps. Reconnecting them can create a direct new hopping bond between
their endpoint ports, allowing a matter transition at order time^2.
Matrix verification uses the exact full N=4 joint Hamiltonian.
This is a failure of *microscopic fixed-graph locality*, NOT a theorem
that a physical large-N continuum violates relativistic causality.
"""
import json
from pathlib import Path
import numpy as np
from full_joint import assemble,fermion_basis,canonical

def check_full(n=4,gamma=.1,t=1.):
    H,ms,fb,idx=assemble(n,relative_sign=-1.,gamma=gamma)
    lookup={canonical(p):i for i,p in enumerate(ms)}
    ring=tuple((2*i+1,2*((i+1)%n)) for i in range(n))
    source,target=0,2*(n//2)
    ia=next(i for i,e in enumerate(ring) if source in e)
    ib=next(i for i,e in enumerate(ring) if target in e)
    u=next(x for x in ring[ia] if x!=source)
    v=next(x for x in ring[ib] if x!=target)
    rewired=list(ring)
    rewired[ia]=(source,target)
    rewired[ib]=(u,v)
    # Different matchings. Fermion changes occupancy at source and target.
    occupied=[source]+[i for i in range(2*n)
                        if i!=source and i!=target][:n-1]
    F=sum(1<<i for i in occupied)
    G=(F^(1<<source))|(1<<target)
    dimF=len(fb)
    initial=lookup[canonical(ring)]*dimF+idx[F]
    final=lookup[canonical(rewired)]*dimF+idx[G]
    first=complex(H[final,initial])
    second=complex((H@H[:,initial]).toarray()[final,0])
    assert abs(first)<1e-12
    assert abs(abs(second)-t*gamma/n)<1e-11
    return dict(N=n,original_vertex_distance=n//2,
        source_port=source,target_port=target,
        Hamiltonian_direct_amplitude=[first.real,first.imag],
        H_squared_amplitude=[second.real,second.imag],
        unitary_t_squared_coefficient=[(-second/2).real,(-second/2).imag])

def analytic_general(n,gamma=.1,t=1.):
    # The matching-sector offdiagonal H_reconnect is -Gamma/N,
    # and the newly formed fermion bond contributes -t times the CAR sign.
    # Their product (H^2)_{final,initial} is +/- Gamma*t/N,
    # independent of the two bonds' ORIGINAL ring separation.
    return dict(N=n,original_graph_vertex_distance=n//2,
        predicted_abs_H_squared=gamma*t/n,
        predicted_abs_time_squared_unitary_amplitude=gamma*t/(2*n))

def main():
    exact=check_full()
    rows=[analytic_general(n) for n in (4,8,16,32,64,128)]
    assert abs(abs(exact['H_squared_amplitude'][0])-
               rows[0]['predicted_abs_H_squared'])<1e-12
    assert rows[-1]['predicted_abs_H_squared']>0
    result=dict(full_joint_matrix_control=exact,
                algebraic_largeN_examples=rows,
                limitation=("No pre-existing graph metric cone at finite N; "
                 "per-remote-channel shortcut suppressed by 1/N. "
                 "Continuum signal and actual experimental clock/locality "
                 "must still be analyzed; not a superluminality claim."))
    Path(__file__).with_name('shortcut_results.json').write_text(
        json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    print('SCOUT_RG02_NONLOCAL_RECONNECTION_SHORTCUT_PASS')

if __name__=='__main__':main()
