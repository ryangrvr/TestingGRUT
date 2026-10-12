#!/usr/bin/env python3
"""Joint graph+fermion relational detector test in full N=3,4 ground states.

Anchor physical source at fermion port 0. A_one measures occupation at the
port actually connected to 0 in each fluctuating matching. A_two follows
a connection, traverses that vertex, then follows its other connection.
Probes mathematically require a relational local readout apparatus that
has NOT been explicitly implemented as a physical quantum pointer.

OPENBLAS_NUM_THREADS=1 python relational_response.py
"""
import json
from pathlib import Path
import numpy as np
from scipy.sparse.linalg import eigsh, expm_multiply
from full_joint import assemble, fermion_basis

def probes(n,ms,fb):
    source=np.tile(np.array([(bit&1)!=0 for bit in fb],float),len(ms))
    one=[];two=[];fixed=[]
    for p in ms:
        mate={a:b for a,b in p}
        mate.update({b:a for a,b in p})
        first=mate[0]
        second=mate[first^1]
        one.extend([(bit>>first)&1 for bit in fb])
        two.extend([(bit>>second)&1 for bit in fb])
        fixed.extend([(bit>>2)&1 for bit in fb])
    return source,[np.array(a,float) for a in (one,two,fixed)]

def initial_slope(H,g,A,B):
    """Exact first time derivative of i<[A(t),B]> for diagonal A,B."""
    C=H.tocoo()
    weights=C.data*(A[C.col]-A[C.row])*(B[C.col]-B[C.row])
    contributions=np.bincount(C.row,weights=weights*g[C.col],
                              minlength=len(g))
    return float(-np.vdot(g,contributions).real)

def test(n):
    H,ms,fb,_=assemble(n,relative_sign=-1.)
    E,V=eigsh(H,k=1,which="SA",tol=1e-11,
              v0=np.ones(H.shape[0]))
    g=V[:,0]
    B,detectors=probes(n,ms,fb)
    t_values=(0.,.02,.04,.08,.1,.2,.4,1.)
    rows=[]
    for t in t_values:
        v=expm_multiply(-1j*t*H,g)
        vb=expm_multiply(-1j*t*H,B*g)
        response=[float(-2*np.vdot(v,A*vb).imag) for A in detectors]
        if t==.4:
            eps=1e-5
            vp=expm_multiply(-1j*t*H,(1+(np.exp(1j*eps)-1)*B)*g)
            vm=expm_multiply(-1j*t*H,(1+(np.exp(-1j*eps)-1)*B)*g)
            impulse=[float((np.vdot(vp,A*vp).real-
                            np.vdot(vm,A*vm).real)/(2*eps))
                     for A in detectors]
            assert max(abs(x-y) for x,y in zip(response,impulse))<1e-8
        rows.append(dict(t=t,one_link=response[0],
                        two_links=response[1],fixed_other_vertex_port=response[2]))
    slopes=[initial_slope(H,g,A,B) for A in detectors]
    if n==4:
        assert abs(rows[-2]["one_link"]-(-.26613337098822465))<2e-8
        assert abs(rows[-2]["two_links"]-(-.00381281686655779))<2e-8
        assert abs(slopes[0]-(-.8526272729841))<2e-8
        assert abs(slopes[1])<1e-7
        # Cubic early behavior for this prepared ground state, not a
        # universal speed limit: t=0.02,0.04,0.08.
        cubic=[r['two_links']/r['t']**3 for r in rows if
               r['t'] in (.02,.04,.08)]
        assert max(cubic)-min(cubic)<5e-5
    return dict(N=n,source="physical port occupancy n0",
                detector="one- or two-bond-traversal conditioned on matching",
                initial_response_slopes=slopes,rows=rows)

if __name__=="__main__":
    result=[test(3),test(4)]
    Path(__file__).with_name('relational_response_results.json').write_text(
       json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    print('SCOUT_RG02_RELATIONAL_RESPONSE_PASS')
