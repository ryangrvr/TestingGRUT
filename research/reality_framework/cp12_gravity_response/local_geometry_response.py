"""Complete CP-12 controls for a definite prepared graph's local back-reaction.

The law is imported unchanged. Both Hamiltonian evolutions and exact finite
matrix Taylor coefficients are computed; no geometry phase is assumed.
"""
from pathlib import Path
import argparse
import itertools
import json
import math
import sys
from fractions import Fraction

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

P = Path(__file__).resolve().parent
sys.path.insert(0, str(P.parents[1] / "constructive_physics/attempt12"))
from model import build, gauss, GAMMA, HOP


def local_bounds(n, q, T):
    M = n * (n-1) // 2
    L = M-q
    g = GAMMA*q*math.sqrt(L)/M
    I = (T*math.sqrt(1+(g*T)**2)+math.asinh(g*T)/g)/2 if g else T
    present_amplitude = GAMMA/M*(math.sqrt(q)*g*T*T/2+math.sqrt(L)*I)
    absent_amplitude = GAMMA/M*(math.sqrt(q)*I+math.sqrt(L)*g*T*T/2)
    return dict(outside_mean_bound=g*g*T*T,
                present_loss_bound=min(1., present_amplitude**2),
                absent_creation_bound=min(1., absent_amplitude**2),
                present_response_bound=min(2., 2*present_amplitude**2),
                absent_response_bound=min(2., 2*absent_amplitude**2))


def response_taylor(H, psi, Bpsi, A, order):
    # <exp(-iHt)psi| A |exp(-iHt) Bpsi>; independent of nested-commutator code.
    hp=[psi]; hb=[Bpsi]
    for _ in range(order):
        hp.append(H @ hp[-1]); hb.append(H @ hb[-1])
    out=[]
    for m in range(order+1):
        coefficient=sum(1j**a*(-1j)**(m-a)*np.vdot(hp[a], A*hb[m-a])/
                        (math.factorial(a)*math.factorial(m-a)) for a in range(m+1))
        out.append(-2*float(coefficient.imag))
    return out


def run(n, host, cutoff, label, times):
    pairs=list(itertools.combinations(range(n), 2)); M=len(pairs)
    present=[pairs.index(tuple(sorted(e))) for e in host]
    mask=sum(1 << e for e in present); q=len(present)
    H, states, lookup, _=build(n,q,cutoff)
    masks=np.array([s[0] for s in states], dtype=np.int64)
    bits=np.array([s[1] for s in states], dtype=np.int64)
    coo=H.tocoo(); change=masks[coo.row]!=masks[coo.col]
    V=coo_matrix((coo.data[change], (coo.row[change], coo.col[change])),shape=H.shape).tocsr()
    HB=H-V
    source_vertex=0; d=sum(source_vertex in e for e in host)
    select=np.array([not cg and bool(int(bits[a]^bits[b]) & (1 << source_vertex))
                     and (int(bits[a]^bits[b]) >> n).bit_count()==1
                     for a,b,cg in zip(coo.row,coo.col,change)])
    # Modulate all existing frozen hoppings at vertex 0; no F-membership query.
    B=coo_matrix((coo.data[select]/math.sqrt(d), (coo.row[select],coo.col[select])),
                 shape=H.shape).tocsr()
    sea=sum(1 << (n+e) for e in range(M))
    psi=np.zeros(len(states)); psi[lookup[(mask,sea,(0,)*(2*M))]]=1
    Bpsi=B @ psi
    assert abs(np.linalg.norm(Bpsi)-HOP)<1e-12
    assert np.all(masks[np.flatnonzero(Bpsi)]==mask)
    rhs=np.column_stack([psi,Bpsi]); rows=[]; taylor=[]
    for e, pair in enumerate(pairs):
        A=(masks >> e & 1).astype(float)
        # Subtract the initially constant readout to avoid cancellation of
        # a large, identically real conserved inner product.
        coeff=response_taylor(H,psi,Bpsi,A-float(e in present),5)
        predicted=-(GAMMA/M)**2*(M-q)*HOP**2/math.sqrt(d) if e in present and 0 in pair else (
                    (GAMMA/M)**2*math.sqrt(d)*HOP**2 if e not in present else 0.)
        taylor.append(dict(edge=list(pair),initial_present=e in present,
                           coefficients=coeff,predicted_cubic=predicted,
                           cubic_error=abs(coeff[3]-predicted)))
        assert abs(coeff[3]-predicted)<3e-11
        assert max(abs(x) for x in coeff[:3])<3e-11
    for T in times:
        evolved=expm_multiply(-1j*T*H,rhs)
        fixed=expm_multiply(-1j*T*HB,rhs)
        p=abs(evolved[:,0])**2
        outside=np.array([(int(m) & ~mask).bit_count() for m in masks])
        bounds=local_bounds(n,q,T)
        assert float(outside @ p)<=bounds["outside_mean_bound"]+2e-10
        for e,pair in enumerate(pairs):
            A=(masks >> e & 1).astype(float)
            R=-2*float(np.imag(np.vdot(evolved[:,0],A*evolved[:,1])))
            Rfixed=-2*float(np.imag(np.vdot(fixed[:,0],A*fixed[:,1])))
            assert abs(Rfixed)<2e-10
            delta=(1-float(A @ p)) if e in present else float(A @ p)
            prob_bound=bounds["present_loss_bound" if e in present else "absent_creation_bound"]
            response_bound=bounds["present_response_bound" if e in present else "absent_response_bound"]
            assert delta<=prob_bound+2e-10
            assert abs(R)<=response_bound+2e-10
            item=dict(T=T,edge=list(pair),initial_present=e in present,
                      occupation_change_probability=delta,probability_bound=prob_bound,
                      response=R,fixed_block_response=Rfixed,response_bound=response_bound,
                      mean_outside_bonds=float(outside @ p),outside_mean_bound=bounds["outside_mean_bound"])
            if T in [.4,2.]:
                eps=1e-4
                if e==0:
                    pulseplus=expm_multiply(-1j*T*H,expm_multiply(1j*eps*B,psi))
                    pulseminus=expm_multiply(-1j*T*H,expm_multiply(-1j*eps*B,psi))
                derivative=float((A @ abs(pulseplus)**2-A @ abs(pulseminus)**2)/(2*eps))
                assert abs(derivative-R)<2e-7
                item.update(independent_impulse_derivative=derivative,impulse_error=abs(derivative-R))
            rows.append(item)
    out=dict(status="CP12_LOCAL_GEOMETRY_RESPONSE_PASS",N=n,Q=q,L=cutoff,
             graph=host,dimension=len(states),source_vertex=0,source_prepared_norm=HOP,
             taylor=taylor,rows=rows,
             scope="Complete finite CP-12 in a Gauss-law electric regulator. Definite supplied host; no spin-2 or continuum inference.")
    (P/f"local_response_{label}_L{cutoff}.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(dict(N=n,Q=q,L=cutoff,dimension=len(states),
                         maximum_cubic_error=max(x['cubic_error'] for x in taylor),
                         first_edges=taylor[:3],maximum_response=max(abs(x['response']) for x in rows)),indent=2),flush=True)
    return out


if __name__=="__main__":
    ap=argparse.ArgumentParser(); ap.add_argument("--case",choices=["path","cycle","star"],default="path")
    ap.add_argument("--cutoff",type=int,default=2); args=ap.parse_args()
    cases={"path":(4,[(0,1),(1,2),(2,3)]),"cycle":(4,[(0,1),(1,2),(2,3),(3,0)]),
           "star":(5,[(0,1),(0,2),(0,3),(0,4)])}
    n,edges=cases[args.case]
    run(n,edges,args.cutoff,args.case,[.01,.1,.4,1.,2.,5.])
