"""Reproduce the completed numerical session calculations.

This wrapper was assembled after the container went offline. It has not been
executed as a packaged script in this session. The underlying commands were.
No approximation to the frozen Hamiltonian is silently promoted to full law.
"""
from pathlib import Path
import argparse, contextlib, io, json, os, sys, time
import numpy as np
from scipy.linalg import eigh

P=Path(__file__).resolve().parent
ROOT=P.parent
sys.path[:0]=[str(ROOT/'cp12_unmasked_formation'),
              str(ROOT/'cp12_flux_loop_gate'),
              str(ROOT/'cp12_distance_history')]
from run import model
from loops import incidence_masks, directed_support
from evolution import Evolution


def target(n,k,states):
    return np.array([
        (s[1]&((1<<n)-1)).bit_count()>=k
        and directed_support(n,*incidence_masks(n,s))==s[0]
        for s in states
    ])


def probability_interval(value,error):
    return [max(0.,float(np.sqrt(value))-error)**2,
            min(1.,(float(np.sqrt(value))+error)**2)]


def spectral():
    sub,H,B,states,sizes,psi=model(4,4,256)
    G=target(4,1,states)
    en,V=eigh(H.toarray())
    c=V.conj().T@psi
    reconstruction=float(np.linalg.norm(psi-V@c))
    internal=np.linalg.norm(H@V-V*en,axis=0)
    boundary=np.linalg.norm(B@V,axis=0)
    full=np.sqrt(internal**2+boundary**2)
    residual=float(abs(c)@full)
    X=V[G]*c[None,:]
    C=X.conj().T@X
    gram=float(abs(C).sum())
    rows=[]
    for T in [1.,2.,5.,10.]:
        value=float(np.linalg.norm(X@np.exp(-1j*en*T))**2)
        error=reconstruction+T*residual
        rows.append(dict(T=T,p=value,full_state_error_bound=error,
            full_probability_interval=probability_interval(value,error),
            full_spectral_envelope_upper=min(1.,(np.sqrt(gram)+error)**2)))
    assert G.sum()==86 and len(states)==577
    return dict(N=4,Q=4,budget=256,retained_dimension=len(states),
        outside_orbits=B.shape[0],target_orbits=int(G.sum()),
        trial_modes=len(en),coefficient_pruning=False,
        initial_reconstruction_norm=reconstruction,
        full_H_weighted_trial_eigenvector_residual=residual,
        boundary_weighted_residual=float(abs(c)@boundary),
        internal_weighted_residual=float(abs(c)@internal),
        absolute_weighted_target_Gram_sum=gram,
        projected_infinite_time_average_non_degenerate=float(np.trace(C).real),
        rows=rows)


def capacity():
    sub,H,B,states,sizes,psi=model(5,10,32)
    G=target(5,2,states)
    nv=np.array([(s[1]&31).bit_count() for s in states])
    rows=[]
    for T in [1.,2.,3.]:
        u,errors,_=Evolution(H,B,5,steps=32)(psi,T)
        value=float(np.linalg.norm(u[G])**2);error=float(errors[0])
        rows.append(dict(T=T,projected_probability=value,
            full_state_error_bound=error,
            full_probability_interval=probability_interval(value,error)))
    assert len(states)==39969 and G.sum()==66
    return dict(N=5,Q=10,budget=32,retained_dimension=len(states),
        outside_orbits=B.shape[0],nnz=H.nnz,boundary_nnz=B.nnz,
        target_orbits_by_vertex_fermions=np.bincount(nv[G],minlength=6).tolist(),
        time_steps=32,rows=rows)


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--allow-cache',action='store_true')
    ap.add_argument('--spectral-only',action='store_true')
    args=ap.parse_args()
    if not args.allow_cache:os.environ.pop('CP12_FORMATION_CACHE',None)
    start=time.monotonic()
    with contextlib.redirect_stdout(io.StringIO()):
        out=dict(spectral=spectral())
        if not args.spectral_only:out['capacity']=capacity()
    out['seconds']=time.monotonic()-start
    print(json.dumps(out,indent=2))
