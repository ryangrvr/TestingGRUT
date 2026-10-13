"""Exact bare neutral-rewiring spectra and conditional vacancy propagation.

This isolates an unchanged term of CP-12. It is not a full-law evolution:
fermion hopping leaves the neutral subspace.
"""
from pathlib import Path
from fractions import Fraction
import itertools
import json
import math

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import expm_multiply

P=Path(__file__).resolve().parent


def spectrum(M,Q):
    assert Q<=M//2
    total=math.comb(M,Q)
    weights=[Fraction(math.comb(M,j)-(math.comb(M,j-1) if j else 0),total) for j in range(Q+1)]
    energies=[Fraction(j*(M+1-j),10*M) for j in range(Q+1)]
    assert sum(weights)==1
    return energies,weights


def small_control(n,Q):
    M=n*(n-1)//2;subsets=list(itertools.combinations(range(M),Q))
    masks=[sum(1 << e for e in es) for es in subsets];index={mask:i for i,mask in enumerate(masks)}
    rows=[];cols=[];data=[]
    for i,mask in enumerate(masks):
        for e in range(M):
            if not(mask >> e & 1):continue
            for f in range(M):
                if mask >> f & 1:continue
                rows.append(index[mask ^ (1 << e) ^ (1 << f)]);cols.append(i);data.append(-.1/M)
        rows.append(i);cols.append(i);data.append(.1/M*Q*(M-Q))
    H=coo_matrix((data,(rows,cols)),shape=(len(masks),)*2).tocsr()
    energies,weights=spectrum(M,Q); initial=np.zeros(len(masks));initial[0]=1
    checks=[]
    for T in [.1,1.,10.,50.,100.]:
        evolved=expm_multiply(-1j*T*H,initial)
        exact=sum(float(w)*np.exp(-1j*float(E)*T) for w,E in zip(weights,energies))
        assert abs(evolved[0]-exact)<2e-11
        checks.append(dict(T=T,survival_probability=float(abs(evolved[0])**2),
                           spectral_survival_probability=float(abs(exact)**2),amplitude_error=float(abs(evolved[0]-exact))))
    return dict(N=n,Q=Q,dimension=len(masks),checks=checks)


def run():
    small=[small_control(4,3),small_control(5,4),small_control(6,6)]
    uniform=[]
    for n in [16,32,128,512,1024,8192,32768]:
        M=n*(n-1)//2;Q=2*n;alpha=Fraction(Q,M-Q+1)
        assert alpha<=Fraction(1,2)
        maximum_departure=4*alpha*(1-alpha)
        uniform.append(dict(N=n,Q=Q,highest_sector_weight=str(1-alpha),
                            exact_all_time_departure_upper_bound=str(maximum_departure),
                            all_time_departure_upper_bound=float(maximum_departure)))
    vacancy=[]
    for L in [4,8,16,32,64,128,256]:
        n=2*L**3;M=n*(n-1)//2;Q=2*n;S=Q+1
        gap=Fraction(S,10*M)
        vacancy.append(dict(supplied_host_linear_size=L,N=n,Q=Q,pool_edges=S,
                 exact_nonuniform_vacancy_gap=str(gap),gap=float(gap),
                 all_nonuniform_momenta_have_same_energy=True,
                 exact_all_time_localized_vacancy_departure_bound=str(Fraction(4*(S-1),S*S)),
                 observation_time=L,accumulated_nonuniform_phase=float(gap*L)))
    # Independently diagonalize the one-vacancy COMPLETE-GRAPH kinetic kernel.
    finite=[]
    for S in [7,13,37]:
        M=4*S;H=.1/M*(S*np.eye(S)-np.ones((S,S)))
        eig=np.linalg.eigvalsh(H)
        assert abs(eig[0])<1e-12 and np.max(abs(eig[1:]-.1*S/M))<1e-12
        initial=np.eye(S)[:,0]
        for T in [.4,10.,50.]:
            evolved=expm_multiply(-1j*T*H,initial)
            amplitude=np.exp(-1j*.1*S/M*T)*initial+(1-np.exp(-1j*.1*S/M*T))/S
            assert np.linalg.norm(evolved-amplitude)<1e-11
            finite.append(dict(S=S,T=T,departure=1-float(abs(evolved[0])**2),
                               uniform_all_time_bound=4*(S-1)/S**2))
    out=dict(status="CP12_NEUTRAL_SPECTRUM_CHECKS_PASS",independent_sparse_controls=small,
             complete_graph_survival_bounds=uniform,conditional_vacancy_spectra=vacancy,
             independently_diagonalized_vacancy_controls=finite,
             scope="Frozen neutral rewiring term and conditional pool compression only. Matter/gauge hopping destroys invariance of this sector; no all-time full-CP12 claim.")
    (P/"neutral_spectrum_results.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(dict(status=out['status'],largest_control=small[-1]['dimension'],
                          large_N_departure_bound=uniform[-1]['all_time_departure_upper_bound']),indent=2))


if __name__=="__main__":run()
