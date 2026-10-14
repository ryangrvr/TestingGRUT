"""Exact finite-matrix certificates for the UNBOUNDED rotor Hamiltonian.

The analytic tail bound is derived in DERIVATION.md. Rational Schur elimination
proves a lower bound; an integer-amplitude trial vector proves an upper bound.
Floating eigensolves only suggest that trial vector, not its certified energy.
The default small controls do not certify the N=5 or N=6 cutoff extrapolation.
"""
from fractions import Fraction as F
from pathlib import Path
import json
import math
import numpy as np
from frozen_model import build, gauss
from exact_destination import transitions

P=Path(__file__).resolve().parent

def tail_lower(n,L):
    # sqrt(3) < 1.732051, checked with integers; no floating rounding premise.
    root_upper=F(1732051,1000000)
    assert root_upper*root_upper>3
    return n*(F(3,2)-root_upper)+F((L+1)**2+L*L,8)

def boundary_data(n,L,states):
    scale=math.lcm(8,10*(n*(n-1)//2));outside={};counts=[0]*len(states)
    for i,s in enumerate(states):
        for target,c in transitions(n,s,scale):
            if max(map(abs,target[2]))<=L:continue
            assert abs(c)==scale and target[0]==s[0]
            assert not np.any(gauss(n,target))
            bad=[j for j,x in enumerate(target[2]) if abs(x)>L]
            assert len(bad)==1 and target[2][bad[0]]==L+1
            # The unique out-of-cutoff leg fixes the inverse fermionic hop.
            assert target not in outside
            outside[target]=(i,c//scale);counts[i]+=1
    return counts,outside

def integer_matrix(n,L,states,scale):
    ix={s:i for i,s in enumerate(states)};A=[[0]*len(states) for _ in states]
    for j,s in enumerate(states):
        for t,c in transitions(n,s,scale):
            if t in ix:A[ix[t]][j]+=c
    assert all(A[i][j]==A[j][i] for i in range(len(A)) for j in range(i))
    return A

def positive_ldl(A):
    """Exact rational Schur pivots; a positive sequence proves positive definite."""
    pivots=[]
    for k in range(len(A)):
        pivot=A[k][k];assert pivot>0,('nonpositive pivot',k,str(pivot))
        pivots.append(pivot)
        for i in range(k+1,len(A)):
            aik=A[i][k]
            if aik:
                for j in range(i,len(A)):
                    A[j][i]-=aik*A[j][k]/pivot
                    if j!=i:A[i][j]=A[j][i]
        for i in range(k+1,len(A)):A[i][k]=A[k][i]=F(0)
    return pivots

def certify(n,L,z):
    H,states,_,_=build(n,n,L);scale=math.lcm(8,10*(n*(n-1)//2))
    A=integer_matrix(n,L,states,scale)
    assert np.max(abs(np.array(A)/scale-H.toarray()))<1e-12
    count,outside=boundary_data(n,L,states)
    # Independently build B and verify B^T B is the claimed diagonal.
    B=np.zeros((len(outside),len(states)),int)
    for row,(i,c) in enumerate(outside.values()):B[row,i]=c
    assert np.array_equal(B.T@B,np.diag(count))
    c=tail_lower(n,L);assert c>z
    S=[[F(a,scale) for a in row] for row in A]
    for i,d in enumerate(count):S[i][i]-=z+F(d)/(c-z)
    pivots=positive_ldl(S)
    _,vec=np.linalg.eigh(H.toarray());trial=[int(round(x*10**8)) for x in vec[:,0]]
    norm=sum(x*x for x in trial)
    numerator=sum(trial[i]*A[i][j]*trial[j] for i in range(len(A)) for j in range(len(A)))
    upper=F(numerator,scale*norm)
    assert upper>=z
    return dict(N=n,Q=n,L=L,dimension=len(states),boundary_states=len(outside),
        tail_lower_rational=str(c),lower_rational=str(z),upper_rational=str(upper),
        lower=float(z),upper=float(upper),interval_width=float(upper-z),
        all_rational_pivots_positive=True,minimum_pivot=float(min(pivots)),
        unique_outside_predecessor=True,
        scope='Unregulated infinite-rotor ground energy at this N; exact rational certificate')

def main():
    rows=[certify(3,2,F(-467300,1000000)),certify(3,3,F(-46729956,100000000))]
    result=dict(status='CP12_UNREGULATED_ENERGY_CERTIFICATES_PASS',certificates=rows,
        limitation='These small exact certificates do not establish the N=5 or N=6 infinite-cutoff result.')
    (P/'cutoff_certificates.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

if __name__=='__main__':main()
