"""Integer-certified global norm bound for any weighted CP-12 hopping.

Physical hop weights must be Hermitian and have magnitude at most one.
Gauss implies that a sector with k vertex fermions has exactly k active
edge holes. A small Jacobi majorant bounds all electric fields at once.
"""
from math import comb


def endpoint_maximum(n,k,edges):
    """Maximum incidences in a k-vertex set among 'edges' distinct pairs."""
    within=min(edges,comb(k,2))
    mixed=min(max(0,edges-within),k*(n-k))
    return 2*within+mixed


def hopping_bound(n,q):
    r=min(n,q)
    squares=[endpoint_maximum(n,n-k,q-k)*endpoint_maximum(n,k+1,k+1)
        for k in range(r)]
    for c in range(1,2*q+1):
        minors=[1,c]
        for a2 in squares:minors.append(c*minors[-1]-a2*minors[-2])
        # Sylvester: c I-J positive definite, exactly in integer arithmetic.
        if all(d>0 for d in minors):return c,squares,minors
    # The usual physical row/column bound remains valid if strict positivity
    # happens not to be certified before that bound.
    return 2*q,squares,None


if __name__=='__main__':
    import json
    print(json.dumps({str(n):dict(zip(['bound','Jacobi_offdiagonal_squares','leading_minors'],hopping_bound(n,q)))
        for n,q in [(4,5),(5,10),(6,12)]},indent=2))
