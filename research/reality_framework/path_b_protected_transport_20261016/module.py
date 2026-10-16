"""Exact two-band walk mediator and its independently checkable bounds.

Graphs here are supplied fixtures, not claimed dynamically formed space.
"""
from __future__ import annotations
import math
import numpy as np
from scipy import sparse
from scipy.linalg import eigh


def cycle(n):
    rows = np.arange(n)
    a = sparse.coo_matrix((np.ones(2*n),
        (np.repeat(rows, 2), np.column_stack(((rows-1)%n, (rows+1)%n)).ravel())),
        shape=(n, n)).tocsr()
    return a


def torus(side):
    n = side*side
    rr, cc = [], []
    for x in range(side):
        for y in range(side):
            i = x*side+y
            for xx, yy in [((x-1)%side,y), ((x+1)%side,y),
                           (x,(y-1)%side), (x,(y+1)%side)]:
                rr.append(i); cc.append(xx*side+yy)
    return sparse.coo_matrix((np.ones(len(rr)),(rr,cc)),shape=(n,n)).tocsr()


def distances(a, source=0):
    return sparse.csgraph.shortest_path(a, directed=False,
        unweighted=True, indices=source).astype(int)


def microscopic(a, gap=8., v=1., g=.25):
    n = a.shape[0]; eye = sparse.eye(n, format='csr')
    return sparse.bmat([[None,-g*eye],[-g*eye,gap*eye-v*a]], format='csr')


def exact_low_band(a, gap=8., v=1., g=.25):
    av, u = eigh(a.toarray())
    d = gap-v*av
    # Stable negative root, avoiding cancellation in sqrt(d*d+4*g*g)-d.
    f = -2*g*g/(d+np.sqrt(d*d+4*g*g))
    x = -f/g
    c = 1/np.sqrt(1+x*x); s = x*c
    fmat = (u*f)@u.T
    w = np.vstack(((u*c)@u.T,(u*s)@u.T))
    return fmat, w, dict(low_min=float(f.min()), low_max=float(f.max()),
        high_min=float(((d+np.sqrt(d*d+4*g*g))/2).min()),
        exact_projector_distance=float(s.max()))


def constants(z=4, gap=8., v=1., g=.25, mu=math.log(1.5)):
    d0 = gap-abs(v)*z
    dm = gap-abs(v)*z*math.exp(mu)
    assert d0 > 2*abs(g) and dm > 2*abs(g)
    b = 2*g*g/(dm+math.sqrt(dm*dm-4*g*g))
    b0 = 2*g*g/(d0+math.sqrt(d0*d0-4*g*g))
    correction = b0-g*g/d0
    x = b/abs(g)
    return dict(z=z, Delta=gap, v=v, g=g, mu=mu, xi=1/mu,
        Q_gap_lower=d0, low_band_bottom_lower=-g*g/d0,
        spectral_gap_lower=d0, projector_distance_upper=abs(g)/d0,
        bare_high_probability_upper=4*g*g/(d0*d0),
        feshbach_E0_C=g*g/d0, feshbach_E0_ratio=abs(v)*z/gap,
        dressed_C=b, dressed_near_lower=g*g*abs(v)/(gap*gap)-correction,
        dressing_L_minus_identity_weighted_upper=1/math.sqrt(1-x*x)-1,
        dressing_H_weighted_upper=x/math.sqrt(1-x*x),
        microscopic_propagation_exponent_rate=z*abs(v)*math.sinh(mu))


def ring_green(n, d, a=8., v=1.):
    root=math.sqrt(a*a-4*v*v)
    r=2*v/(a+root)
    return (r**d+r**(n-d))/(root*(1-r**n))


def latch(a, j=None, gap=8., delta=8., v=1., g=.25, lam=.25):
    """Order: (geometry 0,L), (0,H), (1,L), (1,H). All terms retained."""
    n=a.shape[0]; eye=sparse.eye(n,format='csr')
    j=n//2 if j is None else j
    shortcut=sparse.coo_matrix(([1.,1.],([0,j],[j,0])),shape=(n,n)).tocsr()
    return sparse.bmat([
        [None,-g*eye,-lam*eye,None],
        [-g*eye,gap*eye-v*a,None,-lam*eye],
        [-lam*eye,None,delta*eye,-g*eye],
        [None,-lam*eye,-g*eye,(gap+delta)*eye-v*(a+shortcut)]
    ],format='csr')


def latch_constants(n, gap=8., delta=8., v=1., g=.25, lam=.25):
    # Gershgorin/Schur bound deliberately conservative and size independent.
    m=min(gap,delta)-(3*abs(v)+abs(g)+abs(lam))
    assert m>0
    emin=-(g*g+lam*lam)/m
    floor=(g*g*lam*lam*abs(v)/(gap+delta-emin)**2
        *(1/(gap-emin)+1/(delta-emin))**2)
    return dict(Q_gap_lower=m, full_spectral_gap_lower=m,
        low_bottom_lower=emin, shortcut_floor_all_low_energies=floor,
        bare_P0_total_Q_probability_upper=4*(g*g+lam*lam)/(m*m+4*(g*g+lam*lam)),
        shortcut_five_edge_lower_E0=(g*g*lam*lam*abs(v)/(gap+delta)**2
            *(1/gap+1/delta)**2),
        nearest_lower_all_low_energies=g*g*abs(v)/(gap-emin)**2,
        N=n, distance=n//2)
