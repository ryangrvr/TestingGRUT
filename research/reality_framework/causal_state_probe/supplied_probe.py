"""Independent causal-order -> retarded response -> SJ-state diagnostic.

This is a finite toy using a 1+1-dimensional *massless* scalar Green
kernel G_R(i,j)=C(i,j)/2.  It does not derive that kernel, a physical
quantum state, geometry dynamics, or gravitational observables.

Requirements: numpy.  Run: python grut_causal_response_state_probe.py
"""
import numpy as np

SEED=20261011
N=80
TOL=1e-8

def two_dimensional_diamond(n, seed):
    """Events uniform in null coordinates u,v on [0,1]^2.

    Event j precedes event i if u_j < u_i and v_j < v_i.
    """
    rng=np.random.default_rng(seed)
    uv=rng.uniform(size=(n,2))
    C=(uv[:,None,0]>uv[None,:,0]) & (uv[:,None,1]>uv[None,:,1])
    return C.astype(float)

def two_layer(n):
    assert n%2==0
    a=n//2
    C=np.zeros((n,n))
    C[a:,:a]=1.0
    return C


def four_dimensional_diamond(n,seed):
    """Uniform spacetime volume in the Alexandrov interval (0,0) -> (1,0)."""
    rng=np.random.default_rng(seed)
    t_half=0.5*rng.random(n)**0.25   # half-cone time density proportional to t^3
    t=np.where(rng.random(n)<0.5,t_half,1-t_half)
    R=np.minimum(t,1-t)
    vec=rng.normal(size=(n,3))
    vec/=np.linalg.norm(vec,axis=1)[:,None]
    x=vec*(R*rng.random(n)**(1/3))[:,None]
    dt=t[:,None]-t[None,:]
    dx2=np.sum((x[:,None,:]-x[None,:,:])**2,axis=-1)
    return ((dt>0)&(dt*dt>dx2)).astype(float)

def links(C):
    """Links (covering relations), not transitive closure."""
    return C*( (C@C)==0 )

def adversarial_4d_layers(n,a,b):
    C=np.zeros((n,n))
    C[a:a+b,:a]=1
    return C

def height(C):
    # Longest chain. A strictly causal order is acyclic; use topological sorting
    n=len(C)
    pending=list(np.sum(C,axis=1).astype(int))
    frontier=[i for i,d in enumerate(pending) if d==0]
    h=np.ones(n,dtype=int)
    seen=0
    while frontier:
        i=frontier.pop()
        seen+=1
        for child in np.flatnonzero(C[:,i]):
            h[child]=max(h[child],h[i]+1)
            pending[child]-=1
            if pending[child]==0:frontier.append(int(child))
    assert seen==n
    return int(h.max())

def analyze(C, propagator="2D"):
    n=len(C)
    assert np.all(np.diag(C)==0)
    # Verify transitivity: a relation via an intermediate must be present.
    assert np.all(((C@C)>0)<=C.astype(bool))
    G=0.5*C if propagator=="2D" else links(C)
    delta=G-G.T
    H=1j*delta  # Hermitian Pauli-Jordan matrix
    assert np.linalg.norm(H-H.conj().T)<1e-12
    eigen,U=np.linalg.eigh(H)
    W=(U*np.maximum(eigen,0))@U.conj().T
    # SJ Wightman: W = (|H|+H)/2, symmetrized noise Re(W).
    noise=W.real
    assert eigen.min()<=0 and eigen.max()>=0
    assert np.linalg.eigvalsh(W).min()>-1e-10
    assert np.linalg.eigvalsh(noise).min()>-1e-10
    assert np.linalg.norm(W-W.T-H)<1e-10
    spectral_rank=int(np.sum(np.abs(eigen)>TOL))
    min_els=np.where(np.sum(C,axis=1)==0)[0]
    if len(min_els)>=2:
        # The two minimal events with most divergent future sizes, as an
        # explicitly order-defined probe.  This is a purpose-built witness.
        out=np.sum(C[:,min_els],axis=0)
        i=int(min_els[np.argmin(out)])
        j=int(min_els[np.argmax(out)])
        if i==j:
            j=int(next(x for x in min_els if x!=i))
    else:
        i,j=0,1
    source=np.zeros(n);source[i]=1;source[j]=-1
    response_energy=float(np.linalg.norm(G@source))
    noise_quadratic=float(source@noise@source)
    assert noise_quadratic>=-1e-8
    return dict(pair_fraction=float(np.sum(C)/(n*(n-1)/2)),
                height=height(C), minimum_events=len(min_els),
                nonzero_commutator_eigenvalues=spectral_rank,
                sj_positive_modes=int(np.sum(eigen>TOL)),
                largest_abs_commutator_eigenvalue=float(np.max(np.abs(eigen))),
                test_probe=(i,j),
                future_profile_distance=response_energy,
                sj_difference_noise=noise_quadratic,
                noise_rank=int(np.sum(np.linalg.eigvalsh(noise)>TOL)),
                smallest_eigenvalue_W=float(np.linalg.eigvalsh(W).min()),
                commutator_reconstruction_error=float(np.linalg.norm(W-W.T-H)),
                C=C,G=G,W=W,noise=noise)

def run():
    results={"2D_sprinkling":analyze(two_dimensional_diamond(N,SEED)),
             "2_layer":analyze(two_layer(N))}
    for name,r in results.items():
        print(name)
        for k,v in r.items():
            if k not in {'C','G','W','noise'}:print('  {}: {}'.format(k,v))
    # Relabelling covariance: test an independent permutation, no re-fitting.
    r=results['2D_sprinkling']
    rng=np.random.default_rng(41)
    perm=rng.permutation(N)
    r2=analyze(r['C'][np.ix_(perm,perm)])
    assert abs(r['pair_fraction']-r2['pair_fraction'])<1e-15
    assert r['height']==r2['height']
    assert r['nonzero_commutator_eigenvalues']==r2['nonzero_commutator_eigenvalues']
    assert np.linalg.norm(r['W'][np.ix_(perm,perm)]-r2['W'])<1e-10
    print('permutation covariance: PASS')
    # Independent 4D benchmark: use link-matrix Green prescription, not 2D C/2.
    n4=180
    a4,b4=37,38
    geom4=analyze(four_dimensional_diamond(n4,SEED),"4D-link")
    fake4=analyze(adversarial_4d_layers(n4,a4,b4),"4D-link")
    print('4D_sprinkling_link_kernel')
    for k,v in geom4.items():
        if k not in {'C','G','W','noise'}: print('  {}: {}'.format(k,v))
    print('4D_layers_plus_isolates_link_kernel')
    for k,v in fake4.items():
        if k not in {'C','G','W','noise'}: print('  {}: {}'.format(k,v))
    assert fake4['nonzero_commutator_eigenvalues']==2
    assert fake4['sj_positive_modes']==1
    assert abs(fake4['pair_fraction']-geom4['pair_fraction'])<0.001
    assert geom4['nonzero_commutator_eigenvalues']>50
    print('4D link-kernel adversarial control: PASS')
    # Same commutator but a different positive quasi-free state.  Unlike W+I,
    # (1+lambda)W + lambda W.T preserves the kernel of the field commutator.
    lam=0.2
    alternative=(1+lam)*r['W']+lam*r['W'].T
    assert np.linalg.eigvalsh(alternative).min()>-1e-9
    assert np.linalg.norm(alternative-alternative.T-1j*(r['G']-r['G'].T))<1e-10
    assert np.linalg.norm(alternative-(1+2*lam)*r['noise']-0.5j*(r['G']-r['G'].T))<1e-10
    print('non-uniqueness of positive state at fixed commutator: PASS (W_mix=(1+λ)W+λW^T)')
    # Exact rank-two analytical control for complete bipartite causal relation:
    layer=results['2_layer']
    expected=np.sqrt((N//2)*(N//2))/2
    assert abs(layer['largest_abs_commutator_eigenvalue']-expected)<1e-9
    assert layer['nonzero_commutator_eigenvalues']==2
    assert layer['sj_positive_modes']==1
    print('exact bipartite rank/eigenvalue: PASS')

if __name__=='__main__':run()
