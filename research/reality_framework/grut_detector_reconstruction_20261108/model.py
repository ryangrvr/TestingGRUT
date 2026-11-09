"""B2: conservative quadratic reconstruction, without probability fitting."""
from pathlib import Path
import sys
import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / 'grut_detector_matching_20261107'
B1 = HERE.parent / 'grut_massless_dictionary_20261108'
sys.path.insert(0, str(OLD))
from continuum import all_histories, gaussian_histories, interval_spectrum, key, tv

def smearing(card):
    length, sigma = card['length'], card['sigma']
    norm = sigma*np.sqrt(2*np.pi)*(ndtr(length/(2*sigma))-ndtr(-length/(2*sigma)))
    return lambda z: np.exp(-.5*((z-length/2)/sigma)**2)/norm

def base_weights(card,n):
    length,sigma=card['length'],card['sigma']; a=length/(n+1)
    edges=np.arange(1,n+1)*a
    norm=ndtr(length/(2*sigma))-ndtr(-length/(2*sigma))
    weights=np.diff(ndtr((edges-length/2)/sigma))/norm
    return a,edges,weights,float(1-weights.sum())

def reconstruction(card,n):
    """Integrate local quadratics from exact average/value constraints only.

    First polynomial covers [0,2a]; last is reflected. Other cells use their
    own and two neighboring averages. No physical result enters the rule.
    """
    a,edges,_,_=base_weights(card,n); f=smearing(card); m=n-1
    weights=np.zeros(m); boundary=np.zeros(2)
    edge_matrix=np.array([[1,0,0],[1,1.5,7/3],[1,2.5,19/3]])
    interior_matrix=np.array([[1,r,r*r+1/12] for r in (-1,0,1)])
    quadrature_error=0.
    for side in (0,1):
        moments=[]
        for power in range(3):
            val,err=quad(lambda xi:a*f(a*xi if side==0 else card['length']-a*xi)*xi**power,
                         0,2,epsabs=1e-13,epsrel=1e-13)
            moments.append(val);quadrature_error+=err
        coeff=np.linalg.solve(edge_matrix.T,moments)
        boundary[side]+=coeff[0]
        indices=(0,1) if side==0 else (m-1,m-2)
        weights[list(indices)]+=coeff[1:]
    for j in range(1,m-1):
        center=j+1.5;moments=[]
        for power in range(3):
            val,err=quad(lambda xi:a*f(a*(center+xi))*xi**power,-.5,.5,epsabs=1e-13,epsrel=1e-13)
            moments.append(val);quadrature_error+=err
        weights[j-1:j+2]+=np.linalg.solve(interior_matrix.T,moments)
    return a,edges,weights,boundary,{
        'partition_error':float(abs(weights.sum()+boundary.sum()-1)),
        'reflection_error':float(np.max(abs(weights-weights[::-1]))),
        'weight_quadrature_error_estimate':float(quadrature_error),
        'min_weight':float(weights.min()),'boundary_coefficients':boundary.tolist()}

def spectrum_from_weights(card,n,modes,weights,boundary):
    a,edges,_,_=base_weights(card,n);length=card['length']
    k=np.arange(1,modes+1)*np.pi/length;omega=np.sqrt(k*k+1/np.pi)
    heights=np.asarray(weights)/a
    fn=np.sqrt(2/length)*(heights@(np.cos(edges[:-1,None]*k)-np.cos(edges[1:,None]*k)))/k
    mu=1/np.sqrt(np.pi)
    mean=float(np.sum(boundary))/4+float(heights@(np.sinh(mu*(edges[1:]-length/2))-
        np.sinh(mu*(edges[:-1]-length/2))))/(4*mu*np.cosh(mu*length/2))
    variation=float(abs(heights[0])+np.abs(np.diff(heights)).sum()+abs(heights[-1]))
    return omega,fn,mean,variation

def history_set(card,omega,fn,mean,variation):
    h2=fn*fn/(2*np.pi*omega)
    p,imag=all_histories(card,omega,h2,mean)
    tail=variation**2*card['length']**2/(2*np.pi**4*len(omega)**2)
    tails={key(t):float(2**len(t)*len(t)**2*card['coupling']**2*tail) for t in card['protocols']}
    return {'histories':p,'mean':float(mean),'variance':float(h2.sum()),
        'tail_variance_bound':float(tail),'history_tail_bounds':tails,
        'max_imaginary_residual':float(imag),'smearing_total_variation':variation}

def references(card,n,modes,diagnostic=False):
    omega,h2,fn,errors,meta=interval_spectrum(card,modes,card['quadrature_epsabs'])
    target=history_set(card,omega,fn,meta['mean'],meta['integration_by_parts_B'])
    target['covariance_quadrature_estimate']=float(meta['covariance_quad_error_estimate'])
    a,edges,v,b=base_weights(card,n)
    _,cellfn,cellmean,cellvar=spectrum_from_weights(card,n,modes,v,[b])
    out={'target':target,'B1':history_set(card,omega,cellfn,cellmean,cellvar)}
    if diagnostic:
        # Four references: neither approximation, interior only, edge only, both.
        f=smearing(card);k=np.arange(1,modes+1)*np.pi/card['length']
        edgefn=np.zeros(modes);error=np.zeros(modes)
        for i in range(0,modes,2):
            val,err=quad(f,0,a,weight='sin',wvar=k[i],epsabs=1e-13,epsrel=1e-13,limit=150)
            edgefn[i]=2*np.sqrt(2/card['length'])*val
            error[i]=2*np.sqrt(2/card['length'])*err
        mu=1/np.sqrt(np.pi);length=card['length']
        edgemean,meanerr=quad(lambda z:2*f(z)*np.cosh(mu*(z-length/2))/(4*np.cosh(mu*length/2)),
                             0,a,epsabs=1e-13,epsrel=1e-13)
        heights=v/a
        interiorvar=2*f(a)+2*abs(heights[0]-f(a))+float(abs(np.diff(heights)).sum())
        out['interior_only']=history_set(card,omega,cellfn+edgefn,cellmean-b/4+edgemean,interiorvar)
        out['edge_only']=history_set(card,omega,fn-edgefn,meta['mean']-edgemean+b/4,meta['integration_by_parts_B'])
        out['diagnostic_quad_estimates']={'max_edge_coefficient':float(error.max()),'edge_mean':float(meanerr)}
    else:
        _,_,w,boundary,checks=reconstruction(card,n)
        _,rfn,rmean,rvar=spectrum_from_weights(card,n,modes,w,boundary)
        out['B2']=history_set(card,omega,rfn,rmean,rvar)
        out['reconstruction']={'weights':w.tolist(),'boundary':boundary.tolist(),'checks':checks}
    return out

def distances(p,q):
    if set(p)!=set(q):raise ValueError('Protocol dictionaries differ')
    return {k:tv(p[k],q[k]) for k in p}

def signed_difference(p,q):
    return {k:{s:p[k][s]-q[k][s] for s in p[k]} for k in p}

def half_l1(diff):
    return {k:.5*sum(abs(v) for v in p.values()) for k,p in diff.items()}

def algebra_checks(card):
    polynomial_error=0.;quadrature_error=0.;partition=0.;reflection=0.
    for n in card['N']:
        a,edges,w,boundary,checks=reconstruction(card,n)
        partition=max(partition,checks['partition_error']);reflection=max(reflection,checks['reflection_error'])
        f=smearing(card);length=card['length']
        for p in range(3):
            averages=(edges[1:]**(p+1)-edges[:-1]**(p+1))/((p+1)*a*length**p)
            reconstructed=w@averages+boundary@np.array([float(p==0),1.])
            exact=quad(lambda z:f(z)*(z/length)**p,0,length,epsabs=1e-13,epsrel=1e-13)[0]
            polynomial_error=max(polynomial_error,abs(reconstructed-exact))
        omega,fn,_,_=spectrum_from_weights(card,n,16,w,boundary)
        for j in (1,3,8,15):
            k=j*np.pi/length
            direct=sum(wi/a*quad(lambda z:np.sqrt(2/length)*np.sin(k*z),lo,hi,epsabs=1e-13)[0]
                       for wi,lo,hi in zip(w,edges[:-1],edges[1:]))
            quadrature_error=max(quadrature_error,abs(direct-fn[j-1]))
    return {'quadratic_reproduction_error':float(polynomial_error),'partition_error':float(partition),
        'reflection_error':float(reflection),'sine_coefficient_quadrature_error':float(quadrature_error)}
