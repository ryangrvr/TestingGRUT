"""Frozen B1 controls. Run with single-threaded BLAS; <=60s, <=1GiB RSS.

Exact rational scalar checks are separate from ordinary double matrix tests.
No finite calculation establishes the full-geometry or many-size claims.
"""
import hashlib, json, math, os, pathlib, resource, time
from fractions import Fraction as F
import numpy as np
from scipy.linalg import eigh, norm
from scipy.special import jv
from model import PathModel, JP, JM, SELF, JPAIR, OMEGA

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def rational_checks():
    tau=F(5,96); gap=F(5,4); a=F(13,46080); j=F(1,3840)
    k=2*a+2*j; r=tau/gap; u=r**3+tau*k/gap**2
    p=tau*u; q=p+r*r*k+u*k; rho2=p*p+q*q
    assert rho2<F(1,105000)**2
    e0=j*tau/gap**2
    assert e0==F(1,115200)
    assert e0+F(15360,105000)<F(147,1000)
    assert e0+F(15360,105000)+F(8,1000)+F(2,10**18)<F(31,200)
    intervals=[]; unarrived=F(0)
    for n in range(1,6):
        total=sum((F((-1)**k*4**(2*k+n),math.factorial(k)*math.factorial(k+n))
                   for k in range(81)),F(0))
        nxt=F(4**(162+n),math.factorial(81)*math.factorial(81+n))
        lo,hi=total-nxt,total+nxt
        assert nxt<F(1,10**70)
        unarrived+=F(n,4)**2*max(lo*lo,hi*hi)
        intervals.append({'order':n,'lo':float(lo),'hi':float(hi),'tail_upper':float(nxt)})
    assert unarrived<F(1,8)
    c=F(31,200); b=F(7,8); target=F(3,5)
    u=b+c*c-target
    assert u>0 and u*u>4*c*c*b
    return {'rho_upper':1/105000,'rho_formula':math.sqrt(float(rho2)),
            'initial_error_upper':float(e0),'time_error_plus_parameter_end_upper':
            float(e0+F(15360,105000)+F(8,1000)+F(2,10**18)),
            'bessel_intervals':intervals,'unarrived_upper':float(unarrived),
            'matter_uniform_lower':float(1-(F(1,2)+c)**2),
            'arrival_lower':(math.sqrt(float(b))-float(c))**2,
            'arithmetic':'Fraction integer/rational tests; displayed decimals are not interval endpoints'}


def one_resource(length):
    model=PathModel(length); v=model.initial(); h=model.H.toarray()
    e,u=eigh(h); coeff=u.T@v
    data=model.intertwiner(); k=data['K']; z=data['Z']; ids=data['ids']
    ek,uk=eigh(2*np.eye(len(ids))+k); ck=uk.T@np.eye(len(ids))[:,0]
    tau=5/96; gap=5/4; kn=2*SELF+2*JPAIR
    uu=(tau/gap)**3+tau*kn/gap**2; pp=tau*uu
    qq=pp+(tau/gap)**2*kn+uu*kn
    rho=math.hypot(pp,qq); e0=JPAIR*tau/gap**2
    expected=-2*SELF*np.eye(len(ids))
    for x in range(len(ids)-1): expected[x,x+1]=expected[x+1,x]=-JPAIR
    assert norm(k-expected,2)<1e-14
    assert norm(data['residual']-data['formula'],2)<1e-14
    assert norm(data['parity_zero'],2)<1e-14
    assert norm(data['residual'],2)<=rho+1e-14
    points=[]
    initial_energy=float(np.vdot(v,h@v).real)
    for theta in [0,1,2,3,4]:
        t=theta/JPAIR
        psi=u@(np.exp(-1j*e*t)*coeff)
        phi=uk@(np.exp(-1j*ek*t)*ck)
        err=float(norm(psi-z@phi)); budget=e0+rho*t
        assert err<=budget+1e-9
        fuel=float(np.sum(abs(psi[model.fuels==1])**2))
        arrived=[ix for ix,(p,m,f) in enumerate(model.states) if len(p)==1 and p==m and p[0]>=4]
        arrival=float(np.sum(abs(psi[arrived])**2))
        edrift=abs(float(np.vdot(psi,h@psi).real)-initial_energy)
        assert edrift<1e-10 and abs(norm(psi)-1)<1e-10
        points.append({'theta':theta,'time':t,'matter_probability':1-fuel,
                       'co_pair_arrival_distance_ge4':arrival,'full_state_error_to_Zwalk':err,
                       'analytic_error_budget':budget,'energy_drift':edrift})
    grid=None
    if length==16:
        vals=[]
        for theta in np.linspace(2,4,41):
            psi=u@(np.exp(-1j*e*(theta/JPAIR))*coeff)
            vals.append(float(np.sum(abs(psi[model.fuels==0])**2)))
        grid={'points':41,'minimum_matter_probability':min(vals),
              'scope':'finite control only, not continuous-time proof'}
        # Color-blind physical matrix is H tensor I2. Check coherent input
        # and a local arrival color readout without conditional renormalizing.
        t=4/JPAIR; psi=u@(np.exp(-1j*e*t)*coeff)
        spin=np.array([math.sqrt(.3),1j*math.sqrt(.7)])
        full=np.kron(psi,spin)
        hc=np.kron(h,np.eye(2))
        assert norm(hc@full-np.kron(h@psi,spin))<1e-12
        success=np.array([int(len(p)==1 and p==m and p[0]>=4) for p,m,f in model.states])
        out=full.reshape(model.dim,2)[success.astype(bool)]
        output=out.T@out.conj()
        pa=float(np.sum(abs(psi[success.astype(bool)])**2))
        assert norm(output-pa*np.outer(spin,spin.conj()))<1e-12
    return {'length':length,'dimension':model.dim,'points':points,'grid':grid,
            'residual_norm':float(norm(data['residual'],2)),
            'residual_formula_error':float(norm(data['residual']-data['formula'],2)),
            'K_identity_error':float(norm(k-expected,2))}


def two_cells():
    model=PathModel(8,sites=(0,7),resource=2)
    h=model.H.toarray(); e,u=eigh(h); v=model.initial(); coeff=u.T@v
    initial=float(np.vdot(v,h@v).real)
    observations=[]
    for theta in [0,1,2,3,4]:
        psi=u@(np.exp(-1j*e*(theta/JPAIR))*coeff); prob=abs(psi)**2
        drift=abs(float(np.vdot(psi,h@psi).real)-initial)
        assert drift<1e-10 and abs(np.sum(prob)-1)<1e-10
        observations.append({'theta':theta,'time':theta/JPAIR,
            'at_least_one_pair_probability':float(np.sum(prob[model.pairs>=1])),
            'two_pairs_probability':float(np.sum(prob[model.pairs==2])),
            'expected_pairs':float(prob@model.pairs),'energy_drift':drift})
    # A fixed exact state already exhibits collision energy, invalidating
    # independent-cell tensorization of the one-resource reduction.
    near=model.index[((2,3),(2,3),0)]; far=model.index[((1,5),(1,5),0)]
    collision=float(model.D[near,near]-model.D[far,far])
    assert abs(collision-.5)<1e-13
    return {'length':8,'sites':[0,7],'dimension':model.dim,'points':observations,
            'adjacent_minus_separated_diagonal_energy':collision,
            'scope':'exact finite common-path two-resource diagnostic; not a many-cell theorem'}


def main():
    started=time.perf_counter()
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())['sources']
    unchanged={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==sha for name,sha in pins.items()}
    assert all(unchanged.values())
    result={'scalar':rational_checks(),'one_resource':[one_resource(n) for n in [2,8,16]],
            'two_cells':two_cells(),'baseline_pins_unchanged':unchanged}
    result['elapsed_seconds']=time.perf_counter()-started
    result['peak_rss_MiB']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
    assert result['elapsed_seconds']<60 and result['peak_rss_MiB']<1024
    result['thread_settings']={key:os.environ.get(key) for key in
                              ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']}
    print(json.dumps(result,indent=2))


if __name__=='__main__': main()
