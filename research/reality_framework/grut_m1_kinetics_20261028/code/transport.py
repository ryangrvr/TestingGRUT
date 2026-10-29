"""Conservative isotropic fixed-core kinetic diagnostic. No energy sink."""
from setup import *
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_ivp
from scipy import sparse
from numpy.polynomial.legendre import leggauss
import argparse,time

def grid(n,emax,kind='quadratic',shift_bins=64):
    if kind=='translated':
        # A graded low-energy block repeated at exact 2*Omega intervals.
        # Forward/reverse exchange is a permutation on supported nodes;
        # the quadratic first block also resolves the nonrelativistic endpoint.
        block=2*OMEGA*(np.arange(shift_bins)/shift_bins)**2
        e=np.concatenate([1.+1e-6+block+2*OMEGA*j for j in range(int(np.ceil((emax-1)/(2*OMEGA)))+1)])
        return e[e<=emax+1e-12]
    if kind=='aligned':
        step=2*OMEGA/shift_bins
        return 1.+1e-6+step*np.arange(int(np.ceil((emax-1.-1e-6)/step))+1)
    return 1.+1e-6+(emax-1.-1e-6)*np.linspace(0,1,n)**2

def deposit(e,values,weights):
    """Linear number/energy preserving projection; caller handles support."""
    values=np.atleast_1d(values);weights=np.broadcast_to(weights,values.shape)
    if np.any(values<e[0]-1e-11) or np.any(values>e[-1]+1e-11):
        raise ValueError('Unsupported daughter energy')
    values=np.clip(values,e[0],e[-1])  # roundoff at an already-checked boundary only
    j=np.clip(np.searchsorted(e,values)-1,0,len(e)-2)
    f=(values-e[j])/(e[j+1]-e[j]);out=np.zeros(len(e))
    np.add.at(out,j,weights*(1-f));np.add.at(out,j+1,weights*f)
    return out

def uniform(e,lo,hi):
    """Integral of linear nodal basis over a uniform two-body energy interval."""
    if hi-lo<1e-12:return deposit(e,[(lo+hi)/2],[1.])
    if lo<e[0]-1e-11 or hi>e[-1]+1e-11: raise ValueError('Uniform source outside support')
    left=np.maximum(e[:-1],lo);right=np.minimum(e[1:],hi)
    width=np.maximum(right-left,0);mid=(right+left)/2
    out=np.zeros(len(e));out[:-1]+=width*(e[1:]-mid)/np.diff(e)/(hi-lo)
    out[1:]+=width*(mid-e[:-1])/np.diff(e)/(hi-lo)
    return out

def daughters(e,eh):
    half=np.sqrt(max(eh*eh-6,0))*np.sqrt(1/3)/2
    lo,hi=eh/2-half,eh/2+half
    lower=max(lo,e[0],eh-e[-1]);upper=min(hi,e[-1],eh-e[0])
    if upper<=lower:return np.zeros(len(e)),0.
    fraction=(upper-lower)/(hi-lo) if hi>lo else 1.
    return uniform(e,lower,upper),fraction

def initial(e):
    z=np.load(ROOT/'SOURCE.npz');plus=np.zeros(len(e));minus=plus.copy();outside={}
    for s,target in [('plus',plus),('minus',minus)]:
        v,w=z[s+'_E'],z[s+'_N'];ok=(v>=e[0])&(v<=e[-1])
        target+=deposit(e,v[ok],w[ok])
        outside[s]={'number':float(w[~ok].sum()),'energy':float(w[~ok]@v[~ok])}
    nmiss=emiss=0.
    for eh,w in zip(z['neutral_E'],z['neutral_N']):
        if w==0:continue
        d,fraction=daughters(e,eh)
        nmiss+=w*(1-fraction);emiss+=w*eh*(1-fraction)
        plus+=w*fraction*d;minus+=w*fraction*d
    outside['neutral']={'number':nmiss,'energy':emiss}
    return np.r_[plus,minus],outside

def core_operator(e):
    z=np.load(ROOT/'CORE_KERNELS.npz');ks=z['k'];w=z['W']
    interp=PchipInterpolator(ks,w,axis=0,extrapolate=False)
    n=len(e);L=np.zeros((2*n,2*n));omitted=np.zeros(2*n)
    rates=np.zeros((4,2*n))
    def get(k):
        if k<ks[0]:return interp(ks[0])*k/ks[0]
        if k>ks[-1]:raise ValueError('Core kernel support exceeded')
        return interp(k)
    for species in (0,1):
        for i,energy in enumerate(e):
            col=species*n+i;p=np.sqrt(energy*energy-1)
            ea=energy if species==1 else energy-2*OMEGA
            if ea<=1:continue
            k=np.sqrt(ea*ea-1);weights=get(k)
            # Symmetrize only roundoff reciprocity error, not physical rates.
            weights=(weights+weights.T)/2
            exchange=np.pi*weights[0,1]/(p*energy)
            eh=energy+(OMEGA if species==1 else -OMEGA)
            neutral=np.pi*weights[2,species]/(p*energy) if eh>MH else 0.
            eo=energy+(2*OMEGA if species==1 else -2*OMEGA)
            if e[0]-1e-11<=eo<=e[-1]+1e-11:
                L[(1-species)*n:(2-species)*n,col]+=exchange*deposit(e,[eo],[1.])
                L[col,col]-=exchange;rates[species,col]=exchange
            else:omitted[col]+=exchange
            if neutral:
                d,fraction=daughters(e,eh);omitted[col]+=neutral*(1-fraction)
                neutral*=fraction
                L[:n,col]+=neutral*d;L[n:,col]+=neutral*d;L[col,col]-=neutral
                rates[2+species,col]=neutral
    invariant=np.r_[e-OMEGA,e+OMEGA]
    err=float(np.max(abs(invariant@L)))
    return L,omitted,rates,err

def resonant_operator(e):
    n=len(e);p=np.sqrt(e*e-1);ii=[];jj=[];cols=[];ks=[];omit=[]
    for i in range(n):
        for j in range(n):
            mu=(e[i]*e[j]-2)/(p[i]*p[j])
            if abs(mu)>=1:continue
            # Integral over isotropic incident directions, using Gamma(h->XXbar).
            k=np.pi*GH**2/(8*e[i]*e[j]*p[i]*p[j])
            total=e[i]+e[j];half=np.sqrt(total*total-6)*np.sqrt(1/3)/2
            lo,hi=total/2-half,total/2+half
            lower=max(lo,e[0],total-e[-1]);upper=min(hi,e[-1],total-e[0])
            fraction=max(0,upper-lower)/(hi-lo)
            if fraction==0:continue
            ii.append(i);jj.append(j);ks.append(k*fraction);omit.append(k*(1-fraction))
            cols.append(uniform(e,lower,upper))
    i=np.array(ii);j=np.array(jj);K=np.array(ks);Gmat=np.array(cols).T
    Dp=Gmat.copy();Dm=Gmat.copy();Dp[i,np.arange(len(i))]-=1;Dm[j,np.arange(len(i))]-=1
    D=np.vstack((Dp,Dm))
    energy=np.r_[e,e]@D
    return {'i':i,'j':j,'K':K,'D':D,'omit':np.array(omit),
            'energy_error':float(np.max(abs(energy))),
            'number_error':float(max(np.max(abs(Dp.sum(axis=0))),np.max(abs(Dm.sum(axis=0)))))}

def elastic_operator(e,nangle):
    """Identical charged bosons: quartic plus t,u neutral exchange.

    Integrate initial relative direction, CM scattering angle and azimuth.
    Initial unordered-pair counting is separate from identical final-state 1/2.
    Events with either daughter beyond support are omitted and tallied.
    """
    n=len(e);mu,wm=leggauss(nangle);ct,wt=leggauss(nangle)
    phi=2*np.pi*(np.arange(nangle)+.5)/nangle
    cosphi=np.cos(phi);p=np.sqrt(e*e-1);ii,jj=np.triu_indices(n)
    D=np.zeros((n,len(ii)));omit=np.zeros(len(ii));totalrate=np.zeros(len(ii))
    start=time.monotonic()
    for col,(i,j) in enumerate(zip(ii,jj)):
        s=2+2*(e[i]*e[j]-p[i]*p[j]*mu)
        p2=np.maximum(s/4-1,0);beta=np.sqrt(np.maximum(1-4/s,0))
        momentum=np.sqrt(np.maximum((e[i]+e[j])**2-s,0))
        boost=momentum*beta
        ca=np.divide(e[i]-e[j],boost,out=np.zeros_like(boost),where=boost>1e-14)
        ca=np.clip(ca,-1,1);sa=np.sqrt(1-ca*ca)
        t=-2*p2[:,None]*(1-ct[None,:]);u=-2*p2[:,None]*(1+ct[None,:])
        amp=4*.001+GH**2*(1/(t-6)+1/(u-6))
        vm=np.sqrt(np.maximum(s*(s-4),0))/(2*e[i]*e[j])
        weight=(wm[:,None]/2)*wt[None,:]*vm[:,None]*amp**2/(64*np.pi*s[:,None])
        cosine=ca[:,None,None]*ct[None,:,None]+sa[:,None,None]*np.sqrt(1-ct**2)[None,:,None]*cosphi[None,None,:]
        e3=(e[i]+e[j])/2+boost[:,None,None]*cosine/2;e4=e[i]+e[j]-e3
        weights=np.broadcast_to(weight[:,:,None]/nangle,e3.shape).ravel();e3=e3.ravel();e4=e4.ravel()
        ok=(e3>=e[0])&(e3<=e[-1])&(e4>=e[0])&(e4<=e[-1]);rate=weights[ok].sum()
        D[:,col]=deposit(e,e3[ok],weights[ok])+deposit(e,e4[ok],weights[ok])
        D[i,col]-=rate;D[j,col]-=rate;omit[col]=weights[~ok].sum();totalrate[col]=weights.sum()
    return {'i':ii,'j':jj,'D':D,'omit':omit,'K':totalrate,
            'energy_error':float(np.max(abs(e@D))),
            'number_error':float(np.max(abs(D.sum(axis=0)))),
            'build_seconds':time.monotonic()-start}

class Model:
    def __init__(self,e,mode,nangle=8):
        self.e=e;self.n=len(e);self.mode=mode
        self.L,self.omit_core,self.core_rates,self.core_error=core_operator(e)
        self.res=resonant_operator(e) if mode!='core' else None
        self.el=elastic_operator(e,nangle) if mode=='full' else None
        self.dim=2*self.n
    def rates(self,y,derivative=False):
        n=self.n;v=y[:2*n];plus=v[:n];minus=v[n:];dy=self.L@v
        lost=float(self.omit_core@v);J=self.L.copy() if derivative else None
        if self.res is not None:
            r=self.res;i,j=r['i'],r['j'];rate=r['K']*plus[i]*minus[j]
            dy+=r['D']@rate;lost+=float(r['omit']@(plus[i]*minus[j]))
            if derivative:
                ar=np.arange(len(i));d=sparse.csr_matrix((np.r_[r['K']*minus[j],r['K']*plus[i]],
                   (np.r_[ar,ar],np.r_[i,n+j])),shape=(len(i),2*n))
                J+=r['D']@d
        if self.el is not None:
            r=self.el;i,j=r['i'],r['j'];factor=np.where(i==j,.5,1.)
            for offset,x in ((0,plus),(n,minus)):
                rate=factor*x[i]*x[j];dy[offset:offset+n]+=r['D']@rate
                lost+=float(r['omit']@rate)
                if derivative:
                    ar=np.arange(len(i));d=sparse.csr_matrix((np.r_[factor*x[j],factor*x[i]],
                        (np.r_[ar,ar],np.r_[i,j])),shape=(len(i),n))
                    J[offset:offset+n,offset:offset+n]+=r['D']@d
        if derivative:return J
        return dy,lost
    def rhs(self,t,y):
        dy,lost=self.rates(y);return np.r_[dy,lost]
    def jac(self,t,y):
        j=np.zeros((self.dim+1,self.dim+1));j[:-1,:-1]=self.rates(y,True)
        # Omission counter does not feed populations; a zero last row affects
        # Newton efficiency only and is separately checked to be negligible.
        return j

def diagnostics(e,y):
    n=len(e);p=y[:n];a=y[n:2*n];np_=float(p.sum());nm=float(a.sum());q=QREF-np_+nm
    eg=float(e@(p+a));kin=eg-np_-nm
    total_tangent=EREF+OMEGA*(q-QREF)+eg
    cold_tangent=EREF+OMEGA*(q-QREF)+np_+nm
    reverse=e>=1+2*OMEGA
    return {'Nplus':np_,'Nminus':nm,'Qcore':q,'gas_energy_mX':eg,'gas_kinetic_mX':kin,
        'gas_invariant_E_minus_Omega_Q':eg-OMEGA*(np_-nm),
        'total_tangent_energy_mX':total_tangent,'core_charge_fraction':q/QREF,
        'positive_number_above_reverse_threshold':float(p[reverse].sum()),
        'tangent_cold_mass_excess_fraction':cold_tangent/EREF-1,
        'min_bin_number':float(min(p.min(),a.min())),
        'omitted_event_integral':float(y[-1])}

def equilibrium_checks(model):
    e=model.e;n=len(e);width=np.empty(n);width[1:-1]=(e[2:]-e[:-2])/2
    width[0]=(e[1]-e[0])/2;width[-1]=(e[-1]-e[-2])/2
    out=[]
    for T in (.2,.5):
        x=e*np.sqrt(e*e-1)*np.exp(-(e-1)/T)*width;x*=1000/x.sum()
        if model.res is not None:
            r=model.res;rr=r['K']*x[r['i']]*x[r['j']];der=r['D']@rr
            out.append({'channel':'resonant','T':T,'L1_residual_per_event':float(np.sum(abs(der))/(4*rr.sum()))})
        if model.el is not None:
            r=model.el;i,j=r['i'],r['j'];rr=np.where(i==j,.5,1.)*x[i]*x[j];der=r['D']@rr
            out.append({'channel':'identical_elastic','T':T,
                        'L1_residual_per_event':float(np.sum(abs(der))/(4*(r['K']@rr)))})
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=96);ap.add_argument('--emax',type=float,default=8)
    ap.add_argument('--mode',choices=['core','resonant','full'],default='full');ap.add_argument('--angles',type=int,default=8)
    ap.add_argument('--tmax',type=float,default=1000);ap.add_argument('--rtol',type=float,default=2e-7)
    ap.add_argument('--grid',choices=['quadratic','aligned','translated'],default='quadratic')
    ap.add_argument('--shift-bins',type=int,default=64)
    ap.add_argument('--tag',default=None);args=ap.parse_args()
    tag=args.tag or f'{args.mode}_n{args.n}_e{args.emax:g}_a{args.angles}'
    start=time.monotonic();e=grid(args.n,args.emax,args.grid,args.shift_bins);y0,outside=initial(e);model=Model(e,args.mode,args.angles)
    print(tag,'built',round(time.monotonic()-start,2),'s',flush=True)
    init=np.r_[y0,0.]
    result=solve_ivp(model.rhs,(0,args.tmax),init,method='BDF',jac=model.jac,
                     rtol=args.rtol,atol=1e-9,dense_output=True)
    if not result.success:raise RuntimeError(result.message)
    times=np.unique(np.r_[0,np.geomspace(.001,args.tmax,181),[.1,1,10,100,args.tmax]])
    times=times[times<=args.tmax];ys=result.sol(times);rows=[]
    for t,y in zip(times,ys.T):rows.append({'tau':float(t),**diagnostics(e,y)})
    energy=np.array([r['gas_invariant_E_minus_Omega_Q'] for r in rows]);minpop=min(r['min_bin_number'] for r in rows)
    checks={'core_column_invariant_error':model.core_error,
        'invariant_relative_drift':float(np.max(abs(energy-energy[0]))/abs(energy[0])),
        'min_bin_population':minpop,'initial_unsupported_spectrum':outside,
        'equilibrium':equilibrium_checks(model)}
    for name,op in [('resonant',model.res),('elastic',model.el)]:
        if op is not None:checks[name]={k:v for k,v in op.items() if k.endswith('_error') or k=='build_seconds'}
    np.savez_compressed(ROOT/(tag+'.npz'),E=e,tau=times,Y=ys,ode_times=result.t)
    out={'configuration':vars(args),'tag':tag,'solver':{'success':result.success,'nfev':result.nfev,'njev':result.njev,
        'nlu':result.nlu,'accepted_steps':len(result.t),'seconds':time.monotonic()-start},'checks':checks,
        'snapshots':[rows[int(np.argmin(abs(times-t)))] for t in [0,.1,1,10,100,args.tmax]]}
    write_json(tag+'.json',out);print(json.dumps(out,indent=2),flush=True)

if __name__=='__main__':main()
