"""Original vectorized Bose thermodynamics, dimensionless m_X=1 integrals.

This reimplements, rather than imports, the review's supplied equilibrium
proxy. Fixed Gauss quadrature is checked separately against adaptive quad.
"""
import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

M=.09240896373753621


class BoseGas:
    def __init__(self,nodes=96):
        x,w=roots_legendre(nodes)
        edges=[0.,.25,1.,4.,12.,40.,80.]
        self.y=np.concatenate([a+(b-a)*(x+1)/2 for a,b in zip(edges[:-1],edges[1:])])
        self.w=np.concatenate([(b-a)*w/2 for a,b in zip(edges[:-1],edges[1:])])
        self.nodes=nodes

    def thermo(self,t,u):
        if not t>0 or not 0<=u<=1:raise ValueError((t,u))
        scale=t if t>=1 else np.sqrt(t)
        p=scale*self.y;ee=np.sqrt(1+p*p)
        xm=(p*p/(ee+1)+(1-u))/t;xp=xm+2*u/t
        fm=np.where(xm<700,1/np.expm1(np.minimum(xm,700)),0.)
        fp=np.where(xp<700,1/np.expm1(np.minimum(xp,700)),0.)
        q=self.w*p*p*scale/(2*np.pi**2)
        np_,nm=float(q@fm),float(q@fp)
        # Algebraic odd-in-mu difference: exactly zero at mu=0, even when
        # the net charge is much smaller than either thermal population.
        df=np.exp(-np.minimum(xm,745))*(-np.expm1(-2*u/t))/((-np.expm1(-xm))*(-np.expm1(-xp)))
        nq=float(q@df);rho=float(q@(ee*(fm+fp)))
        pressure=float(q@(p*p/(3*ee)*(fm+fp)))
        s=(rho+pressure-u*nq)/t
        return dict(npl=np_,nmi=nm,nq=nq,rho=rho,P=pressure,s=s)

    def state(self,nQ,uc,sigma=None,energy=None,ec=None):
        n=nQ/M**3;en=None if energy is None else energy/M**4
        assert (sigma is None)!=(energy is None)
        lo,hi=np.log(1e-10),np.log(1e4)
        lb=brentq(lambda l:self.thermo(np.exp(l),uc)['nq']-n,lo,hi,xtol=2e-11)
        tb=np.exp(lb);thb=self.thermo(tb,uc)
        def residual(t,th):
            return th['s']-sigma*n if sigma is not None else th['rho']+ec*(n-th['nq'])-en
        if residual(tb,thb)>=0:
            lr=brentq(lambda l:residual(np.exp(l),self.thermo(np.exp(l),uc)),lo,lb,xtol=2e-11)
            t=np.exp(lr);u=uc;th=self.thermo(t,u);nc=max(n-th['nq'],0.)
        else:
            def chem(t):
                if self.thermo(t,uc)['nq']<=n:return uc
                return brentq(lambda u:self.thermo(t,u)['nq']-n,0.,uc,xtol=2e-13)
            def target(l):
                t=np.exp(l);u=chem(t);th=self.thermo(t,u)
                return th['s']-sigma*n if sigma is not None else th['rho']-en
            lr=brentq(target,lb,hi,xtol=2e-11)
            t=np.exp(lr);u=chem(t);th=self.thermo(t,u);nc=0.
        out={k:v*M**(4 if k in ['rho','P'] else 3) for k,v in th.items()}
        out.update(T=t*M,mu=u*M,nc=nc*M**3)
        return out


def adaptive(t,u):
    from scipy.integrate import quad
    def f(y,k,sgn):
        p=y*t;ee=np.sqrt(1+p*p);x=(p*p/(ee+1)+(1-sgn*u))/t
        if x>700:return 0.
        weight=1 if k==0 else ee if k==1 else p*p/(3*ee)
        return y*y*weight/np.expm1(x)
    scale=np.sqrt(1+2/t);edges=np.array([0,.1,1,5,20,80])*scale
    a=[]
    for sgn in [1,-1]:
        a.append([sum(quad(lambda x:f(x,k,sgn),x,y,epsabs=0,epsrel=1e-10,limit=200)[0]
                      for x,y in zip(edges[:-1],edges[1:]))*t**3/(2*np.pi**2) for k in [0,1,2]])
    npl,rpl,ppl=a[0];nmi,rmi,pmi=a[1];rho=rpl+rmi;pressure=ppl+pmi
    return dict(npl=npl,nmi=nmi,nq=npl-nmi,rho=rho,P=pressure,s=(rho+pressure-u*(npl-nmi))/t)
