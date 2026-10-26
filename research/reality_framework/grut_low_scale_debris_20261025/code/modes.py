"""Real-matrix particle projection and independent vacuum-width dictionary.

Physical canonical complex fields in H_i units; comoving numbers in H_i^3.
The two-real-mode convention follows the f1aa7df review. Evolution is not
changed or copied from the reviewer.
"""
import numpy as np
from model import parameters


def project(H,a,state):
    p=parameters(H);v=np.asarray(state,float)
    u=v[:6:2]+1j*v[1:6:2];pi=v[6::2]+1j*v[7::2]
    z=u/a;zd=(pi-u/a)/a**2;X,Y,S=z;Xd,Yd,Sd=zd
    A=p['M']**2;B=2*p['mu']*S
    mass=np.array([[A+B.real,B.imag],[B.imag,A-B.real]])
    J=1e-3*X**3;vec=np.array([J.real,J.imag])
    yr=np.linalg.solve(mass,-vec);yp=yr[0]+1j*yr[1]
    Bd=2*p['mu']*Sd;md=np.array([[Bd.real,Bd.imag],[Bd.imag,-Bd.real]])
    Jd=3e-3*X**2*Xd
    ypr=np.linalg.solve(mass,-np.array([Jd.real,Jd.imag])-md@yr)
    ypd=ypr[0]+1j*ypr[1]
    d=np.array([(Y-yp).real,(Y-yp).imag]);dd=np.array([(Yd-ypd).real,(Yd-ypd).imag])
    phase_dot=float(np.imag(np.conj(S)*Sd)/abs(S)**2)
    rotation=np.array([[0.,-1.],[1.,0.]])
    co_velocity=dd-.5*phase_dot*rotation@d
    eig,basis=np.linalg.eigh(mass)
    assert eig[0]>0
    coord=basis.T@d;velocity=basis.T@co_velocity;omega=np.sqrt(eig)
    number=a**3*(velocity**2+eig*coord**2)/omega
    r=abs(S);rd=np.real(np.conj(S)*Sd)/r
    Erad=a**4*(rd**2+(r*r-p['h']**2/2)**2)
    ns=Erad/(a*np.sqrt(2)*p['h'])
    EX=a**4*(abs(Xd)**2+1e-3*abs(X)**4)
    N=-2*a**3*np.imag(np.conj(z)*zd)
    nmax=2*np.sqrt(2e-3)*(EX/.003)**.75
    JY=a**3*(np.dot(dd,dd)+d@mass@d)/p['M']
    # Independent finite difference of the forced-response derivative.
    step=1e-6
    def forced(ZZ):
        xx,_,ss=ZZ;bb=2*p['mu']*ss
        mat=np.array([[A+bb.real,bb.imag],[bb.imag,A-bb.real]])
        jj=1e-3*xx**3;vv=np.linalg.solve(mat,-np.array([jj.real,jj.imag]))
        return vv[0]+1j*vv[1]
    fd=(forced(z+step*zd)-forced(z-step*zd))/(2*step)
    return dict(NX=float(N[0]),EX=float(EX),eps=float(N[0]/nmax),
                NYL=float(number[0]),NYH=float(number[1]),NSIG=float(ns),JY=float(JY),
                ESIG=float(Erad),ESANG=float(a**4*(abs(Sd)**2-rd**2)),
                Sigma_phase_dot=float(phase_dot),Y_split_ratio=float(abs(B)/A),
                Y_self_force_ratio=float(2*abs(Y)**2/A),
                forced_derivative_FD_relative_error=float(abs(fd-ypd)/max(abs(ypd),1e-30)),
                Gauss=float(N@np.array([1.,3.,6.])))


def vacuum():
    M,h,mu,e=1e12,1e11,5e11,1e-7;mx=.09240896373753621
    mh=np.sqrt(M*M+np.sqrt(2)*mu*h);ml=np.sqrt(M*M-np.sqrt(2)*mu*h)
    ms=np.sqrt(2)*h;ma=6*e*h
    lam=(mh*mh-(ml+ma)**2)*(mh*mh-(ml-ma)**2)
    ps=np.sqrt(lam)/(2*mh)
    gHL=9*e*e*lam**1.5/(16*np.pi*mh**3*ma**2)
    r=ma*ma/(ms*ms)
    gs=ms**3/(32*np.pi*h*h)*np.sqrt(1-4*r)*(1-4*r+12*r*r)
    ga=e*e*ma/(48*np.pi)*(1-4*mx*mx/(ma*ma))**1.5
    return dict(mYH=mh,mYL=ml,mSigma=ms,mA=ma,mX=mx,
                GammaYH3=3e-6*mh/(256*np.pi**3),GammaYL3=3e-6*ml/(256*np.pi**3),
                GammaYHLA=gHL,GammaSigmaAA=gs,GammaAXX=ga,
                pA_from_YH=ps,pA_from_Sigma=np.sqrt(ms*ms/4-ma*ma))
