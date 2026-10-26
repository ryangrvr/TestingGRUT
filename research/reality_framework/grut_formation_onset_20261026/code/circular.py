"""Exact finite-response dispersion about the supplied circular stationary family.

tau=m_X t, k in m_X; canonical (u,v,binding perturbation).
The circular-state family is a diagnostic, not the matched elliptic release.
"""
import numpy as np
from scipy.optimize import minimize_scalar,brentq
from common import *

def state(a):
    r=RC*a;z=NU*np.sqrt(max(1-a*a,0))
    w2=G*G*z*z+LX*r*r;w=np.sqrt(w2)
    A=2*LX*r*r
    # Use the exact stationary identity: subtracting the equal terms at
    # a=1 creates a spurious sqrt(machine-epsilon) tachyon.
    D=2*LP*z*z if a<=1 else LP*NU*NU*(a*a-1)
    B=2*G*G*r*z
    return dict(a=a,r=r,phi=z,omega=w,A=A,D=D,B=B,nQ=w*r*r)

def matrix(s,k):
    A,D,B,w=s['A'],s['D'],s['B'],s['omega']
    M=np.zeros((6,6));M[:3,3:]=np.eye(3)
    M[3,0]=-k*k-A;M[3,2]=-B;M[3,4]=2*w
    M[4,1]=-k*k;M[4,3]=-2*w
    M[5,0]=-B;M[5,2]=-k*k-D
    return M

def growth(s,k):
    return max(0.,float(np.max(np.linalg.eigvals(matrix(s,k)).real)))

def polynomial(s,k):
    z=k*k;A,D,B,w=s['A'],s['D'],s['B'],s['omega']
    p=np.polymul([1,z+D],[1,2*z+A+4*w*w,z*(z+A)])
    p[-2]-=B*B;p[-1]-=B*B*z
    return p

def cutoff(s):
    A,D,B=s['A'],s['D'],s['B']
    return np.sqrt(max(0.,(-(A+D)+np.sqrt((A-D)**2+4*B*B))/2))

def old_growth(s,k):
    A=s['A']-4*G**4/LP*(s['r']**2/2);w=s['omega'];z=k*k
    y=np.roots([1,2*z+A+4*w*w,z*(z+A)])
    return max(0.,float(max(np.sqrt(y.astype(complex)).real)))

def maximum(fun,hi):
    if hi<=1e-10:return 0.,None
    ks=np.linspace(0,hi,251);gs=np.array([fun(k) for k in ks]);i=int(gs.argmax())
    lo=ks[max(i-1,0)];up=ks[min(i+1,len(ks)-1)]
    opt=minimize_scalar(lambda k:-fun(k),bounds=(lo,up),method='bounded',options={'xatol':1e-12})
    candidates=[(gs[i],ks[i]),(-opt.fun,opt.x),(fun(0),0.)]
    g,k=max(candidates)
    return (float(g),float(k)) if g>1e-9 else (0.,None)

def row(a):
    s=state(a);kc=cutoff(s);g,k=maximum(lambda k:growth(s,k),kc)
    oldcut=np.sqrt(max(0.,4*(G**4/LP-LX)*(s['r']**2/2)))
    go,ko=maximum(lambda k:old_growth(s,k),oldcut) if a<=1 else (None,None)
    if a>1:oldcut=None
    q=lambda kk,shape:shape*s['nQ']*(2*np.pi/kk)**3 if kk is not None and kk>1e-10 else None
    return dict(**s,full_cutoff=kc,full_growth=g,full_peak_k=k,adiabatic_cutoff=oldcut,adiabatic_growth=go,
        adiabatic_peak_k=ko,adiabatic_old_edge_Q_shape1=q(oldcut,1.),
        adiabatic_peak_Q_by_shape={str(f):q(ko,f) for f in [.1,.3,1.]},
        full_peak_Q_by_shape={str(f):q(k,f) for f in [.1,.3,1.]},
        response_parameter_at_peak=(g*g+k*k)/s['D'] if s['D']>1e-12 and k is not None else None)

def main():
    aa=np.r_[np.linspace(.05,.97,65),.98,.99,.995,.999,.9999,.99999,1.,1.05,1.2]
    rows=[row(float(a)) for a in aa]
    am=brentq(lambda a:state(a)['nQ']-Q0,1e-5,.95,xtol=1e-14);s=state(am)
    E=.5*s['omega']**2*s['r']**2+potential([s['r'],0,s['phi']])
    matched=dict(**s,energy_per_charge_m=E/Q0,release_energy_per_charge_m=(EX+V0)/Q0,
        fraction_of_release_energy_in_this_circular_state=E/(EX+V0),
        statement='Matches charge only; the missing energy must be in other degrees or excitations. No identification with the release.')
    endpoint=state(1.)
    endpoint_energy=.5*endpoint['omega']**2*endpoint['r']**2+potential([endpoint['r'],0.,0.])
    save('CIRCULAR_RESULTS.json',dict(rows=rows,examples=[row(a) for a in [.5,.8,.99,1.]],charge_matched_circular=matched,
        endpoint_charge_over_release_charge=endpoint['nQ']/Q0,
        endpoint_energy_over_release_energy=endpoint_energy/(EX+V0),
        old_onset_edge_Q=row(1.)['adiabatic_old_edge_Q_shape1'],
        claim='At the exact endpoint phi=0 the full circular vacuum system has no exponentially growing finite-k mode; the adiabatic negative quartic there is a singular dynamical elimination.',
        scope='Vacuum, circular, no plasma, gauge/heavy modes or expansion. Shape-based charge volumes are not final objects.'))
    print('CIRCULAR',rows[-3],matched,flush=True)

if __name__=='__main__':main()
