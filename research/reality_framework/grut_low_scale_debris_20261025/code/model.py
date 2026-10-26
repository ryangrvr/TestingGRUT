"""Fixed classical point; full mediator potential, with no one-field truncation."""
import numpy as np
from scipy.optimize import brentq

MRED = 2.435e18
YQ = 4.8294009076988175e-9
GSTAR = 106.75
LX, LY, LS, KAPPA, E = 1e-3, 1., 1., 1e-3, 1e-7
M_GEV, HIGGS_GEV, MU_GEV = 1e12, 1e11, 5e11
CHARGES = np.array([1., 3., 6.])


def parameters(H):
    return dict(lx=LX, ly=LY, ls=LS, k=KAPPA, M=M_GEV/H,
                h=HIGGS_GEV/H, mu=MU_GEV/H, e=E)


def potential(z, p, cH=0.):
    x, y, s = z
    return float(LX*abs(x)**4+LY*abs(y)**4+LS*(abs(s)**2-p['h']**2/2)**2
                 + p['M']**2*abs(y)**2
                 + 2*np.real(KAPPA*np.conj(y)*x**3+p['mu']*np.conj(s)*y*y)
                 - cH*abs(x)**2)


def force(z, p, cH=0.):
    x, y, s = z
    return np.array([2*LX*abs(x)**2*x+3*KAPPA*y*np.conj(x)**2-cH*x,
                     p['M']**2*y+2*LY*abs(y)**2*y+KAPPA*x**3+2*p['mu']*s*np.conj(y),
                     2*LS*(abs(s)**2-p['h']**2/2)*s+p['mu']*y*y])


def minimum(H, cH):
    """Aligned radial stationary point reduced to one scalar root.

    At fixed y, the positive x and s stationary roots are unique minima.
    The returned y root is also checked by its full radial Hessian.
    Multiple-start minimization is a separate numerical control.
    """
    p = parameters(H)
    def x_of_y(y):
        return (3*KAPPA*y+np.sqrt(9*KAPPA**2*y*y+8*LX*cH))/(4*LX)
    def s_of_y(y):
        hv=p['h']/np.sqrt(2)
        if y == 0: return hv
        f=lambda s: 2*LS*s*(s*s-hv*hv)-p['mu']*y*y
        upper=hv+(p['mu']*y*y/(2*LS))**(1/3)+1e-12
        return brentq(f, hv, upper, xtol=1e-14, rtol=1e-14)
    def condition(y):
        x, s=x_of_y(y), s_of_y(y)
        return 2*LY*y**3+p['M']**2*y-KAPPA*x**3-2*p['mu']*s*y
    upper=max(1., np.sqrt(cH))
    while condition(upper)<0: upper*=2
    y=brentq(condition, 0., upper, xtol=1e-13, rtol=1e-14)
    x, s=x_of_y(y), s_of_y(y)
    z=np.array([x,-y,-s], complex)
    radial=np.array([[12*LX*x*x-12*KAPPA*y*x-2*cH, -6*KAPPA*x*x,0.],
                    [-6*KAPPA*x*x,12*LY*y*y+2*p['M']**2-4*p['mu']*s,-4*p['mu']*y],
                    [0.,-4*p['mu']*y,12*LS*s*s-2*LS*p['h']**2]])/2
    v1=np.array([3.,-1.,0.]);v2=np.array([0.,2.,-1.])
    angular=2*KAPPA*y*x**3*np.outer(v1,v1)+2*p['mu']*s*y*y*np.outer(v2,v2)
    norm=np.diag(1/(np.sqrt(2)*np.array([x,y,s])))
    mass2=np.linalg.eigvalsh(norm@angular@norm)
    assert min(np.linalg.eigvalsh(radial))>0
    return dict(H_GeV=H,cH=cH,amplitudes_over_H=[x,y,s],
                gradient_norm=float(np.linalg.norm(2*force(z,p,cH))),
                V_inflation_over_H4=potential(z,p,cH),
                radial_mass2_over_H2=np.linalg.eigvalsh(radial).tolist(),
                angular_mass2_over_H2=mass2.tolist(),
                quadratic_Y_force_ratio=2*LY*y*y/p['M']**2,
                kappa_X2_over_M2=KAPPA*x*x/p['M']**2,
                vector_mass_over_H=float(np.sqrt(2)*E*np.sqrt(x*x+9*y*y+36*s*s)))


def initial_z(point, delta):
    x,y,s=point['amplitudes_over_H']
    return np.array([x,-y,-s*np.exp(1j*delta)],complex)


def entropy(H):
    T=(90*MRED*MRED*H*H/(np.pi**2*GSTAR))**.25
    return 2*np.pi**2/45*GSTAR*T**3, T


def yield_factor(H):
    return H**3/entropy(H)[0]
