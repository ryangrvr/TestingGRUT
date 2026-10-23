#!/usr/bin/env python3
"""Three-channel linear scattering at the fixed C1-HP-03/v1 point.

The amplitude vector is (eta_plus, eta_minus_conjugate, sqrt(2)*s_plus).
Its real symmetric radial potential gives a positive number-flux Wronskian.
The negative-frequency charged component has negative physical U(1) charge.
"""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy.integrate import solve_bvp, simpson
from scipy.interpolate import CubicSpline
from scipy.special import spherical_jn, spherical_yn, spherical_kn

ROOT = Path(__file__).resolve().parents[1]
PARENT = ROOT.parent / 'grut_portal_debris_window_20261023'
G, V, A, B = .1, .924089637375362, .1, 3.
MX = G*V


class Background:
    def __init__(self, refined=True):
        tag = '_refined' if refined else ''
        path = PARENT/'data'/f'profile_Q200000_lambda0.03{tag}.csv'
        arr = np.loadtxt(path, delimiter=',', skiprows=1)
        self.path = path
        self.hash = hashlib.sha256(path.read_bytes()).hexdigest()
        self.f = CubicSpline(arr[:, 0], arr[:, 1])
        self.s = CubicSpline(arr[:, 0], arr[:, 3])
        data = json.loads((PARENT/'FIXED_YIELD_RESULTS.json').read_text())
        row = data['refinement'] if refined else data['new_profile']
        self.omega = row['Omega']

    def fields(self, r):
        inside = r <= self.f.x[-1]
        f = np.where(inside, self.f(np.minimum(r, self.f.x[-1])), 0.)
        s = np.where(inside, self.s(np.minimum(r, self.s.x[-1])), 1.)
        return f, s

    def potential(self, r, scale=1.):
        f, s = self.fields(r)
        n = len(np.atleast_1d(r))
        out = np.zeros((3, 3, n))
        out[0, 0] = out[1, 1] = s*s-1+2*A*f*f
        out[2, 2] = B*(3*s*s-1)+f*f-6
        out[0, 1] = out[1, 0] = A*f*f
        out[0, 2] = out[2, 0] = np.sqrt(2)*s*f
        out[1, 2] = out[2, 1] = np.sqrt(2)*s*f
        return scale*out


def boundary_wave(ell, k2, radius):
    if k2 > 0:
        k = np.sqrt(k2); z = k*radius
        h = spherical_jn(ell, z)+1j*spherical_yn(ell, z)
        hp = spherical_jn(ell, z, derivative=True)+1j*spherical_yn(ell, z, derivative=True)
        H = z*h
        return dict(open=True, k=k, H=H, L=1/radius+k*hp/h)
    kappa = np.sqrt(-k2); z = kappa*radius
    h = spherical_kn(ell, z)
    hp = spherical_kn(ell, z, derivative=True)
    return dict(open=False, k=kappa, H=0j, L=complex(1/radius+kappa*hp/h))


def scatter(bg, energy, ell=0, radius=40., tol=2e-7, mesh=801,
            scale=1., incoming=(1,), max_nodes=30000, return_solutions=False,
            basis='radial', rmin=1e-4, absorption=0.):
    """energy is the incident antiparticle energy divided by m_X.

    Boundary f_i(R)=1 for a unit incoming partial wave; exact Riccati-Hankel
    derivatives set the outgoing condition. Regular amplitudes are factored
    as radial_u=r**(ell+1)*f, removing the centrifugal singularity.
    """
    nu = energy+bg.omega
    freq = np.array([nu+bg.omega, nu-bg.omega, nu])
    k2 = freq*freq-np.array([1., 1., 6.])
    waves = [boundary_wave(ell, z, radius) for z in k2]
    active = [i for i, w in enumerate(waves) if w['open']]
    L = np.array([w['L']-((ell+1)/radius if basis=='factored' else 0) for w in waves])
    r = np.linspace(0. if basis=='factored' else rmin, radius, mesh)
    singular = np.zeros((6, 6)); singular[3:, 3:] = -2*(ell+1)*np.eye(3)
    origin_matrix=bg.potential(np.array([0.]),scale)[:,:,0]-np.diag(k2)
    origin_L=(ell+1)/rmin*np.eye(3)+rmin*origin_matrix/(2*ell+3)

    def rhs(x, y):
        mat = bg.potential(x, scale).astype(complex)
        mat[np.arange(3), np.arange(3)] -= k2[:, None]
        if basis=='radial':mat[np.arange(3),np.arange(3)]+=ell*(ell+1)/x[None,:]**2
        mat[2,2]-=1j*absorption*np.exp(-(x/12)**8)
        return np.vstack([y[3:], np.einsum('ijn,jn->in', mat, y[:3])])

    columns=[]; diagnostics=[]; solutions=[]
    for inc in incoming:
        if inc not in active:
            raise ValueError('Incoming channel is closed')
        source = np.zeros(3, complex)
        source[inc] = waves[inc]['L'].conjugate()-waves[inc]['L']
        def bc(left, right):
            left_res=left[3:] if basis=='factored' else left[3:]-origin_L@left[:3]
            return np.r_[left_res, right[3:]-L*right[:3]-source]
        # Analytic free regular spherical wave is a useful linear guess.
        guess = np.zeros((6, len(r)), complex)
        k = np.sqrt(k2[inc]); z = k*r
        small = r==0
        rr = np.where(small, 1e-6, r)
        H = waves[inc]['H'].conjugate()
        if basis=='factored':
            fac = 2*k*radius**(ell+1)/H
            guess[inc] = fac*spherical_jn(ell, k*rr)/rr**ell
            guess[inc+3] = fac*(k*spherical_jn(ell, k*rr, derivative=True)/rr**ell
                                   -ell*spherical_jn(ell, k*rr)/rr**(ell+1))
            guess[inc+3, 0] = 0.
        else:
            guess[inc]=2*k*rr*spherical_jn(ell,k*rr)/H
            guess[inc+3]=2*k*(spherical_jn(ell,k*rr)+k*rr*spherical_jn(ell,k*rr,derivative=True))/H
        # Real/imaginary blocks avoid complex-Newton refinement pathology
        # in scipy's singular-endpoint solver for the high-energy column.
        def rhs_real(x, y):
            val=rhs(x, y[:6]+1j*y[6:])
            return np.vstack([val.real, val.imag])
        def bc_real(left, right):
            val=bc(left[:6]+1j*left[6:], right[:6]+1j*right[6:])
            return np.r_[val.real, val.imag]
        sr=np.zeros((12,12));sr[:6,:6]=singular;sr[6:,6:]=singular
        sol = solve_bvp(rhs_real, bc_real, r, np.vstack([guess.real,guess.imag]),
                        S=sr if basis=='factored' else None, tol=tol, bc_tol=tol*.05, max_nodes=max_nodes)
        if not sol.success:
            raise RuntimeError(sol.message)
        right = sol.y[:3, -1]+1j*sol.y[6:9,-1]
        col = np.zeros(3, complex)
        for j in active:
            col[j] = (right[j]-(j==inc))*np.sqrt(waves[j]['k']/waves[inc]['k'])*H/waves[j]['H']
        columns.append(col)
        solutions.append(sol)
        p = abs(col)**2
        diagnostics.append(dict(incoming=inc, probabilities=p.tolist(),
            number_flux_residual=float(sum(p)-1),
            max_BVP_residual=float(max(sol.rms_residuals)), nodes=len(sol.x),
            core_field_max=float(np.max(abs(sol.y[:3])))))
    smat = np.column_stack(columns)
    result=dict(energy_over_mX=float(energy), incident_K_MeV=float((energy-1)*MX*1000),basis=basis,rmin=rmin,absorption_fixture=absorption,
        nu=float(nu), ell=ell, radius=radius, tolerance=tol, initial_mesh=mesh,
        active_channels=active, frequencies=freq.tolist(), k_squared=k2.tolist(),
        S_real=smat.real.tolist(), S_imag=smat.imag.tolist(), columns=diagnostics)
    if return_solutions:return result,solutions,waves
    return result


def visible_spectral_probability(bg, result, sol, waves, inc=1, nr=4001, nq=401,
                                 photon_coefficient=0.):
    """First-order inclusive SM loss for a CLOSED scalar channel.

    Retain the spatial form factor and the invariant timelike SM spectral
    kernel. No frequency-local absorption or trapped-state lifetime is used.
    Lepton cuts are exact tree level; photons are an explicit optional
    low-energy coefficient sensitivity, not an exact hadronic computation.
    """
    if 2 in result['active_channels']:
        raise ValueError('Open scalar decay must use a finite-time outer receiver')
    radius=result['radius'];ell=result['ell'];nu=result['nu']
    r=np.linspace(0,radius,nr);y=sol.sol(r);z=y[2]+1j*y[8]
    # Divide the physical radial amplitude by the incident Hankel coefficient.
    uz=waves[inc]['H'].conjugate()*z
    if result['basis']=='factored':uz*=(r/radius)**(ell+1)
    else:uz[0]=0.
    _,s=bg.fields(r)
    q=np.linspace(0,nu,nq)
    form=simpson(s[None,:]*uz[None,:]*r[None,:]*spherical_jn(ell,q[:,None]*r[None,:]),x=r,axis=1)
    invariant=np.maximum((nu*nu-q*q)*MX*MX,0.)
    mh,vh,kappa=125.,246.,9.006401966916295e-10
    # Im Pi_h(s) = sqrt(s) Gamma_unit(s). The s/h mixing insertion is
    # kappa*v*vh; factors S(r) on each end are in the form factor above.
    im_h=np.zeros_like(q)
    for ml in [.00051099895,.1056583755]:
        beta=np.sqrt(np.maximum(1-4*ml*ml/np.maximum(invariant,1e-200),0.))
        im_h+=ml*ml*invariant/(8*np.pi*vh*vh)*beta**3
    alpha=1/137.035999084
    im_h+=photon_coefficient**2*alpha**2*invariant**2/(256*np.pi**3*vh*vh)
    kernel=(kappa*V*vh/(mh*mh-invariant))**2*im_h/(MX*MX)
    probability=2/np.pi*simpson(q*q*abs(form)**2*kernel,x=q)/waves[inc]['k']
    return dict(probability=float(probability),radial_points=nr,momentum_points=nq,
                photon_coefficient=photon_coefficient,
                W_max_GeV=float(nu*MX),invariant_hadronic_threshold_open=bool(nu*MX>.268))


def main():
    p=argparse.ArgumentParser();p.add_argument('--energy',type=float,default=1.001)
    p.add_argument('--ell',type=int,default=0);p.add_argument('--radius',type=float,default=40.)
    args=p.parse_args();bg=Background()
    start=time.monotonic()
    r=scatter(bg,args.energy,args.ell,args.radius,incoming=(0,1))
    print(json.dumps(dict(omega=bg.omega,result=r,seconds=time.monotonic()-start),indent=2))


if __name__=='__main__':main()
