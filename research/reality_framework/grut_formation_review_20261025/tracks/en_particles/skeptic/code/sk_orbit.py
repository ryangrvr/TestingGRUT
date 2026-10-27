"""en_particles SKEPTIC (own code): exact lab-frame +/- decomposition of the homogeneous tau=20 orbit.

Method (no time window): for the central quartic potential V = (lambda_X/4)|psi|^4 the orbit psi = r e^{i theta}
is quasi-periodic. Over one radial period T_r, theta advances by D_theta. With Omega = D_theta/T_r,
f(s) = psi(s) e^{-i Omega s} is exactly T_r-periodic, so psi = sum_n c_n exp(i(Omega + n w_r) s), w_r = 2 pi/T_r.
Then n_+ = sum_{w_n>0} w_n |c_n|^2, n_- = sum_{w_n<0} |w_n| |c_n|^2 (per unit volume), exactly.
The radial period and D_theta come from quadrature of the radial equation between turning points.
Also reports the circular-orbit energy at the same charge and the eccentricity measures.
Usage: python3 -I sk_orbit.py <state_s0000000.bin> <out.json>
"""
import sys, json, struct
import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.optimize import brentq

LX = 1e-3


def main(path, out):
    with open(path, 'rb') as f:
        N, L, s, bits = struct.unpack('<iddi', f.read(24))
    a = np.fromfile(path, dtype='<f8', offset=24).reshape(6, -1)
    m = a.mean(axis=1)
    Qtot = float(np.mean(a[0]*a[4]-a[1]*a[3]))
    x, y, px, py = m[0], m[1], m[3], m[4]
    l = x*py-y*px                      # angular momentum = zero-mode charge density
    E = .5*(px*px+py*py)+LX/4*(x*x+y*y)**2
    Veff = lambda r: l*l/(2*r*r)+LX/4*r**4
    # turning points: E = Veff(r)
    rc = (l*l/LX)**(1/6)               # minimum of Veff
    rmin = brentq(lambda r: Veff(r)-E, 1e-8*rc, rc); rmax = brentq(lambda r: Veff(r)-E, rc, 1e3*rc)
    # T_r/2 = int dr / sqrt(2(E-Veff)); D_theta/2 = int l/r^2 dr / sqrt(...): substitution r = (rmin+rmax)/2 + (rmax-rmin)/2 sin u
    c, h = (rmin+rmax)/2, (rmax-rmin)/2
    def integrand(u, kind):
        r = c+h*np.sin(u); dr = h*np.cos(u); den = np.sqrt(max(2*(E-Veff(r)), 1e-300))
        return dr/den if kind == 0 else l/r**2*dr/den
    Tr = 2*quad(integrand, -np.pi/2, np.pi/2, args=(0,), limit=400, epsabs=0, epsrel=1e-13)[0]
    Dth = 2*quad(integrand, -np.pi/2, np.pi/2, args=(1,), limit=400, epsabs=0, epsrel=1e-13)[0]
    Om = Dth/Tr; wr = 2*np.pi/Tr
    # integrate one radial period with a dense grid and take the exact Fourier series of f(s)
    M = 1 << 14; t = np.linspace(0, Tr, M, endpoint=False)
    rhs = lambda tt, Y: [Y[2], Y[3], -LX*(Y[0]**2+Y[1]**2)*Y[0], -LX*(Y[0]**2+Y[1]**2)*Y[1]]
    sol = solve_ivp(rhs, (0, Tr), [x, y, px, py], t_eval=t, rtol=1e-13, atol=1e-11, method='DOP853')
    psi = sol.y[0]+1j*sol.y[1]
    f = psi*np.exp(-1j*Om*t); cn = np.fft.fft(f)/M; n = np.fft.fftfreq(M)*M
    wn = Om+n*wr; P = np.abs(cn)**2
    npl = float((wn*P)[wn > 0].sum()); nmi = float((-wn*P)[wn < 0].sum())
    # end-of-period closure check: psi(Tr) should equal psi(0) e^{i Dth}
    end = solve_ivp(rhs, (0, Tr), [x, y, px, py], rtol=1e-13, atol=1e-11, method='DOP853').y[:, -1]
    closure = abs((end[0]+1j*end[1])-(x+1j*y)*np.exp(1j*Dth))/abs(x+1j*y)
    top = np.argsort(-P)[:5]
    Ecirc = .75*LX**(1/3)*abs(l)**(4/3)
    res = dict(source='en_particles_skeptic own code (exact quasi-periodic Fourier series, no window)',
               zero_mode_charge=l, total_charge=Qtot, energy_zero_mode=E, r_min=rmin, r_max=rmax, r_circ=rc,
               radial_period=Tr, delta_theta=Dth, Omega=Om, omega_r=wr, closure_error=float(closure),
               n_plus_per_Q0=npl/l, n_minus_per_Q0=nmi/l, n_minus_per_Qtotal=nmi/Qtot, charge_check=(npl-nmi)/l,
               lines=[[float(wn[i]), float(P[i]/P.sum())] for i in top],
               E_circular_over_E0=Ecirc/E)
    json.dump(res, open(out, 'w'), indent=1); print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main(*sys.argv[1:])
