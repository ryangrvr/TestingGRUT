"""en_particles: lab-frame frequency decomposition of the homogeneous orbit at the matching time (own code) [C/E].

Takes the spatial mean (chi1,chi2,pi1,pi2) of the regenerated tau=20 state and integrates the homogeneous orbit of
H = |pi|^2/2 + (lambda_X/4)|chi|^4 (phi mean ~0 and neglected; fluctuations neglected) for many periods; then
n_pm = sum_{w>0 / w<0} |w| P(w) per net charge, the same estimator as spectral_window.py.
Usage: python3 -I zero_mode_orbit.py <state_s0000000.bin> <out.json>
"""
import sys, json, struct
import numpy as np
from scipy.integrate import solve_ivp

LX = 1e-3


def main(path, out):
    with open(path, 'rb') as f:
        N, L, s, bits = struct.unpack('<iddi', f.read(24))
    a = np.fromfile(path, dtype='<f8', offset=24).reshape(6, N, N, N)
    m = a.reshape(6, -1).mean(axis=1)
    nQ_total = float(np.mean(a[0]*a[4]-a[1]*a[3]))
    y0 = np.array([m[0], m[1], m[3], m[4]]); Q0 = m[0]*m[4]-m[1]*m[3]
    def rhs(t, y):
        r2 = y[0]**2+y[1]**2
        return [y[2], y[3], -LX*r2*y[0], -LX*r2*y[1]]
    T = 200.; ds = .01
    t = np.arange(0, T, ds)
    sol = solve_ivp(rhs, (0, T), y0, t_eval=t, rtol=1e-12, atol=1e-10, method='DOP853')
    psi = sol.y[0]+1j*sol.y[1]
    M = len(t); h = np.hanning(M); h2 = (h**2).sum()
    P = np.abs(np.fft.fft(h*psi))**2/(M*h2); w = 2*np.pi*np.fft.fftfreq(M, ds)
    npl = float((np.abs(w)*P)[w > 0].sum()); nmi = float((np.abs(w)*P)[w < 0].sum())
    r = np.abs(psi)
    order = np.argsort(-P)[:6]
    res = dict(source='own code', mean_state=dict(chi=[m[0], m[1]], pi=[m[3], m[4]], phi=m[2], pi_phi=m[5]),
               Q_zero=Q0, Q_total=nQ_total, n_plus_per_Q0=npl/Q0, n_minus_per_Q0=nmi/Q0,
               n_minus_per_Qtotal=nmi/nQ_total, charge_check=(npl-nmi)/Q0,
               r_min=float(r.min()), r_max=float(r.max()), top_lines=[[float(w[i]), float(P[i]/P.sum())] for i in order],
               energy=float(.5*(m[3]**2+m[4]**2)+LX/4*(m[0]**2+m[1]**2)**2))
    json.dump(res, open(out, 'w'), indent=1); print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main(*sys.argv[1:])
