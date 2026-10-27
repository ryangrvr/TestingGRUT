"""en_particles SKEPTIC (own code): bias test of the en_particles response transform (cos^2 taper, T=5, ds=0.05,
rho=-2 Im G~, peak on a 0.01 grid) on (a) free fields with known m^2 and (b) an oscillator whose mass is modulated
like the measured isotropic Hartree mass, m^2(s) = 20.1 + 8.3 cos(2*3.49 s + phase), averaged over 8 kick phases.
Usage: python3 -I sk_synth_response.py <out.json>"""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp
ds = .05; D = np.arange(0, 5.0001, ds); T = D[-1]; taper = np.cos(np.pi*D/(2*T))**2
trap = np.full(len(D), ds); trap[[0, -1]] = ds/2
wg = np.linspace(-20, 20, 4001); dw = wg[1]-wg[0]; E = np.exp(-1j*np.outer(wg, D))
def transform(G):
    rho = -2*(E@(G*taper*trap)).imag; ip = np.argmax(np.where(wg > 0, rho, -np.inf))
    return float(wg[ip]), float((wg*rho)[wg > 0].sum()*dw/(2*np.pi)), float((wg*rho)[wg < 0].sum()*dw/(2*np.pi))
res = dict(free=[], modulated=[])
for k, m2 in [(0.52, 20.1), (1.05, 20.1), (2.09, 20.1), (4.19, 20.1), (6.28, 20.1), (4.19, 12.), (1.05, 12.)]:
    w = np.sqrt(k*k+m2); pk, Zp, Zm = transform(np.sin(w*D)/w)
    res['free'].append(dict(k=k, m2_true=m2, omega_true=w, rho_peak=pk, m2_recovered=pk*pk-k*k, Z_plus=Zp, Z_minus=Zm))
for k in [1.05, 2.09, 4.19, 5.24, 6.28]:
    pk = []
    for ph in np.linspace(0, np.pi, 8, endpoint=False):
        f = lambda t, y: [y[1], -(k*k+20.1+8.3*np.cos(2*3.49*t+ph))*y[0]]
        G = solve_ivp(f, (0, T), [0, 1], t_eval=D, rtol=1e-10, atol=1e-12).y[0]; pk.append(transform(G)[0])
    res['modulated'].append(dict(k=k, rho_peak_mean=float(np.mean(pk)), rho_peak_std=float(np.std(pk)),
                                 m2_from_peak=float(np.mean(pk)**2-k*k), omega_mean_mass=float(np.sqrt(k*k+20.1))))
json.dump(res, open(sys.argv[1], 'w'), indent=1); print(json.dumps(res, indent=1))
