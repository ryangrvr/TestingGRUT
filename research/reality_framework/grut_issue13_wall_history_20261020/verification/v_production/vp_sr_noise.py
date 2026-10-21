#!/usr/bin/env python3
"""Sampling noise of the 'sampled incl. slow-roll' eps coefficient vs the exact equilibrium average.
Same rejection sampler as vp_spectator.py (seed 20261021): mean over N samples of E_c/3 with chi' = phi0 - gradU/3."""
import json, math, os
import numpy as np
lphi = 0.01; gX = 0.1; lX = 0.001; hD = 0.01; lD = 0.001; g2 = gX**2; h2 = hD**2
K8 = 8 * math.pi**2 / 3
def U4(p, r, d): return lphi/4*p**4 + g2*p**2*r**2/2 + lX/4*r**4 + h2*p**2*d**2/2 + lD/4*d**4
def gradU(p, r, d): return (lphi*p**3 + g2*p*r**2 + h2*p*d**2, g2*p**2*r + lX*r**3, h2*p**2*d + lD*d**3)
rng = np.random.default_rng(7)
def sample(n):
    cp = K8*lphi/4; cd = K8*lD/4; sig2 = 3/(4*math.pi**2*lX); P=[];R=[];D=[]
    while sum(len(x) for x in P) < n:
        m = 4*n
        p = (rng.gamma(0.25, size=m)/cp)**0.25*rng.choice([-1,1], size=m)
        d = (rng.gamma(0.25, size=m)/cd)**0.25*rng.choice([-1,1], size=m)
        r = np.sqrt(np.abs(rng.normal(0, math.sqrt(sig2), size=m)))*rng.choice([-1,1], size=m)
        a = rng.random(m) < np.exp(-K8*(g2*p**2*r**2/2 + h2*p**2*d**2/2))
        P.append(p[a]); R.append(r[a]); D.append(d[a])
    return np.concatenate(P)[:n], np.concatenate(R)[:n], np.concatenate(D)[:n]
p, r, d = sample(2_000_000)
g = gradU(p, r, d)
E_sr = 0.5*((p-g[0]/3)**2 + (r-g[1]/3)**2 + (d-g[2]/3)**2) + U4(p, r, d)
E_0 = 0.5*(p**2 + r**2 + d**2) + U4(p, r, d)
out = dict(N=len(p), c_eps_SR_mean=float(E_sr.mean()/3), c_eps_noSR_mean=float(E_0.mean()/3),
           per_sample_std_of_Ec_over_3=float(E_sr.std()/3),
           SEM_for_1500=float(E_sr.std()/3/math.sqrt(1500)), SEM_for_2000=float(E_sr.std()/3/math.sqrt(2000)))
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "vp_sr_noise.json"), "w"), indent=1)
print(out)
