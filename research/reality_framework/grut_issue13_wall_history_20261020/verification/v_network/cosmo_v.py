#!/usr/bin/env python3
"""Closed-form checks of the module's unit conversions (verifier's own arithmetic, numpy only).
- A definition: rho_w = A sigma/t; comoving S_c/V_c; a ~ tau; t = int a dtau = a tau / 2.
- t_w = 3 M_P^2/(4 A sigma) (reduced Planck mass; rho_rad = 3 M_P^2 H^2, H = 1/(2t)).
- H^-1/delta at T = 1.5476 v, g* = 106.75; module-unit relation H^-1/delta = tau^2/2.
"""
import json
import numpy as np

MP = 2.435323e18          # reduced Planck mass, GeV (BRIEF)
hbar = 6.582119569e-25    # GeV s
v, lam = 1000.0, 0.01
sigma = 2 * np.sqrt(2 * lam) / 3 * v ** 3
m_phi = np.sqrt(2 * lam) * v
delta = 2 / m_phi          # phi = v tanh(z/delta), delta = sqrt(2/lam)/v
out = {"sigma_GeV3": sigma, "sigma_cuberoot": sigma ** (1 / 3), "m_phi": m_phi, "delta_GeV^-1": delta}
# numerical check of a = c tau, t = int_0^tau a dtau'
c = 0.37
tau = 5.3
t = c * tau ** 2 / 2
out["t_equals_a_tau_over_2_check"] = t / (c * tau * tau / 2)
# in module units a = tau (tau in 1/m_phi): H = a'/a^2 = 1/tau^2 ; H^-1/delta = tau^2/(2) since delta = 2/m_phi
out["module_H^-1/delta_at_tau=32"] = (32.0 ** 2) / 2
out["module_H^-1/delta_at_tau=12"] = (12.0 ** 2) / 2
for A in (0.8, 0.5, 0.652, 0.777):
    tw = 3 * MP ** 2 / (4 * A * sigma)  # GeV^-1
    out[f"t_w_s_A={A}"] = tw * hbar
T = 1.547646465 * v
gs = 106.75
H = np.sqrt(np.pi ** 2 * gs / 90) * T ** 2 / MP
out["H^-1/delta_at_T=1.5476v_g106.75"] = (1 / H) / delta
json.dump(out, open("results/cosmo_v.json", "w"), indent=1)
print(json.dumps(out, indent=1))
