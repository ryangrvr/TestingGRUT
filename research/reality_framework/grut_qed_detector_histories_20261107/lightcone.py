"""Shared light cone: Wilson fermion (r=1, massless limit) vs plaquette photon on the same cubic lattice."""
import json, numpy as np
rng = np.random.default_rng(7)
def v_photon(k):
    s, c = np.sin(k/2), np.cos(k/2); return np.linalg.norm(s*c/np.sqrt(np.sum(s**2)))
def v_wilson(k, r=1.0):
    E = np.sqrt(np.sum(np.sin(k)**2) + (r*np.sum(1-np.cos(k)))**2)
    grad = (np.sin(k)*np.cos(k) + r*np.sum(1-np.cos(k))*r*np.sin(k))/E
    return np.linalg.norm(grad)
ka = 2e-2
dirs = rng.normal(size=(20000, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
Ce = np.array([(v_wilson(ka*n)-1)/ka**2 for n in dirs])
Cg = np.array([(v_photon(ka*n)-1)/ka**2 for n in dirs])
S = np.sum(dirs**4, 1)
res = {"check_formula_Ce=3/8-S/2_maxerr": float(np.max(np.abs(Ce-(3/8-S/2)))),
       "check_formula_Cg=-S/8_maxerr": float(np.max(np.abs(Cg+S/8))),
       "fraction_dirs_electron_faster_than_photon": float(np.mean(Ce > Cg)),
       "fraction_dirs_electron_superluminal": float(np.mean(Ce > 0)),
       "Ce_max_body_diag": 5/24, "Ce_axis": -1/8}
# Crab: electrons >= 2.3 PeV (LHAASO, SSC model) do not vacuum-Cherenkov radiate.
# soft threshold: C_e (E a)^2 = m^2/(2E^2)  ->  a_max = m / (sqrt(2 C_e) E^2)
hbarc, me, Ee = 1.973269804e-16, 0.51099895e-3, 2.3e6   # GeV m, GeV, GeV
def amax(C): return me/(np.sqrt(2*C)*Ee**2)*hbarc
res["a_max_m_body_diagonal"] = amax(5/24)
res["a_max_m_median_superluminal_dir"] = amax(float(np.median(Ce[Ce > 0])))
res["photon_only_LHAASO_GRB_ceiling_m"] = 1.6e-27
json.dump(res, open("lightcone_results.json", "w"), indent=1); print(json.dumps(res, indent=1))
