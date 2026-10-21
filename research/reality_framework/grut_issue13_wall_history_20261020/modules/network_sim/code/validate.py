"""Validation of the wall-area estimator, periodic labelling, IC embedding and the integrator energy budget.
Writes results/validation.json."""
import json, os, sys, math, time
import numpy as np
from scipy import ndimage
sys.path.insert(0, os.path.dirname(__file__))
import wallsim as ws

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")
out = {}

# 1. PRS area estimator on circles / spheres / tilted planes (tanh profile of width DELTA, dx varied)
tests = []
for d in (2, 3):
    for dx in (1.0, 0.5):
        N = 128 if d == 3 else 512
        L = N * dx
        x = (np.arange(N) - N / 2 + 0.5) * dx
        grids = np.meshgrid(*([x] * d), indexing="ij")
        r = np.sqrt(sum(g * g for g in grids))
        for R in (8.0, 16.0, 24.0):
            if R > 0.4 * L:
                continue
            phi = np.tanh((r - R) / ws.DELTA).astype(np.float32)
            pos, links = ws.wall_links(phi)
            a_prs, a_iso, _ = ws.prs_area(phi, dx, links)
            exact = 2 * math.pi * R if d == 2 else 4 * math.pi * R * R
            tests.append({"shape": "circle" if d == 2 else "sphere", "d": d, "dx": dx, "R": R,
                          "prs_over_exact": a_prs / exact, "iso_over_exact": a_iso / exact})
        # tilted periodic planes: phi = tanh(c sin(2 pi k.x / L)); zero set area = 2 |k| L^(d-1)
        xs = np.arange(N) * dx
        grids = np.meshgrid(*([xs] * d), indexing="ij")
        for k in ([1, 0, 0], [1, 1, 0], [1, 2, 0], [1, 2, 3]):
            k = k[:d]
            if all(v == 0 for v in k[1:]) and d == 3 and k != [1, 0, 0]:
                continue
            arg = sum(kk * g for kk, g in zip(k, grids)) * 2 * math.pi / L
            phi = np.tanh(3 * np.sin(arg)).astype(np.float32)
            pos, links = ws.wall_links(phi)
            a_prs, a_iso, _ = ws.prs_area(phi, dx, links)
            exact = 2 * math.sqrt(sum(kk * kk for kk in k)) * L ** (d - 1)
            tests.append({"shape": "planes", "d": d, "dx": dx, "k": k,
                          "prs_over_exact": a_prs / exact, "iso_over_exact": a_iso / exact})
out["area_estimator"] = tests
prs_err = [abs(t["prs_over_exact"] - 1) for t in tests]
out["area_estimator_max_abs_err_prs"] = max(prs_err)
out["area_estimator_max_abs_err_iso"] = max(abs(t["iso_over_exact"] - 1) for t in tests)

# 2. periodic labelling
face2 = ndimage.generate_binary_structure(2, 1)
m = np.zeros((16, 16), bool)
m[3, :] = True          # horizontal stripe wrapping in x: 1 component, spans
m[8:10, 0:2] = True     # blob crossing boundary with
m[8:10, 14:16] = True   # its periodic partner -> 1 component
lab, n = ws.plabel(m, face2)
lab_info = {"stripe_plus_wrapped_blob_n": n, "expected": 2,
            "stripe_spans": bool(ws.spans(lab, lab[3, 0])), "blob_spans": bool(ws.spans(lab, lab[8, 0]))}
m3 = np.zeros((12, 12, 12), bool)
m3[0, 0, 0] = True; m3[11, 11, 11] = True  # corner cells: face-disconnected, full-connected through corner
lab_info["corner_face_n"] = ws.plabel(m3, ndimage.generate_binary_structure(3, 1))[1]
lab_info["corner_full_n"] = ws.plabel(m3, ndimage.generate_binary_structure(3, 3))[1]
lab_info["pass"] = (n == 2 and lab_info["stripe_spans"] and not lab_info["blob_spans"]
                    and lab_info["corner_face_n"] == 2 and lab_info["corner_full_n"] == 1)
out["labelling"] = lab_info

# 3. IC embedding: same continuum field on N and 2N grids
f1 = ws.gaussian_ic(2, 256, 64.0, 0.75, 0.1, 7, N_base=256, dtype=np.float64)
f2 = ws.gaussian_ic(2, 512, 64.0, 0.75, 0.1, 7, N_base=256, dtype=np.float64)
out["ic_embedding_max_abs_diff_on_common_points"] = float(np.abs(f2[::2, ::2] - f1).max())
out["ic_std"] = [float(f1.std()), float(f2.std())]

# 4. integrator energy budget: E(tau) - E(tau_i) - int (dE/dtau)_friction/expansion
ech = {}
for name, cfg in {
    "phys2d_N256": dict(d=2, N=256, L=16.0, alpha=2, beta=2, tau_a=1.0, tau_i=2.0, tau_end=8.0, courant=1 / 3.5,
                        s_c=0.25, amp=0.1, seed=11, n_out=4, dtype="f64", clusters=False, energy_check=True),
    "prs2d_N256": dict(d=2, N=256, L=256.0, alpha=2, beta=0, tau_i=2.0, tau_end=128.0, courant=1 / 3,
                       s_c=1.5, amp=0.1, seed=11, n_out=4, dtype="f64", clusters=False, energy_check=True),
    "prs3d_N48": dict(d=3, N=48, L=48.0, alpha=3, beta=0, tau_i=2.0, tau_end=24.0, courant=1 / 3,
                      s_c=1.5, amp=0.1, seed=11, n_out=4, dtype="f64", clusters=False, energy_check=True),
}.items():
    cfg["name"] = "echeck_" + name
    for dtfac in (1.0, 0.5):
        c = dict(cfg); c["courant"] = cfg["courant"] * dtfac
        r = ws.evolve(c, outdir=None, verbose=False)
        ech[f"{name}_courant{c['courant']:.4f}"] = {"max_abs_resid_over_E_tau": r["energy_max_abs_rel_resid"],
                                                    "max_abs_resid_over_maxE": r["energy_max_abs_resid_over_maxE"]}
out["energy_budget_max_abs_rel_resid"] = ech
with open(os.path.join(RES, "validation.json"), "w") as f:
    json.dump(out, f, indent=1)
print(json.dumps({k: v for k, v in out.items() if k != "area_estimator"}, indent=1))
for t in tests:
    print(t)
