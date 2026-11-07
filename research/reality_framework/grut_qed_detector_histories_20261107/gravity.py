"""Gravity-sector probe in the style of the GRUT photon bridge.

G1  Lattice Newton law: the static exchange between two masses on the cubic
    lattice is the lattice Laplacian Green function. Does it approach 1/(4 pi r)?
G2  Lattice spin-2 counting: gauge redundancy of the Fierz-Pauli operator with
    lattice momenta, and the number/sign of propagating polarizations seen by
    conserved sources (massless vs massive spin-2).
G3  Calibration vs holdout: tune each candidate force law to reproduce Newton
    for slow masses, then predict light bending. Compare with Cassini.
G4  Observational ceiling on the regulator spacing a from the photon lattice's
    own dispersion and LHAASO GRB 221009A time-of-flight limits.
"""
import json
import numpy as np
from scipy.integrate import quad
from scipy.special import ive

out = {}

# ---------------- G1 ----------------
def G_lat(x, y, z):
    # G(r) = int_0^inf prod_i e^{-2t} I_{r_i}(2t) dt  for -Delta_lat G = delta
    f = lambda t: ive(abs(x), 2*t) * ive(abs(y), 2*t) * ive(abs(z), 2*t)
    v1, _ = quad(f, 0, 50, limit=400)
    v2, _ = quad(f, 50, np.inf, limit=400)
    return v1 + v2

g1 = {"G(0)": G_lat(0, 0, 0), "watson_reference_G(0)": 0.252731}
rows = []
for n in [1, 2, 4, 8, 16]:
    for name, v in [("axis", (n, 0, 0)), ("diagonal", (n, n, n))]:
        r = np.linalg.norm(v)
        rows.append({"dir": name, "r_over_a": round(float(r), 3),
                     "4pi_r_G": round(4*np.pi*r*G_lat(*v), 6)})
g1["ratio_to_newton"] = rows
out["G1_lattice_newton"] = g1

# ---------------- G2 ----------------
eta = np.diag([-1.0, 1, 1, 1])
pairs = [(m, n) for m in range(4) for n in range(m, 4)]          # 10 components

def sym_basis():
    B = []
    for (m, n) in pairs:
        h = np.zeros((4, 4)); h[m, n] = h[n, m] = 1.0
        B.append(h)
    return B
SB = sym_basis()

def lower(h):
    return eta @ h @ eta

def FP_operator(k):
    """Linearized Einstein operator (Lichnerowicz/Fierz-Pauli) acting on h_{mu nu}
    (lower indices), with k_mu the (lattice) wave covector."""
    ku = eta @ k          # k^mu
    k2 = k @ ku
    def K(h):              # h lower, returns lower
        hu = eta @ h @ eta
        trh = np.trace(eta @ h)
        kh = h @ ku        # k^a h_{a nu} -> vector (lower nu)
        kkh = ku @ h @ ku
        return (k2*h - np.outer(k, kh) - np.outer(kh, k)
                + np.outer(k, k)*trh - eta*(k2*trh - kkh))
    M = np.array([[np.sum(SB[i]*K(SB[j])) for j in range(10)] for i in range(10)])
    return M, K

def lat_mom(kvec, a=1.0):
    return (2/a)*np.sin(np.asarray(kvec)*a/2)

g2 = []
for label, kvec in [("axis", [0.9, 0, 0]), ("off_axis", [0.7, 0.5, 0.3]),
                    ("near_zone_edge", [2.6, 1.9, 0.4])]:
    kh = lat_mom(kvec)
    for shell, w in [("off_shell", 0.37*np.linalg.norm(kh)), ("on_shell", np.linalg.norm(kh))]:
        k = np.array([-w, *kh])            # k_mu = (-omega, khat)
        M, K = FP_operator(k)
        # gauge: h = k_mu xi_nu + k_nu xi_mu  must be annihilated
        gauge_res = max(np.abs(K(np.outer(k, xi) + np.outer(xi, k))).max()
                        for xi in np.eye(4))
        sv = np.linalg.svd(M, compute_uv=False)
        g2.append({"k": label, "shell": shell,
                   "gauge_residual": float(gauge_res),
                   "nullity": int(np.sum(sv < 1e-9*sv.max()))})

def source_residue_rank(kvec, lam, mass=0.0):
    kh = lat_mom(kvec)
    w = np.sqrt(kh @ kh + mass**2)
    ku = np.array([w, *kh])                       # k^mu
    # conserved symmetric sources T^{mu nu}: k_mu T^{mu nu} = 0
    kl = eta @ ku
    C = np.array([[ (SB[j] @ kl)[nu] for j in range(10)] for nu in range(4)])
    _, s, Vt = np.linalg.svd(C)
    null = Vt[np.sum(s > 1e-10):]                 # basis of conserved sources (6)
    def bil(Ta, Tb):                              # Ta^{mn} P_{mnrs} Tb^{rs}
        A, B = lower(Ta), lower(Tb)
        return 0.5*(np.sum(Ta*B) + np.sum(Ta*B.T)) - lam*np.trace(eta@Ta)*np.trace(eta@Tb)
    Ts = [sum(c*SB[j] for j, c in enumerate(v)) for v in null]
    R = np.array([[bil(a, b) for b in Ts] for a in Ts])
    ev = np.linalg.eigvalsh(R)
    tol = 1e-9*np.abs(ev).max()
    return {"conserved_source_dim": len(Ts), "positive": int(np.sum(ev > tol)),
            "negative": int(np.sum(ev < -tol))}

out["G2_spin2_counting"] = {
    "gauge_and_nullity": g2,
    "massless_spin2_lambda_1/2": source_residue_rank([0.7, 0.5, 0.3], 0.5),
    "massive_spin2_lambda_1/3": source_residue_rank([0.7, 0.5, 0.3], 1/3, mass=0.4),
}

# ---------------- G3 ----------------
def amp(T1, T2, lam):
    A, B = lower(T1), lower(T2)
    return 0.5*(np.sum(T1*B) + np.sum(T1*B.T)) - lam*np.trace(eta@T1)*np.trace(eta@T2)
Tsun = np.diag([1.0, 0, 0, 0])                    # static mass, T^{00}=M
Tslow = np.diag([1.0, 0, 0, 0])                   # slow test mass
n = np.array([1.0, 1, 0, 0]); Tlight = np.outer(n, n)  # photon along x
cands = {}
ref = amp(Tsun, Tlight, 0.5) / amp(Tsun, Tslow, 0.5)
for name, lam in [("massless_spin2 (GR)", 0.5), ("massive_spin2_m->0 (Fierz-Pauli)", 1/3)]:
    bend = (amp(Tsun, Tlight, lam)/amp(Tsun, Tslow, lam)) / ref
    gamma = 2*bend - 1
    cands[name] = {"light_bending_vs_GR": round(float(bend), 6), "implied_gamma": round(float(gamma), 6)}
# scalar coupled to trace T
tr = lambda T: np.trace(eta@T)
bend_s = (tr(Tsun)*tr(Tlight)/(tr(Tsun)*tr(Tslow))) / ref
cands["scalar_coupled_to_trace"] = {"light_bending_vs_GR": float(bend_s), "implied_gamma": float(2*bend_s-1)}
cassini = {"gamma_minus_1": 2.1e-5, "sigma": 2.3e-5}
for c in cands.values():
    c["sigmas_from_cassini"] = round(abs(c["implied_gamma"] - 1 - cassini["gamma_minus_1"])/cassini["sigma"], 1)
out["G3_newton_calibrated_light_bending"] = {"cassini": cassini, "candidates": cands}

# ---------------- G4 ----------------
# plaquette-only photon: omega = (2c/a) sqrt(sum sin^2(k_i a/2))
def group_velocity(kvec, a=1.0):
    k = np.asarray(kvec, float); s = np.sin(k*a/2); c = np.cos(k*a/2)
    w = (2/a)*np.sqrt(np.sum(s**2))
    grad = (s*c)/np.sqrt(np.sum(s**2))   # d omega / dk_i  (units c=1)
    return np.linalg.norm(grad)
coef = {}
for name, u in [("axis", [1, 0, 0]), ("face_diag", [1, 1, 0]), ("body_diag", [1, 1, 1])]:
    u = np.array(u)/np.linalg.norm(u); ka = 1e-2
    coef[name] = (1 - group_velocity(ka*u))/ka**2   # dv/c = -C (k a)^2
hbar_c = 1.973269804e-16       # GeV m
E_pl = 1.22089e19              # GeV
E_qg2 = 6e-8*E_pl              # LHAASO GRB 221009A, quadratic, 95% CL
bounds = {}
for name, C in coef.items():
    # LIV convention: dv/c = -(3/2)(E/E_QG2)^2 ; lattice: dv/c = -C (E a / hbar c)^2
    a_max = np.sqrt(1.5/C)*hbar_c/E_qg2
    bounds[name] = {"C_dv_over_c_per_(ka)^2": round(float(C), 6), "a_max_m": float(a_max)}
out["G4_spacing_ceiling"] = {"E_QG2_GeV": E_qg2, "by_direction": bounds,
                             "planck_length_m": 1.616255e-35}

json.dump(out, open("gravity_results.json", "w"), indent=2)
print(json.dumps(out, indent=2))
