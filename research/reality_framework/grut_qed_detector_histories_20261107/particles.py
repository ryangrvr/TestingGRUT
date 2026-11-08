"""Particle-sector probe in the style of the GRUT photon bridge.

Question: can the photon-bridge recipe (local, Hermitian, translation-invariant
lattice Hamiltonian on a spatial torus, continuous time) also host the
Standard Model's matter particles?

P1  Naive lattice Weyl fermion: count low-energy species and their chirality.
P2  Try to beat it: random local perturbations; measure net chirality with a
    Berry-flux (Chern number) sweep that does not assume where nodes are.
P3  Wilson fermion: one light Dirac particle, doublers at the cutoff, and the
    price (explicit chiral-symmetry breaking).
"""
import itertools, json
import numpy as np

rng = np.random.default_rng(20261106)
sx = np.array([[0, 1], [1, 0]], complex)
sy = np.array([[0, -1j], [1j, 0]], complex)
sz = np.array([[1, 0], [0, -1]], complex)
S = [sx, sy, sz]
out = {}

# ---------------- P1: naive Weyl nodes ----------------
nodes = []
for n in itertools.product([0, 1], repeat=3):
    k = np.pi * np.array(n)
    J = np.diag(np.cos(k))            # d_i = sin k_i  ->  dd_i/dk_j
    nodes.append({"k/pi": list(n), "chirality": int(np.sign(np.linalg.det(J)))})
out["P1_naive_weyl"] = {
    "species": len(nodes),
    "right_handed": sum(c["chirality"] > 0 for c in nodes),
    "left_handed": sum(c["chirality"] < 0 for c in nodes),
    "net_chirality": sum(c["chirality"] for c in nodes),
}

# ---------------- P2: random local deformations ----------------
# d_i(k) = sin k_i + eps * sum over hoppings |n|_inf <= 2 of (a cos(n.k) + b sin(n.k))
hops = [np.array(n) for n in itertools.product(range(-2, 3), repeat=3) if any(n)]


def make_d(eps, seed):
    r = np.random.default_rng(seed)
    A = r.normal(size=(3, len(hops))) * eps
    B = r.normal(size=(3, len(hops))) * eps
    H0 = r.normal(size=len(hops)) * eps     # identity part (tilts nodes, harmless)

    def d(K):  # K shape (...,3)
        ph = np.tensordot(K, np.array(hops).T, axes=1)  # (...,nh)
        base = np.sin(K)
        pert = np.stack([np.cos(ph) @ A[i] + np.sin(ph) @ B[i] for i in range(3)], -1)
        return base + pert
    return d


def lower_band(dv):
    # 2x2 H = d.sigma ; lower-band eigenvector, vectorised
    H = sum(dv[..., i, None, None] * S[i] for i in range(3))
    w, v = np.linalg.eigh(H)
    return v[..., :, 0]


def chern_slices(d, N=48, Nz=97):
    ks = 2 * np.pi * (np.arange(N) + 0.5) / N
    kzs = 2 * np.pi * (np.arange(Nz) + 0.37) / Nz
    KX, KY = np.meshgrid(ks, ks, indexing="ij")
    C = []
    for kz in kzs:
        K = np.stack([KX, KY, np.full_like(KX, kz)], -1)
        u = lower_band(d(K))
        def link(a, b):
            z = np.sum(np.conj(a) * b, -1)
            return z / np.abs(z)
        U1 = link(u, np.roll(u, -1, 0))
        U2 = link(u, np.roll(u, -1, 1))
        F = np.angle(U1 * np.roll(U2, -1, 0) * np.conj(np.roll(U1, -1, 1)) * np.conj(U2))
        C.append(int(np.rint(F.sum() / (2 * np.pi))))
    return np.array(C)


from scipy.optimize import root

def find_nodes(d, starts=4000, seed=0):
    r = np.random.default_rng(seed)
    found = []
    def jac(k, h=1e-6):
        return np.array([(d(k + h*e) - d(k - h*e)) / (2*h) for e in np.eye(3)]).T
    for k0 in r.uniform(0, 2*np.pi, size=(starts, 3)):
        sol = root(lambda k: d(k), k0, jac=jac, tol=1e-12)
        if not sol.success or np.linalg.norm(d(sol.x)) > 1e-9:
            continue
        k = np.mod(sol.x, 2*np.pi)
        if any(np.linalg.norm(np.angle(np.exp(1j*(k - q)))) < 1e-5 for q, _ in found):
            continue
        found.append((k, int(np.sign(np.linalg.det(jac(k))))))
    return found

trials = []
for seed, eps in [(2, 0.03), (3, 0.08), (4, 0.15), (5, 0.25)]:
    d = make_d(eps / np.sqrt(len(hops)), seed)
    nd = find_nodes(d, seed=seed)
    C = chern_slices(d)
    jumps = np.diff(C)
    trials.append({
        "eps_rms_scale": eps,
        "nodes_found": len(nd),
        "right_handed": sum(c > 0 for _, c in nd),
        "left_handed": sum(c < 0 for _, c in nd),
        "net_chirality_from_nodes": int(sum(c for _, c in nd)),
        "berry_flux_slice_chern_range": [int(C.min()), int(C.max())],
        "sum_abs_slice_jumps": int(np.abs(jumps).sum()),
    })
out["P2_random_local_deformations"] = trials

# ---------------- P3: Wilson fermion ----------------
m, r = 0.05, 1.0   # lattice units, a = 1
species = []
for n in itertools.product([0, 1], repeat=3):
    k = np.pi * np.array(n)
    M = m + r * np.sum(1 - np.cos(k))
    species.append(M)
species = sorted(species)
# small-k dispersion: E^2 = sum sin^2 k + (m + r sum(1-cos k))^2 vs k^2 + m^2
def E2(k):
    k = np.asarray(k, float)
    return np.sum(np.sin(k) ** 2) + (m + r * np.sum(1 - np.cos(k))) ** 2
disp = {}
for name, dirn in [("axis", [1, 0, 0]), ("diagonal", [1, 1, 1])]:
    u = np.array(dirn) / np.linalg.norm(dirn)
    rel = []
    for kk in [0.02, 0.05, 0.1]:
        rel.append(float(E2(kk * u) / (kk ** 2 + m ** 2) - 1))
    disp[name] = rel
# chiral breaking: gamma5 anticommutator with Wilson term at k=0 and at a doubler
# Explicit chiral basis: alpha_i = diag(s_i,-s_i), beta = offdiag(I), gamma5 = diag(I,-I)
Z = np.zeros((2, 2))
alpha = [np.block([[s, Z], [Z, -s]]) for s in S]
beta = np.block([[Z, np.eye(2)], [np.eye(2), Z]])
g5 = np.block([[np.eye(2), Z], [Z, -np.eye(2)]])
def HW(k):
    return sum(np.sin(k[i]) * alpha[i] for i in range(3)) + (m + r * np.sum(1 - np.cos(k))) * beta
def HW0(k):  # massless naive part only
    return sum(np.sin(k[i]) * alpha[i] for i in range(3))
kk = np.array([0.3, 0.2, 0.1])
out["P3_wilson"] = {
    "m_lattice": m, "r": r,
    "species_masses_lattice_units": [round(x, 4) for x in species],
    "light_species": int(sum(x < 0.5 for x in species)),
    "dispersion_rel_error_E2_vs_k2+m2": disp,
    "naive_part_commutes_with_g5_norm": float(np.linalg.norm(HW0(kk) @ g5 - g5 @ HW0(kk))),
    "wilson_term_breaks_chirality_norm_at_doubler": float(np.linalg.norm(
        (HW(np.array([np.pi, 0, 0])) - HW0(np.array([np.pi, 0, 0]))) @ g5
        - g5 @ (HW(np.array([np.pi, 0, 0])) - HW0(np.array([np.pi, 0, 0]))))),
}

json.dump(out, open("particles_results.json", "w"), indent=2)
print(json.dumps(out, indent=2))
