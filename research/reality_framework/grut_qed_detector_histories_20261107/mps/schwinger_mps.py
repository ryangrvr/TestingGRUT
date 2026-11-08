"""Lattice Schwinger model as an MPO (bond dim 6), Gauss law eliminated exactly.
W = x sum(S+S- + S-S+) + sum_n [mu/2 + mu(-1)^n Sz_n] + sum_{n=0}^{N-1} Q_n^2,
Q_n = sum_{k<=n} q_k, q_k = Sz_k + (-1)^k/2.  (Q_{N-1}=0 in the Sz=0 sector, so
including it matches the ED convention of links 0..N-2.)  Conventions identical to
../qed/schwinger.py: sz=2*Sz, hopping +x, H/g = W/(2 sqrt x).
MPO channels: 0=Id-left, 1=Q, 2=Q^2, 3=S+ pending, 4=S- pending, 5=done.
"""
import numpy as np
from tenpy.networks.site import SpinHalfSite
from tenpy.networks.mpo import MPO
from tenpy.networks.mps import MPS
from tenpy.models.model import MPOModel
from tenpy.models.lattice import Chain
from tenpy.algorithms import dmrg

def build(N, x, m_over_g):
    site = SpinHalfSite(conserve="Sz")
    mu = 2 * np.sqrt(x) * m_over_g
    Id, Sz, Sp, Sm = (site.get_op(o) for o in ["Id", "Sz", "Sp", "Sm"])
    grids = []
    for n in range(N):
        s = (-1) ** n
        q = Sz + 0.5 * s * Id
        q2 = 0.5 * Id + s * Sz
        onsite = q2 + (mu / 2) * Id + (mu * s) * Sz
        G = [[None] * 6 for _ in range(6)]
        G[0][0] = Id; G[0][1] = q; G[0][2] = q2; G[0][3] = x * Sp; G[0][4] = x * Sm; G[0][5] = onsite
        G[1][1] = Id; G[1][2] = 2 * q; G[1][5] = 2 * q
        G[2][2] = Id; G[2][5] = Id
        G[3][5] = Sm; G[4][5] = Sp
        G[5][5] = Id
        if n == 0:
            G = [G[0]]
        if n == N - 1:
            G = [[row[5]] for row in G]
        grids.append(G)
    sites = [site] * N
    H = MPO.from_grids(sites, grids, bc="finite", IdL=0, IdR=-1)
    lat = Chain(N, site, bc="open", bc_MPS="finite")
    return MPOModel(lat, H), sites

def gap(N, x, m_over_g, chi=(40, 80, 120), tol=1e-10, max_sweeps=40, return_info=False):
    M, sites = build(N, x, m_over_g)
    params = {"trunc_params": {"chi_max": chi[-1], "svd_min": 1e-12},
              "chi_list": {i * 4: c for i, c in enumerate(chi)},
              "max_E_err": tol, "max_S_err": 1e-6, "max_sweeps": max_sweeps,
              "mixer": True, "combine": True}
    neel = ["up", "down"] * (N // 2)
    psi0 = MPS.from_product_state(sites, neel, bc="finite")
    E0 = dmrg.run(psi0, M, params)["E"]
    psi1 = MPS.from_product_state(sites, ["down", "up"] * (N // 2), bc="finite")
    E1 = dmrg.run(psi1, M, dict(params), orthogonal_to=[psi0])["E"]
    om = (E1 - E0) / (2 * np.sqrt(x))
    if return_info:
        return om, {"E0": E0, "E1": E1, "chi0": max(psi0.chi), "chi1": max(psi1.chi),
                    "maxS0": float(max(psi0.entanglement_entropy()))}
    return om
