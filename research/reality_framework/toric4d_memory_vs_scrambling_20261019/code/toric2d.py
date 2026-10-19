"""2D toric code control: Z-error strings on the edges of an L x L torus, Metropolis dynamics,
minimum-weight-perfect-matching readout (pymatching).  Same observables as toric4d.py.
"""
import json, sys, time
import numpy as np
import numba as nb
import scipy.sparse as sp
import pymatching


def build_tables(L):
    nq, nv = 2 * L * L, L * L
    e2v = np.zeros((nq, 2), np.int64)
    rows, cols = [], []
    for x in range(L):
        for y in range(L):
            for d in range(2):
                e = (x * L + y) * 2 + d
                v0 = x * L + y
                v1 = (((x + 1) % L) * L + y) if d == 0 else (x * L + (y + 1) % L)
                e2v[e] = (v0, v1)
                rows += [v0, v1]; cols += [e, e]
    H = sp.csc_matrix((np.ones(len(rows), np.uint8), (rows, cols)), shape=(nv, nq))
    # logical X-bar support: horizontal edges in column x=0 (crossed by a vertical dual loop)
    MX = np.array([(0 * L + y) * 2 + 0 for y in range(L)], np.int64)
    matcher = pymatching.Matching(H)
    return e2v, H, MX, matcher, nq, nv


@nb.njit(cache=True)
def _seed(s):
    np.random.seed(s)


@nb.njit(cache=True)
def metropolis(spins, synd, e2v, beta, nsweeps):
    nq = spins.shape[0]
    for _ in range(nsweeps):
        for _k in range(nq):
            e = np.random.randint(nq)
            dE = 0
            for j in range(2):
                dE += -2 if synd[e2v[e, j]] else 2
            if dE <= 0 or np.random.random() < np.exp(-beta * dE):
                spins[e] ^= 1
                synd[e2v[e, 0]] ^= 1; synd[e2v[e, 1]] ^= 1


@nb.njit(cache=True)
def metropolis_twin(sA, yA, sB, yB, e2v, beta, nsweeps):
    nq = sA.shape[0]
    for _ in range(nsweeps):
        for _k in range(nq):
            e = np.random.randint(nq); u = np.random.random()
            dA = (-2 if yA[e2v[e, 0]] else 2) + (-2 if yA[e2v[e, 1]] else 2)
            dB = (-2 if yB[e2v[e, 0]] else 2) + (-2 if yB[e2v[e, 1]] else 2)
            if dA <= 0 or u < np.exp(-beta * dA):
                sA[e] ^= 1; yA[e2v[e, 0]] ^= 1; yA[e2v[e, 1]] ^= 1
            if dB <= 0 or u < np.exp(-beta * dB):
                sB[e] ^= 1; yB[e2v[e, 0]] ^= 1; yB[e2v[e, 1]] ^= 1


def decoded_logical(spins, synd, matcher, MX):
    corr = matcher.decode(synd.astype(np.uint8))
    tot = (spins.astype(np.uint8) ^ corr.astype(np.uint8))
    return 1 - 2 * (int(tot[MX].sum()) % 2)


def memory_trajectory(L, beta, t_grid, seed, e2v, matcher, MX, nq, nv):
    _seed(seed)
    spins = np.zeros(nq, np.uint8); synd = np.zeros(nv, np.uint8)
    out = []; t_prev = 0
    for t in t_grid:
        metropolis(spins, synd, e2v, beta, int(t - t_prev)); t_prev = t
        out.append(decoded_logical(spins, synd, matcher, MX))
    return out, 2 * int(synd.sum()) / nq


def local_autocorr(L, beta, seed, e2v, nq, nv, n_eq=400, n_obs=400, n_sites=64):
    _seed(seed)
    spins = np.zeros(nq, np.uint8); synd = np.zeros(nv, np.uint8)
    metropolis(spins, synd, e2v, beta, n_eq)
    rs = np.random.RandomState(seed)
    n_sites = min(n_sites, nq)
    qs = rs.choice(nq, n_sites, replace=False); vs = rs.choice(nv, min(n_sites, nv), replace=False)
    S = np.zeros((n_obs, len(qs))); Y = np.zeros((n_obs, len(vs)))
    for t in range(n_obs):
        metropolis(spins, synd, e2v, beta, 1)
        S[t] = spins[qs]; Y[t] = synd[vs]
    def acf(X, maxlag):
        X = X - X.mean(0); var = (X * X).mean()
        if var <= 0:
            return [1.0] + [0.0] * maxlag
        return [float((X[:len(X) - k] * X[k:]).mean() / var) for k in range(maxlag + 1)]
    return acf(S, 30), acf(Y, 30)


def damage_spreading(L, beta, seed, e2v, nq, nv, n_eq=400, t_grid=(1, 2, 4, 8, 16, 32, 64, 128)):
    _seed(seed)
    sA = np.zeros(nq, np.uint8); yA = np.zeros(nv, np.uint8)
    metropolis(sA, yA, e2v, beta, n_eq)
    sB = sA.copy(); yB = yA.copy()
    e0 = nq // 2
    sB[e0] ^= 1; yB[e2v[e0, 0]] ^= 1; yB[e2v[e0, 1]] ^= 1
    out = []; t_prev = 0
    for t in t_grid:
        metropolis_twin(sA, yA, sB, yB, e2v, beta, int(t - t_prev)); t_prev = t
        out.append(int((sA != sB).sum()))
    return out


def fit_tau(t, P):
    t = np.asarray(t, float); P = np.asarray(P, float)
    if P.min() > 0.9:
        return {"tau": None, "lower_bound": float(t.max()), "note": "no decay within t_max"}
    m = (P > 0.05) & (P < 0.95)
    if m.sum() < 2:
        m = P < 0.95
    if m.sum() < 1:
        return {"tau": None, "note": "insufficient points"}
    slope = np.sum(-np.log(P[m]) * t[m]) / np.sum(t[m] ** 2)
    return {"tau": float(1.0 / slope) if slope > 0 else None, "n_points": int(m.sum())}


def main(L, betas, ntraj, tmax, seed):
    e2v, H, MX, matcher, nq, nv = build_tables(L)
    t_grid = np.unique(np.round(np.logspace(0, np.log10(tmax), 18)).astype(int))
    out = {"L": L, "n_qubits": nq, "t_grid": t_grid.tolist(), "betas": {}}
    t0 = time.time()
    for beta in betas:
        bits, eng = [], []
        for k in range(ntraj):
            b, e = memory_trajectory(L, beta, t_grid, seed + 1000 * k + int(beta * 100), e2v, matcher, MX, nq, nv)
            bits.append(b); eng.append(e)
        bits = np.array(bits, float); P = bits.mean(0)
        acf_s, acf_y = local_autocorr(L, beta, seed + 7, e2v, nq, nv)
        dmg = damage_spreading(L, beta, seed + 11, e2v, nq, nv)
        out["betas"][f"{beta:.2f}"] = {
            "P_log": P.tolist(), "P_log_se": (bits.std(0) / np.sqrt(ntraj)).tolist(),
            "tau_log": fit_tau(t_grid, P), "energy_per_qubit": float(np.mean(eng)),
            "acf_edge": acf_s, "acf_vertex": acf_y,
            "damage": {"t": [1, 2, 4, 8, 16, 32, 64, 128], "hamming": dmg},
        }
        print(f"2D L={L} beta={beta:.2f} P_log={np.round(P,3).tolist()} tau={out['betas'][f'{beta:.2f}']['tau_log']} ({time.time()-t0:.0f}s)", flush=True)
    return out


if __name__ == "__main__":
    L = int(sys.argv[1]); betas = [float(b) for b in sys.argv[2].split(",")]
    ntraj = int(sys.argv[3]); tmax = int(sys.argv[4]); seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    res = main(L, betas, ntraj, tmax, seed)
    outp = sys.argv[6] if len(sys.argv) > 6 else f"toric2d_L{L}.json"
    json.dump(res, open(outp, "w"), indent=1)
