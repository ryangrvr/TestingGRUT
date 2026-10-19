"""4D toric code, Z-error sector: Metropolis thermal dynamics, cooling decoder, memory time,
local autocorrelation and damage spreading.  See ../README.md for the pre-registered design.

Faces (qubits) of the L^4 torus, 6 planes per site.  Edge stabilizers; E = 2 * (#violated edges).
"""
import json, sys, time
import numpy as np
import numba as nb

PLANES = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def build_tables(L):
    N = L ** 4
    nf, ne = 6 * N, 4 * N
    f2e = np.zeros((nf, 4), np.int64)

    def site(x):
        return (((x[0] % L) * L + x[1] % L) * L + x[2] % L) * L + x[3] % L

    for i, x in enumerate(np.ndindex(L, L, L, L)):
        x = np.array(x)
        for p, (a, b) in enumerate(PLANES):
            xa = x.copy(); xa[a] += 1
            xb = x.copy(); xb[b] += 1
            f = i * 6 + p
            f2e[f, 0] = site(x) * 4 + a
            f2e[f, 1] = site(xb) * 4 + a
            f2e[f, 2] = site(x) * 4 + b
            f2e[f, 3] = site(xa) * 4 + b
    # logical X-membrane for plane (01): all (01)-faces at x0 = x1 = 0
    M = np.array([((0 * L + 0) * L + x2) * L * 6 + x3 * 6 + 0 for x2 in range(L) for x3 in range(L)], np.int64)
    # check: each edge is in exactly 6 faces
    cnt = np.bincount(f2e.ravel(), minlength=ne)
    assert cnt.min() == 6 and cnt.max() == 6
    return f2e, M, nf, ne


@nb.njit(cache=True)
def _seed(s):
    np.random.seed(s)


@nb.njit(cache=True)
def metropolis(spins, synd, f2e, beta, nsweeps):
    nf = spins.shape[0]
    for _ in range(nsweeps):
        for _k in range(nf):
            f = np.random.randint(nf)
            dE = 0
            for j in range(4):
                dE += -2 if synd[f2e[f, j]] else 2
            if dE <= 0 or np.random.random() < np.exp(-beta * dE):
                spins[f] ^= 1
                for j in range(4):
                    synd[f2e[f, j]] ^= 1


@nb.njit(cache=True)
def metropolis_twin(sA, synA, sB, synB, f2e, beta, nsweeps):
    """Two copies driven by the identical random stream (damage spreading)."""
    nf = sA.shape[0]
    for _ in range(nsweeps):
        for _k in range(nf):
            f = np.random.randint(nf)
            u = np.random.random()
            dA = 0; dB = 0
            for j in range(4):
                e = f2e[f, j]
                dA += -2 if synA[e] else 2
                dB += -2 if synB[e] else 2
            if dA <= 0 or u < np.exp(-beta * dA):
                sA[f] ^= 1
                for j in range(4):
                    synA[f2e[f, j]] ^= 1
            if dB <= 0 or u < np.exp(-beta * dB):
                sB[f] ^= 1
                for j in range(4):
                    synB[f2e[f, j]] ^= 1


@nb.njit(cache=True)
def cool(spins, synd, f2e, nsweeps):
    """T=0 local decoder: accept dE<0; accept dE=0 with prob 1/2 (erodes contractible membranes)."""
    nf = spins.shape[0]
    for _ in range(nsweeps):
        for _k in range(nf):
            f = np.random.randint(nf)
            dE = 0
            for j in range(4):
                dE += -2 if synd[f2e[f, j]] else 2
            if dE < 0 or (dE == 0 and np.random.random() < 0.5):
                spins[f] ^= 1
                for j in range(4):
                    synd[f2e[f, j]] ^= 1


def energy(synd):
    return 2 * int(synd.sum())


def logical_bit(spins, M):
    return 1 - 2 * (int(spins[M].sum()) % 2)


def memory_trajectory(L, beta, t_grid, seed, f2e, M, nf, ne):
    """Start in the zero-error code state; at each readout time, decode a copy and record the bit."""
    _seed(seed)
    spins = np.zeros(nf, np.uint8); synd = np.zeros(ne, np.uint8)
    out = []
    t_prev = 0
    for t in t_grid:
        metropolis(spins, synd, f2e, beta, int(t - t_prev)); t_prev = t
        s2 = spins.copy(); y2 = synd.copy()
        cool(s2, y2, f2e, L * L)
        out.append(logical_bit(s2, M))
    return out, energy(synd) / nf


def local_autocorr(L, beta, seed, f2e, nf, ne, n_eq=400, n_obs=400, n_faces=64):
    """Equilibrium autocorrelation of face spins and edge syndromes (connected, normalised)."""
    _seed(seed)
    spins = np.zeros(nf, np.uint8); synd = np.zeros(ne, np.uint8)
    metropolis(spins, synd, f2e, beta, n_eq)
    rs = np.random.RandomState(seed)
    faces = rs.choice(nf, n_faces, replace=False); edges = rs.choice(ne, n_faces, replace=False)
    S = np.zeros((n_obs, n_faces)); Y = np.zeros((n_obs, n_faces))
    for t in range(n_obs):
        metropolis(spins, synd, f2e, beta, 1)
        S[t] = spins[faces]; Y[t] = synd[edges]
    def acf(X, maxlag):
        X = X - X.mean(0)
        var = (X * X).mean()
        if var <= 0:
            return [1.0] + [0.0] * maxlag
        return [float((X[:len(X) - k] * X[k:]).mean() / var) for k in range(maxlag + 1)]
    return acf(S, 30), acf(Y, 30)


def damage_spreading(L, beta, seed, f2e, nf, ne, n_eq=400, t_grid=(1, 2, 4, 8, 16, 32, 64, 128)):
    _seed(seed)
    sA = np.zeros(nf, np.uint8); yA = np.zeros(ne, np.uint8)
    metropolis(sA, yA, f2e, beta, n_eq)
    sB = sA.copy(); yB = yA.copy()
    f0 = nf // 2
    sB[f0] ^= 1
    for j in range(4):
        yB[f2e[f0, j]] ^= 1
    out = []; t_prev = 0
    for t in t_grid:
        metropolis_twin(sA, yA, sB, yB, f2e, beta, int(t - t_prev)); t_prev = t
        out.append(int((sA != sB).sum()))
    return out


def energy_scan(L, betas, seed, f2e, nf, ne, n_eq=300, n_obs=200):
    res = {}
    for direction, bl in (("heating", sorted(betas, reverse=True)), ("cooling", sorted(betas))):
        _seed(seed)
        spins = np.zeros(nf, np.uint8); synd = np.zeros(ne, np.uint8)
        if direction == "heating":
            pass  # start ordered (zero errors), go to low beta
        else:
            metropolis(spins, synd, f2e, 0.05, n_eq)  # start disordered
        for b in bl:
            metropolis(spins, synd, f2e, b, n_eq)
            acc = 0.0
            for _ in range(n_obs):
                metropolis(spins, synd, f2e, b, 1); acc += energy(synd)
            res.setdefault(f"{b:.3f}", {})[direction] = acc / n_obs / nf
    return res


def fit_tau(t, P):
    """Robust memory time: the half-life t_half at which P_log first falls to 0.5 (log-interpolated).
    Reported instead of an exponential fit, which is dominated by late-time noise once P ~ 0.
    tau_log = t_half / ln 2 under P = exp(-t/tau).  If P never reaches 0.5, t_max is a lower bound."""
    t = np.asarray(t, float); P = np.asarray(P, float)
    idx = np.where(P <= 0.5)[0]
    if len(idx) == 0:
        return {"t_half": None, "tau": None, "lower_bound_t_half": float(t.max()), "note": "P_log > 0.5 at t_max"}
    i = int(idx[0])
    if i == 0:
        return {"t_half": float(t[0]), "tau": float(t[0] / np.log(2)), "note": "decayed by first readout"}
    lt = np.log(t[i - 1]) + (0.5 - P[i - 1]) * (np.log(t[i]) - np.log(t[i - 1])) / (P[i] - P[i - 1])
    th = float(np.exp(lt))
    return {"t_half": th, "tau": th / np.log(2), "P_before": float(P[i - 1]), "P_after": float(P[i])}


def main(L, betas, ntraj, tmax, seed):
    f2e, M, nf, ne = build_tables(L)
    t_grid = np.unique(np.round(np.logspace(0, np.log10(tmax), 18)).astype(int))
    out = {"L": L, "n_faces": nf, "t_grid": t_grid.tolist(), "betas": {}}
    t0 = time.time()
    for beta in betas:
        bits = []; eng = []
        for k in range(ntraj):
            b, e = memory_trajectory(L, beta, t_grid, seed + 1000 * k + int(beta * 100), f2e, M, nf, ne)
            bits.append(b); eng.append(e)
        bits = np.array(bits, float)
        P = bits.mean(0)
        acf_s, acf_y = local_autocorr(L, beta, seed + 7, f2e, nf, ne)
        dmg = damage_spreading(L, beta, seed + 11, f2e, nf, ne)
        out["betas"][f"{beta:.2f}"] = {
            "P_log": P.tolist(), "P_log_se": (bits.std(0) / np.sqrt(ntraj)).tolist(),
            "tau_log": fit_tau(t_grid, P), "energy_per_face": float(np.mean(eng)),
            "acf_face": acf_s, "acf_edge": acf_y,
            "damage": {"t": [1, 2, 4, 8, 16, 32, 64, 128], "hamming": dmg},
        }
        print(f"L={L} beta={beta:.2f} P_log={np.round(P,3).tolist()} tau={out['betas'][f'{beta:.2f}']['tau_log']} ({time.time()-t0:.0f}s)", flush=True)
    return out


if __name__ == "__main__":
    L = int(sys.argv[1]); betas = [float(b) for b in sys.argv[2].split(",")]
    ntraj = int(sys.argv[3]); tmax = int(sys.argv[4]); seed = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    res = main(L, betas, ntraj, tmax, seed)
    outp = sys.argv[6] if len(sys.argv) > 6 else f"toric4d_L{L}.json"
    json.dump(res, open(outp, "w"), indent=1)
