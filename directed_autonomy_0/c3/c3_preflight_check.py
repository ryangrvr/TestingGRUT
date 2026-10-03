"""DA0 / C3 analytic preflight -- SECONDARY sanity check only (the results are the analytic derivations in
C3_ANALYTIC_PREFLIGHT.md).  Tiny finite graphs, exact enumeration, mu = 1 (T_lin = 1).  Not a C3 simulation.

Checks: (1) envelope inequality F(J) >= F(|J|) with equality iff unfrustrated on the support; (2) the local-minimum
frustration bound: a frustrated cycle of length l can lie in the support of a local minimum only if T >= 2/(mu*l);
(3) small-K series of the uniform branch e(K)/K on the single cube.
"""
import itertools
import numpy as np
from scipy.optimize import minimize

MU = 1.0


def graph_plaquette():
    return 4, [(0, 1), (1, 2), (2, 3), (3, 0)], [[0, 1, 2, 3]]


def graph_cube():
    V = list(itertools.product([0, 1], repeat=3)); idx = {v: i for i, v in enumerate(V)}
    E = [(idx[a], idx[b]) for a in V for b in V if idx[a] < idx[b] and sum(x != y for x, y in zip(a, b)) == 1]
    faces = []
    for axis in range(3):
        for val in (0, 1):
            vs = [v for v in V if v[axis] == val]; o = [w for w in range(3) if w != axis]
            cyc = [vs[0]]
            for _ in range(3):
                cur = cyc[-1]
                nxt = [w for w in vs if sum(x != y for x, y in zip(cur, w)) == 1 and w not in cyc][0]; cyc.append(nxt)
            bonds = []
            for k in range(4):
                a, b = idx[cyc[k]], idx[cyc[(k + 1) % 4]]; bonds.append(E.index((min(a, b), max(a, b))))
            faces.append(bonds)
    return len(V), E, faces


def make(nV, E):
    S = np.array(list(itertools.product([-1, 1], repeat=nV)))
    s = np.array([S[:, i] * S[:, j] for i, j in E]).T            # configurations x bonds
    def stats(J, T):
        w = s @ J / T; w = np.exp(w - w.max()); p = w / w.sum()
        m = p @ s; C = (s - m).T @ (p[:, None] * (s - m))
        return m, C, np.log(np.sum(np.exp(s @ J / T)))
    def F(J, T):
        return 0.5 * MU * J @ J - T * np.log(np.sum(np.exp(s @ J / T)))
    def grad(J, T):
        return MU * J - stats(J, T)[0]
    return stats, F, grad


ZERO = 1e-6     # computational tolerance: |J_b| below this is treated as an absent bond (round-off of the J = 0 minimum)


def frustrated(J, faces):
    return [f for f in faces if np.all(np.abs(J[f]) > ZERO) and np.prod(np.sign(J[f])) < 0]


def check_envelope(name, nV, E, faces):
    stats, F, grad = make(nV, E); rng = np.random.default_rng(1); worst = np.inf; eqv = 0; strict = 0
    for _ in range(2000):
        J = rng.normal(size=len(E)); T = rng.uniform(0.2, 2.0)
        d = F(J, T) - F(np.abs(J), T); worst = min(worst, d)
        if frustrated(J, faces):
            strict += d > 1e-12
        else:
            eqv += 1
    print(f"  {name}: min[F(J) - F(|J|)] over 2000 random (J,T) = {worst:.2e} (>= 0 required); "
          f"frustrated samples with strict inequality: {strict}")


def local_minima(name, nV, E, faces, Ts, starts=400):
    stats, F, grad = make(nV, E); rng = np.random.default_rng(2)
    shortest = 4
    for T in Ts:
        mins = []
        for _ in range(starts):
            J0 = rng.normal(scale=1.5, size=len(E))
            r = minimize(F, J0, args=(T,), jac=grad, method='BFGS', options=dict(gtol=1e-10, maxiter=4000))
            J = r.x; m, C, _ = stats(J, T); H = MU * np.eye(len(E)) - C / T
            if np.linalg.eigvalsh(H).min() <= 1e-9 or np.linalg.norm(grad(J, T)) > 1e-7:
                continue
            key = (round(F(J, T), 8),)
            if not any(abs(F(J, T) - F(J2, T)) < 1e-8 and np.allclose(np.abs(J), np.abs(J2), atol=1e-6) for J2 in mins):
                mins.append(J)
        nfr = sum(1 for J in mins if frustrated(J, faces))
        diag_ok = all(np.all(1 - stats(J, T)[0] ** 2 <= MU * T + 1e-9) for J in mins)
        desc = sorted(set(f"|J|~{np.round(np.abs(J).mean(), 3)}" for J in mins))
        print(f"  {name} T={T:.2f} (bound: frustrated {shortest}-cycle needs T >= {2 / (MU * shortest):.2f}): minima {desc}; "
              f"distinct stable minima {len(mins)}, of which frustrated {nfr}; all satisfy Var(s_b) <= mu*T: {diag_ok}")


def uniform_branch():
    nV, E, faces = graph_cube(); stats, F, grad = make(nV, E)
    print("  single cube, uniform J = K*T: e(K)/K where e = <s_b> (bond-averaged); small-K series predicts 1 + c K^2 with c > 0")
    for K in [0.05, 0.1, 0.2, 0.4, 0.8, 1.6]:
        m, C, _ = stats(np.full(len(E), K), 1.0); e = m.mean()
        print(f"    K={K}: e/K = {e / K:.4f}")


if __name__ == "__main__":
    print("== (1) envelope inequality ==")
    check_envelope("plaquette", *graph_plaquette()); check_envelope("cube (8 spins, 12 bonds, 6 faces)", *graph_cube())
    print("== (2) local minima vs frustration bound (mu = 1) ==")
    Ts = [0.3, 0.45, 0.55, 0.7, 0.9, 1.1, 1.3]
    local_minima("plaquette", *graph_plaquette(), Ts)
    local_minima("cube", *graph_cube(), Ts, starts=300)
    print("== (3) uniform branch on the cube (finite-graph illustration only) ==")
    uniform_branch()
