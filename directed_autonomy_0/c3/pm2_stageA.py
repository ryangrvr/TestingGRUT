"""DA0 / C3 / PM2 Stage A (deterministic) -- adaptive conserved-flow network.  NUMERICAL ILLUSTRATION; independent code
path, not independent reviewer.  Law, drive, constants, sizes, seeds and statistics are frozen by C3_PM2_CHARTER.md +
REPAIR 01 (sec. R).  Implementation (numerical) parameters are declared in IMPL and are not model parameters.
"""
import sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from collections import deque

GAMMA, NU, KAPPA, DELTA = 0.5, 1.0, 1.0, 0.01
IMPL = dict(max_rel_change=0.0125, prune_floor=1e-12, stat_tol=1e-6, stable_steps=2000, step_cap=400000, y_ref=1e-3)


def lattice(L):
    idx = lambda x, y: x * L + y
    E = []
    for x in range(L):
        for y in range(L):
            if x + 1 < L: E.append((idx(x, y), idx(x + 1, y)))
            if y + 1 < L: E.append((idx(x, y), idx(x, y + 1)))
    return L * L, np.array(E)


def solve_flow(N, E, C, h):
    """Kirchhoff: Q_e = C_e (p_i - p_j); grounded at outlet node 0.  Active edges only."""
    act = C > 0; Ea, Ca = E[act], C[act]
    i, j = Ea[:, 0], Ea[:, 1]
    Lap = sp.coo_matrix((np.concatenate([Ca, Ca, -Ca, -Ca]), (np.concatenate([i, j, i, j]), np.concatenate([i, j, j, i]))),
                        shape=(N, N)).tocsr()
    p = np.zeros(N); p[1:] = spl.spsolve(Lap[1:, 1:].tocsc(), h[1:])
    Q = np.zeros(len(C)); Q[act] = Ca * (p[i] - p[j]); return Q


def energy(C, Q):
    act = C > 0
    return float(np.sum(Q[act] ** 2 / C[act]) + NU * np.sum(C[act] ** GAMMA))


def adapt(N, E, h, C0, early_stop=True):
    C = C0.copy(); stable = 0; nact_prev = -1; steps = 0; early = False
    while steps < IMPL['step_cap']:
        Q = solve_flow(N, E, C, h); act = C > 0
        rate = np.zeros_like(C); rate[act] = KAPPA * (Q[act] ** 2 / C[act] - NU * GAMMA * C[act] ** GAMMA)
        rel = np.abs(rate[act]) / C[act]
        resid = np.max(np.abs(Q[act] ** 2 / C[act] - NU * GAMMA * C[act] ** GAMMA) / (NU * GAMMA * C[act] ** GAMMA))
        nact = int(act.sum())
        stable = stable + 1 if nact == nact_prev else 0; nact_prev = nact
        if resid < IMPL['stat_tol'] and stable >= IMPL['stable_steps']:
            break
        if early_stop and nact == N - 1 and tree_from_edges(N, E, act)[2]:
            # EXACT early stop (implementation refinement, ledger PM2-I1): on a connected spanning-tree support the
            # Kirchhoff flows are Q_e = +-S_e != 0 independent of C, so each active edge obeys an autonomous 1-D ODE with
            # unique attracting root C* = (S_e^2/(nu*gamma))^(1/(gamma+1)) > 0; no active edge can be pruned, and pruned
            # edges never revive (C = 0 absorbing).  Final topology is therefore fixed; set C to its exact limit.
            parent, order, _ = tree_from_edges(N, E, act)
            C = stationary_C_tree(N, parent, order, E)[0]; resid = 0.0; early = True
            break
        dt = IMPL['max_rel_change'] / max(rel.max(), 1e-300)
        Cn = C + dt * rate
        floor = IMPL['prune_floor'] * Cn.max()
        Cn[(Cn < floor)] = 0.0
        C = Cn; steps += 1
    Q = solve_flow(N, E, C, h)
    return C, Q, steps, resid


def adapt_y(N, E, h, C0, early_stop=True):
    """Same ODE, integrated in y = sqrt(C) (gamma = 1/2 only).  Using Q_e = C_e dp_e:
        dy_e/dt = (kappa/2)(y_e dp_e^2 - nu*gamma)      [identical to dC/dt = kappa(Q^2/C - nu*gamma*C^gamma)]
    which is non-stiff as C -> 0 (finite-time extinction at speed kappa*nu*gamma/2).  Explicit Euler in y; dt so that
    max_e |dy_e| / max(y_e, y_ref * max y) <= max_rel_change; an edge whose y crosses 0 (or C < prune_floor*max C) is
    pruned (exact finite-time extinction of the ODE).  Implementation refinement PM2-I2; same stop rules as adapt."""
    assert GAMMA == 0.5
    C = C0.copy(); stable = 0; nact_prev = -1; steps = 0
    while steps < IMPL['step_cap']:
        Q = solve_flow(N, E, C, h); act = C > 0; y = np.sqrt(C[act])
        dp = Q[act] / C[act]
        resid = np.max(np.abs(Q[act] ** 2 / C[act] - NU * GAMMA * C[act] ** GAMMA) / (NU * GAMMA * C[act] ** GAMMA))
        nact = int(act.sum()); stable = stable + 1 if nact == nact_prev else 0; nact_prev = nact
        if resid < IMPL['stat_tol'] and stable >= IMPL['stable_steps']:
            break
        if early_stop and nact == N - 1 and tree_from_edges(N, E, act)[2]:
            parent, order, _ = tree_from_edges(N, E, act)
            C = stationary_C_tree(N, parent, order, E)[0]; resid = 0.0
            break
        ydot = 0.5 * KAPPA * (y * dp ** 2 - NU * GAMMA)
        dt = IMPL['max_rel_change'] / max(np.max(np.abs(ydot) / np.maximum(y, IMPL['y_ref'] * y.max())), 1e-300)
        yn = y + dt * ydot; Cn = np.where(yn > 0, yn, 0.0) ** 2
        Cn[Cn < IMPL['prune_floor'] * Cn.max()] = 0.0
        C = C.copy(); C[act] = Cn; steps += 1
    Q = solve_flow(N, E, C, h)
    return C, Q, steps, resid


# ------------------------------------------------------------------ tree statistics
def tree_from_edges(N, E, mask):
    adj = [[] for _ in range(N)]
    for (a, b) in E[mask]:
        adj[a].append(b); adj[b].append(a)
    parent = -np.ones(N, int); order = []; seen = np.zeros(N, bool); seen[0] = True; dq = deque([0])
    while dq:
        u = dq.popleft(); order.append(u)
        for v in sorted(adj[u]):
            if not seen[v]:
                seen[v] = True; parent[v] = u; dq.append(v)
    return parent, order, bool(seen.all())


def subtree_stats(N, parent, order):
    S = np.ones(N, int); H = np.zeros(N, int); depth = np.zeros(N, int)
    for u in order[1:]:
        depth[u] = depth[parent[u]] + 1
    for u in reversed(order[1:]):
        p = parent[u]; S[p] += S[u]; H[p] = max(H[p], H[u] + 1)
    return S, H, depth


def strahler(N, parent, order):
    kids = [[] for _ in range(N)]
    for u in order[1:]:
        kids[parent[u]].append(u)
    o = np.ones(N, int)
    for u in reversed(order):
        if kids[u]:
            ks = [o[k] for k in kids[u]]; m = max(ks); o[u] = m + 1 if ks.count(m) >= 2 else m
    return int(o[0])


def p1_p2(N, parent, order):
    S, H, _ = subtree_stats(N, parent, order)
    Se = S[order[1:]]; Lb = 1 + H[order[1:]]
    lo, hi = np.sqrt(N / 10), np.sqrt(10 * N)
    sk = np.unique(np.round(np.logspace(np.log10(lo), np.log10(hi), 20)).astype(int))
    ccdf = np.array([np.mean(Se >= s) for s in sk])
    ok = ccdf > 0
    tau = -np.polyfit(np.log(sk[ok]), np.log(ccdf[ok]), 1)[0] if ok.sum() >= 3 else np.nan
    w = (Se >= lo) & (Se <= hi)
    eta = np.polyfit(np.log(Se[w]), np.log(Lb[w]), 1)[0] if w.sum() >= 3 and len(np.unique(Se[w])) >= 3 else np.nan
    return tau, eta


def p3(N, E, parent, order):
    S, _, depth = subtree_stats(N, parent, order)
    tree = set()
    for u in order[1:]:
        tree.add((min(u, parent[u]), max(u, parent[u])))
    vals = []
    for (a, b) in E:
        if (min(a, b), max(a, b)) in tree:
            continue
        x, y = a, b
        while x != y:                       # walk to LCA; each tree edge (x -> parent) on the cycle reroutes S[x]
            if depth[x] >= depth[y]:
                vals.append(S[x]); x = parent[x]
            else:
                vals.append(S[y]); y = parent[y]
    v = np.array(vals); return float(np.median(v)), float(np.mean(v))


def stationary_C_tree(N, parent, order, E):
    S, _, _ = subtree_stats(N, parent, order); C = np.zeros(len(E)); eid = {(min(a, b), max(a, b)): k for k, (a, b) in enumerate(E)}
    for u in order[1:]:
        C[eid[(min(u, parent[u]), max(u, parent[u]))]] = (S[u] ** 2 / (NU * GAMMA)) ** (1 / (GAMMA + 1))
    return C, eid


def static_barriers(N, E, h, parent, order, rng, nswap=10):
    """PM2R-06 diagnostic: straight-path max of E between stationary tree T_a and a random fundamental swap T_b."""
    S, _, depth = subtree_stats(N, parent, order); Ca, eid = stationary_C_tree(N, parent, order, E)
    nontree = [k for k in range(len(E)) if Ca[k] == 0]; out = []
    for _ in range(nswap):
        f = nontree[rng.integers(len(nontree))]; a, b = E[f]; x, y = a, b; cyc = []
        while x != y:
            if depth[x] >= depth[y]: cyc.append(x); x = parent[x]
            else: cyc.append(y); y = parent[y]
        u = cyc[rng.integers(len(cyc))]
        newpar = parent.copy()
        # re-root: remove edge (u, parent[u]); attach B_u through f
        side = a if _in_subtree(u, a, parent) else b; other = b if side == a else a
        path = []; z = side
        while z != u: path.append(z); z = parent[z]
        path.append(u)
        prev = other
        for z in path:
            nxt = newpar[z]; newpar[z] = prev; prev = z
        order_b = _bfs_order(N, newpar)
        Cb, _ = stationary_C_tree(N, newpar, order_b, E)
        Ea_ = energy(Ca, solve_flow(N, E, Ca, h)); Es = []
        for t in np.linspace(0, 1, 11):
            Ct = (1 - t) * Ca + t * Cb; Es.append(energy(Ct, solve_flow(N, E, Ct, h)))
        out.append(max(Es) - Ea_)
    return np.array(out)


def _in_subtree(u, w, parent):
    while w != -1:
        if w == u: return True
        w = parent[w]
    return False


def _bfs_order(N, parent):
    kids = [[] for _ in range(N)]
    for v in range(1, N):
        kids[parent[v]].append(v)
    order = [0]; i = 0
    while i < len(order):
        order.extend(kids[order[i]]); i += 1
    return order


# ------------------------------------------------------------------ comparators
def wilson_ust(N, E, rng):
    adj = [[] for _ in range(N)]
    for (a, b) in E:
        adj[a].append(b); adj[b].append(a)
    intree = np.zeros(N, bool); intree[0] = True; parent = -np.ones(N, int); nxt = -np.ones(N, int)
    for s in range(N):
        u = s
        while not intree[u]:
            nxt[u] = adj[u][rng.integers(len(adj[u]))]; u = nxt[u]
        u = s
        while not intree[u]:
            parent[u] = nxt[u]; intree[u] = True; u = nxt[u]
    return parent, _bfs_order(N, parent)


def sp_tree(N, E):
    mask = np.ones(len(E), bool); parent, order, _ = tree_from_edges(N, E, mask); return parent, order


def bisection_tree(L):
    idx = lambda x, y: x * L + y; parent = -np.ones(L * L, int)
    def build(x0, x1, y0, y1, r):                     # rectangle [x0,x1) x [y0,y1), exit cell r inside
        if (x1 - x0) * (y1 - y0) == 1:
            return
        rx, ry = divmod(r, L)
        if x1 - x0 >= y1 - y0:
            xm = (x0 + x1) // 2; A = (x0, xm, y0, y1) if rx < xm else (xm, x1, y0, y1)
            B = (xm, x1, y0, y1) if rx < xm else (x0, xm, y0, y1); ym = (y0 + y1) // 2
            bcell = idx(xm, ym) if rx < xm else idx(xm - 1, ym); acell = idx(xm - 1, ym) if rx < xm else idx(xm, ym)
        else:
            ym = (y0 + y1) // 2; A = (x0, x1, y0, ym) if ry < ym else (x0, x1, ym, y1)
            B = (x0, x1, ym, y1) if ry < ym else (x0, x1, y0, ym); xm = (x0 + x1) // 2
            bcell = idx(xm, ym) if ry < ym else idx(xm, ym - 1); acell = idx(xm, ym - 1) if ry < ym else idx(xm, ym)
        build(*A, r); build(*B, bcell); parent[bcell] = acell
    build(0, L, 0, L, 0)
    return parent, _bfs_order(L * L, parent)


def tree_summary(N, E, parent, order):
    tau, eta = p1_p2(N, parent, order); med, mean = p3(N, E, parent, order); return tau, eta, med, mean, strahler(N, parent, order)


def tint(x):
    x = np.asarray([v for v in x if np.isfinite(v)]); n = len(x)
    from scipy.stats import t as tdist
    h = tdist.ppf(0.975, n - 1) * x.std(ddof=1) / np.sqrt(n); return x.mean(), x.mean() - h, x.mean() + h


def pm2_run(args):
    """One PM2 run (seed 0 = uniform start).  Returns tree flag, summary, active-edge set, energies, static barriers."""
    L, seed = args; N, E = lattice(L); h = np.ones(N); h[0] = -(N - 1); t0 = time.time()
    C0 = np.ones(len(E)) if seed == 0 else 1 + DELTA * np.random.default_rng(seed).uniform(-1, 1, len(E))
    C, Q, steps, resid = INTEGRATOR(N, E, h, C0); act = C > 0
    parent, order, conn = tree_from_edges(N, E, act); is_tree = bool(conn and act.sum() == N - 1)
    out = dict(L=L, seed=seed, steps=steps, resid=float(resid), nact=int(act.sum()), tree=is_tree, E=energy(C, Q),
               act=tuple(int(k) for k in np.nonzero(act)[0]), secs=time.time() - t0)
    if is_tree:
        out['summary'] = tree_summary(N, E, parent, order)
        out['bars'] = static_barriers(N, E, h, parent, order, np.random.default_rng(1000 + seed)) if seed > 0 else None
    return out


INTEGRATOR = adapt_y          # PM2-I2 (validated against adapt, which is validated against full Euler; ledger)

if __name__ == "__main__":
    from multiprocessing import Pool
    argv = sys.argv[1:]
    if argv and argv[0].startswith("--rel="):        # PM2-I4 robustness pass (y-Euler at 0.05)
        IMPL['max_rel_change'] = float(argv.pop(0)[6:])
    Ls = [int(a) for a in argv] or [16, 24, 32, 48, 64]
    names = ["tau", "eta_H", "P3_median", "P3_mean", "Strahler"]
    print("IMPL (implementation parameters, not model parameters):", IMPL, "integrator:", INTEGRATOR.__name__, flush=True)
    for L in Ls:
        N, E = lattice(L); h = np.ones(N); h[0] = -(N - 1); t0 = time.time()
        print(f"=== L={L} N={N} edges={len(E)} ===", flush=True)
        with Pool(4, initializer=IMPL.update, initargs=(dict(IMPL),)) as pool:
            res = pool.map(pm2_run, [(L, s) for s in range(0, 21)], chunksize=1)
        u = res[0]
        print(f"  PM2 uniform start: steps {u['steps']}, resid {u['resid']:.1e}, active {u['nact']} (tree needs {N - 1}), "
              f"tree {u['tree']}, E {u['E']:.6e}, Ebar {u['E'] / (N - 1) ** (2 / 3):.6e}", flush=True)
        if u['tree']:
            print("    uniform-start tree stats (tau, eta_H, P3 median, P3 mean, Strahler):", tuple(np.round(u['summary'], 4)))
        pr = res[1:]; good = [r for r in pr if r['tree']]
        for r in pr:
            print(f"    seed {r['seed']:2d}: steps {r['steps']}, active {r['nact']}, tree {r['tree']}, Ebar "
                  f"{r['E'] / (N - 1) ** (2 / 3):.6e}, " + (", ".join(f"{n} {v:.4f}" for n, v in zip(names, r['summary']))
                                                            if r['tree'] else "") + f", {r['secs']:.0f}s")
        rows = np.array([r['summary'] for r in good]); Eb = np.array([r['E'] for r in good]) / (N - 1) ** (2 / 3)
        print(f"  PM2 perturbed: tree runs {len(good)}/20 (non-tree {20 - len(good)}); distinct trees (diagnostic only) "
              f"{len(set(r['act'] for r in good))}; uniform-start tree among them: {u['act'] in set(r['act'] for r in good)}")
        print(f"    Ebar over seeds: min {Eb.min():.6e}  median {np.median(Eb):.6e}  max {Eb.max():.6e}  (max-min)/median "
              f"{(Eb.max() - Eb.min()) / np.median(Eb):.2e}")
        for k, name in enumerate(names):
            m, lo, hi = tint(rows[:, k]); print(f"    PM2 {name}: mean {m:.4f}  95% [{lo:.4f}, {hi:.4f}]")
        bars = np.concatenate([r['bars'] for r in good])
        print(f"    static straight-path barrier diagnostic ({len(bars)} swaps): median raw {np.median(bars):.4e}, median reduced "
              f"{np.median(bars) / (N - 1) ** (2 / 3):.4e}, max raw {bars.max():.4e}, fraction > 0 {np.mean(bars > 0):.3f}")
        rr = np.array([tree_summary(N, E, *wilson_ust(N, E, np.random.default_rng(s))) for s in range(1, 21)])
        for k, name in enumerate(names):
            m, lo, hi = tint(rr[:, k]); print(f"    R-TREE {name}: mean {m:.4f}  95% [{lo:.4f}, {hi:.4f}]")
        sp_ = tree_summary(N, E, *sp_tree(N, E)); a2 = tree_summary(N, E, *bisection_tree(L))
        print(f"    SP-TREE (tau, eta_H, P3 med, P3 mean, Strahler): {tuple(np.round(sp_, 4))}")
        print(f"    A2 bisection (supplied hierarchy): {tuple(np.round(a2, 4))}")
        Q0 = solve_flow(N, E, np.ones(len(E)), h)
        print(f"    A0 frozen C=1: active edges {len(E)} (no topology selection); max |Q| {np.abs(Q0).max():.1f}; "
              f"E {energy(np.ones(len(E)), Q0):.6e}, Ebar {energy(np.ones(len(E)), Q0) / (N - 1) ** (2 / 3):.6e}")
        CA1 = (Q0 ** 2 / (NU * GAMMA)) ** (1 / (GAMMA + 1))
        print(f"    A1 flow-decoupled (stationary C under fixed Q0): active edges {int((CA1 > IMPL['prune_floor'] * CA1.max()).sum())}"
              f" of {len(E)}; spanning tree? {int((CA1 > IMPL['prune_floor'] * CA1.max()).sum()) == N - 1}")
        print(f"  [L={L} done in {time.time() - t0:.0f}s]", flush=True)
