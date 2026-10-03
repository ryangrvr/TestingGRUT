"""C3-D0 control computation (sec. R, D0-O8): N0, N1, P0, P1, P0-eps, P1-eps on exactly G2 = {16,32,64,128} and
G3 = {27,81,243}, seeds 1-20.  Primary Phi + preregistered Q1-Q6 algebra (sec. R, D0-O4) + secondary diagnostics
(cannot rescue Phi).  PM2 trees are NOT scored.  Independent code path, not independent reviewer."""
import sys, os, json, time
import numpy as np
from multiprocessing import Pool
from scipy.stats import t as tdist, ks_2samp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d0_detector import (phi, exact_classes, light_depth, light_fraction_top, horton_tokunaga, ternary_tree,
                         bisection_eps, icbrt_floor)
from pm2_stageA import lattice, sp_tree, bisection_tree, wilson_ust, p1_p2, strahler, _bfs_order

G2, G3, SEEDS = [16, 32, 64, 128], [27, 81, 243], list(range(1, 21))
TQ = tdist.ppf(0.975, 19)


def make(kind, L, seed):
    N, E = lattice(L)
    if kind == "N0": return np.asarray(wilson_ust(N, E, np.random.default_rng(seed))[0])
    if kind == "N1": return np.asarray(sp_tree(N, E)[0])
    if kind == "P0": return np.asarray(bisection_tree(L)[0])
    if kind == "P1": return ternary_tree(L)
    if kind == "P0e": return bisection_eps(L, np.random.default_rng(seed), icbrt_floor(L * L))
    if kind == "P1e": return ternary_tree(L, np.random.default_rng(seed), icbrt_floor(L * L))
    raise ValueError(kind)


def score(args):
    kind, L, seed = args; t0 = time.time(); par = make(kind, L, seed); N = L * L
    v, nU, nc = phi(par, True); order = _bfs_order(N, par)
    tau, eta = p1_p2(N, par, order); RB, T1, T2, c, K = horton_tokunaga(par)
    return dict(kind=kind, L=L, seed=seed, phi=v, nU=nU, ncls=nc, D=exact_classes(par), lam=light_depth(par),
                tau=float(tau), eta=float(eta), strahler=strahler(N, par, order), RB=RB, T1=T1, T2=T2, c=c,
                q=light_fraction_top(par).tolist(), secs=time.time() - t0)


def ens(x):
    x = np.asarray(x, float); m = x.mean(); s = x.std(ddof=1)
    h = TQ * s / np.sqrt(len(x)) if s > 0 else 0.0
    return m, h


if __name__ == "__main__":
    tasks = []
    for L in G2 + G3:
        tasks += [("N0", L, s) for s in SEEDS] + [("N1", L, 0)]
    tasks += [("P0", L, 0) for L in G2] + [("P0e", L, s) for L in G2 for s in SEEDS]
    tasks += [("P1", L, 0) for L in G3] + [("P1e", L, s) for L in G3 for s in SEEDS]
    tasks.sort(key=lambda a: -a[1])                       # largest first for load balance
    t0 = time.time()
    with Pool(4) as pool:
        res = pool.map(score, tasks, chunksize=1)
    print(f"computed {len(res)} trees in {time.time() - t0:.0f}s", flush=True)
    R = {}
    for r in res:
        R.setdefault((r['kind'], r['L']), []).append(r)
    for k in R: R[k].sort(key=lambda r: r['seed'])
    json.dump({f"{k[0]}_{k[1]}": [{kk: vv for kk, vv in r.items() if kk != 'q'} for r in v] for k, v in R.items()},
              open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "d0_results.json"), "w"), indent=0)

    def row(kind, L):
        rs = R[(kind, L)]; ph = [r['phi'] for r in rs]
        if len(rs) == 1:
            r = rs[0]; return dict(m=r['phi'], h=0.0, nU=r['nU'], ncls=r['ncls'], D=r['D'], lam=r['lam'], tau=r['tau'],
                                   eta=r['eta'], st=r['strahler'], RB=r['RB'], T1=r['T1'], T2=r['T2'], c=r['c'])
        m, h = ens(ph); g = lambda key: float(np.nanmean([r[key] for r in rs]))
        return dict(m=m, h=h, nU=g('nU'), ncls=g('ncls'), D=g('D'), lam=g('lam'), tau=g('tau'), eta=g('eta'),
                    st=g('strahler'), RB=g('RB'), T1=g('T1'), T2=g('T2'), c=g('c'))

    print("\n=== PRIMARY Phi per family and size (mean, 95% half-width; |U| and #classes are means for ensembles) ===")
    for grid, fam, fame in [(G2, "P0", "P0e"), (G3, "P1", "P1e")]:
        for L in grid:
            for kind in ["N0", "N1", fam, fame]:
                x = row(kind, L)
                print(f"  L={L:3d} {kind:4s}: Phi {x['m']:.4f} +- {x['h']:.4f} | |U| {x['nU']:.1f} | classes {x['ncls']:.1f}")

    print("\n=== QUALIFICATION (exact algebra of sec. R, D0-O4) ===")
    verdict = {}
    for grid, fam, fame in [(G2, "P0", "P0e"), (G3, "P1", "P1e")]:
        Rr = [row("N0", L) for L in grid]; N1 = [row("N1", L)['m'] for L in grid]
        F = [row(fam, L)['m'] for L in grid]; Fe = [row(fame, L) for L in grid]
        upR = [r['m'] + r['h'] for r in Rr]
        Q1 = all(F[j] > upR[j] for j in range(len(grid)))
        Q2 = all(F[j] > N1[j] for j in range(len(grid)))
        D = [F[j] - Rr[j]['m'] for j in range(len(grid))]
        Q3 = Q1 and all(D[j + 1] >= D[j] - (Rr[j]['h'] + Rr[j + 1]['h']) for j in range(len(grid) - 1))
        lo = [x['m'] - x['h'] for x in Fe]
        Q61 = all(lo[j] - upR[j] > 0 for j in range(len(grid)))
        Q62 = all(lo[j] > N1[j] for j in range(len(grid)))
        De = [Fe[j]['m'] - Rr[j]['m'] for j in range(len(grid))]
        Q63 = all(De[j + 1] >= De[j] - (Fe[j]['h'] + Fe[j + 1]['h'] + Rr[j]['h'] + Rr[j + 1]['h']) for j in range(len(grid) - 1))
        Q4 = True    # detector reads the parent array only: code inspection + relabelling tests (d0_tests.log)
        Q5 = all(np.isfinite(row("N0", L)['m']) and np.isfinite(row("N1", L)['m']) for L in grid)
        verdict[fam] = dict(Q1=Q1, Q2=Q2, Q3=Q3, Q4=Q4, Q5=Q5); verdict[fame] = dict(Q61=Q61, Q62=Q62, Q63=Q63, Q6=Q61 and Q62 and Q63)
        print(f"  {fam}: Delta_j {[round(d, 4) for d in D]}; up_R {[round(u, 4) for u in upR]}; h_R {[round(r['h'], 4) for r in Rr]}")
        print(f"  {fam}: Q1 {Q1}  Q2 {Q2}  Q3 {Q3}  Q4 {Q4}  Q5 {Q5}")
        print(f"  {fame}: Delta_low_j {[round(lo[j] - upR[j], 4) for j in range(len(grid))]}; Delta_eps_j {[round(d, 4) for d in De]}; "
              f"h_eps {[round(x['h'], 4) for x in Fe]}")
        print(f"  {fame}: Q6.1 {Q61}  Q6.2 {Q62}  Q6.3 {Q63}  => Q6 {Q61 and Q62 and Q63}")
    det = all(all(verdict[f].values()) for f in ("P0", "P1")); eps = verdict["P0e"]["Q6"] and verdict["P1e"]["Q6"]
    grade = "D0-A-N" if det and eps else ("D0-A-T" if det else "D0-B")
    print(f"\n  GRADE (mechanical): {grade}")

    print("\n=== SECONDARY DIAGNOSTICS (reported; cannot rescue Phi) ===")
    for grid, fam, fame in [(G2, "P0", "P0e"), (G3, "P1", "P1e")]:
        for kind in ["N0", "N1", fam, fame]:
            xs = [row(kind, L) for L in grid]; Ns = [L * L for L in grid]
            sl = np.polyfit(np.log(Ns), np.log([x['D'] for x in xs]), 1)[0]
            print(f"  {kind:4s} D(T) {[round(x['D'], 1) for x in xs]} slope(log D vs log N) {sl:.3f}")
            for L, x in zip(grid, xs):
                print(f"     L={L:3d}: lambda {x['lam']:.4f}  tau {x['tau']:.4f}  eta_H {x['eta']:.4f}  Strahler {x['st']:.2f}  "
                      f"R_B {x['RB']:.3f}  T1 {x['T1']:.3f}  T2 {x['T2']:.3f}  c {x['c']:.3f}")
        for L in grid:
            qR = np.concatenate([r['q'] for r in R[("N0", L)]])
            for kind in [fam, fame, "N1"]:
                qF = np.concatenate([r['q'] for r in R[(kind, L)]])
                ks = ks_2samp(qF, qR).statistic if len(qF) and len(qR) else float('nan')
                print(f"     KS(light fraction, top window) {kind} vs N0 at L={L}: {ks:.4f}  (n_F {len(qF)}, n_R {len(qR)})")
