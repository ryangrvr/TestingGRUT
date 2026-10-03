"""PM2 Stage-A aggregation: criteria across L (PM2R-04/05/06, Revised Stage-A positive, PM2-I4).  Parses per-seed PM2
rows from the Stage-A logs; recomputes the (cheap, deterministic) R-TREE / SP-TREE / A2 statistics with the same code.
Independent code path, not independent reviewer.  Usage: python pm2_stageA_summary.py pm2_stageA.log [pm2_stageA_rel005.log]

Declared reading of item 4 (charter gives no numeric rule; fixed here before aggregation, conservative):
  "Strahler growing beyond the controls" := the PM2 ensemble Strahler slope d(Strahler)/d(log2 N) has a 95% interval
  entirely above the R-TREE slope interval AND above the SP-TREE slope (OLS on its 5 points), AND at L = 48 and 64 the PM2
  Strahler interval lies strictly above the R-TREE interval and above the SP-TREE value.
P3: rho from pooled OLS of log(P3 median) on log N over all PM2 tree runs (100 points), 95% t-interval; rho > 0 passes iff
  the interval lower bound > 0.  The same fit is reported for R-TREE (pooled) and SP-TREE (5 points) as context.
"""
import re, sys
import numpy as np
from scipy.stats import t as tdist
from pm2_stageA import lattice, tree_summary, wilson_ust, sp_tree, bisection_tree, tint

NAMES = ["tau", "eta_H", "P3_median", "P3_mean", "Strahler"]


def parse(path):
    out = {}; L = None
    for line in open(path):
        m = re.match(r"=== L=(\d+)", line)
        if m: L = int(m.group(1)); out[L] = []; continue
        m = re.match(r"\s+seed\s+(\d+): .*tree True.*Ebar ([\d.e+-]+), tau ([\d.na-]+), eta_H ([\d.na-]+), P3_median ([\d.]+), "
                     r"P3_mean ([\d.]+), Strahler ([\d.]+)", line)
        if m and L: out[L].append([float(m.group(k)) for k in range(3, 8)])
    return {L: np.array(v) for L, v in out.items() if len(v)}


def slope_ci(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float); n = len(x); X = np.vstack([x, np.ones(n)]).T
    b, res, *_ = np.linalg.lstsq(X, y, rcond=None); r = y - X @ b; s2 = r @ r / (n - 2)
    se = np.sqrt(s2 / np.sum((x - x.mean()) ** 2)); h = tdist.ppf(0.975, n - 2) * se; return b[0], b[0] - h, b[0] + h


def controls(Ls):
    c = {}
    for L in Ls:
        N, E = lattice(L)
        c[L] = dict(R=np.array([tree_summary(N, E, *wilson_ust(N, E, np.random.default_rng(s))) for s in range(1, 21)]),
                    SP=np.array(tree_summary(N, E, *sp_tree(N, E))), A2=np.array(tree_summary(N, E, *bisection_tree(L))))
    return c


def verdicts(pm, c):
    Ls = sorted(pm); v = {}; big = [L for L in Ls if L in (48, 64)]
    for k, name in [(0, "P1"), (1, "P2")]:
        ok = []
        for L in big:
            m, lo, hi = tint(pm[L][:, k]); rm, rlo, rhi = tint(c[L]['R'][:, k]); sp = c[L]['SP'][k]
            ok.append((hi < rlo or lo > rhi) and not (lo <= sp <= hi))
        v[name] = bool(len(big) == 2 and all(ok))
    x = np.concatenate([[np.log(L * L)] * len(pm[L]) for L in Ls]); y = np.concatenate([np.log(pm[L][:, 2]) for L in Ls])
    rho = slope_ci(x, y); v['rho'] = rho; v['P3'] = bool(rho[1] > 0)
    xs = np.concatenate([[np.log2(L * L)] * len(pm[L]) for L in Ls]); ys = np.concatenate([pm[L][:, 4] for L in Ls])
    st = slope_ci(xs, ys)
    xr = np.concatenate([[np.log2(L * L)] * 20 for L in Ls]); yr = np.concatenate([c[L]['R'][:, 4] for L in Ls])
    sr = slope_ci(xr, yr); ssp = np.polyfit([np.log2(L * L) for L in Ls], [c[L]['SP'][4] for L in Ls], 1)[0]
    above = []
    for L in big:
        m, lo, hi = tint(pm[L][:, 4]) if np.std(pm[L][:, 4]) > 0 else (pm[L][0, 4],) * 3
        rm, rlo, rhi = tint(c[L]['R'][:, 4]); above.append(lo > rhi and lo > c[L]['SP'][4])
    v['strahler_slope'] = st; v['R_strahler_slope'] = sr; v['SP_strahler_slope'] = ssp
    v['Strahler'] = bool(st[1] > sr[2] and st[1] > ssp and len(big) == 2 and all(above))
    return v


if __name__ == "__main__":
    runs = [(p, parse(p)) for p in sys.argv[1:]]
    Ls = sorted(runs[0][1]); c = controls(Ls)
    for L in Ls:
        print(f"--- L={L} N={L * L}  controls ---")
        for k, n in enumerate(NAMES):
            m, lo, hi = tint(c[L]['R'][:, k])
            print(f"  {n:10s} R-TREE {m:9.4f} [{lo:9.4f},{hi:9.4f}]  SP {c[L]['SP'][k]:9.4f}  A2 {c[L]['A2'][k]:9.4f}")
    for k, nm in [(2, "P3 median"), (4, "Strahler")]:
        xr = np.concatenate([[np.log(L * L) if k == 2 else np.log2(L * L)] * 20 for L in Ls])
        yr = np.concatenate([np.log(c[L]['R'][:, k]) if k == 2 else c[L]['R'][:, k] for L in Ls])
        xs = [np.log(L * L) if k == 2 else np.log2(L * L) for L in Ls]
        ys = [np.log(c[L]['SP'][k]) if k == 2 else c[L]['SP'][k] for L in Ls]
        ya = [np.log(c[L]['A2'][k]) if k == 2 else c[L]['A2'][k] for L in Ls]
        print(f"controls {nm} slope vs {'log N' if k == 2 else 'log2 N'}: R-TREE {tuple(np.round(slope_ci(xr, yr), 4))}, "
              f"SP {np.polyfit(xs, ys, 1)[0]:.4f}, A2 {np.polyfit(xs, ya, 1)[0]:.4f}")
    for p, pm in runs:
        print(f"=== {p} ===")
        for L in sorted(pm):
            print(f"  L={L}: tree runs {len(pm[L])}/20; " + "; ".join(
                f"{n} {tint(pm[L][:, k])[0]:.4f} [{tint(pm[L][:, k])[1]:.4f},{tint(pm[L][:, k])[2]:.4f}]" for k, n in enumerate(NAMES)))
        v = verdicts(pm, c)
        print(f"  P1 pass {v['P1']} | P2 pass {v['P2']} | rho {tuple(np.round(v['rho'], 4))} P3 pass {v['P3']} | Strahler slope "
              f"{tuple(np.round(v['strahler_slope'], 4))} vs R {tuple(np.round(v['R_strahler_slope'], 4))} SP "
              f"{v['SP_strahler_slope']:.4f}; Strahler-beyond-controls {v['Strahler']}")
