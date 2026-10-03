"""BRI1-PF4Q frozen-tau run: 20 K_q[a,b,c] (P1, P2 x 10 unique components) by the amended pre-declared evidence method
(BRI1_PF4Q.md sec. 3 + 3.1).  NUMERICAL EVIDENCE ONLY -- not certified.  No Monte Carlo."""
import json, numpy as np, pf4q_core as C
res = {}
for h in (0.08, 0.06):
    X, P, W = C.rule_trap(h); sol = C.variational(X, P, list(C.TAU))
    for pr in (1, 2): res[f"A_h{h}_P{pr}"] = C.K_from_variational(W, sol, pr)
X, P, W, _ = C.rule_gh(192); sol = C.variational(X, P, list(C.TAU))
for pr in (1, 2): res[f"GH192_P{pr}"] = C.K_from_variational(W, sol, pr)
for h in (0.12, 0.08):
    X, P, W = C.rule_trap(h)
    for pr in (1, 2):
        k1 = C.K_from_fd(W, X, P, pr, 1e-3); k2 = C.K_from_fd(W, X, P, pr, 5e-4)
        res[f"B_h{h}_P{pr}"] = {k: (4 * k2[k] - k1[k]) / 3 for k in k1}
names = ["pi", "3pi/2", "2pi"]
summary = {}
print("component        proto   K_A(h=.06)        K_A(h=.08)        K_B_R(h=.08)      K_B_R(h=.12)      GH192        "
      "  criterion  sign")
for pr in (1, 2):
    for k in C.IDX:
        a6, a8 = res[f"A_h0.06_P{pr}"][k], res[f"A_h0.08_P{pr}"][k]
        b8, b12 = res[f"B_h0.08_P{pr}"][k], res[f"B_h0.12_P{pr}"][k]; g = res[f"GH192_P{pr}"][k]
        spread = max(abs(a6 - a8), abs(a6 - b8), abs(b8 - b12))
        ok = abs(a6 - b8) < 1e-6 and abs(a6) > 100 * spread
        lab = "(" + ",".join(names[i] for i in k) + ")"
        summary[f"P{pr}{lab}"] = dict(K_A06=a6, K_A08=a8, K_B08=b8, K_B12=b12, GH192=g, spread=spread, evidence_nonzero=bool(ok))
        print(f"{lab:16s} P{pr}  {a6: .10e} {a8: .10e} {b8: .10e} {b12: .10e} {g: .3e}  {str(ok):5s}  {'+' if a6 > 0 else '-'}")
json.dump(summary, open("pf4q_results.json", "w"), indent=1)
n = sum(v["evidence_nonzero"] for v in summary.values())
print(f"\ncomponents meeting the pre-declared evidence criterion: {n} / 20  (EVIDENCE ONLY -- NOT CERTIFIED)")
