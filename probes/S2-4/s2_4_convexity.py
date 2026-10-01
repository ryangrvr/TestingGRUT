"""SCOUT-2 S2-4: convexity / mixtures without supplied randomness?

Operational setup (deterministic throughout):
  run n: a hidden coin c_n in {0,1} selects a deterministic preparation P_{c_n} (a microstate);
  the tester chooses a setting s_n in {0,1}; the outcome is deterministic: o = O[c_n][s_n].
  'State' of the procedure = the run-frequencies of outcomes per setting (equivalence class over tests).
Coins:
  (a) endogenous, uniquely ergodic: c_n = [n*phi mod 1 < lam]  (Sturmian; no measure needed for its frequency)
  (b) the same coin with an adaptive tester who can see the run index.
Tester settings:
  (i)   s_n = [n*sqrt2 mod 1 < mu]         rationally independent of phi  -> joint equidistribution (Weyl)
  (ii)  s_n = [2 n*phi mod 1 < mu]         rationally DEPENDENT on the coin's rotation
  (iii) s_n = c_n                          tester reads the coin (access to the hidden variable)
Question: is the observed statistic the convex mixture lam*O[1][s] + (1-lam)*O[0][s]?
"""
import numpy as np

phi = (np.sqrt(5) - 1) / 2
O = np.array([[0, 1],     # preparation 0: outcome 0 under setting 0, outcome 1 under setting 1
              [1, 1]])    # preparation 1: outcome 1 under both settings
N = 2_000_000
n = np.arange(N)

for lam in (0.3, phi / 2):
    c = ((n * phi) % 1.0 < lam).astype(int)
    print(f"=== hidden uniquely ergodic coin, lam = {lam:.6f}: coin frequency = {c.mean():.6f} ===")
    mu = 0.4
    testers = {
        "(i)   s = [n sqrt2 mod 1 < mu]  (rationally independent)": ((n * np.sqrt(2)) % 1.0 < mu).astype(int),
        "(ii)  s = [2 n phi mod 1 < mu]  (rationally dependent)": ((2 * n * phi) % 1.0 < mu).astype(int),
        "(iii) s = c                     (tester reads the coin)": c.copy(),
    }
    for name, s in testers.items():
        joint = np.mean((c == 1) & (s == 1)); prod = c.mean() * s.mean()
        out = O[c, s]
        line = []
        for setting in (0, 1):
            sel = s == setting
            if sel.sum() == 0:
                line.append("n/a"); continue
            f_obs = out[sel].mean()
            f_mix = lam * O[1][setting] + (1 - lam) * O[0][setting]
            line.append(f"setting {setting}: observed {f_obs:.4f} vs convex mixture {f_mix:.4f}")
        print(f"  {name}: P(c=1,s=1) = {joint:.4f} vs product {prod:.4f}")
        print("       " + "; ".join(line))
print("\n-> mixtures (convexity) emerge from a deterministic hidden coin with EARNED weights, iff the tester's settings are")
print("   statistically independent of the coin: earned by rational independence of the rotation numbers (Weyl), lost")
print("   for rationally dependent or coin-reading testers. Prices: coin inaccessibility (access) + independence condition.")
