#!/usr/bin/env python3
"""SCOUT_0 W1 P-06: the CM/non-CM boundary of the S5-1 memory kernel.

Question (charter P-06, frozen): is the completely-monotone (CM) boundary of the
S5-1 conservative-parent memory kernel a stable surface in a controlled
deformation family, or is the S5-1 parent CM-isolated?

Objects (all from the canonical S5-1 record, not invented here):
  - scaled-limit spectral density: Lorentzian  S(x) = kappa/(x^2 + kappa^2),
    kappa = 1/(2*sqrt(2.3))  [S5_OWNER_RULING_03.md:49]
  - controlled deformation family: Lorentzian exponent p,
    S_p(x) = (x^2 + kappa^2)^(-p)   (p=1 recovers the recorded Lorentzian)
  - inverse Fourier transform:  f_p(t) ~ t^(p-1/2) K_(p-1/2)(kappa t)
    (Gradshteyn-Ryzhik 3.723.4; p=1 gives the pure exponential e^{-kappa t},
    which is the known CM anchor).

CM test (Bernstein): f is completely monotone on t>0 iff
  (-1)^n f^(n)(t) >= 0 for all n >= 0, t > 0.
We test n = 0..N_MAX on a log grid with mpmath high-precision derivatives.

Self-tests:
  S-1: p=1 reduces to e^{-kappa t} (up to normalization) — must PASS CM.
  S-2: mutation check — a known NON-CM function (e^{-t} cos t, whose derivative
       sign pattern violates alternation) must FAIL the tester.
  S-3: p=1/2 gives K_0(kappa t), which has the exact Laplace representation
       K_0(t) = int_0^inf e^{-t cosh s} ds — positive measure, must PASS.

Mutation policy: S-2 failing means the tester is live. If S-2 passes (tester
blind), the whole probe result is VOID.
"""
import mpmath as mp
import json, sys

mp.mp.dps = 22
KAPPA = 1 / (2 * mp.sqrt(mp.mpf("2.3")))
N_MAX = 6          # derivative order tested
T_GRID = [mp.mpf("0.08") * mp.mpf("1.8")**i for i in range(12)]  # 0.08 .. ~12

def f_p(t, p):
    """Unnormalized kernel of the deformed family (normalization is positive
    and irrelevant to CM)."""
    t = mp.mpf(t); p = mp.mpf(p)
    nu = p - mp.mpf("0.5")
    return t**nu * mp.besselk(nu, KAPPA * t)

def cm_test(f, t_grid, n_max=N_MAX, tol="1e-25"):
    """Return (ok, worst evidence). ok=True iff (-1)^n f^(n)(t) >= -tol
    for all tested n, t."""
    worst_n, worst_t, worst_val = None, None, mp.mpf("inf")
    for n in range(n_max + 1):
        for t in t_grid:
            d = mp.diff(f, t, n)
            val = (-1)**n * d
            if val < mp.mpf(tol) and val < worst_val:
                worst_n, worst_t, worst_val = n, t, val
    return (worst_n is None), (None if worst_n is None else
                               {"n": worst_n, "t": float(worst_t), "val": float(worst_val)})

def main():
    out = {"probe": "P-06", "kappa": float(KAPPA), "n_max": N_MAX,
           "self_tests": {}, "family_sweep": [], "verdict": None}

    # S-1: p=1 is the exponential
    f1 = lambda t: f_p(t, 1)
    ok, ev = cm_test(f1, T_GRID)
    out["self_tests"]["S1_p1_is_exponential"] = {"ok": ok, "evidence": ev}
    # numeric identity check
    dev = max(abs(f1(t) / f1(mp.mpf("0.1")) - mp.e**(-KAPPA * (t - mp.mpf("0.1"))))
              for t in T_GRID)
    out["self_tests"]["S1_exponential_identity_dev"] = float(dev)

    # S-2: mutation check with a known non-CM function
    bad = lambda t: mp.e**(-t) * mp.cos(t)
    ok_bad, ev_bad = cm_test(bad, T_GRID[:12], n_max=4)
    out["self_tests"]["S2_mutation_nonCM_detected"] = {"tester_live": (not ok_bad),
                                                       "evidence": ev_bad}
    if ok_bad:
        out["verdict"] = "VOID — tester failed the mutation check"
        print(json.dumps(out, indent=1)); return 1

    # S-3: p=1/2 is K_0
    f05 = lambda t: f_p(t, mp.mpf("0.5"))
    ok, ev = cm_test(f05, T_GRID)
    out["self_tests"]["S3_p05_K0"] = {"ok": ok, "evidence": ev}

    # family sweep
    ps = [mp.mpf("0.05"), "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.8",
          "1.0", "1.2", "1.5", "2.0"]
    for p in ps:
        fp = lambda t, p=p: f_p(t, p)
        ok, ev = cm_test(fp, T_GRID)
        out["family_sweep"].append({"p": float(p), "CM": ok, "evidence": ev})

    # verdict: locate the boundary from the sweep
    cm_flags = [r["CM"] for r in out["family_sweep"]]
    if all(cm_flags):
        out["verdict"] = "CM-OPEN: whole tested family completely monotone; " \
                         "S5-1 parent not CM-isolated on this family"
    elif not any(cm_flags):
        out["verdict"] = "CM-ISOLATED-INDICATED: no tested p is CM"
    else:
        first_bad = next(r["p"] for r in out["family_sweep"] if not r["CM"])
        last_ok = max(r["p"] for r in out["family_sweep"]
                      if r["CM"] and r["p"] < first_bad) \
                  if any(r["CM"] and r["p"] < first_bad for r in out["family_sweep"]) else None
        out["verdict"] = {"boundary_bracketed": [last_ok, first_bad]}

    print(json.dumps(out, indent=1))
    return 0

if __name__ == "__main__":
    sys.exit(main())
