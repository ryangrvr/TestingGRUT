"""CARD #1 v1R2 shared logic: result pools and the owner-ruled LOWER-MINIMUM REPRODUCIBILITY criterion.

A target is NUMERICALLY REPRODUCED iff
  (1) at least two genuinely distinct, cycle-converged starts reach chi2 within TOL (0.2) of each other, and
  (2) those agreeing minima constitute the lowest basin found, and
  (3) no executed start (any arm, any campaign v1R/v1R2, converged or not) is more than TOL lower.
Operationally: L = min over ALL executed starts; reproduced iff #(cycle-converged starts with chi2 <= L + TOL) >= 2.
Higher local minima never cause failure by themselves; they are preserved in the record.
"""
import glob, json, os

TOL = 0.2


def key_of(r):
    if r["model"] == "cpl":
        return "cpl"
    if r["model"] == "free":
        return "free"
    return f"card:{round(float(r['eps_target']), 3)}"


def pool(*dirs):
    P = {}
    for d in dirs:
        for f in glob.glob(os.path.join(d, "*.json")):
            if f.endswith(".tmp.json"):
                continue
            r = json.load(open(f))
            if "chi2_eff" not in r or "model" not in r:
                continue
            r["_file"] = os.path.basename(f)
            P.setdefault(key_of(r), []).append(r)
    return P


def reproduced(rows):
    if not rows:
        return False, None, "no starts"
    L = min(r["chi2_eff"] for r in rows)
    agree = [r for r in rows if r.get("converged_cycles") and r["chi2_eff"] <= L + TOL]
    ok = len(agree) >= 2
    return ok, L, f"lowest {L:.3f}; {len(agree)} converged start(s) within {TOL} of lowest; {len(rows)} executed"


def best(rows):
    return min(rows, key=lambda r: r["chi2_eff"])
