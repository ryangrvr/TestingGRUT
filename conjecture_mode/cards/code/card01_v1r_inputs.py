#!/usr/bin/env python3
"""CARD-01 v1R numerical inputs from the checksum-verified OFFICIAL DESI DR2 chains (numerical aids only).

Produces (in OUTDIR):
  v1r_covmat_lcdm.txt   proposal covariance for Card/LCDM targets, from cobaya/base/<primary>
  v1r_covmat_w0wa.txt   proposal covariance for the CPL comparator, from cobaya/base_w_wa/<primary>
  v1r_cpl_startB.json   CPL alternate start: the highest-likelihood (lowest minuslogpost) sample of the official
                        base_w_wa primary chains.  Owner ruling v1R s.5: a starting location only, never a fit result.
Covariances use the post-burn-in (0.3) weighted samples, restricted to the parameters sampled in v1R.
usage: card01_v1r_inputs.py CHAINS_ROOT OUTDIR
"""
import json, os, sys
import numpy as np
from getdist import loadMCSamples

PRIMARY = "desi-bao-all_planck2018-lowl-TT-clik_planck2018-lowl-EE-clik_planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing"
P_LCDM = ["logA", "ns", "H0", "ombh2", "omch2", "tau", "A_planck", "amp_143", "amp_217", "amp_143x217",
          "n_143", "n_217", "n_143x217", "calTE", "calEE"]
P_W0WA = P_LCDM + ["w", "wa"]


def covmat(root, model, params, out):
    s = loadMCSamples(os.path.join(root, "cobaya", model, PRIMARY, "chain"), settings={"ignore_rows": 0.3})
    X = np.vstack([s.getParams().__dict__[p] for p in params])
    C = np.cov(X, aweights=s.weights)
    with open(out, "w") as f:
        f.write("# " + " ".join(params) + "\n")
        np.savetxt(f, C)
    return s


def main():
    root, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    covmat(root, "base", P_LCDM, os.path.join(outdir, "v1r_covmat_lcdm.txt"))
    covmat(root, "base_w_wa", P_W0WA, os.path.join(outdir, "v1r_covmat_w0wa.txt"))
    best = None
    for k in range(1, 5):
        fn = os.path.join(root, "cobaya", "base_w_wa", PRIMARY, f"chain.{k}.txt")
        hdr = open(fn).readline().lstrip("#").split()
        A = np.loadtxt(fn)
        i = int(np.argmin(A[:, hdr.index("minuslogpost")]))
        if best is None or A[i, hdr.index("minuslogpost")] < best[0]:
            best = (A[i, hdr.index("minuslogpost")], {p: float(A[i, hdr.index(p)]) for p in P_W0WA}, fn, i)
    json.dump({"source": os.path.relpath(best[2], root), "row": best[3],
               "official_minuslogpost_in_DESI_conventions": best[0], "start": best[1],
               "note": "numerical starting location only (owner ruling v1R s.5); not a fit result"},
              open(os.path.join(outdir, "v1r_cpl_startB.json"), "w"), indent=1)
    print("CPL alt start:", best[1])


if __name__ == "__main__":
    main()
