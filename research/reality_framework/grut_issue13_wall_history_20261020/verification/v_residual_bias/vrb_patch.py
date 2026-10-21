"""Verifier v_residual_bias, part 3b: route B patch realization spread (module section 4.4), spectral method.
rms of m(y*, N) over equilibrium draws y* ~ P_eq, and the sd of m over sub-patches one e-fold after a fresh start.
F is represented by exact cell averages around each node (thresholds from vrb_trapping)."""
import json
import math
import os

import numpy as np

import vrb_trapping as T


def cell_avg_F(y, h, th):
    def G(x):          # Int_0^x s(a) da for x >= 0
        pts = np.concatenate([[0.0], th])
        out = np.zeros_like(x)
        for k in range(len(pts)):
            r = pts[k + 1] if k + 1 < len(pts) else np.inf
            out += (-1) ** k * np.clip(np.minimum(x, r) - pts[k], 0, None)
        return out
    I = lambda x: G(np.abs(x))       # Int_0^x F (F odd -> even primitive)
    return (I(y + h / 2) - I(y - h / 2)) / h


def m_of_y0(fp, F, N):
    """m(y0, N) on all nodes: e^{s(y0)/2} sum_n u_n(y0) e^{-Lam N} (u_n . e^{-s/2} F)."""
    c = (fp.U.T @ (np.exp(-fp.s / 2) * F)) * np.exp(-fp.Lam * N)
    return np.exp(fp.s / 2) * (fp.U @ c)


def main():
    out = []
    for H in (5000.0, 1e4, 2e4):
        fp = T.FP(H, h_target=0.002, Y=6.0)
        th, _ = T.thresholds(fp.u, 6.0)
        F = cell_avg_F(fp.y, fp.h, th)
        peq = np.exp(-(fp.s - fp.s.min()))
        peq /= peq.sum()
        r = dict(H=H)
        for N in (50.0,):
            m = m_of_y0(fp, F, N)
            r["m_fresh_check"] = float(m[fp.i0])
            r["eq_rms"] = float(np.sqrt((peq * m * m).sum()))
            r["eq_mean"] = float((peq * m).sum())
            m1 = m_of_y0(fp, F, N - 1)
            P1 = fp.P_fresh(1.0) * fp.h
            P1 = np.clip(P1, 0, None)
            P1 /= P1.sum()
            mu = float((P1 * m1).sum())
            r["subpatch_mean"] = mu
            r["subpatch_sd"] = float(np.sqrt(max((P1 * m1 * m1).sum() - mu * mu, 0)))
        out.append(r)
        print(r, flush=True)
    with open(os.path.join(T.OUT, "patch.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    main()
