#!/usr/bin/env python3
"""CARD-01 v1 D1/D2 diagnostics — DIAGNOSTIC — NOT THE KILL TEST.
Uses ONLY the official DESI DR2 base_w_wa chains (checksum-verified against DESI's published sha256sum).
D1: posterior geometry in (w0, wa): means, covariance, 68/95% HPD levels (S8), quadrant mass P(w0>-1, wa<0).
D2: overlay of the frozen Card #1 CPL line wa = 1.55590467*(1+w0) (registered unweighted LSQ projection):
    primary ray eps<0 -> w0<-1; control ray eps>0 -> w0>-1 (CONTROL - NOT A GRUT CLAIM).
    Reported: the smallest HPD credible level whose region the ray reaches (density-ranked, KDE on samples).
usage: card01_d1d2_diagnostic.py CHAINS_ROOT OUT_JSON
"""
import json, os, sys
import numpy as np
from getdist import loadMCSamples
from scipy.stats import gaussian_kde

SLOPE = 1.55590467046443
C0 = 0.433262795583067
COMBOS = {
    "PRIMARY_DESI_CMB": "desi-bao-all_planck2018-lowl-TT-clik_planck2018-lowl-EE-clik_planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing",
    "EXT_PANTHEONPLUS": "desi-bao-all_pantheonplus_planck2018-lowl-TT-clik_planck2018-lowl-EE-clik_planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing",
    "EXT_UNION3": "desi-bao-all_union3_planck2018-lowl-TT-clik_planck2018-lowl-EE-clik_planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing",
    "EXT_DESY5": "desi-bao-all_desy5sn_planck2018-lowl-TT-clik_planck2018-lowl-EE-clik_planck-NPIPE-highl-CamSpec-TTTEEE_planck-act-dr6-lensing",
}


def hpd_level_of(dens_query, dens_samples, w):
    """credible level of the HPD region whose boundary passes through density dens_query."""
    return float(np.sum(w[dens_samples >= dens_query]) / np.sum(w))


def main():
    root, out = sys.argv[1], sys.argv[2]
    res = {"label": "DIAGNOSTIC - NOT THE KILL TEST", "ignore_rows": 0.3}
    rng = np.random.default_rng(1)
    for key, d in COMBOS.items():
        s = loadMCSamples(os.path.join(root, "cobaya/base_w_wa", d, "chain"), settings={"ignore_rows": 0.3})
        w0 = s.getParams().w; wa = s.getParams().wa; wt = s.weights
        m = np.array([np.average(w0, weights=wt), np.average(wa, weights=wt)])
        cov = np.cov(np.vstack([w0, wa]), aweights=wt)
        quad = float(np.sum(wt[(w0 > -1) & (wa < 0)]) / np.sum(wt))
        # KDE on a weighted subsample for density ranking
        n = min(20000, len(w0)); idx = rng.choice(len(w0), size=n, replace=True, p=wt / wt.sum())
        kde = gaussian_kde(np.vstack([w0[idx], wa[idx]]))
        dens_s = kde(np.vstack([w0[idx], wa[idx]])); ws = np.ones(n)
        lev68 = np.quantile(dens_s, 1 - 0.68); lev95 = np.quantile(dens_s, 1 - 0.95)
        eps = np.linspace(-1, 0, 401)
        ray_p = np.vstack([-1 + C0 * eps, SLOPE * C0 * eps]); dp = kde(ray_p)
        ray_c = np.vstack([-1 + C0 * (-eps), SLOPE * C0 * (-eps)]); dc = kde(ray_c)
        ip, ic = int(np.argmax(dp)), int(np.argmax(dc))
        res[key] = {
            "n_samples": int(len(w0)), "mean_w0_wa": m.tolist(), "cov_w0_wa": cov.tolist(),
            "P_quadrant_w0>-1_wa<0": quad,
            "primary_ray_eps<=0": {"max_density_eps": float(eps[ip]),
                                   "point": ray_p[:, ip].tolist(),
                                   "smallest_HPD_level_reached": hpd_level_of(dp[ip], dens_s, ws),
                                   "inside_95": bool(dp[ip] >= lev95), "inside_68": bool(dp[ip] >= lev68)},
            "control_ray_eps>=0 (CONTROL - NOT A GRUT CLAIM)": {
                "max_density_eps": float(-eps[ic]), "point": ray_c[:, ic].tolist(),
                "smallest_HPD_level_reached": hpd_level_of(dc[ic], dens_s, ws),
                "inside_95": bool(dc[ic] >= lev95), "inside_68": bool(dc[ic] >= lev68)},
            "LCDM_point_HPD_level": hpd_level_of(float(kde(np.array([[-1.0], [0.0]]))[0]), dens_s, ws),
        }
        print(key, json.dumps(res[key]["primary_ray_eps<=0"]), "quadrant", round(quad, 4),
              "LCDM level", round(res[key]["LCDM_point_HPD_level"], 4))
    json.dump(res, open(out, "w"), indent=1)


if __name__ == "__main__":
    main()
