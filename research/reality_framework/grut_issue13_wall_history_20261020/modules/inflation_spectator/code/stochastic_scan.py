"""Part 2: stochastic spectator phi during an (outside-action) inflationary era, X = D = 0 along the path.

For each H_inf (GeV) on a log grid, with couplings and xi run to mu = H_inf (xi_running.py):
  * Lambda_1 (smallest phi-odd FP eigenvalue; sign-flip relaxation rate per e-fold): exact birth-death
    MFPT identity, cross-checked against the tridiagonal eigenvalue where the latter is resolvable;
  * delta start at phi = +v exactly ("fresh" start: the patch is in the +v vacuum with no fluctuations
    when the last N_* e-folds begin): p_-(n) for 0 <= n <= 60 (relative precision at any barrier);
  * relaxed start (quasi-stationary + well; long prior inflation conditioned on the + sign):
    p_- = (1 - e^{-Lambda_1 N})/2;
  * expected number of opposite-sign regions created inside the observable patch,
        N_reg = int_0^{N_*} (dp_-/dn) e^{3n} dn,
    and the number of those whose comoving size re-enters after T_X (1 MeV; T_dom where one tree wall per
    Hubble patch carries sigma/(M_P^2 H) = 1), using n'_X = ln(a_X H_X / a_0 H_0);
  * X and D entropic factors R_X = Z_X(saddle)/Z_X(well), R_D likewise (multi-field Langer estimate,
    checked against the 2D (phi,|X|) Fokker-Planck in fp2d_checks.py);
  * Hawking-Moss exponent B = 8 pi^2 DeltaV/(3H^4) and the Basu-Guth-Vilenkin wall-bubble exponent
    2 pi^2 sigma/H^3 (derived: tree wall on the equatorial S^3 of Euclidean de Sitter).
Thresholds H_* are obtained by interpolating log10(quantity) in log10 H.
"""
import time
import numpy as np
from scipy.integrate import quad
from common import Spectator, V, BETA_FP, SIGMA_TREE, PC_SITE_SC, dump, couplings
import fp1d
import cosmo

NSTARS = (50.0, 60.0)


def entropic_factors(sp):
    """R_X = Z_X(y=0)/Z_X(y=ymin), Z_X(y) = int r dr exp(-s_X), s_X = beta[c y^2 r^2/2 + lx r^4/4 + b2 r^2/2];
    R_D likewise for D with m_D^2(phi) = m_D0^2 + 12 xi_D H^2 + h_D^2 phi^2 and lam_D D^4/4."""
    p = sp.p
    ym = sp.ymin()

    def ZX(y):
        f = lambda r: r * np.exp(-BETA_FP * (0.5 * sp.c * y * y * r * r + 0.25 * sp.lx * r ** 4 + 0.5 * sp.b2 * r * r))
        rs = 3.0 / max(1e-9, np.sqrt(BETA_FP * max(sp.c * y * y + sp.b2, 1e-12)))
        rq = (40.0 / (BETA_FP * 0.25 * sp.lx)) ** 0.25
        return quad(f, 0, min(rs * 6, rq), limit=400)[0] + quad(f, min(rs * 6, rq), rq * 1.5, limit=400)[0]

    def ZD(y):
        mD2 = (p["md"] + 12 * p["xd"] * sp.H ** 2) / sp.H ** 2 + p["d"] * y * y
        f = lambda D: np.exp(-BETA_FP * (0.5 * mD2 * D * D + 0.25 * p["ld"] * D ** 4))
        Dq = (60.0 / (BETA_FP * 0.25 * p["ld"])) ** 0.25
        if mD2 > 0:
            Dq = min(Dq, 12.0 / np.sqrt(BETA_FP * mD2))
        return 2 * quad(f, 0, Dq, limit=400)[0]

    if ym == 0:
        return 1.0, 1.0
    return ZX(0.0) / ZX(ym), ZD(0.0) / ZD(ym)


def kramers(sp):
    """Overdamped Kramers rate (per e-fold, one direction) for U = w/3, D = 1/(8 pi^2); valid for B >> 1."""
    if sp.a2 >= 0:
        return np.nan
    um = 2 * abs(sp.a2) / 3.0
    ut = abs(sp.a2) / 3.0
    return np.sqrt(um * ut) / (2 * np.pi) * np.exp(-sp.barrier_B())


def n_cross_integral(ns, logp, nmax):
    """int_0^nmax (dp/dn) e^{3n} dn with piecewise-linear p; returns log10 (or -inf)."""
    m = ns <= nmax + 1e-12
    n = ns[m]
    p = np.exp(np.where(np.isfinite(logp[m]), logp[m], -np.inf))
    dp = np.diff(p)
    w = (np.exp(3 * n[1:]) - np.exp(3 * n[:-1])) / (3 * np.diff(n))
    tot = np.sum(np.clip(dp, 0, None) * w)
    return np.log10(tot) if tot > 0 else -np.inf


def log10_nreg_qsd(L1, nmax):
    # (Lambda_1/2) int_0^nmax e^{-Lambda_1 n} e^{3n} dn ~ (Lambda_1/2)(e^{3 nmax}-1)/3
    return np.log10(L1 / 2.0) + np.log10(np.expm1(3 * nmax) / 3.0)


def scan_point(H, eq_draws=False):
    sp = Spectator(H)
    yy, h = fp1d.make_grid(sp)
    y0 = V / H
    m0 = int(np.argmin(np.abs(yy - y0)))
    L1, logE = fp1d.lambda1_mfpt(sp, yy, h)
    B = sp.barrier_B()
    row = {"H_GeV": H, "v_over_H": y0, "xi_phi": sp.xp, "lam_phi": sp.lp, "B_HM": B,
           "B_BGV_wall": 2 * np.pi ** 2 * SIGMA_TREE / H ** 3, "m_phi_over_H": np.sqrt(2 * 0.01) * V / H,
           "grid_h": h, "grid_cells_half": len(yy), "y_start_cell": yy[m0], "y_min": sp.ymin(),
           "Lambda1_mfpt": L1, "log10_Lambda1": -logE / np.log(10)}
    s, lam, vec = fp1d.odd_spectrum(sp, yy, h)
    spectral_ok = lam[0] > 1e-8
    row["Lambda1_eig"] = lam[0]
    row["Lambda_odd2_eig"] = lam[1]
    row["spectral_used"] = bool(spectral_ok)
    kr = kramers(sp)
    row["kramers_rate_one_way"] = kr
    row["log10_kramers_x2_over_Lambda1"] = np.log10(2 * kr / L1) if np.isfinite(kr) and kr > 0 else None
    if B >= 3.0:
        RX, RD = entropic_factors(sp)
    else:
        RX, RD = None, None          # no metastable well: Langer/entropic estimate not applicable
    row["R_X_entropic"] = RX
    row["R_D_entropic"] = RD
    ns, logp, floor = fp1d.p_minus_delta_curve(sp, yy, h, m0)
    if y0 < 2 * h:
        # start inside the first cell: use oddness, S(y) ~ S'(0) y, p_- = 1/2 - (1/2 - p(h/2)) y0/(h/2)
        ns1, logp1, _ = fp1d.p_minus_delta_curve(sp, yy, h, 0)
        p1 = np.exp(logp1)
        logp = np.log(0.5 - (0.5 - p1) * y0 / (h / 2))
        row["delta_start_note"] = "y0 = v/H below grid resolution: linear (odd) interpolation between y=0 and the first cell"
    # relaxed start: equilibrium restricted to phi > 0 (exact spectral when resolvable, else quasi-stationary)
    A = vec.T @ np.exp(-(s - s.min()) / 2)
    A2 = A ** 2 / np.sum(A ** 2)
    if spectral_ok:
        logp_rel = np.log(0.5 * np.array([np.sum(A2 * (-np.expm1(-lam * n))) for n in ns]) + 1e-300)
    else:
        logp_rel = np.log(0.5 * (-np.expm1(-L1 * ns)) + 1e-300)
    nX = {"1MeV": cosmo.efolds_obs_to_reentry(1e-3),
          "Tdom": cosmo.efolds_obs_to_reentry(cosmo.T_wall_dominates_hubble_patch(SIGMA_TREE))}
    for N in NSTARS:
        tag = "N%d" % N
        i = int(round(N / 0.25))
        lp_delta = logp[i] / np.log(10)
        lp_rel = logp_rel[i] / np.log(10)
        row[tag] = {
            "delta_start": {"log10_p_minus": lp_delta, "p_minus": 10 ** lp_delta,
                            "log10_Nreg_all_sizes": n_cross_integral(ns, logp, N),
                            "log10_Nreg_reenter_after_1MeV": n_cross_integral(ns, logp, nX["1MeV"]),
                            "log10_Nreg_reenter_after_Tdom": n_cross_integral(ns, logp, nX["Tdom"]),
                            "solver_floor_at_nX_1MeV": bool(floor[int(nX["1MeV"] / 0.25)])},
            "relaxed_start": {"p_minus": 10 ** lp_rel, "log10_p_minus": lp_rel,
                              "log10_Nreg_all_sizes": (n_cross_integral(ns, logp_rel, N) if spectral_ok
                                                       else log10_nreg_qsd(L1, N)),
                              "log10_Nreg_reenter_after_1MeV": (n_cross_integral(ns, logp_rel, nX["1MeV"]) if spectral_ok
                                                                else log10_nreg_qsd(L1, nX["1MeV"])),
                              "log10_Nreg_reenter_after_Tdom": (n_cross_integral(ns, logp_rel, nX["Tdom"]) if spectral_ok
                                                                else log10_nreg_qsd(L1, nX["Tdom"]))},
        }
        # region counting assumes rare, independent crossings; it is meaningless once p_- is O(1)
        for key in ("delta_start", "relaxed_start"):
            if row[tag][key]["p_minus"] > 0.05:
                for f in ("log10_Nreg_all_sizes", "log10_Nreg_reenter_after_1MeV", "log10_Nreg_reenter_after_Tdom"):
                    row[tag][key][f] = None
        if RX is not None:
            Lx = L1 * RX * RD
            row[tag]["relaxed_start_with_XD_entropic"] = {
                "log10_p_minus": np.log10(0.5 * (-np.expm1(-Lx * N))),
                "log10_Nreg_all_sizes": log10_nreg_qsd(Lx, N),
                "log10_Nreg_reenter_after_1MeV": log10_nreg_qsd(Lx, nX["1MeV"])}
        else:
            row[tag]["relaxed_start_with_XD_entropic"] = None
        if eq_draws and spectral_ok:
            # minority fraction for every start z on the + side: F(z) = sum_n phi_n(z) e^{s_z/2} A_n (1-e^{-lam N})
            Fz = (vec * (A * (-np.expm1(-lam * N)))[None, :]).sum(axis=1) * np.exp((s - s.min()) / 2)
            mfrac = np.clip(0.5 * Fz, 0.0, 0.5)
            pi = np.exp(-(s - s.min()))
            pi /= pi.sum()
            order = np.argsort(mfrac)
            cdf = np.cumsum(pi[order])
            q = lambda x: float(mfrac[order][min(np.searchsorted(cdf, x), len(cdf) - 1)])
            row[tag]["equilibrium_draws"] = {
                "mean_p_minus_over_full_equilibrium": 0.5,
                "mean_minority_fraction": float(np.sum(pi * mfrac)),
                "minority_fraction_quantiles_10_50_90": [q(0.1), q(0.5), q(0.9)],
                "prob_minority_below_pc": float(np.sum(pi[mfrac < PC_SITE_SC])),
                "prob_minority_below_1e-3": float(np.sum(pi[mfrac < 1e-3]))}
    if eq_draws:
        row["equilibrium"] = fp1d.equilibrium_stats(sp, yy, h)
    return row


def threshold(Hs, vals, target):
    """Largest H below which log10(quantity) < log10(target): interpolate on the increasing branch."""
    Hs = np.asarray(Hs)
    vals = np.asarray(vals, dtype=float)
    lt = np.log10(target)
    ok = np.isfinite(vals)
    x = np.log10(Hs[ok])
    y = vals[ok]
    idx = np.where((y[:-1] < lt) & (y[1:] >= lt))[0]
    if idx.size == 0:
        return None
    i = idx[0]
    return float(10 ** (x[i] + (lt - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i])))


def main():
    t0 = time.time()
    Hs = np.unique(np.concatenate([np.logspace(2.0, 3.0, 41), np.logspace(3.05, 4.0, 12),
                                   np.logspace(4.5, 13.0, 18)]))
    rows = []
    for H in Hs:
        rows.append(scan_point(H, eq_draws=True))
    out = {"module": "inflation_spectator", "part": "2 stochastic spectator (1D, X = D = 0 along the path)",
           "labels": "computed (SqRA Fokker-Planck; Starobinsky-Yokoyama astro-ph/9407016); thresholds interpolated; "
                     "inflation, H_inf, N_* and the start condition are OUTSIDE-ACTION preparation data",
           "validation_pure_quartic": None, "rows": rows}
    # validation vs literature for V = lam phi^4/4
    sp = Spectator(1e13)
    sp.a2 = 0.0
    sp.lp = 0.01
    yy, h = fp1d.make_grid(sp)
    s, lam, vec = fp1d.odd_spectrum(sp, yy, h)
    s_, up, dn, kill, k = fp1d.odd_chain(sp, yy, h)
    from scipy.linalg import eigh_tridiagonal
    lame, _ = eigh_tridiagonal(up + dn, -k * np.ones(len(s_) - 1))
    out["validation_pure_quartic"] = {
        "Lambda1_over_sqrt_lam": lam[0] / 0.1, "literature_SY_quoted_from_memory": 0.08892,
        "two_Lambda2_over_sqrt_lam": 2 * lame[1] / 0.1, "literature_MRST_1904.11917_abstract": 0.579}
    # convergence: halve h at two points
    conv = []
    for H in (250.0, 150.0, 1e13):
        sp = Spectator(H)
        yy, h = fp1d.make_grid(sp)
        yy2, h2 = fp1d.make_grid(sp, h=h / 2)
        a, _ = fp1d.lambda1_mfpt(sp, yy, h)
        b, _ = fp1d.lambda1_mfpt(sp, yy2, h2)
        m0 = int(np.argmin(np.abs(yy - V / H)))
        m02 = int(np.argmin(np.abs(yy2 - V / H)))
        pa = fp1d.p_minus_delta(sp, yy, h, m0, [50.0])[50.0]["log10_p_minus"]
        pb = fp1d.p_minus_delta(sp, yy2, h2, m02, [50.0])[50.0]["log10_p_minus"]
        conv.append({"H": H, "Lambda1_ratio_h_over_h2": a / b, "log10_p50_delta_h_minus_h2": pa - pb})
    out["grid_convergence"] = conv
    # thresholds
    thr = {}
    for N in NSTARS:
        tag = "N%d" % N
        d = {}
        for start in ("delta_start", "relaxed_start", "relaxed_start_with_XD_entropic"):
            get = lambda r, key: (r[tag][start][key] if (r[tag][start] is not None and r[tag][start][key] is not None)
                                  else np.nan)
            lp_ = [get(r, "log10_p_minus") for r in rows]
            nr = [get(r, "log10_Nreg_all_sizes") for r in rows]
            n1 = [get(r, "log10_Nreg_reenter_after_1MeV") for r in rows]
            e = {"H_star_p_minus_1e-3": threshold(Hs, lp_, 1e-3),
                 "H_star_p_minus_1e-10": threshold(Hs, lp_, 1e-10),
                 "H_perc_p_minus_eq_pc": threshold(Hs, lp_, PC_SITE_SC),
                 "H_star_Nreg_all_1e-3": threshold(Hs, nr, 1e-3),
                 "H_star_Nreg_all_1": threshold(Hs, nr, 1.0),
                 "H_star_Nreg_after_1MeV_1e-3": threshold(Hs, n1, 1e-3),
                 "H_star_Nreg_after_1MeV_1": threshold(Hs, n1, 1.0)}
            if start != "relaxed_start_with_XD_entropic":
                nd = [get(r, "log10_Nreg_reenter_after_Tdom") for r in rows]
                e["H_star_Nreg_after_Tdom_1e-3"] = threshold(Hs, nd, 1e-3)
            d[start] = e
        thr[tag] = d
    out["thresholds_GeV"] = thr
    out["reentry_efolds"] = {"n_obs_to_1MeV": cosmo.efolds_obs_to_reentry(1e-3),
                             "T_dom_GeV": cosmo.T_wall_dominates_hubble_patch(SIGMA_TREE),
                             "n_obs_to_Tdom": cosmo.efolds_obs_to_reentry(cosmo.T_wall_dominates_hubble_patch(SIGMA_TREE))}
    out["Nstar_instant_reheating"] = {str(H): cosmo.Nstar_instant(H) for H in (150.0, 300.0, 1e3, 1e8, 1e13)}
    out["runtime_s"] = time.time() - t0
    dump("stochastic_1d.json", out)
    return out


if __name__ == "__main__":
    o = main()
    print("runtime", o["runtime_s"])
    print(o["validation_pure_quartic"])
    print(o["grid_convergence"])
    import json
    print(json.dumps(o["thresholds_GeV"], indent=1))
    for r in o["rows"]:
        ed = r["N50"].get("equilibrium_draws", {})
        print("%10.4g B=%8.3f L1=%.3e RX=%s p50d=%.3e p50r=%.3e eqq=%s" % (
            r["H_GeV"], r["B_HM"], r["Lambda1_mfpt"], r["R_X_entropic"], r["N50"]["delta_start"]["p_minus"],
            r["N50"]["relaxed_start"]["p_minus"], ed.get("minority_fraction_quantiles_10_50_90")))
