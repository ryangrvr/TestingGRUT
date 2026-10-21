"""Aggregate results/runs/*.json into results/summary.json (per-run metrics, comparisons, conservation,
resolution check, sign statistics, physical-scale numbers).  Missing runs are skipped."""
import glob
import json
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "..", "results")


def load(name):
    p = os.path.join(RES, "runs", name + ".json")
    if not os.path.exists(p):
        return None
    with open(p) as f:
        return json.load(f)


def arr(series, key, default=np.nan):
    return np.array([(s.get(key) if s.get(key) is not None else default) for s in series], dtype=float)


def interp_at(x, xs, ys):
    return float(np.interp(x, xs, ys))


def metrics(run):
    cfg = run["config"]
    s = run["series"]
    tau = arr(s, "tau"); a = arr(s, "a")
    fp = arr(s, "f_plus"); fm = arr(s, "f_minus"); fc = arr(s, "f_core")
    Lw = arr(s, "L_wall"); Lb = arr(s, "L_bare"); Lall = arr(s, "L_all"); Lint = arr(s, "L_internal")
    Lv = arr(s, "L_vor")                       # vacuum-basin interface length (NaN for runs without it)
    has_vor = bool(np.isfinite(Lv).any())
    d = cfg["d"]
    Lbox = cfg["N"] * cfg["dx"]
    out = dict(name=cfg["name"], kind=cfg["kind"], kappa=cfg.get("kappa"), tau_i=cfg["tau_i"], seed=cfg["seed"],
               d=d, N=cfg["N"], dx=cfg["dx"], L_box_comoving=Lbox, a_final=float(a[-1]), tau_final=float(tau[-1]),
               a_freeze=cfg.get("a_freeze"), wall_s=run.get("wall_s"))
    out["init"] = run.get("info")
    vac = fp + fm
    iform = np.nonzero(vac >= 0.5)[0]
    if cfg["kind"].startswith("memory"):
        # sign-memory controls: formation only counts after a restoration (core volume >= 1/2)
        irest = np.nonzero(fc >= 0.5)[0]
        out["restored"] = bool(irest.size)
        out["max_f_core"] = float(fc.max())
        out["max_f_minus"] = float(fm.max())
        out["min_f_plus"] = float(fp.min())
        out["sign_changes_ever"] = bool((Lall > 0).any())
        if irest.size:
            out["a_restored"] = float(a[irest[0]])
            mS = arr(s, "mean_S")
            out["mean_S_min"] = float(mS.min()); out["mean_S_max_after_restore"] = float(mS[irest[0]:].max())
            # two-sign formation: first output after restoration with both vacua >= 5% of the volume
            i2 = np.nonzero((np.minimum(fp, fm) >= 0.05) & (np.arange(len(fp)) > irest[0]))[0]
            out["a_two_sign"] = float(a[i2[0]]) if i2.size else None
            iform = iform[iform > (i2[0] if i2.size else len(fp))]
            # first output after restoration at which the -v phase exceeds 10% of the volume
            im = np.nonzero((fm >= 0.1) & (np.arange(len(fm)) > irest[0]))[0]
            out["a_first_fminus_gt_0.1"] = float(a[im[0]]) if im.size else None
        else:
            iform = iform[:0]
    if iform.size:
        i0 = int(iform[0])
        out["a_form"] = float(a[i0]); out["tau_form"] = float(tau[i0])
        out["p_plus_form"] = float(fp[i0] / (fp[i0] + fm[i0]))
        # first time both vacua have >= 10% volume (or neither)
        post = slice(i0, None)
        both = np.array([s[k]["wraps_plus"] and s[k]["wraps_minus"] for k in range(i0, len(s))])
        anyw = np.array([s[k]["wraps_plus"] or s[k]["wraps_minus"] for k in range(i0, len(s))])
        out["frac_outputs_both_signs_wrap"] = float(both.mean())
        out["frac_outputs_some_sign_wraps"] = float(anyw.mean())
        out["walls_at_every_output_after_form"] = bool((Lw[post] > 0).all())
        out["min_f_minus_after_form"] = float(fm[post].min()); out["max_f_minus_after_form"] = float(fm[post].max())
        out["min_f_plus_after_form"] = float(fp[post].min()); out["max_f_plus_after_form"] = float(fp[post].max())
        # coarsening rate d xi_c / d tau (comoving domain scale xi_c = 1/L_wall) in the box-safe window
        # xi_c <= L/8, starting 5 conformal units after formation
        xic = np.where(Lw > 0, 1.0 / np.maximum(Lw, 1e-300), np.nan)
        afr0 = cfg.get("a_freeze")
        tend0 = tau[-1] if not afr0 else afr0 * cfg["tau_i"]
        selc = (tau >= tau[i0] + 5) & (tau <= tend0) & (xic <= Lbox / 8)
        if selc.sum() >= 4:
            out["dxi_dtau_boxsafe"] = float(np.polyfit(tau[selc], xic[selc], 1)[0])
        if has_vor:
            xiv = np.where(Lv > 0, 1.0 / np.maximum(Lv, 1e-300), np.nan)
            t_start = tau[i0] + (5 if not cfg.get("init_from") else 0)
            selv = (tau >= t_start) & (tau <= tend0) & (xiv <= Lbox / 8)
            if selv.sum() >= 4:
                out["dxi_vor_dtau_boxsafe"] = float(np.polyfit(tau[selv], xiv[selv], 1)[0])
            if afr0:
                selfz = tau >= afr0 * cfg["tau_i"]
                if selfz.sum() >= 4:
                    out["dxi_vor_dtau_frozen"] = float(np.polyfit(tau[selfz], xiv[selfz], 1)[0])
                    out["xi_vor_frozen_start_end"] = [float(xiv[selfz][0]), float(xiv[selfz][-1])]
            out["L_vor_over_L_wall_end"] = float(Lv[-1] / Lw[-1]) if Lw[-1] > 0 else None
            out["boxsafe_window_tau"] = [float(tau[selc][0]), float(tau[selc][-1])]
            out["boxsafe_window_xi"] = [float(xic[selc][0]), float(xic[selc][-1])]
        # is the coarsening slowing down (approach to freezing)?  rate in the first and second half of the
        # post-formation window that ends at xi_c = L/4 or at the end of the expanding phase
        selw = (tau >= tau[i0] + 5) & (tau <= tend0) & (xic <= Lbox / 4)
        if selw.sum() >= 8:
            tt = tau[selw]; xx = xic[selw]
            h = len(tt) // 2
            out["dxi_dtau_first_half"] = float(np.polyfit(tt[:h], xx[:h], 1)[0])
            out["dxi_dtau_second_half"] = float(np.polyfit(tt[h:], xx[h:], 1)[0])
            out["halves_window_tau"] = [float(tt[0]), float(tt[h]), float(tt[-1])]
        # charge carried by walls per unit physical wall area, expressed as an equivalent Q-matter slab thickness
        # (charge per area / n_*); 2D: per unit z-length, physical
        if cfg["kind"] != "neutral" and run.get("Q0"):
            qw = arr(s, "qfrac_core_walltype")
            Qc = run["Q0"]
            Vc = Lbox ** d
            with np.errstate(divide="ignore", invalid="ignore"):
                sig_q = qw * Qc / Vc / np.where(Lw > 0, Lw, np.nan) / a ** 2     # physical charge per physical area
            thick = sig_q / 1.77828
            out["wall_charge_slab_thickness_ell_form"] = float(thick[i0]) if np.isfinite(thick[i0]) else None
            out["wall_charge_slab_thickness_ell_end"] = float(thick[-1]) if np.isfinite(thick[-1]) else None
        # topology class after formation (2D torus): 'stripe' = both signs wrap (necessarily along the same axis),
        # 'cross' = one sign wraps both axes, the other none, 'single-axis' / 'none' otherwise
        cls = []
        for k in range(i0, len(s)):
            wp, wm = s[k]["wraps_plus"], s[k]["wraps_minus"]
            bp, bm = s[k]["wraps_both_axes_plus"], s[k]["wraps_both_axes_minus"]
            if wp and wm:
                cls.append("both_wrap")
            elif bp or bm:
                cls.append("one_sign_wraps_both_axes")
            elif wp or wm:
                cls.append("one_sign_wraps_one_axis")
            else:
                cls.append("no_wrap")
        out["topology_at_form"] = cls[0]
        out["topology_end"] = cls[-1]
        # slope of the wall length density over the later 3/4 of the post-formation interval (expanding part)
        afr = cfg.get("a_freeze")
        t_end_exp = tau[-1] if not afr else afr * cfg["tau_i"]
        t0 = tau[i0] + 0.25 * (t_end_exp - tau[i0])
        sel = (tau >= t0) & (tau <= t_end_exp) & (Lw > 0)
        if sel.sum() >= 4:
            out["mu_Lwall_vs_tau_late"] = float(np.polyfit(np.log(tau[sel]), np.log(Lw[sel]), 1)[0])
        if afr:
            ifr = np.nonzero(a >= afr - 1e-9)[0]
            if ifr.size:
                j0 = int(ifr[0])
                tphys = (tau - tau[j0]) * afr
                xi = arr(s, "xi_wall_physical")
                out["freeze"] = dict(tau_freeze=float(tau[j0]), t_phys_end=float(tphys[-1]),
                                     xi_phys_at_freeze=float(xi[j0]), xi_phys_end=float(xi[-1]),
                                     L_wall_at_freeze=float(Lw[j0]), L_wall_end=float(Lw[-1]))
                # growth exponent of xi_phys vs t_phys measured from the formation time
                tf = (tau - tau[i0]) * np.where(a >= afr, afr, a)
                selq = (tau > tau[j0] + 0.2 * (tau[-1] - tau[j0])) & np.isfinite(xi)
                if selq.sum() >= 4:
                    tt = (tau[selq] - tau[i0]) * afr
                    out["freeze"]["beta_xi_vs_t_since_form"] = float(np.polyfit(np.log(tt), np.log(xi[selq]), 1)[0])
                # flat-space wall velocity proxy: d xi_comoving / d tau over the frozen interval
                xic = 1.0 / Lw
                sel2 = tau >= tau[j0]
                if sel2.sum() >= 4:
                    out["freeze"]["dxi_comoving_dtau"] = float(np.polyfit(tau[sel2], xic[sel2], 1)[0])
    out["end"] = {k: s[-1].get(k) for k in (
        "f_plus", "f_minus", "f_core", "mean_S", "L_wall", "L_bare", "L_all", "L_internal", "A_wall",
        "xi_wall_comoving", "xi_wall_physical", "qfrac_core", "qfrac_vac", "qfrac_F2gt1", "qfrac_core_walltype",
        "qfrac_core_walltype_F2lt1", "qfrac_core_single", "qfrac_core_restored", "qfrac_lumps", "qfrac_lumps_attached",
        "n_lumps", "n_lumps_attached", "lump_Q_quantiles", "lump_R_phys_median", "lump_mean_core_frac",
        "wraps_plus", "wraps_minus", "spans_plus", "spans_minus", "n_comp_plus", "n_comp_minus",
        "largest_frac_plus", "largest_frac_minus", "wrapping_volume_frac_plus", "wrapping_volume_frac_minus",
        "n_core_walltype", "n_core_single", "mean_n_phys")}
    if out["end"]["xi_wall_comoving"]:
        out["end"]["xi_comoving_over_Lbox"] = out["end"]["xi_wall_comoving"] / Lbox
    if out["end"]["f_plus"] + out["end"]["f_minus"] > 0:
        out["end"]["p_plus"] = out["end"]["f_plus"] / (out["end"]["f_plus"] + out["end"]["f_minus"])
    # max over the run of the charge on walls and in lumps
    for k in ("qfrac_lumps", "qfrac_lumps_attached", "qfrac_core_walltype", "n_lumps"):
        v = arr(s, k)
        if np.isfinite(v).any():
            out[f"max_{k}"] = float(np.nanmax(v))
    # Q-lump charge lists (selected outputs): lumps above the 2D thin-wall binding size Q > 17
    # (tube radius > 1.75 ell; charge per ell of z in units sqrt(lam) v^3 ell^3 = 100 X quanta)
    ll = []
    for m in s:
        if "lump_Q_list" in m:
            q = np.array(m["lump_Q_list"], float)
            att = np.array(m["lump_attached_list"], bool)
            zf = m["a"] if d == 2 else 1.0
            Qt = m["Q_comoving"] / zf
            big = q > 17.0
            ll.append(dict(a=m["a"], n_lumps=int(q.size), n_Qgt17=int(big.sum()),
                           frac_Qgt17_attached=float(att[big].mean()) if big.any() else None,
                           qfrac_Qgt17=float(q[big].sum() / Qt) if Qt else 0.0,
                           qfrac_Qgt17_attached=float(q[big & att].sum() / Qt) if Qt else 0.0,
                           Q_max=float(q.max()) if q.size else 0.0,
                           Q_median_Qgt17=float(np.median(q[big])) if big.any() else None))
    out["lump_lists"] = ll
    # conservation
    dQ = arr(s, "dQ_rel")
    out["max_abs_dQ_rel"] = float(np.nanmax(np.abs(dQ)))
    out["max_abs_Q_drift_comoving"] = float(np.nanmax(np.abs(arr(s, "Q_abs_drift", 0.0))))
    out["Q0"] = run.get("Q0")
    kg = arr(s, "H_kin") + arr(s, "H_grad")
    res = arr(s, "budget_resid", 0.0)
    flow = np.abs(arr(s, "budget_flow", 0.0))
    out["max_budget_resid_over_kin_grad"] = float(np.nanmax(np.abs(res[1:]) / kg[1:]))
    out["final_budget_resid_over_kin_grad"] = float(res[-1] / kg[-1])
    out["max_budget_resid_over_H"] = float(np.nanmax(np.abs(res) / np.abs(arr(s, "H"))))
    out["final_budget_resid_over_flow"] = float(abs(res[-1]) / flow[-1]) if flow[-1] > 0 else None
    # resolution of the wall at the end: cells across 2 delta (core |phi| < tanh 1)
    out["cells_across_wall_core_end"] = 2 * math.sqrt(2) / (a[-1] * cfg["dx"])
    out["H_over_mX_at_start"] = 1.0 / cfg["tau_i"]
    if cfg["kind"] == "charged_exit":
        ac = cfg["kappa"] ** (1 / 3)
        out["a_cross_nc"] = ac
        out["H_over_mX_at_nc_crossing"] = 1.0 / (cfg["tau_i"] * ac * ac)
        out["n_at_form_over_nc"] = (cfg["kappa"] / out["a_form"] ** 3) if "a_form" in out else None
    return out, dict(tau=tau, a=a, Lw=Lw, fp=fp, fm=fm, Lv=Lv, has_vor=has_vor)


def timeseries_table(run, keys, every=1):
    s = run["series"]
    rows = []
    for k in range(0, len(s), every):
        rows.append([s[k].get(x) for x in keys])
    if (len(s) - 1) % every:
        rows.append([s[-1].get(x) for x in keys])
    return dict(columns=keys, rows=rows)


def main():
    names = sorted(os.path.basename(p)[:-5] for p in glob.glob(os.path.join(RES, "runs", "*.json")))
    runs = {n: load(n) for n in names}
    per = {}
    ts = {}
    for n, r in runs.items():
        per[n], ts[n] = metrics(r)
    summ = dict(per_run=per)
    # charged vs neutral at matched tau
    comp = {}
    for n in names:
        if not n.startswith(("A_", "F_k")):
            continue
        ti = int(per[n]["tau_i"])
        partner = f"B_t{ti}_s{per[n]['seed']}" if n.startswith("A_") else "F_neutral_t150_freeze1.6"
        if partner not in ts:
            partner = f"B_t{ti}_s1"
        if partner not in ts:
            continue
        A, B = ts[n], ts[partner]
        tmax = min(A["tau"][-1], B["tau"][-1])
        rat = []
        for frac in (0.5, 0.75, 1.0):
            t = per[n].get("tau_form", A["tau"][0]) + frac * (tmax - per[n].get("tau_form", A["tau"][0]))
            la = interp_at(t, A["tau"], A["Lw"]); lb = interp_at(t, B["tau"], B["Lw"])
            rat.append(dict(tau=t, a=t / per[n]["tau_i"], L_charged=la, L_neutral=lb,
                            ratio=la / lb if lb > 0 else None))
        comp[n] = dict(partner=partner, matched=rat, mu_charged=per[n].get("mu_Lwall_vs_tau_late"),
                       mu_neutral=per[partner].get("mu_Lwall_vs_tau_late"))
    summ["charged_vs_neutral"] = comp
    dc = {}
    for ch, dcn in (("CC_k2_t150_s1_from1.6", "DC_k2_t150_s1_from1.6"), ("CC_k2_t40_s1_from1.8", "DC_k2_t40_s1_from1.8"),
                    ("F_k2_t150_freeze1.6", "F_dc_k2_t150_freeze1.6"), ("A_k2_t150_s1", "DC_k2_t150_s1_from1.6"),
                    ("A_k2_t40_s1", "DC_k2_t40_s1_from1.8")):
        if ch in ts and dcn in ts:
            A, B = dict(ts[ch]), dict(ts[dcn])
            use_vor = A["has_vor"] and B["has_vor"]
            if use_vor:
                A["Lw"], B["Lw"] = A["Lv"], B["Lv"]
            t0 = B["tau"][0]
            rows = []
            for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
                t = t0 + frac * (min(A["tau"][-1], B["tau"][-1]) - t0)
                la = interp_at(t, A["tau"], A["Lw"]); lb = interp_at(t, B["tau"], B["Lw"])
                rows.append(dict(tau=t, L_charged=la, L_decharged=lb, ratio=la / lb if lb > 0 else None,
                                 xi_c_charged=1 / la if la > 0 else None, xi_c_decharged=1 / lb if lb > 0 else None,
                                 fplus_ch=interp_at(t, A["tau"], A["fp"]), fplus_dc=interp_at(t, B["tau"], B["fp"])))
            sel = A["tau"] >= t0
            selb = B["tau"] >= t0
            fa = np.polyfit(A["tau"][sel], 1 / np.maximum(A["Lw"][sel], 1e-300), 1)[0]
            fb = np.polyfit(B["tau"][selb], 1 / np.maximum(B["Lw"][selb], 1e-300), 1)[0]
            dc[ch] = dict(decharged=dcn, length_measure="L_vor" if use_vor else "L_wall", rows=rows,
                          dxi_dtau_charged=float(fa), dxi_dtau_decharged=float(fb),
                          rate_ratio_charged_over_decharged=float(fa / fb))
    summ["decharged_controls"] = dc
    if "S_k2_t150_s1_to1.6" in runs and "A_k2_t150_s1" in runs:
        sa = runs["S_k2_t150_s1_to1.6"]["series"]; sb = runs["A_k2_t150_s1"]["series"]
        n = min(len(sa), len(sb))
        dif = max(abs(sa[k]["L_wall"] - sb[k]["L_wall"]) + abs(sa[k]["f_plus"] - sb[k]["f_plus"]) for k in range(n))
        summ["determinism_check"] = dict(outputs_compared=n, max_abs_diff_Lwall_plus_fplus=dif)
    # sign statistics
    sg = {}
    for grp, pref in (("charged_all", "A_"), ("neutral_all", "B_"), ("charged_main_k2_t40", "A_k2_t40")):
        vals_end = [per[n]["end"].get("p_plus") for n in names if n.startswith(pref) and per[n]["end"].get("p_plus")]
        vals_form = [per[n].get("p_plus_form") for n in names if n.startswith(pref) and per[n].get("p_plus_form")]
        if vals_end:
            sg[grp] = dict(n=len(vals_end), p_plus_end=vals_end, mean_end=float(np.mean(vals_end)),
                           std_end=float(np.std(vals_end, ddof=1)) if len(vals_end) > 1 else None,
                           p_plus_form=vals_form, mean_form=float(np.mean(vals_form)) if vals_form else None)
    summ["sign_statistics"] = sg
    # resolution check
    if "R_lo_k2_t40" in runs and "R_hi_k2_t40" in runs:
        lo, hi = runs["R_lo_k2_t40"]["series"], runs["R_hi_k2_t40"]["series"]
        tl = arr(lo, "tau"); th = arr(hi, "tau")
        rc = {}
        i0 = int(np.nonzero(arr(lo, "f_plus") + arr(lo, "f_minus") >= 0.5)[0][0])
        for k in ("L_wall", "f_plus", "f_minus", "qfrac_core_walltype", "qfrac_lumps", "qfrac_core", "n_lumps"):
            yl = arr(lo, k)[i0:]
            yh = np.interp(tl[i0:], th, arr(hi, k))
            rel = np.abs(yl - yh) / np.maximum(np.abs(yh), 1e-12)
            rc[k] = dict(mean_rel_diff=float(np.mean(rel)), max_rel_diff=float(np.max(rel)),
                         late_mean_rel_diff=float(np.mean(rel[len(rel) // 2:])),
                         end_lo=float(yl[-1]), end_hi=float(yh[-1]))
        rc["a_form_lo"] = per["R_lo_k2_t40"].get("a_form"); rc["a_form_hi"] = per["R_hi_k2_t40"].get("a_form")
        summ["resolution_check"] = rc
    # physical-scale numbers (derived)
    MP = 2.435323e18
    gs = 106.75
    phys = {}
    for T in (2.0e5, 1.982e5):
        H = math.sqrt(math.pi ** 2 * gs / 90) * T * T / MP
        phys[f"T_SM_{T:.4g}_GeV"] = dict(H_GeV=H, H_over_mX=H / 100.0, H_over_mphi=H / 141.42136,
                                         tau_Q_over_ell=1 / (2 * H / 100.0),
                                         xi_KZ_ell_estimate=(1 / (2 * H / 100.0)) ** (1 / 3))
    phys["simulated_H_over_mX_at_nc_crossing"] = {n: per[n].get("H_over_mX_at_nc_crossing") for n in names
                                                 if per[n].get("H_over_mX_at_nc_crossing")}
    summ["physical_scales"] = phys
    # compact time-series tables for NOTES (every ~10th output)
    tab = {}
    keys = ["a", "tau", "f_plus", "f_minus", "f_core", "L_wall", "A_wall", "xi_wall_physical", "qfrac_vac",
            "qfrac_core_walltype", "qfrac_core_single", "qfrac_lumps", "qfrac_lumps_attached", "n_lumps",
            "n_lumps_attached", "wraps_plus", "wraps_minus"]
    for n in names:
        ev = max(1, len(runs[n]["series"]) // 12)
        tab[n] = timeseries_table(runs[n], keys, ev)
    summ["tables"] = tab
    with open(os.path.join(RES, "summary.json"), "w") as f:
        json.dump(summ, f, indent=1, default=float)
    print(json.dumps({n: dict(a_form=per[n].get("a_form"), p_form=per[n].get("p_plus_form"),
                              p_end=per[n]["end"].get("p_plus"), both_wrap=per[n].get("frac_outputs_both_signs_wrap"),
                              mu=per[n].get("mu_Lwall_vs_tau_late"), Aend=per[n]["end"].get("A_wall"),
                              qlumps=per[n]["end"].get("qfrac_lumps"), qwall=per[n]["end"].get("qfrac_core_walltype"),
                              dQ=per[n]["max_abs_dQ_rel"], res=per[n]["max_budget_resid_over_kin_grad"])
                      for n in names}, indent=0, default=float))
    print(json.dumps(comp, indent=0, default=float))
    print(json.dumps(sg, indent=0, default=float))


if __name__ == "__main__":
    main()
