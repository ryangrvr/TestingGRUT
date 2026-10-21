"""Combine the per-e-fold minority production (production.py) with re-entry, collapse energy and the relic bound.

Bookkeeping [derived, NOTES §2-4]. A minority domain created n e-folds after the start of the last-N_* window, with
physical radius c_R/H_inf at creation, has comoving wavenumber k = a_n H_inf/c_R and re-enters (a R = 1/H) at T_re with

    a(T_re) H(T_re) = a0 H0 exp(n'),   n' = n - DeltaN + s_off,   DeltaN = N_* - N_0,   s_off = -ln c_R,
    N_0 = ln(a_end H_inf/(a0 H0))  (e-folds before the end at which today's Hubble scale exits).

Instant reheating [P]: a_end = a(T_RH), 3 M_P^2 H_inf^2 = (pi^2/30) g_rho T_RH^4, so N_0 = N_inst(H_inf) ~ 48.1-49.4
and in the radiation era T_re = T_RH exp(-(N_* - n) - s_off ... ) x [g^(1/2) g_s^(-1/3)](T_RH)/[same](T_re).
Alternative [P]: "exit" convention DeltaN = 0 (fresh start when today's Hubble scale exits; N_* only marks the end).

Counts in today's Hubble volume (comoving radius 1/(a0 H0)):
    N(T_re < T_X) = int_{0}^{n'_X} (dp/dn) e^{3 n'} dn'.
Super-horizon minority walls crossing today's Hubble volume (n' < 0): N_wall,0 = 6 int_{n'<0} (dp/dn) e^{n'} dn'.
Collapse relic of stable hidden quanta (rest-mass fraction f, energy factor C_E = 3 kappa_E c_A/(4 pi)):
    rho/s = f C_E sigma int_{T_re >= T_cut} (dp/dn) H(T_re)/s(T_re) dn'   <= 4.373e-10 GeV  (Omega h^2 <= 0.12)
Late domains (T_re < T_dom): upper estimate of a black-hole/inhomogeneity mass density, all enclosed energy plus the
wall energy (capped at the enclosed energy, i.e. at ~ the horizon mass) collapsing:
    rho_late/s <= int_{T_re < T_dom} (dp/dn) [rho_tot/s](T_re) [1 + min(kappa_E sigma/(4 pi M_P^2 H_re), 1)] dn'.
"""
import math
import json
import os
import numpy as np
import mr_common as C

DN = 0.05
# precomputed background tables on n'
_NP = np.arange(-30.0, 75.0 + 1e-9, 0.01)
_T = C.T_of_nprime(np.clip(_NP, C.nprime_of_T(C.T0 * 2.2) + 1e-9, None))
_valid = np.isfinite(_T)
_lnHs = np.full_like(_NP, np.nan)
_lnRs = np.full_like(_NP, np.nan)
for i in np.where(_valid)[0]:
    T = _T[i]
    _lnHs[i] = math.log(C.H_of_T(T) / C.s_of_T(T))
    _lnRs[i] = math.log(3 * C.M_P ** 2 * C.H_of_T(T) ** 2 / C.s_of_T(T))


def tab(x, arr):
    return np.interp(x, _NP[_valid], arr[_valid], left=np.nan, right=np.nan)


KAPPA_E_RD = None
C_E_SPHERE = None


def load_collapse():
    global KAPPA_E_RD, C_E_SPHERE
    d = json.load(open(os.path.join(C.RESULTS, "collapse.json")))
    KAPPA_E_RD = d["adopted"]["kappa_E_RD"]
    C_E_SPHERE = 3 * KAPPA_E_RD / (4 * math.pi)


def cum_log(x, logf):
    """log of cumulative trapezoid int_{x0}^{x} exp(logf), robust for exponentially varying integrands."""
    f = np.exp(logf - np.nanmax(logf))
    f = np.where(np.isfinite(f), f, 0.0)
    c = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))])
    with np.errstate(divide="ignore"):
        return np.log(c) + np.nanmax(logf)


def interp_log(xq, x, logc):
    ok = np.isfinite(logc)
    if not np.any(ok) or xq < x[ok][0]:
        return -np.inf
    if xq >= x[ok][-1]:
        return logc[ok][-1]
    return float(np.interp(xq, x[ok], logc[ok]))


def quantities(ns, ldpdn, N_star, DeltaN, s_off=0.0, FX=1.0, T_cuts=None):
    """All criteria for one (H, start, N_*, DeltaN, s_off, F_X). ldpdn = ln(dp/dn) on ns."""
    m = ns <= N_star + 1e-9
    n = ns[m]
    lf = ldpdn[m] + math.log(FX)
    npr = n - DeltaN + s_off
    out = {}
    # counts after T_X (n' in [0, n'_X])
    pos = npr >= 0
    x = npr[pos]
    if x.size > 1:
        lc = cum_log(x, lf[pos] + 3 * x)
        for key, TX in T_cuts["counts"].items():
            nX = C.nprime_of_T(TX)
            out["log10_N_after_" + key] = interp_log(nX, x, lc) / math.log(10)
        out["log10_N_all_sizes_in_Hubble_volume"] = lc[-1] / math.log(10) if np.isfinite(lc[-1]) else -np.inf
    # super-horizon walls crossing today's Hubble volume
    neg = npr < 0
    if neg.sum() > 1:
        xn = npr[neg]
        ln6 = cum_log(xn, lf[neg] + xn)
        out["log10_N_wall_superhorizon_today"] = (math.log(6.0) + ln6[-1]) / math.log(10)
    else:
        out["log10_N_wall_superhorizon_today"] = -np.inf
    # relic: reverse cumulative from the end of the window (largest n') down to n'_cut
    lHs = tab(npr, _lnHs)
    ok = np.isfinite(lHs)
    xr = npr[ok][::-1]
    lr = cum_log(-xr, (lf[ok] + lHs[ok])[::-1])       # integrate from max n' downward
    for key, Tc in T_cuts["relic"].items():
        ncut = C.nprime_of_T(Tc)
        val = interp_log(-ncut, -xr, lr)
        out["log10_rho_over_s_unit_" + key] = (val + math.log(C.SIGMA)) / math.log(10)   # times f C_E
    # late-domain enclosed energy (upper estimate) for T_re < T_dom,iso, n' >= 0
    nd = C.nprime_of_T(C.T_DOM_ISO)
    lm = (npr >= 0) & (npr <= nd)
    if lm.sum() > 1:
        lR = tab(npr[lm], _lnRs)
        # per domain: enclosed energy (horizon mass) + min(wall energy, enclosed energy)  [upper estimate]
        wall_over_encl = 3 * KAPPA_E_RD * C.SIGMA * np.exp(tab(npr[lm], _lnHs) - lR) / (4 * math.pi)   # = kappa_E sigma/(4 pi M_P^2 H)
        lR = lR + np.log1p(np.minimum(wall_over_encl, 1.0))
        lcl = cum_log(npr[lm], lf[lm] + lR)
        out["log10_rho_late_over_s_max"] = lcl[-1] / math.log(10)
    else:
        out["log10_rho_late_over_s_max"] = -np.inf
    # where the relic integrand peaks (dominant re-entry temperature)
    if ok.any():
        w = lf[ok] + lHs[ok]
        sel = npr[ok] <= C.nprime_of_T(T_cuts["relic"]["Tdom_iso"])
        if sel.any():
            i = np.argmax(np.where(sel, w, -np.inf))
            out["relic_dominant_nprime"] = float(npr[ok][i])
            out["relic_dominant_T_re_GeV"] = float(C.T_of_nprime(npr[ok][i]))
    return out


T_CUTS = {"counts": {"Tdom_iso": C.T_DOM_ISO, "Tdom_net": C.T_DOM_NET, "BBN_1MeV": C.T_BBN, "BBN_70keV": C.T_BBN_D,
                     "Teq": None, "Trec": C.T_REC},
          "relic": {"Tdom_iso": C.T_DOM_ISO, "Tdom_net": C.T_DOM_NET, "BBN_1MeV": C.T_BBN, "late_0.03eV": 3e-11}}


def T_eq():
    from scipy.optimize import brentq
    import core
    f = lambda lT: math.log(core.rho_rad(math.exp(lT))) - math.log(core.rho_matter(math.exp(lT)))
    return math.exp(brentq(f, math.log(1e-10), math.log(1e-8)))


T_CUTS["counts"]["Teq"] = T_eq()


def threshold(Hs, q, target):
    """Lowest H at which q(H) (log10) crosses log10(target) from below; log-log interpolation."""
    Hs = np.asarray(Hs, float)
    q = np.asarray(q, float)
    lt = math.log10(target)
    ok = np.isfinite(q)
    if ok.sum() < 2:
        return None
    x = np.log10(Hs[ok])
    y = q[ok]
    if y[0] >= lt:
        return {"H_GeV": None, "note": "above target at the lowest H on the grid (%.1f GeV)" % Hs[ok][0]}
    idx = np.where((y[:-1] < lt) & (y[1:] >= lt))[0]
    if idx.size == 0:
        return {"H_GeV": None, "note": "not reached on the grid (max %.1f GeV)" % Hs[ok][-1]}
    i = idx[0]
    Hstar = 10 ** (x[i] + (lt - y[i]) * (x[i + 1] - x[i]) / (y[i + 1] - y[i]))
    return {"H_GeV": float(Hstar), "monotone_below": bool(np.all(np.diff(y[: i + 1]) >= -1e-9))}
