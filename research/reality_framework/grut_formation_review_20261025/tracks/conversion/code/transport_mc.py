#!/usr/bin/env python3
"""Independent conversion-product transport (Claude review, track `conversion`).

Own methods, different from Codex 59a1a71 (grid in delta with remapping):
  1. first-event survival: adaptive quadrature of the optical depth along the
     never-converted trajectory u(a) = u0 a_f / a;
  2. population: grid-free Monte Carlo of individual particles in (species,
     comoving momentum), exponential waiting times inside short ln(a) steps,
     several events per step allowed, exact analytic threshold crossing
     a_cross = q / u_esc for positive particles (E_+ = 1 + 2 Omega);
     late anti-X tail by an analytic remaining-survival weight;
  3. stiff 20 MeV case: local-equilibrium drift ODE in nu (Gamma/H ~ 1e6),
     cross-checked by a two-state Monte Carlo at the same epoch;
  4. horizons: own quadrature of the collisionless comoving path after escape,
     with a pure-radiation asinh control and a g*(T) background sensitivity.
Inputs quoted [Q]: rate table sigma_v(p) of bb211d1 RATE_RESULTS.json (also my own
independent dense scan, see scatter.py), Omega, Y_Q, s0, Q=200000, cosmology.
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/conversion')
RATE_JSON = SCR / 'pbb/research/reality_framework/grut_core_conversion_20261023/RATE_RESULTS.json'
RES = Path(__file__).resolve().parents[1] / 'results'

MPC_CM = 3.085677581491367e24
H0 = 67.4e5 / MPC_CM                     # s^-1
OR = 9.20917977964354e-5                 # quoted cosmology of 59a1a71
OM = 0.315
OL = 1 - OM - OR
C_MPC_S = 2.99792458e10 / MPC_CM
S0 = 2891.2                              # cm^-3
YQ = 4.8294009076988175e-9
QOBJ = 2.0e5
NOBJ0 = S0 * YQ / QOBJ                   # objects per cm^3 today (comoving)
OMEGA = 0.8611951308217225               # confirmed by profile.py
EB = 1 + 2 * OMEGA                       # escape energy / m_X
UESC = np.sqrt(EB ** 2 - 1)
E0 = np.sqrt(1.5); U0 = np.sqrt(0.5)
T0_GEV = 2.7255 * 8.617333262e-14
GS0 = 2 + 7 / 8 * 6 * (4 / 11)
EPOCHS = {'T20MeV': T0_GEV / 0.02 * (GS0 / 10.75) ** (1 / 3),
          'T1keV': T0_GEV / 1e-6,
          'z1e6': 1 / 1000001.}


def hubble(a):
    return H0 * np.sqrt(OR / a ** 4 + OM / a ** 3 + OL)


# --- rate tables -----------------------------------------------------------
class Rate:
    """sigma_v(p) for anti-X; interpolation 'loglin' (my choice) or 'pchip'."""

    def __init__(self, p, sv, kind='loglin'):
        o = np.argsort(p)
        self.lp = np.log(np.asarray(p)[o]); self.ls = np.log(np.asarray(sv)[o])
        self.kind = kind
        self.pchip = PchipInterpolator(self.lp, self.ls, extrapolate=False)

    def __call__(self, p):
        x = np.log(np.clip(p, np.exp(self.lp[0]), np.exp(self.lp[-1])))
        if self.kind == 'pchip':
            return np.exp(self.pchip(x))
        return np.exp(np.interp(x, self.lp, self.ls))


def codex_table(kind='loglin'):
    rows = json.loads(RATE_JSON.read_text())['rows']
    return Rate([r['p_over_m'] for r in rows], [r['sigma_v_conversion_cm3_s'] for r in rows], kind)


def my_table(kind='loglin'):
    rows = []
    for tag in 'abc':
        f = SCR / f'scan_{tag}.json'
        if f.exists():
            rows += json.loads(f.read_text())
    for extra in sorted(SCR.glob('scan_refine*.json')):
        rows += json.loads(extra.read_text())
    spot = json.loads((RES / 'scatter_spot.json').read_text())['rows']
    p = [r['p_over_m'] for r in rows] + [r['p_over_m'] for r in spot]
    sv = [r['sigma_v_cm3_s'] for r in rows] + [r['sigma_v_cm3_s'] for r in spot]
    p, i = np.unique(p, return_index=True)
    return Rate(p, np.asarray(sv)[i], kind)


def lam_minus(rate, a, p):
    """anti-X conversion rate per e-fold of a."""
    return NOBJ0 / a ** 3 * rate(p) / hubble(a)


def lam_plus(rate, a, p_plus):
    """inverse (X -> anti-X) rate per e-fold from detailed balance."""
    Ep = np.sqrt(1 + p_plus ** 2)
    Em = Ep - 2 * OMEGA
    ok = Em > 1
    pm = np.sqrt(np.where(ok, Em * Em - 1, 0.0))
    sv = np.where(ok, rate(np.maximum(pm, 1e-12)) * pm * Em / (p_plus * Ep), 0.0)
    return NOBJ0 / a ** 3 * sv / hubble(a)


# --- 1. first-event survival -----------------------------------------------
def first_event(rate, af):
    def integrand(y):
        a = af * np.exp(y)
        return lam_minus(rate, a, U0 * np.exp(-y))
    span = np.log(1 / af)
    # break points where the never-converted momentum crosses table knots
    knots = np.sort(np.log(U0) - rate.lp); knots = knots[(knots > 0) & (knots < span)]

    def depth(y):
        pts = np.r_[0.0, knots[knots < y], y]
        return sum(quad(integrand, x0, x1, epsabs=0, epsrel=1e-12, limit=200)[0]
                   for x0, x1 in zip(pts[:-1], pts[1:]) if x1 > x0)
    total = depth(span)
    out = dict(a_form=af, z_form=1 / af - 1, Gamma_over_H_initial=float(integrand(0.0)),
               optical_depth_to_today=total, never_converted_today=float(np.exp(-total)))
    for frac in [0.1, 0.5, 0.9, 0.99]:
        tgt = -np.log1p(-frac)
        if total > tgt:
            y = brentq(lambda y: depth(y) - tgt, 0, span, xtol=1e-14)
            out[f'z_first_{frac}'] = float(1 / (af * np.exp(y)) - 1)
            out[f'dlna_first_{frac}'] = float(y)
    return out


# --- 2. Monte Carlo population ------------------------------------------------
def monte_carlo(rate, af, n=200000, ds_fac=0.004, seed=1, p_cold=0.05, max_dlna=6.0):
    """Grid-free MC.  Step ln(a) by ds = ds_fac / (largest current rate per e-fold)."""
    rng = np.random.default_rng(seed)
    sp = -np.ones(n, dtype=np.int8)          # -1 anti, +1 active X, 0 escaped
    qc = np.full(n, af * U0)                 # comoving momentum q = a p (units m_X)
    never = np.ones(n, dtype=bool)           # never converted (first-event survival)
    a_esc = np.full(n, np.nan)
    nfwd = 0; ninv = 0
    lna = np.log(af); lna_end = lna + max_dlna
    hist = []; first_hist = []
    while lna < lna_end:
        act = np.nonzero(sp != 0)[0]
        if act.size == 0:
            break
        a = np.exp(lna)
        npos = np.count_nonzero(sp == 1)
        anti = sp == -1
        if npos == 0 and (not anti.any() or np.max(qc[anti]) / a < p_cold):
            break
        s_ = sp[act]; p_ = qc[act] / a
        lam0 = np.where(s_ == -1, lam_minus(rate, a, p_), lam_plus(rate, a, p_))
        ds = min(ds_fac / max(lam0.max(), 1e-3), 0.01)
        rem = np.full(act.size, ds)
        t0 = np.full(act.size, lna)
        idx = act.copy()
        while idx.size:
            s = sp[idx]
            amid = np.exp(t0 + rem / 2)
            p = qc[idx] / amid
            lamv = np.where(s == -1, lam_minus(rate, amid, p), lam_plus(rate, amid, p))
            wait = rng.exponential(1.0, idx.size) / np.maximum(lamv, 1e-300)
            lcross = np.where(s == 1, np.log(qc[idx] / UESC), np.inf)
            tev = t0 + wait
            tend = t0 + rem
            esc = (s == 1) & (lcross <= np.minimum(tev, tend))
            a_esc[idx[esc]] = np.exp(lcross[esc]); sp[idx[esc]] = 0
            ev = (~esc) & (tev < tend)
            j = idx[ev]
            if j.size:
                ae = np.exp(tev[ev])
                pe = qc[j] / ae
                E = np.sqrt(1 + pe * pe)
                fwd = sp[j] == -1
                Enew = np.where(fwd, E + 2 * OMEGA, E - 2 * OMEGA)
                pnew = np.sqrt(np.maximum(Enew * Enew - 1, 0.0))
                qc[j] = ae * pnew
                sp[j] = np.where(fwd, 1, -1).astype(np.int8)
                never[j] = False
                nfwd += int(np.count_nonzero(fwd)); ninv += int(np.count_nonzero(~fwd))
            keep = ev
            rem = (tend - tev)[keep]; t0 = tev[keep]; idx = idx[keep]
        lna += ds
        if len(hist) == 0 or (lna - np.log(af)) - hist[-1][0] > 0.005:
            hist.append([lna - np.log(af), float(np.count_nonzero(sp == -1) / n), float(np.count_nonzero(sp == 1) / n),
                         float(np.count_nonzero(sp == 0) / n), float(1 - never.mean()), nfwd / n, ninv / n])
    anti = sp == -1
    return dict(a_esc=a_esc, anti_mask=anti, qc=qc, a_end=float(np.exp(lna)), forward=nfwd / n,
                inverse=ninv / n, history=hist, n=n, never_converted_at_stop=float(never.mean()))


def cold_tail_survival(rate, a_start, q):
    """Probability that an anti-X with comoving momentum q at a_start never
    converts again until today (single-event survival; re-conversion of the
    resulting X back to anti-X is closed once p_- is small, see NOTES)."""
    f = lambda y: lam_minus(rate, np.exp(y), q / np.exp(y))
    tau = quad(f, np.log(a_start), 0.0, epsrel=1e-10, limit=400)[0]
    return np.exp(-tau)


# --- 3. stiff local-equilibrium drift ----------------------------------------
def equilibrium_drift(rate, af):
    """nu = E/m - q Omega; fast forward/inverse equilibrium at fixed nu.
    N_+/N_- = lam_-/lam_+ = p_+E_+/(p_-E_-).  d nu/dln a = -<e - 1/e>."""
    nu0 = E0 + OMEGA
    nuth = 1 + OMEGA

    def rhs(y, s):
        nu = s[0]
        em = max(nu - OMEGA, 1.0 + 1e-15); ep = nu + OMEGA
        pm = np.sqrt(em * em - 1); pp = np.sqrt(ep * ep - 1)
        fm = pm * em / (pm * em + pp * ep)             # anti-X fraction
        dnu = -(fm * (em - 1 / em) + (1 - fm) * (ep - 1 / ep))
        a = af * np.exp(y)
        fwd = lam_minus(rate, a, pm) * fm
        return [dnu, fwd]

    ev = lambda y, s: s[0] - nuth
    ev.terminal = True
    sol = solve_ivp(rhs, (0, 5), [nu0, 0.0], method='LSODA', events=ev, rtol=1e-11, atol=1e-14)
    y = sol.t_events[0][0]
    return dict(dlna_escape=float(y), a_escape=float(af * np.exp(y)), z_escape=float(1 / (af * np.exp(y)) - 1),
                forward_events_per_particle=float(sol.y_events[0][0][1]),
                Gamma_over_H_initial=float(lam_minus(rate, af, U0)))


# --- 4. horizons ---------------------------------------------------------------
def horizon(a1, q, a2=1.0, hub=hubble):
    if a1 >= a2:
        return 0.0
    f = lambda y: C_MPC_S * (q / np.exp(y)) / np.sqrt(1 + (q / np.exp(y)) ** 2) / (np.exp(y) * hub(np.exp(y)))
    return quad(f, np.log(a1), np.log(a2), epsabs=0, epsrel=1e-11, limit=400)[0]


def horizon_rad(a1, q, a2):
    return C_MPC_S * q / (H0 * np.sqrt(OR)) * (np.arcsinh(a2 / q) - np.arcsinh(a1 / q))


def light_prefix(af, ae):
    return quad(lambda y: C_MPC_S / (np.exp(y) * hubble(np.exp(y))), np.log(af), np.log(ae), epsrel=1e-11)[0]


def gstar_hubble():
    """Sensitivity background: rho_r(a) with g*(T), g*s(T) through e+e- annihilation
    (instantaneous neutrino decoupling, Tnu/T from entropy), normalised to OR today."""
    from scipy.integrate import quad as q_

    def fe(x, sign=1):  # energy and entropy integrals of a massive fermion, per dof, relative to massless
        if x == 0:
            return 1.0, 1.0
        rho = q_(lambda y: y * y * np.sqrt(y * y + x * x) / (np.exp(np.sqrt(y * y + x * x)) + 1), 0, 60 + x)[0]
        p = q_(lambda y: y ** 4 / np.sqrt(y * y + x * x) / 3 / (np.exp(np.sqrt(y * y + x * x)) + 1), 0, 60 + x)[0]
        rho0 = 7 / 8 * np.pi ** 4 / 15; s0 = 7 / 8 * 4 * np.pi ** 4 / 45
        return rho / rho0, (rho + p) / s0
    me = 0.51099895e-3
    Tg = np.geomspace(1e-7, 0.05, 400)     # photon temperature GeV
    out = []
    for T in Tg:
        r_e, s_e = fe(me / T)
        gs_eg = 2 + 7 / 8 * 4 * s_e              # photons + e+-: entropy shared after nu decoupling
        out.append((T, r_e, gs_eg))
    out = np.array(out)
    # entropy in (gamma, e) sector: gs_eg T^3 a^3 = const = (2) T0^3 (late)
    aT = (2.0 / out[:, 2]) ** (1 / 3) * T0_GEV / out[:, 0]     # a as a function of T_gamma
    Tnu = T0_GEV * (4 / 11) ** (1 / 3) / aT
    rho = (np.pi ** 2 / 30) * (2 * out[:, 0] ** 4 + 7 / 8 * 4 * out[:, 1] * out[:, 0] ** 4 + 7 / 8 * 6 * Tnu ** 4)
    rho_late = (np.pi ** 2 / 30) * (2 + 7 / 8 * 6 * (4 / 11) ** (4 / 3)) * (T0_GEV / aT) ** 4
    ratio = rho / rho_late                       # >1 not expected; equals 1 late
    o = np.argsort(aT)
    la, lr = np.log(aT[o]), ratio[o]
    # The quoted OR uses N_eff = 3.046; keep its normalisation and apply the relative change.

    def hub(a):
        r = np.interp(np.log(a), la, lr, left=lr[0], right=1.0)
        return H0 * np.sqrt(OR * r / a ** 4 + OM / a ** 3 + OL)
    return hub, dict(ratio_at_20MeV=float(np.interp(np.log(EPOCHS['T20MeV']), la, lr)))


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else 'all'
    table = sys.argv[2] if len(sys.argv) > 2 else 'codex'
    quick = len(sys.argv) > 3 and sys.argv[3] == 'quick'
    rate = codex_table() if table == 'codex' else my_table()
    rate_p = codex_table('pchip') if table == 'codex' else my_table('pchip')
    out = dict(table=table)
    # first-event survival
    fe = {}
    for k, af in EPOCHS.items():
        fe[k] = first_event(rate, af)
        fe[k + '_pchip'] = first_event(rate_p, af)
        print(k, fe[k], flush=True)
    out['first_event'] = fe
    # stiff equilibrium at 20 MeV (and, as a diagnostic, at the other epochs)
    eq = {k: equilibrium_drift(rate, af) for k, af in EPOCHS.items()}
    out['equilibrium_drift'] = eq
    print('eq', eq, flush=True)
    hub_g, ginfo = gstar_hubble()
    out['gstar_background'] = ginfo
    # Monte Carlo for the non-stiff epochs
    mc = {}
    for k in ['z1e6', 'T1keV']:
        af = EPOCHS[k]
        runs = []
        settings = [(0.004, 100000, 3)] if quick else [(0.004, 200000, 1), (0.001, 100000, 2)]
        for ds_fac, n, seed in settings:
            r = monte_carlo(rate, af, n=n, ds_fac=ds_fac, seed=seed)
            ae = r['a_esc'][np.isfinite(r['a_esc'])]
            anti = r['anti_mask']
            qa = r['qc'][anti]
            surv = np.array([cold_tail_survival(rate, r['a_end'], q) for q in qa])
            anti_today = float(surv.sum() / r['n'])
            S0_stop = float(np.exp(-quad(lambda y: lam_minus(rate, af * np.exp(y), U0 * np.exp(-y)), 0,
                                          np.log(r['a_end'] / af), limit=400)[0]))
            qs = UESC * ae
            grid = np.geomspace(ae.min() * 0.999, ae.max() * 1.001, 240)
            ht = np.array([horizon(x, UESC * x, 1.0) for x in grid])
            hz5 = np.array([horizon(x, UESC * x, 1 / 6) for x in grid])
            lp = np.array([light_prefix(af, x) for x in grid])
            hg = np.array([horizon(x, UESC * x, 1.0, hub_g) for x in grid[::8]])
            runs.append(dict(ds_fac=ds_fac, n=n, seed=seed, escaped_fraction=float(ae.size / r['n']),
                             anti_at_stop=float(anti.mean()), a_stop_over_af=r['a_end'] / af,
                             never_converted_at_stop_MC=r['never_converted_at_stop'],
                             never_converted_at_stop_analytic=S0_stop,
                             anti_today_estimate=anti_today, forward_per_particle=r['forward'],
                             inverse_per_particle=r['inverse'],
                             z_escape_quantiles={str(qq): float(1 / np.quantile(ae, qq) - 1) for qq in [0.1, 0.5, 0.9, 0.99]},
                             q_over_m_mean=float(qs.mean()), q_over_m_rel_spread=float(qs.std() / qs.mean()),
                             mean_horizon_today_Mpc=float(np.mean(np.interp(ae, grid, ht))),
                             mean_horizon_z5_Mpc=float(np.mean(np.interp(ae, grid, hz5))),
                             mean_lightcone_prefix_Mpc=float(np.mean(np.interp(ae, grid, lp))),
                             mean_horizon_today_gstar_background_Mpc=float(np.mean(np.interp(ae, grid[::8], hg))),
                             history=r['history']))
            print(k, {kk: vv for kk, vv in runs[-1].items() if kk != 'history'}, flush=True)
            np.save(SCR / f'mc_aesc_{k}_{table}_{seed}.npy', ae)
        mc[k] = runs
    out['monte_carlo'] = mc
    # 20 MeV: deterministic escape
    e20 = eq['T20MeV']; ae = e20['a_escape']
    out['T20MeV_escape'] = dict(a_escape=ae, q_over_m=UESC * ae,
                                horizon_today_Mpc=horizon(ae, UESC * ae, 1.0),
                                horizon_z5_Mpc=horizon(ae, UESC * ae, 1 / 6),
                                horizon_today_rad_only_Mpc=horizon_rad(ae, UESC * ae, 1.0),
                                lightcone_prefix_Mpc=light_prefix(EPOCHS['T20MeV'], ae),
                                horizon_today_gstar_background_Mpc=horizon(ae, UESC * ae, 1.0, hub_g))
    print('20MeV', out['T20MeV_escape'], flush=True)
    # original positive daughters (q = a_f u0) horizons
    out['original_daughter_horizon_today_Mpc'] = {k: horizon(af, af * U0, 1.0) for k, af in EPOCHS.items()}
    # radiation-only control of my quadrature
    out['quadrature_control'] = dict(rad_numeric=horizon(1e-6, 2.5e-6, 1e-5, hub=lambda a: H0 * np.sqrt(OR) / a ** 2),
                                     rad_analytic=horizon_rad(1e-6, 2.5e-6, 1e-5))
    if table != 'codex' and table != 'mine':
        pass
    (RES / f'transport_{table}.json').write_text(json.dumps(out, indent=1, default=float) + '\n')


if __name__ == '__main__':
    main()
