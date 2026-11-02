"""fx_formation SKEPTIC (Claude, 1 Nov 2026): independent recomputation of the vapour's charge-to-entropy ratio at
coexistence, the survival threshold Y_crit, T_ev / T_s, and the evaporation-window undersaturation.

Method differs from the track's thermo.py (momentum quadrature): ideal relativistic Bose gas by the Bessel series
(Kolb-Turner), n = m^2 T/(2 pi^2) sum_k e^{k mu/T} K2(k m/T)/k, etc., summed to convergence with scaled Bessel
functions. Coexistence: mu_c = Omega_b (T = 0 bulk Q-matter), with the vapour-pressure shift added (P_liq = n_B (mu -
Omega_b) = p_vap). Units m_X = 1.  Vapour species: X (m = 1, +mu), Xbar (m = 1, -mu), phi (m = sqrt 6, mu = 0).
Optional: thermal-mass shift dm^2 = c T^2 on X and phi (sensitivity), and a phonon-gas entropy for the liquid.
"""
import json, math, sys
import numpy as np
from scipy.special import kve
from scipy.optimize import brentq, minimize_scalar

OMB = 0.7758647074595397
NB = 449.31336788454803
MPHI = math.sqrt(6.0)


def bose(m, mu, T, kmax=None):
    """n, P, rho, s for one bosonic dof by the Bessel series."""
    if T <= 0:
        return 0.0, 0.0, 0.0, 0.0
    d = (m - mu) / T
    if d > 700:
        return 0.0, 0.0, 0.0, 0.0
    # terms fall like exp(-k d) k^{-3/2}; sum until < 1e-16 relative
    K = kmax or int(min(200000, max(50, 40.0 / max(d, 1e-6))))
    k = np.arange(1, K + 1, dtype=float)
    x = k * m / T
    w = np.exp(-k * d)                    # e^{k mu/T} e^{-k m/T}
    k2 = kve(2, x) * w
    k1 = kve(1, x) * w
    n = m * m * T / (2 * math.pi ** 2) * np.sum(k2 / k)
    P = m * m * T * T / (2 * math.pi ** 2) * np.sum(k2 / k ** 2)
    rho = m * m * T / (2 * math.pi ** 2) * np.sum(3 * T * k2 / k ** 2 + m * k1 / k)
    s = (rho + P - mu * n) / T
    return n, P, rho, s


def vapour(mu, T, cth=0.0):
    mX = math.sqrt(1.0 + cth * T * T)
    mF = math.sqrt(6.0 + cth * T * T)
    nx, px, ex, sx = bose(mX, mu, T)
    nb, pb, eb, sb = bose(mX, -mu, T)
    nf, pf, ef, sf = bose(mF, 0.0, T)
    return dict(n=nx - nb, nX=nx, nXb=nb, nphi=nf, P=px + pb + pf, rho=ex + eb + ef, s=sx + sb + sf)


_TRK = None


def mu_coex(T, cth=0.0):
    """Own coexistence (bulk liquid at T = 0 plus the vapour-pressure shift) for T <= 1.6; above that the liquid's
    thermal pressure matters and the track's tabulated mu_c(T) is interpolated (only the e = 1e-7 crossing uses it)."""
    if T <= 1.6:
        g = lambda mu: NB * (mu - OMB) - vapour(mu, T, cth)['P']
        return brentq(g, OMB - 0.05, 0.99, xtol=1e-14)
    global _TRK
    if _TRK is None:
        d = json.load(open(sys.argv[2]))['vapour_ratio_at_coexistence']
        _TRK = (np.log([r['T'] for r in d]), np.array([r['mu_c'] for r in d]))
    return float(np.interp(math.log(T), *_TRK))


def ratio_at(mu, T, cth=0.0):
    v = vapour(mu, T, cth)
    return v['n'] / v['s']


def ycrit(mu_fixed=None, cth=0.0):
    f = lambda lt: -ratio_at(mu_fixed if mu_fixed is not None else mu_coex(math.exp(lt), cth), math.exp(lt), cth)
    r = minimize_scalar(f, bounds=(math.log(0.1), math.log(1.0)), method='bounded', options=dict(xatol=1e-9))
    return -r.fun, math.exp(r.x)


def crossings(Yh, cth=0.0):
    """T where vapour n/s at coexistence equals Y_h (T_ev above the maximum, T_s below)."""
    ym, Tm = ycrit(cth=cth)
    if Yh >= ym:
        return None, None
    g = lambda lt: ratio_at(mu_coex(math.exp(lt), cth), math.exp(lt), cth) - Yh
    Tev = math.exp(brentq(g, math.log(Tm), math.log(30.0), xtol=1e-12))
    Ts = math.exp(brentq(g, math.log(1e-3), math.log(Tm), xtol=1e-12))
    return Tev, Ts


def vapour_only_mu(Yh, T):
    g = lambda mu: vapour(mu, T)['n'] - Yh * vapour(mu, T)['s']
    return brentq(g, 1e-9, 0.9999999, xtol=1e-14)


def main(out):
    res = {}
    # 0. checks of the series
    n, P, rho, s = bose(1e-8, 0.0, 1.0)
    res['check_massless'] = dict(n_over_zeta3_pi2=n / (1.2020569031595942 / math.pi ** 2), P_over_pi2_90=P / (math.pi ** 2 / 90),
                                 s_over_2pi2_45=s / (2 * math.pi ** 2 / 45))
    # NR Boltzmann limit: n -> (mT/2pi)^{3/2} e^{-(m-mu)/T} (1 + 15T/8m) for small fugacity
    T = 0.01; mu = 0.5
    n, P, rho, s = bose(1.0, mu, T)
    nnr = (T / (2 * math.pi)) ** 1.5 * math.exp(-(1 - mu) / T) * (1 + 15 * T / 8)
    res['check_NR'] = dict(ratio=n / nnr)
    # 1. compare with the track's tabulated ratio
    trk = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else None
    cmp_rows = []
    if trk:
        for r in trk['vapour_ratio_at_coexistence'][::5]:
            T = r['T']
            vt = vapour(r['mu_c'], T)          # same mu as the track: tests the gas thermodynamics alone
            row = dict(T=T, mu_c_track=r['mu_c'], ratio_track=r['ratio'], ratio_own_at_track_mu=vt['n'] / vt['s'],
                       rel_gas=(vt['n'] / vt['s']) / r['ratio'] - 1)
            if T <= 1.6:                       # own coexistence without liquid thermal pressure (valid at low T)
                mc = mu_coex(T)
                v = vapour(mc, T)
                row.update(mu_c_own=mc, ratio_own=v['n'] / v['s'])
            cmp_rows.append(row)
    res['comparison_with_track'] = cmp_rows
    # 2. Y_crit, bulk and finite objects (mu = Omega(Q)); thermal-mass sensitivity
    yb, Tb = ycrit()
    res['Y_crit_bulk'] = dict(Y_crit=yb, T_at_max=Tb, S_per_Q_max=1 / yb)
    res['Y_crit_finite_Q'] = {}
    for Q, Om in ((2e5, 0.8817204919524271), (1e6, 0.83429668), (1e8, 0.78777238), (1e11, OMB + 2 / 3 * 8.1429 * 1e11 ** (-1 / 3))):
        y, Tq = ycrit(mu_fixed=Om)
        res['Y_crit_finite_Q']['%g' % Q] = dict(Omega=Om, Y_crit=y, T_at_max=Tq)
    res['Y_crit_thermal_mass_sensitivity'] = {}
    for c in (1e-3, 1e-2, 0.1):
        y, Tq = ycrit(cth=c)
        res['Y_crit_thermal_mass_sensitivity']['dm2_over_T2=%g' % c] = dict(Y_crit=y, T_at_max=Tq)
    # mu_c sensitivity: shift by 1e-3
    y1, _ = ycrit(mu_fixed=OMB + 1e-3)
    res['dYcrit_per_1e-3_mu'] = y1 / yb - 1
    # 3. the Y_h of each label: conserved value (onset) vs the track's 'spinodal' value
    Ys = dict(heavy_e1e6_onset=0.2020822595873346, heavy_e1e6_spin_track=0.18538485573344238,
              pilot_e1e6=0.27680265413087224, heavy_R0=0.809575848826171, heavy_med=0.4839026445487382,
              heavy_e1e7=0.045831181928068294)
    res['crossings'] = {}
    for k, y in Ys.items():
        Tev, Ts = crossings(y)
        res['crossings'][k] = dict(Y_h=y, S_per_Q=1 / y, T_ev=Tev, T_s=Ts)
    # 4. undersaturation in the pilot window: with eta = 0 the vapour holds n/s = Y_h at mu < mu_c; the net
    #    evaporation drive is the fugacity deficit 1 - exp((mu - mu_c)/T)
    win = []
    Tev, Ts = res['crossings']['pilot_e1e6']['T_ev'], res['crossings']['pilot_e1e6']['T_s']
    for T in np.geomspace(Ts * 1.0001, Tev * 0.9999, 9):
        mu = vapour_only_mu(0.27680265413087224, T)
        mc = mu_coex(T)
        win.append(dict(T=float(T), mu=mu, mu_c=mc, deficit=1 - math.exp((mu - mc) / T)))
    res['pilot_window_undersaturation'] = win
    winh = []
    Tev, Ts = res['crossings']['heavy_e1e6_onset']['T_ev'], res['crossings']['heavy_e1e6_onset']['T_s']
    for T in np.geomspace(Ts * 1.0001, Tev * 0.9999, 9):
        mu = vapour_only_mu(0.2020822595873346, T)
        mc = mu_coex(T)
        winh.append(dict(T=float(T), mu=mu, mu_c=mc, deficit=1 - math.exp((mu - mc) / T)))
    res['heavy_window_undersaturation'] = winh
    # 5. vapour at the decisive-3D point (T = 0.5, mu = Omega(2e5) = 0.8817): densities and the box volume needed
    T = 0.5; mu = 0.8817204919524271
    v = vapour(mu, T)
    box = 24.0 ** 3
    vobj = 4 / 3 * math.pi * 8.4 ** 3
    res['decisive_box'] = dict(T=T, mu=mu, n_net=v['n'], nX=v['nX'], s_vap=v['s'], n_over_s=v['n'] / v['s'],
                               box_volume=box, vapour_volume=box - vobj,
                               vapour_charge_in_box=v['n'] * (box - vobj), vapour_entropy_in_box=v['s'] * (box - vobj),
                               S_over_Q_in_box_for_Q2e5=v['s'] * (box - vobj) / 2e5,
                               volume_needed_for_SQ5=5 * 2e5 / v['s'], L_needed_for_SQ5=(5 * 2e5 / v['s']) ** (1 / 3),
                               volume_needed_for_SQ3=3 * 2e5 / v['s'], L_needed_for_SQ3=(3 * 2e5 / v['s']) ** (1 / 3))
    json.dump(res, open(out, 'w'), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
