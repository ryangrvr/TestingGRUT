"""Observables of a background history with walls: H0 given T_CMB, age, acceleration intervals, wall fractions.

Used by history_benchmark.py and scale_family.py. Deterministic, numpy/scipy only.
"""
import math
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
import core as C

T_BBN = (1e-3, 1e-4)                     # GeV: 1 MeV and 0.1 MeV
T_REC = C.T0 * (1.0 + C.Z_STAR)          # photon temperature at z_* (standard-history value, used as a temperature)

def H_kmsmpc(H_gev):
    return H_gev / C.kmsmpc_to_gev(1.0)

def rho_total(h, T):
    return C.rho_rad(T) + C.rho_matter(T) + C.wall_density_at(h, T) + h['rho_L_GeV4']

def rho_plus_3p(h, T):
    """rho + 3p of the total; walls: phase 1 w = 1/(3 H t) - 1 (from rho_w = A sigma / t), phase 2 w = n/3 - 1."""
    rw = C.wall_density_at(h, T)
    if h['switched'] and T <= h['T_switch_GeV']:
        one_plus_3w = h['n_post'] - 2.0
    else:
        t = C.time_at(h, T)
        H = math.sqrt(rho_total(h, T) / (3 * C.M_P**2))
        one_plus_3w = 1.0 / (H * t) - 2.0
    return 2 * C.rho_rad(T) + C.rho_matter(T) + one_plus_3w * rw - 2 * h['rho_L_GeV4']

def z_of_T(T):
    return 1.0 / C.a_ratio(T) - 1.0

def acceleration_intervals(h, npts=3000):
    """Redshift intervals [z_hi, z_lo] in which a-double-dot > 0 between T_start and T0."""
    Ts = np.exp(np.linspace(math.log(h['T_start_GeV']) - 1e-9, math.log(C.T0) + 1e-9, npts))
    if h['switched']:
        # make sure the switch point is bracketed cleanly
        Ts = np.sort(np.concatenate([Ts, [h['T_switch_GeV'] * (1 + 1e-9), h['T_switch_GeV'] * (1 - 1e-9)]]))[::-1]
    f = np.array([rho_plus_3p(h, T) for T in Ts])
    acc = f < 0
    intervals = []
    i = 0
    while i < len(Ts):
        if acc[i]:
            j = i
            while j + 1 < len(Ts) and acc[j + 1]:
                j += 1
            T_on = Ts[i] if i == 0 else _refine(h, Ts[i - 1], Ts[i])
            T_off = Ts[j] if j == len(Ts) - 1 else _refine(h, Ts[j], Ts[j + 1])
            intervals.append({'T_begin_GeV': T_on, 'T_end_GeV': T_off, 'z_begin': z_of_T(T_on),
                              'z_end': max(z_of_T(T_off), 0.0)})
            i = j + 1
        else:
            i += 1
    return intervals

def _refine(h, Ta, Tb):
    g = lambda lx: rho_plus_3p(h, math.exp(lx))
    try:
        return math.exp(brentq(g, math.log(Tb), math.log(Ta), xtol=1e-12))
    except ValueError:
        return math.sqrt(Ta * Tb)

def summarize(h, H0_meas=(67.0, 70.0, 74.0), with_accel=True):
    T0 = C.T0
    rho0 = rho_total(h, T0)
    H0 = math.sqrt(rho0 / (3 * C.M_P**2))
    H0k = H_kmsmpc(H0)
    age_s = C.time_at(h, T0) * C.HBAR_GEV_S
    out = {
        'switched': h['switched'],
        'T_switch_GeV': h['T_switch_GeV'],
        'z_switch': z_of_T(h['T_switch_GeV']) if h['switched'] else None,
        't_switch_s': h['t_switch_GeV-1'] * C.HBAR_GEV_S if h['switched'] else None,
        'rho_w0_GeV4': C.wall_density_at(h, T0),
        'Omega_w0': C.wall_density_at(h, T0) / rho0,
        'H0_pred_kmsMpc': H0k,
        'R_H': {('%g' % m): H0k / m for m in H0_meas},
        'age_pred_s': age_s,
        'age_pred_yr': age_s / C.JULIAN_YEAR_S,
        'q0': rho_plus_3p(h, T0) / (2 * rho0),
    }
    fr = {}
    for T in list(T_BBN) + [T_REC]:
        rw = C.wall_density_at(h, T)
        rt = rho_total(h, T)
        key = 'T=%.4g_GeV' % T
        fr[key] = {'rho_w/rho_tot': rw / rt, '1-Omega_w': (rt - rw) / rt, 'rho_w/rho_rad': rw / C.rho_rad(T),
                   'DeltaNeff_equiv': rw / C.rho_one_nu(T), 'Tnu_over_T': C.Tnu_over_T(T)}
    out['fractions'] = fr
    if with_accel:
        out['acceleration_intervals'] = acceleration_intervals(h)
    return out

def tuned_negative_vacuum(h, H0_target=70.0):
    """Choose rho_L < 0 so that H(T0) = H0_target exactly; return rho_L, q0, age. Analytic in phase 2 to avoid
    the 1-part-in-R_H^2 cancellation (computed with expm1)."""
    assert h['switched']
    n = h['n_post']
    rhoc = 3 * C.M_P**2 * C.kmsmpc_to_gev(H0_target)**2
    rw0 = h['rho_w2'](C.T0)
    rr0, rm0 = C.RHO_R0, C.RHO_M0
    rhoL = rhoc - rw0 - rr0 - rm0
    q0 = (n * rw0 + 4 * rr0 + 3 * rm0 - 2 * rhoc) / (2 * rhoc)
    # Age: t_switch + int_{a_sw}^{a0} dln a / H, with
    # 3 M_P^2 H^2 = rw0 expm1(n L) + rr0 expm1(4 L) + rm0 expm1(3 L) + rhoc, L = ln(a0/a) (g_s const below keV).
    Lsw = -math.log(h['a_switch_over_a0'])
    def invH(L):
        val = rw0 * math.expm1(n * L) + rr0 * math.expm1(4 * L) + rm0 * math.expm1(3 * L) + rhoc
        return 1.0 / math.sqrt(val / (3 * C.M_P**2))
    # split integral for accuracy near L = 0 where H is smallest
    pts = [0.0, 1e-20, 1e-16, 1e-12, 1e-8, 1e-4, 1e-2, 1.0, Lsw]
    pts = [p for p in pts if p <= Lsw]
    tot = 0.0
    for a, b in zip(pts[:-1], pts[1:]):
        val, _ = quad(invH, a, b, limit=400, epsrel=1e-10)
        tot += val
    age = (h['t_switch_GeV-1'] + tot) * C.HBAR_GEV_S
    return {'H0_target_kmsMpc': H0_target, 'rho_L_GeV4': rhoL, 'rho_L_over_rho_w0': rhoL / rw0,
            'fine_tuning_1_minus_abs_ratio': (rhoc - rr0 - rm0) / rw0, 'q0': q0, 'age_s': age,
            'age_yr': age / C.JULIAN_YEAR_S}
