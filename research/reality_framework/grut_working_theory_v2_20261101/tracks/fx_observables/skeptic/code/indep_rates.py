"""fx_observables SKEPTIC (Claude, 1 Nov 2026): independent direct-detection rates and the late-channel arithmetic.
[C] own code.

(1) Rates with a Monte-Carlo velocity integral (1e6 SHM samples, Galactic frame -> lab boost, v_esc cut), Helm form
    factor coded from Lewin & Smith (1996) [R], the object form factor from the skeptic's own solve_bvp profile
    (indep_obj_profile.npz).  Same detector toys as the track (LZ WS2024 sigmoid at 5.4 keV, 4.2 t yr; calibration to
    2.2e-48 cm^2 at 40 GeV).  Checks: N90-equivalent, point limit at 17.3 TeV, composite/point ratio, b limit, XLZD
    coupling for 10 events; light-WIMP mimic with a Poisson (Asimov) likelihood ratio instead of Pearson chi^2.
(2) Merger mass conversion per t_H: per volume, collisions (1/2) n^2 sigma v each convert f_rel 2M -> fraction per t_H =
    n sigma v t_H f_rel (the track divides by a further 2).
(3) Hidden cost of X-bar burn-out: each X-bar absorbed on an object releases 2 Omega m_X into the hidden bath,
    Delta(S/Q) = 2 Omega p m_X / T_h.
(4) h -> phi phi freeze-in yield (analytic, K1-weighted number rate) at b = 5e-11 vs the track's 2.237e-8.
"""
import json, math, sys
import numpy as np
from scipy.optimize import minimize_scalar

rng = np.random.default_rng(12345)
C = 299792.458
V0, VESC, VE = 238.0, 544.0, 250.5
RHO, AMU, MN = 0.3, 0.9314941, 0.9389
MX = 0.0924089637375362
GEVKG = 5.60958865e26


def vel_samples(n=2_000_000):
    v = rng.normal(0, V0 / math.sqrt(2), size=(n, 3))
    v = v[np.linalg.norm(v, axis=1) < VESC]
    lab = v - np.array([0, VE, 0])
    return np.linalg.norm(lab, axis=1)


VS = vel_samples()
VS_sorted = np.sort(VS)
INV = 1 / VS_sorted
CUM = np.concatenate([[0], np.cumsum(INV[::-1])])[::-1] / len(VS_sorted)   # CUM[i] = mean of 1/v over v >= VS_sorted[i]


def eta(vmin):
    idx = np.searchsorted(VS_sorted, vmin)
    return CUM[np.minimum(idx, len(VS_sorted))]          # s/km


def helm(q, A):
    hbarc = 0.1973269804
    s = 0.9 / hbarc; a = 0.52 / hbarc; c = (1.23 * A**(1 / 3) - 0.60) / hbarc
    rn = math.sqrt(c * c + 7 / 3 * math.pi**2 * a * a - 5 * s * s)
    x = np.maximum(q * rn, 1e-9)
    j1 = (np.sin(x) - x * np.cos(x)) / x**2
    return 3 * j1 / x * np.exp(-(q * s)**2 / 2)


PROF = None


def objff(q):
    r, S = PROF['r'], PROF['S']
    A = S * S - 1
    w = 4 * math.pi * r * r * A
    F0 = np.trapezoid(w, r)
    qq = np.atleast_1d(q) / MX
    return np.array([np.trapezoid(w * np.sinc(x * r / math.pi), r) for x in qq]) / F0


def dRdE(E_keV, M, A, comp=False):
    mA = A * AMU
    mu = M * mA / (M + mA); mun = M * MN / (M + MN)
    E = E_keV * 1e-6
    q = np.sqrt(2 * mA * E)
    vmin = np.sqrt(mA * E / 2) / mu * C
    sig = A * A * (mu / mun)**2
    F2 = helm(q, A)**2 * (objff(q)**2 if comp else 1.0)
    # events / (kg day keV) per unit sigma_n [cm^2]
    return RHO * sig * F2 * eta(vmin) * C / (2 * M * mu * mu) * (C * 1e5) * GEVKG * 86400 * 1e-6


def eff_lz(E):
    return 0.9 / (1 + np.exp(-(E - 5.4) / 0.9)) / (1 + np.exp((E - 55.0) / 3.0))


def events(M, A, expo, comp=False, eff=eff_lz, Emax=300):
    E = np.linspace(0.5, Emax, 3000)
    return float(np.trapezoid(dRdE(E, M, A, comp) * eff(E), E) * expo * 1000 * 365.25)


def main(prof_fn, out_fn):
    global PROF
    PROF = np.load(prof_fn)
    out = {}
    n40 = events(40.0, 131, 4.2) * 2.2e-48
    M2 = 17325.03
    pt = n40 / events(M2, 131, 4.2)
    co = n40 / events(M2, 131, 4.2, comp=True)
    sig_b = 1.3979e-26
    b_lim = math.sqrt(co / sig_b)
    nX = events(M2, 131, 200.0, comp=True) * sig_b
    out['dd'] = dict(N90_equiv=n40, point_limit_17p3TeV=pt, composite_limit=co, ratio=co / pt, b_LZ=b_lim,
                     b_XLZD_10ev=math.sqrt(10 / nX), track=dict(N90=1.864, point=4.22e-46, comp=5.98e-45, b_LZ=6.54e-10, b_XLZD=2.195e-10))
    print(out['dd'], flush=True)
    # light-WIMP mimic, Xe 4-60 keV in 2 keV bins, Poisson Asimov: data = composite, model = point(m) x free norm
    Eb = np.arange(4, 62, 2.0)
    def binned(M, comp):
        e = np.linspace(4, 62, 581)
        r = dRdE(e, M, 131, comp) * eff_lz(e)
        cs = np.concatenate([[0], np.cumsum(0.5 * (r[1:] + r[:-1]) * np.diff(e))])
        edges = np.interp(np.r_[Eb, 62.0], e, cs)
        return np.diff(edges)
    o = binned(M2, True); o = o / o.sum()
    def q_per_event(lnm):
        p = binned(math.exp(lnm), False); p = p / p.sum()       # free norm -> best norm = total (Poisson)
        return float(2 * np.sum(o * np.log(np.maximum(o, 1e-300) / np.maximum(p, 1e-300)) - o + p))
    best = minimize_scalar(q_per_event, bounds=(math.log(5), math.log(500)), method='bounded', options=dict(xatol=1e-4))
    out['mimic'] = dict(best_point_mass_GeV=math.exp(best.x), q_per_event=best.fun, events_for_3sigma=9 / best.fun,
                        track=dict(mass=28.3, events=30.8))
    print(out['mimic'], flush=True)
    # (2) mergers
    T_H = 4.35e17; cm = 1.973269804e-14 / MX
    sig = math.pi * (2 * 5.943 * cm)**2
    frel = 0.03727
    mer = {}
    for k, (rho, v) in {'MW_solar': (0.4, 220 * math.sqrt(2)), 'MW_inner_1kpc': (8.0, 150 * math.sqrt(2)),
                        'cluster_core': (1.0, 1000 * math.sqrt(2)), 'dwarf_core': (3.8, 15 * math.sqrt(2))}.items():
        G = rho / 17324.83 * sig * v * 1e5 * T_H
        mer[k] = dict(collisions_per_object_per_tH=G, mass_fraction_converted=G * frel, track_value=G * frel / 2)
    out['mergers'] = mer
    print(mer, flush=True)
    # (3) X-bar burn-out entropy
    Om = 0.88172
    out['Xbar_burnout_entropy'] = {('p%g_Th%g' % (p, T)): 2 * Om * p / T for p in (0.0023, 0.005, 0.0143) for T in (0.003, 0.006, 0.01, 0.05)}
    print(out['Xbar_burnout_entropy'])
    # (4) freeze-in yield h -> phi phi
    b = 5e-11; VH, MH, MPL, gs = 246.22, 125.20, 1.22089e19, 106.75
    beta = math.sqrt(1 - 4 * 0.22635**2 / MH**2)
    Gam = b * b * VH * VH * beta / (8 * math.pi * MH)
    Y = 2 * 135 / (8 * math.pi**3 * 1.66 * gs * math.sqrt(gs)) * Gam * MPL / MH**2
    out['freeze_in'] = dict(Gamma_h_phiphi_GeV=Gam, Y_phi=Y, track=2.23656e-8)
    print(out['freeze_in'])
    json.dump(out, open(out_fn, 'w'), indent=1, default=float)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
