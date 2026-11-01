"""fx_observables item 4a: direct detection of M1 objects through a Higgs portal (b |H|^2 phi^2, or the trilinear
A phi |H|^2, or lambda_HX |H|^2 |X|^2).  [C] own code (Claude, 1 Nov 2026).

Physics (derived in NOTES.md section 4):
  the object sources h through J(x) = b v_H (phi^2 - v^2)  (phi^2 portal);  a nucleon has V = (f_N m_N / v_H) h, so
  for q << m_h the object-nucleon amplitude is coherent over the object's scalar source:
      sigma_n(q) = (mu_n^2/pi) (f_N m_N b Ft(q) / m_h^2)^2,   Ft(q) = int d^3x (phi^2 - v^2) e^{iqx}  [GeV^-1]
  (the point-particle limit Ft -> 1/m_s reproduces the standard Higgs-portal scalar formula).  Nucleus:
      dR/dE = rho sigma_n0 A^2 (mu_A/mu_n)^2 |F_Helm(q)|^2 |Ft(q)/Ft(0)|^2 eta(v_min) / (2 M mu_A^2).
  Standard halo (Baxter et al. 2021 conventions, arXiv:2105.00599 [R]): rho = 0.3 GeV/cm^3, v0 = 238, v_esc = 544,
  v_E = 250.5 km/s; Helm form factor (Lewin & Smith 1996 [R]).
Detector models [E]: LZ WS2024 4.2 t yr (arXiv:2410.17036 [Q]: 50 % NR efficiency at 5.4 keV, S1c < 80 phd ROI;
  plateau and the upper roll-off are modelled, E_up bracketed 45/55/65 keV); the limit is converted with N_90 fixed by
  the published 2.2e-48 cm^2 at 40 GeV.  Projections: LZ full 15.3 t yr, XLZD 200 t yr (arXiv:2410.17137 [Q]: 90 %
  sensitivity ~2-3e-49 cm^2 at 40 GeV; scaled the same way); argon 200 t yr with a 30-200 keV window (DarkSide-20k
  class [R]).  Germanium shown for the zero position only.
Outputs results/dd.json.
"""
import json, math, os, sys
import numpy as np
from scipy.special import erf, spherical_jn
from scipy.integrate import quad

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
SP = os.environ.get('SP', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')
PROF = np.load(os.path.join(SP, 'wt2', 'fx_observables', 'profiles.npz'))

MX = 0.0924089637375362
V = 0.924089637375362
VH = 246.22
MH = 125.20
MN = 0.9389
AMU = 0.9314941
GEV2_TO_CM2 = 0.3893794e-27
C_KMS = 299792.458
RHO = 0.3
V0, VESC, VE = 238.0, 544.0, 250.5
SEC_YR = 3.15576e7
CM3_KMS = 1e5                      # km/s -> cm/s
KG_GEV = 5.60958865e26              # GeV per kg
FN_CENTRAL, FN_RANGE = 0.30, (0.26, 0.35)    # Higgs-nucleon coupling [R] (lattice sigma terms; e.g. Hoferichter et al. 2017)


def eta_shm(vmin):
    """mean inverse speed (s/km) for the truncated boosted Maxwellian (Savage, Freese & Gondolo 2006 form)."""
    x, y, z = vmin / V0, VE / V0, VESC / V0
    Nesc = erf(z) - 2 * z * math.exp(-z * z) / math.sqrt(math.pi)
    pre = 1 / (2 * Nesc * V0 * y)
    x = np.asarray(x, float)
    out = np.zeros_like(x)
    a = x < z - y
    out[a] = pre * (erf(x[a] + y) - erf(x[a] - y) - 4 / math.sqrt(math.pi) * y * math.exp(-z * z))
    b = (x >= z - y) & (x < z + y)
    out[b] = pre * (erf(z) - erf(x[b] - y) - 2 / math.sqrt(math.pi) * (y + z - x[b]) * math.exp(-z * z))
    return out


def helm(q_GeV, A):
    """Helm form factor; q in GeV."""
    fm = 1 / 0.1973269804            # GeV^-1 per fm
    s = 0.9 * fm
    c = (1.23 * A**(1 / 3) - 0.60) * fm
    a = 0.52 * fm
    rn = math.sqrt(c * c + 7 / 3 * math.pi**2 * a * a - 5 * s * s)
    x = np.maximum(q_GeV * rn, 1e-9)
    return 3 * spherical_jn(1, x) / x * np.exp(-(q_GeV * s)**2 / 2)


def obj_ff(Qlab, which):
    """normalised object form factor Ft(q)/Ft(0) as a function of q in GeV, from the solved profile (any q)."""
    r, F, S = PROF[Qlab]
    h = r[1] - r[0]
    A = {'phi2': S * S - 1, 'phi1': S - 1, 'X2': F * F, 'point': None}[which]
    F0 = 4 * math.pi * h * np.sum(r * r * A)

    def f(qGeV):
        q = np.atleast_1d(qGeV) / MX
        return np.array([4 * math.pi * h * np.sum(r * r * A * np.sinc(qq * r / math.pi)) for qq in q]) / F0
    return f, F0


def spectrum(E_keV, M, A, ffun):
    """dR/dE per kg per day per keV for sigma_n0 = 1 cm^2."""
    mA = A * AMU
    muA = M * mA / (M + mA)
    mun = M * MN / (M + MN)
    E = E_keV * 1e-6
    q = np.sqrt(2 * mA * E)
    vmin = np.sqrt(mA * E / 2) / muA * C_KMS                      # km/s
    eta = eta_shm(vmin) * C_KMS                                   # dimensionless (units of 1/c)
    sigA = A * A * (muA / mun)**2                                 # cm^2 per unit sigma_n0
    Fo = ffun(q) if ffun is not None else 1.0
    # dR/dE [per GeV per GeV-of-target-mass ...]: rho [GeV/cm^3] sigma [cm^2] / (2 M mu^2 [GeV^3]) eta [1/c] * c
    dRdE = RHO * sigA * helm(q, A)**2 * Fo**2 * eta / (2 * M * muA**2) * (C_KMS * CM3_KMS)   # per s per GeV(target) per GeV(E)
    return dRdE * KG_GEV * 86400 * 1e-6                           # per kg per day per keV


def eff_lz(E, Eup):
    return 0.9 / (1 + np.exp(-(E - 5.4) / 0.9)) / (1 + np.exp((E - Eup) / 3.0))


def events(M, A, ffun, exposure_t_yr, eff, Emax=300.0):
    E = np.linspace(0.5, Emax, 6000)
    s = spectrum(E, M, A, ffun) * eff(E)
    return float(np.trapezoid(s, E) * exposure_t_yr * 1000 * 365.25)


def main():
    out = {'scope': __doc__.strip().splitlines()[0]}
    # --- calibrate LZ: N90 from the published 40 GeV limit
    cal = {}
    for Eup in (45.0, 55.0, 65.0):
        eff = lambda E, Eup=Eup: eff_lz(E, Eup)
        n40 = events(40.0, 131, None, 4.2, eff) * 2.2e-48
        cal[Eup] = dict(N90_equiv=n40)
        # implied point-WIMP limits at higher masses
        for M in (1e3, 1e4, 17324.83):
            cal[Eup]['point_limit_%g' % M] = n40 / events(M, 131, None, 4.2, eff)
    out['lz_calibration'] = cal
    print('LZ calibration', json.dumps(cal, indent=0))
    # --- objects
    labs = {'Q1e+05': 1e5, 'Q2e+05': 2e5, 'Q5e+05': 5e5, 'Q1e+06': 1e6, 'Q1e+07': 1e7, 'Q1e+08': 1e8}
    rows = []
    for lab, Q in labs.items():
        for which in ('phi2', 'phi1', 'X2'):
            ff, F0 = obj_ff(lab, which)
            r, F, S = PROF[lab]
            hh = r[1] - r[0]
            Om_EQ = None
            # mass from the object table
            obj = json.load(open(os.path.join(RES, 'objects.json')))
            row = [z for z in obj['rows'] if abs(z['Q'] - Q) / Q < 1e-6][0]
            M = row['mass_GeV']
            # Ft(0) in GeV^-1 for the three portals (couplings b, A_tri [GeV], lambda_HX):
            if which == 'phi2':
                Ft0 = V * V / MX**3 * F0                      # int (phi^2 - v^2) d^3x
                # sigma_n0 = (mu^2/pi) (fN mN b Ft0 / mh^2)^2
                coup = lambda fN: (MN * M / (M + MN))**2 / math.pi * (fN * MN * Ft0 / MH**2)**2 * GEV2_TO_CM2   # per b^2
                coup_label = 'b'
            elif which == 'phi1':
                Ft0 = V / MX**3 * F0                          # int (phi - v) d^3x  [GeV^-2]
                coup = lambda fN: (MN * M / (M + MN))**2 / math.pi * (fN * MN * Ft0 / (VH * MH**2))**2 * VH**2 / VH**2 * GEV2_TO_CM2
                # A phi |H|^2 -> A v_H h phi : source A v_H (phi - v); same algebra with b v_H -> A v_H / ... (A in GeV)
                coup = lambda fN: (MN * M / (M + MN))**2 / math.pi * (fN * MN * Ft0 / MH**2)**2 * GEV2_TO_CM2   # per (A/GeV)^2
                coup_label = 'A_tri [GeV]'
            else:
                Ft0 = V * V / (2 * MX**3) * F0                # int |X|^2 d^3x = (v^2/2) int F^2 d^3rho / m_X^3
                coup = lambda fN: (MN * M / (M + MN))**2 / math.pi * (fN * MN * Ft0 / MH**2)**2 * GEV2_TO_CM2   # per lambda_HX^2
                coup_label = 'lambda_HX'
            sig_per_c2 = coup(FN_CENTRAL)
            res = dict(Q=Q, portal=which, mass_GeV=M, Ft0_GeVinv=Ft0, coupling=coup_label,
                       sigma_n0_per_coupling2_cm2=sig_per_c2,
                       sigma_n0_per_coupling2_cm2_fN_range=[coup(FN_RANGE[0]), coup(FN_RANGE[1])])
            # LZ-equivalent limits on sigma_n0 and on the coupling
            lim = {}
            for Eup in (45.0, 55.0, 65.0):
                eff = lambda E, Eup=Eup: eff_lz(E, Eup)
                nper = events(M, 131, ff, 4.2, eff)
                s_lim = cal[Eup]['N90_equiv'] / nper
                lim['Eup%g' % Eup] = dict(sigma_n0_limit=s_lim, coupling_limit=math.sqrt(s_lim / sig_per_c2),
                                          ratio_to_point_limit=s_lim / cal[Eup]['point_limit_17324.8'] if Q == 2e5 else None)
            res['LZ2024_equiv'] = lim
            # projections: 10 signal events in the window
            eff55 = lambda E: eff_lz(E, 55.0)
            proj = {}
            for name, A, expo, eff in (('LZ_full_15.3tyr', 131, 15.3, eff55), ('XLZD_200tyr', 131, 200.0, eff55),
                                       ('Ar_200tyr_30-200keV', 40, 200.0, lambda E: 0.7 / (1 + np.exp(-(E - 30) / 2)) / (1 + np.exp((E - 200) / 5)))):
                n1 = events(M, A, ff, expo, eff) * sig_per_c2      # events per coupling^2
                proj[name] = dict(events_per_coupling2=n1, coupling_for_10_events=math.sqrt(10 / n1) if n1 > 0 else None)
            res['projections_10_events'] = proj
            rows.append(res)
            print('%s %-5s M=%.4g Ft0=%.4g sig/c^2=%.3e  LZ lim(55): sig=%.3e coup=%.3e   10ev XLZD coup=%.3e Ar=%.3e' % (
                lab, which, M, Ft0, sig_per_c2, lim['Eup55']['sigma_n0_limit'], lim['Eup55']['coupling_limit'],
                proj['XLZD_200tyr']['coupling_for_10_events'], proj['Ar_200tyr_30-200keV']['coupling_for_10_events'] or -1), flush=True)
    out['rows'] = rows
    # --- spectra for Q = 2e5 (phi2) vs point-like WIMP of the same mass, Xe and Ar; and a Q spread
    E = np.linspace(1, 120, 239)
    ff2, _ = obj_ff('Q2e+05', 'phi2')
    M2 = 17324.83
    spec = dict(E_keV=E.tolist(),
                Xe_obj=(spectrum(E, M2, 131, ff2)).tolist(), Xe_point=(spectrum(E, M2, 131, None)).tolist())
    Ea = np.linspace(5, 300, 296)
    spec['Ar_E_keV'] = Ea.tolist()
    spec['Ar_obj'] = spectrum(Ea, M2, 40, ff2).tolist(); spec['Ar_point'] = spectrum(Ea, M2, 40, None).tolist()
    out['spectra_Q2e5_phi2_per_kg_day_keV_per_cm2'] = spec
    # zero-depth with a lognormal spread of Q (R ~ Q^1/3, q0 ~ 1/R): rate-weighted average of |F|^2 over Q
    # (approximate the shape change by scaling the Q = 2e5 form factor: F_Q(q) = F_2e5(q R(Q)/R(2e5)))
    spread = {}
    for sig_lnQ in (0.0, 0.3, 0.7, 1.5):
        lnQ = np.linspace(-4, 4, 81) * max(sig_lnQ, 1e-9)
        w = np.exp(-lnQ**2 / (2 * max(sig_lnQ, 1e-9)**2)); w /= w.sum()
        qz = 0.08202
        qq = np.linspace(0.3, 2.0, 341) * qz
        Fm = np.zeros_like(qq)
        for lq, ww in zip(lnQ, w):
            sc = math.exp(lq / 3)
            # source scales as volume ~ Q: weight by (Q/2e5)^2 for the coherent rate, and by 1/Q for number density at fixed rho
            Fm += ww * math.exp(lq) * ff2(qq * sc)**2
        Fm /= np.sum(w * np.exp(lnQ))
        spread['sigma_lnQ_%g' % sig_lnQ] = dict(min_over_q_of_F2=float(Fm.min()), q_at_min_GeV=float(qq[np.argmin(Fm)]),
                                               F2_at_half_q0=float(np.interp(0.5 * qz, qq, Fm)),
                                               F2_at_1p5_q0=float(np.interp(1.5 * qz, qq, Fm)))
    out['zero_smearing_Q2e5_phi2'] = spread
    print('zero smearing', json.dumps(spread, indent=0))
    # --- events needed to distinguish the zero (Xe, 4-60 keV, 2 keV bins), Asimov, normalisation free
    Eb = np.arange(4, 62, 2.0)
    eff = eff_lz(Eb + 1, 55.0)
    so = np.array([quad(lambda e: spectrum(np.array([e]), M2, 131, ff2)[0], a, a + 2)[0] for a in Eb]) * eff
    sp_ = np.array([quad(lambda e: spectrum(np.array([e]), M2, 131, None)[0], a, a + 2)[0] for a in Eb]) * eff
    po, pp = so / so.sum(), sp_ / sp_.sum()
    # chi^2 per event for 'point' shape fitted to composite data (free norm): sum (po - pp)^2/pp  -> N needed for 9
    chi1 = float(np.sum((po - pp)**2 / pp))
    out['zero_discrimination'] = dict(chi2_per_event_point_vs_composite=chi1, events_for_3sigma=9 / chi1,
                                      fraction_of_events_in_20to35keV_obj=float(po[(Eb >= 20) & (Eb < 35)].sum()),
                                      fraction_in_20to35keV_point=float(pp[(Eb >= 20) & (Eb < 35)].sum()),
                                      note='background-free Asimov, shape only; a real analysis with ER/neutron backgrounds needs more')
    print('discrimination', out['zero_discrimination'])
    # --- neutron-star kinematics: coherence at q ~ 0.3-1 GeV
    ns = {}
    for lab in ('Q2e+05', 'Q1e+08'):
        f, _ = obj_ff(lab, 'phi2')
        ns[lab] = {('q%.1fGeV' % q): float(f(np.array([q]))[0]**2) for q in (0.1, 0.3, 0.5, 1.0)}
    out['NS_coherence_|F/F0|^2'] = ns
    print('NS', ns)
    with open(os.path.join(RES, 'dd.json'), 'w') as fo:
        json.dump(out, fo, indent=1, default=float)


if __name__ == '__main__' and len(sys.argv) == 1:
    main()


def discrimination_vs_light_wimp():
    """Asimov test: composite Q = 2e5 (phi^2 portal) data fitted by a point-like WIMP of free mass and free
    normalisation, Xe alone and Xe + Ar jointly (same sigma_n per target, i.e. one normalisation).  Background-free;
    bins 2 keV (Xe 4-60 keV, LZ-like efficiency) and 5 keV (Ar 30-200 keV).  Returns events needed (in Xe) for
    Delta chi^2 = 9, with the Ar exposure set to the XLZD/Ar 200 t yr ratio."""
    from scipy.optimize import minimize_scalar
    ff2, _ = obj_ff('Q2e+05', 'phi2')
    M2 = 17324.83
    EbX = np.arange(4, 62, 2.0); EbA = np.arange(30, 200, 5.0)
    def binned(M, A, ffun, Eb, w, eff):
        return np.array([quad(lambda e: spectrum(np.array([e]), M, A, ffun)[0] * eff(np.array([e]))[0], a, a + w, limit=50)[0] for a in Eb])
    effX = lambda E: eff_lz(E, 55.0)
    effA = lambda E: 0.7 / (1 + np.exp(-(E - 30) / 2)) / (1 + np.exp((E - 200) / 5))
    oX = binned(M2, 131, ff2, EbX, 2, effX); oA = binned(M2, 40, ff2, EbA, 5, effA)
    out = {}
    for lab, useA in (('Xe_only', False), ('Xe_plus_Ar', True)):
        def chi2_per_event(lnm):
            m = math.exp(lnm)
            pX = binned(m, 131, None, EbX, 2, effX); pA = binned(m, 40, None, EbA, 5, effA)
            o = np.r_[oX, oA] if useA else oX
            p = np.r_[pX, pA] if useA else pX
            o = o / oX.sum(); p_norm = p / pX.sum()
            # best normalisation for chi^2 = sum (o - c p)^2 / (c p) (Asimov, Neyman-like): c = sqrt(sum o^2/p / sum p)
            c = math.sqrt(np.sum(o**2 / np.maximum(p_norm, 1e-300)) / np.sum(p_norm))
            return float(np.sum((o - c * p_norm)**2 / np.maximum(c * p_norm, 1e-300)))
        best = minimize_scalar(chi2_per_event, bounds=(math.log(5), math.log(2e4)), method='bounded', options=dict(xatol=1e-3))
        out[lab] = dict(best_point_mass_GeV=math.exp(best.x), chi2_per_Xe_event=best.fun, Xe_events_for_3sigma=9 / best.fun if best.fun > 0 else None,
                        Ar_to_Xe_event_ratio=float(oA.sum() / oX.sum()))
    return out


def zero_contrast():
    """Depth of the first zero for a lognormal spread in Q: |F|^2 averaged with weight Q^2 n(Q) (coherent rate at
    fixed mass density ~ Q n(Q)), using the exact solved profiles' scaling F_Q(q) = F_2e5(q R_Q/R_2e5)."""
    ff2, _ = obj_ff('Q2e+05', 'phi2')
    q0 = 0.08202
    qq = np.linspace(0.2, 3.0, 561) * q0
    single = ff2(qq)**2
    i2 = np.argmax(np.where(qq > q0, single, 0))
    res = {'single_secondary_max_over_peak': float(single[i2]), 'q_secondary_max_GeV': float(qq[i2])}
    for s in (0.1, 0.3, 0.7, 1.5):
        lnQ = np.linspace(-4 * s, 4 * s, 161)
        w = np.exp(-lnQ**2 / (2 * s * s)) * np.exp(lnQ)
        w /= w.sum()
        Fm = sum(ww * ff2(qq * math.exp(l / 3))**2 for l, ww in zip(lnQ, w))
        j0 = np.argmin(np.abs(qq - q0))
        res['sigma_lnQ_%g' % s] = dict(F2_at_q0=float(Fm[j0]), F2_secondary_max=float(Fm[i2]),
                                       contrast_secmax_over_q0=float(Fm[i2] / Fm[j0]),
                                       has_local_min=bool(np.any((np.diff(np.sign(np.diff(Fm))) > 0))))
    return res


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == 'extra':
    o = json.load(open(os.path.join(RES, 'dd.json')))
    o['discrimination_vs_light_point_wimp'] = discrimination_vs_light_wimp()
    print(o['discrimination_vs_light_point_wimp'])
    o['zero_contrast'] = zero_contrast()
    print(o['zero_contrast'])
    with open(os.path.join(RES, 'dd.json'), 'w') as fo:
        json.dump(o, fo, indent=1, default=float)
