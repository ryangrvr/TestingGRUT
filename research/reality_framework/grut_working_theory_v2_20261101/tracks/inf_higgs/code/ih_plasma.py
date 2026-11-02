#!/usr/bin/env python3
"""Track inf_higgs: when does a Higgs-coupled plasma appear on the inflaton benchmarks? (own code; model of
inf_reheat's benchmarks.py: plasma_solve, transcribed onto the benchmarks' actual oscillation instead of a matter-era
fit H = 1/(1 + 1.5 tau) and a harmonic amplitude).

Channel S1 of inf_reheat (every element a model choice): -(y/2) phi_I Nbar N, m_N(t) = M_N + y phi_I(t), seesaw
y_nu^2 = mt M_N / v^2 (mt = 0.05 eV, v = 174.1 GeV), N -> L H at Gamma_N = y_nu^2 |m_N(t)| / (8 pi).
  * Greene-Kofman (instant preheating, Felder-Kofman-Linde hep-ph/9812289) at every zero of m_N(t) on the benchmark's
    own anharmonic oscillation: n = (y |phidot_0|)^{3/2} / (4 pi^3) Majorana N per crossing (inf_reheat's formula,
    checked there numerically to 1e-3); Pauli: a crossing only refills up to that occupation (n_N -> max(n_N, n)).
  * Perturbative phi_I -> N N at y^2 m_I beta^3 / (16 pi) times the open fraction (2/pi) arcsin(min(1, (m_I/2)/(y Phi)))
    (Phi from the turning points), N at energy m_I/2.
  * SM sector: number and energy densities of the hard L, H quanta (E = m_N/2 each), radiation redshift; hard-gas
    induced Higgs mass m_ind^2 = c rho_SM / E^2 with E = rho_SM/n_SM and c = 0.0228, normalised so that a thermal SM
    gas (E = 2.70 T, g* = 106.75) gives kappa_T = 0.11 [E; inf_reheat's normalisation].  Composition caveat: a gas of
    L and H alone couples to the Higgs only through g^2, g'^2 (and lam < 0); gauge bosons and tops appear on the
    splitting time [A].
  * Thermalisation, pessimistic (inf_reheat's LPM-type): t_th = alpha^{-2} (E/T_*)^{1/2} / T_*, alpha = 1/40,
    T_* = (30 rho_SM / (pi^2 g*))^{1/4}; Harigaya-Mukaida 1312.3097 [Q] argue it is near-instantaneous.
  * Induced R-odd A phi_I |H|^2 of S1: A_ind ~ y y_nu^2 M_N/(8 pi^2) x [1, ln(M_P^2/M_N^2)] (inf_reheat skeptic) against
    the inflationary bound of ih_inflation.py Aterm.
R-even portal (lam_IH/2) phi_I^2 |H|^2 (m_H^2 = lam phi_I^2/2): lam_min = (20/3)(m_I/M_P)^2 keeps <m_H^2> >= 10 H^2 in
the matter era; q = lam Phi_1^2/(8 m_I^2); condensate annihilation phi phi -> H H* (inf_reheat skeptic N9) as a source.
Writes ../results/plasma.json; the induced-mass curves are then put into the danger function (ih_final.py).
"""
import sys
sys.dont_write_bytecode = True
import json
import math

import numpy as np

import ih_common as C

MNU = 0.05e-9
VEW = 174.1
C_HARD = 0.11 * 2.70 ** 2 / (math.pi ** 2 * C.GSTAR / 30)
ALPHA_G = 1 / 40.0


def y_nat(name, iH=0):
    """radiative flatness y_nat = (0.008 pi^2 V_* / phi_I*^2)^{1/4} (inf_reheat R10), phi_I at CMB exit."""
    B = C.load_bench(name, iH)
    m = B['model']
    bg = C.IC.background(m, Nwant=75.0)
    f, fN = C.IC.bg_at(bg, B['Nstar'])
    fI = f - B['phimin']
    V = B['rec']['V0_MP4'] * m.v(f)
    return (0.008 * math.pi ** 2 * V / fI ** 2) ** 0.25, fI


def solve_full(bg, B, y, MN_frac=0.25, Nf=1, stride=20):
    """S1 plasma on the fine background table of bg (H_i units); output thinned by stride."""
    F = bg.fine()
    t, H, fI, fd = F['tau'], F['H'], F['phiI'], F['phid']
    eta = B['eta']; mI = B['m_I']
    MN = MN_frac * mI
    MNg = MN * B['H_GeV']
    ynu2 = MNU * MNg / VEW ** 2
    tt, pt = F['turn_tau'], np.abs(F['turn_phiI'])
    i0 = int(np.searchsorted(t, bg.tau_endinf - 0.05))
    beta = math.sqrt(max(1 - 4 * MN_frac ** 2, 0.0))
    Phis = np.interp(t, tt, pt) if len(tt) else np.abs(fI)
    nN = rS = nS = 0.0
    mN_prev = MN + y * fI[i0] / eta
    rows = []
    ncross = 0
    first = None
    rN = 0.0
    for i in range(i0 + 1, len(t)):
        dt = t[i] - t[i - 1]
        Hh = H[i]
        mN = MN + y * fI[i] / eta
        if mN * mN_prev < 0:
            ncr = Nf * (y * abs(fd[i]) / eta) ** 1.5 / (4 * math.pi ** 3)
            nN = max(nN, ncr)
            ncross += 1
            first = t[i] if first is None else first
        mN_prev = mN
        amN = abs(mN)
        GN = ynu2 * amN / (8 * math.pi)
        yPhi = y * Phis[i] / eta
        openf = (2 / math.pi) * math.asin(min(1.0, (mI / 2) / yPhi)) if yPhi > 0 else 1.0
        rho_phi = (0.5 * fd[i] ** 2 + bg.V(fI[i] + bg.phimin)) / eta ** 2
        Gp = Nf * y * y * mI * beta ** 3 / (16 * math.pi) * openf
        dec = GN * nN
        nN += (-3 * Hh * nN - dec + 2 * Gp * rho_phi / mI) * dt
        rS += (-4 * Hh * rS + dec * max(amN, 1e-6)) * dt
        nS += (-3 * Hh * nS + 2 * dec) * dt
        rN = nN * max(amN, MN)
        if (i - i0) % stride == 0:
            rows.append((t[i], rS, nS, nN, rho_phi, rN, Gp, openf, yPhi))
    R = np.array(rows)
    tt_, rS_, nS_ = R[:, 0], R[:, 1], R[:, 2]
    E = np.where(nS_ > 0, rS_ / np.maximum(nS_, 1e-300), np.inf)
    m2 = np.where(nS_ > 0, C_HARD * rS_ / E ** 2, 0.0)
    Ts = (30 * np.maximum(rS_, 0) / (math.pi ** 2 * C.GSTAR)) ** 0.25
    tth = np.where(Ts > 0, ALPHA_G ** -2 * np.sqrt(np.maximum(E / np.maximum(Ts, 1e-300), 1.0)) / np.maximum(Ts, 1e-300), np.inf)
    rel = R[:, 1] / np.maximum(R[:, 4], 1e-300)
    return dict(t=tt_, m2=m2, T=Ts, E=E, tth=tth, rS_over_rphi=rel, rN_over_rphi=R[:, 5] / np.maximum(R[:, 4], 1e-300),
                Gp=R[:, 6], openf=R[:, 7], yPhi=R[:, 8], ncross=ncross, first_cross=first, MN=MN, ynu2=ynu2)


def first_time(t, x, thr):
    i = np.where(x >= thr)[0]
    return float(t[i[0]]) if len(i) else None


def thermal_time(r):
    """first time the would-be bath is thermal by the pessimistic estimate: t - t_first >= t_th"""
    t, tth = r['t'], r['tth']
    t0 = r['first_cross'] if r['first_cross'] else t[0]
    i = np.where((t - t0) >= tth)[0]
    return float(t[i[0]]) if len(i) else None


def A_induced(y, ynu2, MN_GeV):
    base = y * ynu2 * MN_GeV / (8 * math.pi ** 2)
    return base, base * math.log(C.MPL ** 2 / MN_GeV ** 2)


def main():
    out = {}
    path = C.RES / 'plasma.json'
    Ainf = json.load(open(C.RES / 'inflation_Aterm.json'))
    for name in C.BENCH:
        B = C.load_bench(name, 0)
        bg = C.Background(B, gam=0.0, Nb0=2.0, tau_end=60.0)
        yn, fIs = y_nat(name)
        o = dict(y_nat=yn, phiI_star=fIs, m_I=B['m_I'], H_GeV=B['H_GeV'], rows={})
        mI = B['m_I']
        # R-even portal window (convention m_H^2 = lam phi^2/2)
        mM = mI * B['eta']
        F = bg.fine()
        Phi1 = float(np.max(np.abs(F['turn_phiI'][:2]))) if len(F['turn_phiI']) else None
        lam_min = (20 / 3) * mM ** 2
        o['lambda_IH'] = dict(lam_min=lam_min, q_at_lam_min=lam_min * Phi1 ** 2 / (8 * mM ** 2), Phi1=Phi1,
                              lam_max_q1=8 * mM ** 2 / Phi1 ** 2,
                              portal_mass2_at_S1_over_Hi2=lam_min * B['rec']['phiI_S1_MP'] ** 2 / (2 * B['eta'] ** 2))
        # condensate annihilation phi phi -> H H* through the portal at lam_min (skeptic N9 formula), hard gas E = m_I
        G_ann = lam_min ** 2 * (Phi1 / B['eta']) ** 2 / (16 * math.pi * mI)
        rho1 = 3 * 0.9 ** 2 / B['eta'] ** 2
        rS1 = G_ann * rho1 * 1.0
        o['lambda_IH']['annihilation'] = dict(Gamma_eff=G_ann, rhoSM_first_Hubble_time=rS1,
                                              m_ind2_first_Hubble_time=C_HARD * rS1 / mI ** 2,
                                              T_star=(30 * rS1 / (math.pi ** 2 * C.GSTAR)) ** 0.25)
        for MNf in (0.25, 0.0625):
            for yv in (1e-5, 3e-5, 1e-4, yn, 1e-3, 3e-3):
                r = solve_full(bg, B, yv, MN_frac=MNf)
                MNg = r['MN'] * B['H_GeV']
                A1, A2 = A_induced(yv, r['ynu2'], MNg)
                Amax = Ainf[name]['172.86']['tree A_max_GeV'] if name in Ainf else None
                row = dict(y=yv, y_over_ynat=yv / yn, MN_over_mI=MNf, MN_GeV=MNg, ynu2=r['ynu2'],
                           Gamma_pert_open=yv ** 2 * mI * (1 - 4 * MNf ** 2) ** 1.5 / (16 * math.pi),
                           yPhi1_over_half_mI=float(r['yPhi'][0] / (mI / 2)), ncross=r['ncross'], first_cross=r['first_cross'],
                           tau_m2_1em3=first_time(r['t'], r['m2'], 1e-3), tau_m2_1em2=first_time(r['t'], r['m2'], 1e-2),
                           tau_m2_1=first_time(r['t'], r['m2'], 1.0), tau_thermal_pessimistic=thermal_time(r),
                           m2_at={f'{x}': float(np.interp(x, r['t'], r['m2'])) for x in (1.5, 2, 3, 5, 8, 12, 20)},
                           T_at={f'{x}': float(np.interp(x, r['t'], r['T'])) for x in (2, 5, 8, 20)},
                           E_at_5=float(np.interp(5.0, r['t'], np.where(np.isfinite(r['E']), r['E'], 0))),
                           rS_over_rphi_at_8=float(np.interp(8.0, r['t'], r['rS_over_rphi'])),
                           max_rN_over_rphi=float(np.max(r['rN_over_rphi'])),
                           tau_unblock=first_time(r['t'], -r['yPhi'] + mI / 2, 0.0),
                           A_ind_GeV=[A1, A2], A_max_inflation_GeV_Mt172_86=Amax,
                           curve=dict(t=r['t'][::5].tolist(), m2=r['m2'][::5].tolist()))
                o['rows'][f'MN={MNf} y={yv:.3g}'] = row
                print(name, f'MN={MNf}', 'y %.2e (y/ynat %.2f)' % (yv, yv / yn), 'ncross', r['ncross'],
                      'tau(m2>=1e-3) %s (1e-2) %s thermal %s unblock %s' % (row['tau_m2_1em3'], row['tau_m2_1em2'],
                                                                         row['tau_thermal_pessimistic'], row['tau_unblock']),
                      'm2@8 %.3g rN/rphi max %.2g A_ind %.3g-%.3g GeV (max %s)' % (row['m2_at']['8'], row['max_rN_over_rphi'],
                                                                              A1, A2, Amax), flush=True)
        out[name] = o
        C.jdump(out, C.SP / 'plasma_full.json')           # full curves (scratchpad; > 2 MB)
        thin = json.loads(json.dumps(out, default=float))
        for nn in thin:
            for rr in thin[nn]['rows'].values():
                t, m2 = np.array(rr['curve']['t']), np.array(rr['curve']['m2'])
                tt = np.unique(np.concatenate([np.linspace(t[0], min(t[0] + 3, t[-1]), 60),
                                               np.geomspace(min(t[0] + 3, t[-1]), t[-1], 60)]))
                rr['curve'] = dict(t=tt.tolist(), m2=np.interp(tt, t, m2).tolist(),
                                   note='thinned to 120 points; full curve in the scratchpad plasma_full.json')
        C.jdump(thin, path)


if __name__ == '__main__':
    main()
