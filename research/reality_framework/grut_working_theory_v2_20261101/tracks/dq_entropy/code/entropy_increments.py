#!/usr/bin/env python3
"""dq_entropy item 3 (Claude, own code): entropy produced after the relaxed pre-onset plasma, per net charge.

(a) phi crossover.  The tilt leaves one phi root at every density (fx_chain [D]); along the homogeneous branch we compute the
    phi curvature m_eff^2 = d^2 U_eff/d phi^2 at phi_eq (T = 0 and with the onset thermal masses scaled as 1/a^2) and the
    adiabaticity H/m_eff.  A quasi-static crossover is isentropic; the coherent excitation is bounded by (H/m_eff)^2.
(b) spinodal decomposition = the latent heat of the supercooled (metastable 449 -> 233 m_X^3) homogeneous condensate, in
    which phi relaxes from phi_eq(n_sp) to ~v between the lumps:  Delta(S/Q) = f_c * [e_hom(n_sp)/n_sp - Omega_b] / T_sp
    (first order in the release; the gas heat capacity ~ 10^2 x the release).  [C] at T = 0 with the m1_vacuum/fx_chain
    homogeneous functional (heat_onset.hom_state, imported read-only); [E] with the frozen thermal masses, where the vacuum
    reference e_vac(T) and the bulk tangent P_rel = n Omega - (e - e_vac) = 0 are recomputed.
    Then the surface (Gibbs-Thomson) free energy of the first lumps, dissipated by ripening/evaporation (fx_formation fit),
    and the kinetic lag of evaporation (net/gross Hertz-Knudsen, Gamma_ex/H >= 1e5 [Q fx_formation]).
(c) Xbar burn-out: absorption of Xbar on lumps in the two-phase era is part of chemical equilibrium (mu_Xbar = -omega =
    -mu_X), hence reversible; an out-of-equilibrium remnant p absorbed at T would cost (m_X + Omega) p / T (fx_observables
    skeptic used 2 Omega).  Evaluated as a hypothetical.
Usage: entropy_increments.py ../results/sq_scan.json ../results/entropy_increments.json
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.optimize import brentq, minimize_scalar
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1].parent / 'fx_chain' / 'code'))
import heat_onset as HO          # read-only: hom_state, hom_at_n, hom_scan, U_phi, constants (m1_vacuum conventions)
import hartree_dq as HQ

LX, G2, LP = 1e-3, 1e-2, 3e-2
EFIT = (0.775865, 8.1429, 125.79, -2847.0)          # fx_formation family fit E/Q = Ob + A Q^-1/3 + B Q^-2/3 + C/Q
MPL = 2.435e18
YQ_FREEX = 4.756e-9                                   # free-X target yield (wt_critic), for the clock only


def thermal_terms(st, mc):
    """onset thermal inputs in m_X units (heat_onset conventions)."""
    tphi = (G2 * st['D'] + 3 * LP * st['P']) / mc ** 2
    tX = (2 * LX * st['D'] + G2 * st['P']) / mc ** 2
    tlin = 3 * HO.C_ * G2 * st['P'] / mc ** 2
    return tphi, tX, tlin


def e_vac(tphi, tlin):
    f = lambda S: HO.U_phi(S) + 0.5 * tphi * S * S - tlin * (S - 1)
    r = minimize_scalar(f, bounds=(0.5, 1.2), method='bounded', options=dict(xatol=1e-13))
    return 100 * float(r.fun), float(r.x)


def thermal_spinodal_bulk(tphi, tX, tlin):
    Fs = np.linspace(0.3, 4.0, 1851)
    rows = [HO.hom_state(F, tphi, tX, tlin, True) for F in Fs]
    Om = np.array([r['Omega'] for r in rows]); i = int(np.argmin(Om))
    r = minimize_scalar(lambda F: HO.hom_state(F, tphi, tX, tlin, True)['Omega'], bracket=(Fs[max(i - 2, 0)], Fs[i], Fs[min(i + 2, len(Fs) - 1)]), tol=1e-12)
    sp = HO.hom_state(r.x, tphi, tX, tlin, True)
    ev, Sv = e_vac(tphi, tlin)
    Prel = lambda F: (lambda h: h['n'] * h['Omega'] - (h['e'] - ev))(HO.hom_state(F, tphi, tX, tlin, True))
    P = np.array([Prel(F) for F in Fs])
    sg = np.where(np.sign(P[:-1]) != np.sign(P[1:]))[0]
    j = sg[-1]
    Fb = brentq(Prel, Fs[j], Fs[j + 1], xtol=1e-14)
    bk = HO.hom_state(Fb, tphi, tX, tlin, True)
    return sp, bk, ev, Sv


def phi_curvature(F, tphi, tlin):
    S = HO.S_of_F2(F * F, tphi, tlin, True)
    return HO.B_ * (3 * S * S - 1) - 6 * HO.C_ * (S - 1) + F * F + tphi, S      # d^2U/dS^2 in m_X^2 (U in g^2 v^4, S = phi/v)


def clock_H_over_mX(nQ_mX3, gstar=96.0):
    """H/m_X when the hidden charge density is nQ (m_X^3), with the free-X yield and g*s ~ g* ~ 96 near 100 GeV [R]."""
    mx = HQ.MX
    s_vis = nQ_mX3 * mx ** 3 / YQ_FREEX
    T = (s_vis * 45 / (2 * math.pi ** 2 * gstar)) ** (1 / 3)
    H = math.sqrt(math.pi ** 2 * gstar / 90) * T * T / MPL
    return H / mx, T


def point(rel, R):
    HQ.set_release(rel['N_X'], rel['E_X'], rel['H_onset_GeV'], HQ.V)
    st, _ = HQ.solve(R, True, 0.0, None)
    a = HQ.onset_a(st); mc = HQ.MX * a / HQ.H
    T_on = st['T'] / mc; n_c_on = st['cond_frac'] * rel['N_X'] / mc ** 3; nQ_on = rel['N_X'] / mc ** 3
    tphi, tX, tlin = thermal_terms(st, mc)
    SQ = st['S'] / rel['N_X']
    out = dict(R=R, S_per_Q=SQ, T_h_on=T_on, n_cond_on=n_c_on, nQ_on=nQ_on, cond_frac=st['cond_frac'], tphi_on=tphi, tX_on=tX, tlin_on=tlin)
    # ---- (b) T = 0 decomposition
    z0 = HO.hom_scan()
    sp0, bk0 = z0['spinodal'], z0['bulk']
    dEQ0 = sp0['e'] / sp0['n'] - bk0['Omega']
    af0 = (n_c_on / sp0['n']) ** (1 / 3); Tsp0 = T_on / af0
    out['T0'] = dict(n_sp=sp0['n'], Omega_sp=sp0['Omega'], phi_eq_sp=sp0['S'], EQ_hom_sp=sp0['e'] / sp0['n'], Omega_b=bk0['Omega'],
                     dE_per_cond_charge=dEQ0, a_sp_over_a_on=af0, T_h_sp=Tsp0, dSQ=st['cond_frac'] * dEQ0 / Tsp0,
                     a_bulk_over_a_on=(n_c_on / bk0['n']) ** (1 / 3))
    # ---- (b) with frozen thermal masses scaled as 1/a^2 [E]: find af where n_c(af) = n_sp(af)
    def gap(af):
        sp, bk, ev, Sv = thermal_spinodal_bulk(tphi / af ** 2, tX / af ** 2, tlin / af ** 2)
        return n_c_on / af ** 3 - sp['n']
    afT = brentq(gap, 0.5, 1.5, xtol=1e-8)
    spT, bkT, evT, SvT = thermal_spinodal_bulk(tphi / afT ** 2, tX / afT ** 2, tlin / afT ** 2)
    dEQT = (spT['e'] - evT) / spT['n'] - bkT['Omega']
    TspT = T_on / afT
    out['thermal'] = dict(a_sp_over_a_on=afT, T_h_sp=TspT, n_sp=spT['n'], Omega_sp=spT['Omega'], phi_eq_sp=spT['S'], n_bulk=bkT['n'],
                          Omega_b=bkT['Omega'], e_vac=evT, phi_vac=SvT, dE_per_cond_charge=dEQT, dSQ=st['cond_frac'] * dEQT / TspT)
    # ---- (a) crossover adiabaticity along the homogeneous branch up to the spinodal (T = 0 and thermal)
    curv = []
    for af in np.linspace(0.45, min(af0, afT), 60):
        n = n_c_on / af ** 3
        h0 = HO.hom_at_n(n); hT = HO.hom_at_n(n, tphi / af ** 2, tX / af ** 2, tlin / af ** 2, True)
        c0, S0 = phi_curvature(h0['F'], 0.0, 0.0); cT, ST = phi_curvature(hT['F'], tphi / af ** 2, tlin / af ** 2)
        curv.append((float(af), float(c0), float(cT), float(S0), float(ST)))
    cmin0 = min(c[1] for c in curv); cminT = min(c[2] for c in curv)
    Hm, Tvis = clock_H_over_mX(nQ_on / min(af0, afT) ** 3)
    out['crossover'] = dict(m_eff2_min_T0=cmin0, m_eff2_min_thermal=cminT, H_over_mX_at_sp=Hm, T_vis_GeV_at_sp=Tvis,
                            adiabaticity_H_over_meff=Hm / math.sqrt(min(cmin0, cminT)),
                            coherent_excitation_bound=(Hm / math.sqrt(min(cmin0, cminT))) ** 2,
                            profile=curv[::6])
    # ---- (b') surface energy of the first lumps and its dissipation
    Qf = (4e11, 1.2e13)
    surf = [EFIT[1] * Q ** (-1 / 3) + EFIT[2] * Q ** (-2 / 3) + EFIT[3] / Q for Q in Qf]
    out['surface'] = dict(Q_frag=Qf, excess_E_per_Q=surf, dSQ_at_T_sp=[st['cond_frac'] * s / TspT for s in surf],
                          dSQ_if_released_at_T_0p8=[st['cond_frac'] * s / 0.8 for s in surf])
    # ---- (b'') evaporation lag (net/gross Hertz-Knudsen): Delta(S/Q) <= eta_max * H/Gamma_ex, Gamma_ex/H >= 1e5 [Q]
    out['evaporation_lag'] = dict(bound=st['cond_frac'] * 1e-5, note='per evaporated charge Delta s = -ln(z/z_sat) ~ H/Gamma_ex')
    return out


def main():
    sq = json.load(open(sys.argv[1]))
    rel = sq['releases']['benchmark_cMY0']
    Rs = dict(R0=0.0, R_crit_bindingQ=rel['R_crit']['Q_bind_71263']['R_crit'], R_crit_Q2e5=rel['R_crit']['Q_2e5']['R_crit'],
              R_crit_bulk=rel['R_crit']['bulk_Q_inf']['R_crit'], R_e1e6_transfer=3.6484, R_e1e6_literal=4.96654, R8=8.0)
    out = dict(points={})
    for k, R in Rs.items():
        o = point(rel, R)
        out['points'][k] = o
        print('%-18s R=%.4f S/Q=%.4f T_on=%.3f | T0: dE/Q=%.5f T_sp=%.3f dSQ=%.5f | thermal: dE/Q=%.5f T_sp=%.3f n_sp=%.1f n_b=%.1f dSQ=%.5f | '
              'm_eff2_min T0 %.3f th %.3f H/m %.2e | surf dSQ %.1e-%.1e' % (
                  k, R, o['S_per_Q'], o['T_h_on'], o['T0']['dE_per_cond_charge'], o['T0']['T_h_sp'], o['T0']['dSQ'],
                  o['thermal']['dE_per_cond_charge'], o['thermal']['T_h_sp'], o['thermal']['n_sp'], o['thermal']['n_bulk'], o['thermal']['dSQ'],
                  o['crossover']['m_eff2_min_T0'], o['crossover']['m_eff2_min_thermal'], o['crossover']['H_over_mX_at_sp'],
                  min(o['surface']['dSQ_at_T_sp']), max(o['surface']['dSQ_if_released_at_T_0p8'])), flush=True)
    # R_crit shift from the decomposition heat: S/Q_final = S/Q_H(R) + dSQ(R); dS/dR from the scan
    scan = rel['scan']
    Rg = np.array([x['R'] for x in scan]); Sg = np.array([x['S_per_Q'] for x in scan])
    slope = lambda R: float(np.interp(R, 0.5 * (Rg[1:] + Rg[:-1]), np.diff(Sg) / np.diff(Rg)))
    shifts = {}
    for thr, key in (('bulk_Q_inf', 'R_crit_bulk'), ('Q_2e5', 'R_crit_Q2e5'), ('Q_bind_71263', 'R_crit_bindingQ')):
        p = out['points'][key]
        for lab in ('T0', 'thermal'):
            d = p[lab]['dSQ'] + max(p['surface']['dSQ_if_released_at_T_0p8'])
            shifts[thr + '_' + lab] = dict(dSQ_total_upper=d, slope_dSQ_dR=slope(p['R']), dR_crit=-d / slope(p['R']))
    out['R_crit_shift'] = shifts
    # (c) burn-out hypothetical
    burn = {}
    for pval in (5.6e-4, 1.3e-3):
        for T in (0.1, 0.05, 0.01):
            burn['p%.1e_T%g' % (pval, T)] = dict(dSQ_mPlusOmega=(1 + EFIT[0]) * pval / T, dSQ_2Omega=2 * EFIT[0] * pval / T)
    out['burnout_hypothetical'] = burn
    out['burnout_M1r'] = ('0 at the equilibrium level in both branches: free-X branch - no objects below T_ev ~ 0.8 m_X, so the gas-frozen '
                          'Xbar (T ~ 0.1) are never absorbed; Q-matter branch - objects present throughout, absorption keeps mu_Xbar = -omega '
                          '(reversible) until capture freeze-out at T_f ~ 0.005-0.015 m_X, where n_Xbar,eq ~ exp(-(m_X+omega)/T) ~ e^-120')
    print(json.dumps(shifts, indent=1))
    json.dump(out, open(sys.argv[2], 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
