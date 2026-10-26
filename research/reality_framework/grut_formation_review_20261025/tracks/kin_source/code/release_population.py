#!/usr/bin/env python3
"""kin_source (Claude, own code): the release's own population (no-source control) and the equilibrium comparators.

A. Onset / case-L coordinate [D]: n_Q = |eps| n_c with n_c = 2 sqrt(2 lambda_X) x_c^3, x_c^2 = lambda_phi v^2/(2 g_X^2);
   a_ref = H_i (|N_X|/n_Q)^(1/3); the released condensate's energy there is E_X H_i^4/a_ref^4 = rho_c = 3 lambda_X x_c^4
   (the two coordinates agree to 2e-5); the phi transition releases V0 = lambda_phi v^4/4 = m_X^2 x_c^2/2.
B. Adiabatic coherent-orbit count [D, own quadrature; method of our gp_leftovers budget.py, re-implemented]:
   planar orbit R(t), theta(t) of the homogeneous X field with L = |Xdot|^2 - V, charge n = 2 R^2 thetadot,
   J_r = (1/2 pi) oint p_R dR, p_R = 2 Rdot.  In the quartic era eps_q = n/n_max(rho) (this is the released eps);
   if a^3 n and a^3 J_r are adiabatic invariants through the transition to the quadratic era, the nonrelativistic
   quanta are N_- = J_r, N_+ = J_r + n per volume, i.e. Xbar per net charge p_0 = j = J_r/n and eps_NR = 1/(1 + 2 j).
C. Energy cap [D]: cold quanta per net charge <= (rho_c + V0)/(m_X n_Q); Xbar per net charge <= (that - 1)/2.
D. Fragmentation template [A]: our tracks/spectrum certified free-mode fractions (prepared spherical object, Q0 = 38969,
   E0/Q0 = 0.99954 m_X) scaled by |N_X|; applicable only as a spectral shape for a charge-saturated object.
E. Equilibrium comparators [C, own Bose thermodynamics]: fully thermalised X, Xbar (and s) gas at fixed energy
   rho_c (1 + R) + V0 and net charge n_Q, condensate at mu_c = omega with linear energy omega n_c (Codex's control
   definition), for R = 0 (no source), 40.5626, 48.5; the reservoir table at zero chemical potential.
Usage: python3 -B release_population.py ../results/release_population.json
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import kve

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kincommon as K

SPEC = K.TRACKS / 'spectrum' / 'results'


# ---------------------------------------------------------------- B. orbit quadrature
def J_r(Vf, rho, n, npts=2000):
    f = lambda R: rho - n * n / (4 * R * R) - Vf(R)
    Rs = np.geomspace(1e-6, 10, 40001)
    fs = np.array([f(R) for R in Rs])
    i = int(np.argmax(fs))
    R2 = brentq(f, Rs[i], 10.0, xtol=1e-15)
    R1 = brentq(f, 1e-9, Rs[i], xtol=1e-15) if n > 0 else 0.0
    # Gauss-Legendre in theta with R = mid - half cos(theta): sqrt(f) dR = sqrt(f) half sin(theta) dtheta
    xg, wg = np.polynomial.legendre.leggauss(npts)
    th = 0.5 * math.pi * (xg + 1); wt = 0.5 * math.pi * wg
    mid, half = 0.5 * (R1 + R2), 0.5 * (R2 - R1)
    R = mid - half * np.cos(th)
    fv = np.maximum([f(r) for r in R], 0.0)
    return (2 / math.pi) * float(np.sum(wt * np.sqrt(fv) * half * np.sin(th))), R1, R2


def adiabatic_count(eps_q):
    nmax = 2 * math.sqrt(2) * (1 / 3)**0.75          # quartic V = R^4, rho = 1
    n = eps_q * nmax
    Jr, R1, R2 = J_r(lambda R: R**4, 1.0, n)
    j = Jr / n
    return dict(eps_q=eps_q, n_over_nmax=eps_q, J_r=Jr, j_Jr_over_n=j, eps_NR=1 / (1 + 2 * j),
                Xbar_per_net_charge=j, quanta_per_net_charge=1 + 2 * j, R1_over_R2=R1 / R2)


# ---------------------------------------------------------------- E. Bose thermodynamics (own)
def bose(m, T, mu=0.0, kmax=None):
    """One real bosonic degree: Bessel series with exponentially scaled K (Kolb-Turner form)."""
    if mu >= m:
        raise ValueError('mu >= m')
    if kmax is None:
        kmax = int(50 * T / (m - mu)) + 200
    k = np.arange(1, kmax + 1, dtype=float)
    z = k * m / T
    w = np.exp(-k * (m - mu) / T)
    k2, k1 = kve(2, z), kve(1, z)
    pref = m * m * T / (2 * math.pi**2)
    n = pref * np.sum(w * k2 / k)
    rho = pref * np.sum(w * (3 * T * k2 / (k * k) + m * k1 / k))
    P = pref * T * np.sum(w * k2 / (k * k))
    return dict(n=float(n), rho=float(rho), P=float(P))


def bose_quad(m, T, mu=0.0):
    """Independent direct momentum quadrature for the cross-check of bose()."""
    def integ(p, mom):
        E = math.sqrt(p * p + m * m)
        x = (E - mu) / T
        f = 0.0 if x > 700 else 1 / math.expm1(x)
        return p * p * f * (1.0 if mom == 0 else E if mom == 1 else p * p / (3 * E))
    out = {}
    for key, mom in (('n', 0), ('rho', 1), ('P', 2)):
        out[key] = quad(integ, 0, np.inf, args=(mom,), epsabs=0, epsrel=1e-11, limit=400)[0] / (2 * math.pi**2)
    return out


def charge(T, mu):
    return bose(K.MX, T, mu)['n'] - bose(K.MX, T, -mu)['n']


def eq_state(T, nQ, include_s, cap=K.OMEGA_OVER_M * K.MX):
    if charge(T, cap) < nQ:
        mu, nc = cap, nQ - charge(T, cap)
    else:
        mu, nc = brentq(lambda u: charge(T, u) - nQ, 0.0, cap, xtol=1e-15), 0.0
    plus, minus = bose(K.MX, T, mu), bose(K.MX, T, -mu)
    s = bose(K.MS, T) if include_s else dict(n=0.0, rho=0.0, P=0.0)
    rho = plus['rho'] + minus['rho'] + s['rho'] + cap * nc
    return dict(T_GeV=T, T_over_m=T / K.MX, mu_over_m=mu / K.MX, condensed_charge=nc, condensed_fraction=nc / nQ,
                plus=plus, minus=minus, neutral=s, rho=rho, pair_inventory=minus['n'] + s['n'],
                Xbar_per_net_charge=minus['n'] / nQ, s_per_net_charge=s['n'] / nQ,
                scalar_energy_fraction=s['rho'] / rho, entropy_per_net_charge=(plus['rho'] + plus['P'] - mu * plus['n'] + minus['rho'] + minus['P'] + mu * minus['n'] + s['rho'] + s['P']) / (T * nQ))


def thermal_state(E, nQ, include_s):
    g = lambda lt: eq_state(math.exp(lt), nQ, include_s)['rho'] - E
    lt = brentq(g, math.log(0.02 * K.MX), math.log(200 * K.MX), xtol=1e-13)
    return eq_state(math.exp(lt), nQ, include_s)


# ---------------------------------------------------------------- D. template from tracks/spectrum
def template():
    f4 = json.loads((SPEC / 'F4long_analysis.json').read_text())
    f3 = json.loads((SPEC / 'F3h_analysis.json').read_text())
    Q0, E0 = f4['Q0'], f4['E0']
    bins = f4['k_bins']
    wX = f4['fixed_radius']['480.0']['windows']['48000']       # X spectrum (certified for p/m_X >= 0.03)
    wXb = f3['fixed_radius']['240.0']['windows']['12000']      # converged Xbar and s (dt = 0.025, h = 0.125)
    NXb, EXb, Ns, Es = 256.07, 433.7, 109.3, 354.0              # certified values quoted from tracks/spectrum NOTES [Q own track]
    w240 = f4['fixed_radius']['240.0']['windows']['48000']    # totals that have passed r = 240 by T = 48000
    NX_esc, EX_esc = w240['NX'], w240['EX']
    scale = abs(K.N_X) / Q0
    xbar_bins = np.array(wXb['Xbar_bins']) * (NXb / wXb['NXbar'])
    x_bins = np.array(wX['X_bins']) * (NX_esc / wX['NX'])     # R=480 shape (certified p/m_X >= 0.03), R=240 total
    cores = Q0 - NX_esc + NXb                                   # charge not escaped by T = 48000 (cores + halo)
    conv = K.H_GEV / (K.A_REF * K.MX)                            # q/H per unit p/m_X at a_ref
    return dict(
        source='tracks/spectrum (own, certified): prepared spherical Gaussian, Q0=%.2f, E0/Q0=%.5f m_X; X bins F4long R=480 T=48000, Xbar bins F3h R=240 T=12000 rescaled to the converged N_Xbar' % (Q0, E0 / Q0),
        fractions=dict(Xbar_over_Q0=NXb / Q0, s_over_Q0=Ns / Q0, X_escaped_over_Q0=NX_esc / Q0, cores_and_halo_over_Q0=cores / Q0,
                       mean_omega_Xbar=EXb / NXb, mean_omega_s=Es / Ns, mean_omega_X=EX_esc / NX_esc, E0_over_Q0_mX=E0 / Q0),
        scaled_to_point=dict(scale_factor=scale, N_Xbar=NXb * scale, N_s=Ns * scale, N_X_free=NX_esc * scale, N_cores=cores * scale,
                             Xbar_per_net_charge=NXb / Q0, s_per_net_charge=Ns / Q0,
                             energy_per_net_charge_mX=E0 / Q0,
                             release_energy_per_net_charge_mX=(K.RHO_C + K.V0) / (K.N_Q * K.MX),
                             energy_fraction_carried=(E0 / Q0) / ((K.RHO_C + K.V0) / (K.N_Q * K.MX))),
        bins_p_over_mX=bins,
        Xbar_dN_per_bin=(xbar_bins * scale).tolist(), X_dN_per_bin=(x_bins * scale).tolist(),
        X_bins_certified_from_p_over_mX=0.03,
        bins_q_over_H_at_a_ref=[b * conv / 1.0 for b in bins],
        s_quanta='N_s = %.1f scaled, <omega> = 3.24 m_X (p ~ 2.1 m_X); each can decay s -> X Xbar (vacuum width 1.73e-5 GeV)' % (Ns * scale),
        assumption='[A] the heavy-point release fragments like a charge-saturated object; it does not: eps_q = 0.146 means ~ (1+2j) quanta per net charge and 5.95 m_X of energy per net charge, against 1.0 m_X in the template. Use only as a spectral shape alternative.')


def main(out):
    H = K.H_GEV
    nQ = K.N_Q
    onset = dict(x_c_GeV=K.XC, n_c_GeV3=K.NC, rho_c_GeV4=K.RHO_C, V0_GeV4=K.V0, V0_over_rho_c=K.V0 / K.RHO_C,
                 mX2_xc2_over_2V0=K.MX**2 * K.XC**2 / (2 * K.V0), n_Q_GeV3=nQ, a_ref_from_charge=K.A_REF, a_ref_from_energy=K.A_REF_E,
                 ratio=K.A_REF / K.A_REF_E, eps_from_window_means=K.EPS,
                 eps_from_NX_EX=abs(K.N_X) / (2 * math.sqrt(2 * K.LAM_X) * (K.E_X / (3 * K.LAM_X))**0.75),
                 T_vis_ref_GeV=K.T_vis_of_a(K.A_REF), T_vis_ref_from_YQ_entropy=brentq(lambda T: K.YQ * K.entropy(T) - nQ, 1.0, 1e3),
                 H_ref_constant_gstar_GeV=H / K.A_REF**2, H_ref_with_gstar_GeV=K.hubble_of_T(K.T_vis_of_a(K.A_REF)),
                 energy_per_net_charge_mX=dict(rho_c=K.RHO_C / (nQ * K.MX), V0=K.V0 / (nQ * K.MX), total=(K.RHO_C + K.V0) / (nQ * K.MX)),
                 condensate_conformal_amplitude_over_H=(K.E_X / K.LAM_X)**0.25,
                 condensate_conformal_frequency_scale_over_H=math.sqrt(2 * K.LAM_X) * (K.E_X / K.LAM_X)**0.25,
                 condensate_momentum_scale_at_ref_over_mX=math.sqrt(2 * K.LAM_X) * (K.E_X / K.LAM_X)**0.25 * H / (K.A_REF * K.MX),
                 phase_note='before the transition phi = 0 (X massless, pure quartic); s exists only after it; the hidden gas temperature at onset may restore the phase (ls_threshold: one-loop T_c = 18-42 m_X) [Q] - open for the solver',
                 V0_release_alternatives=dict(sudden_into_s_at_rest=dict(n_s_GeV3=K.V0 / K.MS, s_per_net_charge=K.V0 / (K.MS * nQ)),
                                              adiabatic_zero_pressure_coexistence='no V0 released to quanta (gp_leftovers)'))
    # B
    ad = adiabatic_count(abs(K.EPS))
    harm = [dict(n=n, Jr=J_r(lambda R: R * R, 1.0, n)[0], exact=(1 - n) / 2) for n in (0.2, 0.5, 0.9)]
    ad_check = adiabatic_count(0.1)
    ad['harmonic_control_max_rel_err'] = max(abs(h['Jr'] / h['exact'] - 1) for h in harm)
    ad['gp_leftovers_check_eps0.1'] = dict(f_adiabatic_here=1.01427 * (1 / ad_check['eps_NR'] - 1), f_adiabatic_gp_leftovers_table=8.09,
                                           f_NR_count_here=1.01427 * (1 / 0.1 - 1))
    ad['frozen_f_if_never_processed'] = 2 * ad['Xbar_per_net_charge'] * K.MX * K.YQ / K.DM_PER_S
    ad['NR_count_law_at_eps_q'] = dict(Xbar_per_net_charge=(1 / abs(K.EPS) - 1) / 2, f=1.01427 * (1 / abs(K.EPS) - 1))
    # C
    cap = dict(quanta_per_net_charge_with_V0=(K.RHO_C + K.V0) / (nQ * K.MX), quanta_per_net_charge_without_V0=K.RHO_C / (nQ * K.MX),
               Xbar_per_net_charge_max_with_V0=((K.RHO_C + K.V0) / (nQ * K.MX) - 1) / 2,
               Xbar_per_net_charge_max_without_V0=(K.RHO_C / (nQ * K.MX) - 1) / 2)
    # D
    tmpl = template()
    # E
    series_vs_quad = []
    for m in (K.MX, K.MS):
        for t in (0.3, 1.0, 12.0, 30.0):
            for mu in (0.0, 0.5 * m):
                a, b = bose(m, t * K.MX, mu), bose_quad(m, t * K.MX, mu)
                series_vs_quad.append(max(abs(a[k] / b[k] - 1) for k in ('n', 'rho', 'P')))
    reservoir = []
    for t in (0.3, 1.0, 3.0, 12.0, 30.0):
        X, s = bose(K.MX, t * K.MX), bose(K.MS, t * K.MX)
        reservoir.append(dict(T_over_m=t, n_s_over_nminus=s['n'] / X['n'], scalar_fraction_of_gas_energy=s['rho'] / (2 * X['rho'] + s['rho'])))
    tri = 2 * K.GX**2 * K.V_PHI
    width = tri**2 / (16 * math.pi * K.MS) * math.sqrt(1 - 4 * K.MX**2 / K.MS**2)
    states = {}
    for R in (0.0, 40.56262288525937, 48.49912):
        for include_s in (False, True):
            for withV0 in (True, False):
                E = K.RHO_C * (1 + R) + (K.V0 if withV0 else 0.0)
                st = thermal_state(E, nQ, include_s)
                states['R%.4f_s%d_V0%d' % (R, include_s, withV0)] = dict(R=R, include_s=include_s, with_V0=withV0, rho_GeV4=E, **st)
    res = dict(method=__doc__, onset=onset, adiabatic_orbit_count=ad, harmonic_control=harm, energy_cap=cap, fragmentation_template=tmpl,
               equilibrium_comparators=dict(bose_series_vs_quadrature_max_rel=max(series_vs_quad), reservoir_zero_mu=reservoir,
                                            ms_GeV=K.MS, trilinear_GeV=tri, vacuum_Gamma_s_to_XXbar_GeV=width, states=states,
                                            scope='fixed vacuum masses, linear condensate energy omega*n_c; comparator only, not the population'))
    K.save(out, res)
    print('onset: a_ref=%.6e (energy %.6e, ratio %.7f) n_Q=%.6e rho_c=%.5e V0=%.5e T_vis=%.2f GeV; E/net charge = %.3f + %.3f = %.3f m_X; condensate p/m at ref %.3f'
          % (K.A_REF, K.A_REF_E, K.A_REF / K.A_REF_E, nQ, K.RHO_C, K.V0, onset['T_vis_ref_GeV'], onset['energy_per_net_charge_mX']['rho_c'],
             onset['energy_per_net_charge_mX']['V0'], onset['energy_per_net_charge_mX']['total'], onset['condensate_momentum_scale_at_ref_over_mX']))
    print('adiabatic: j=%.4f eps_NR=%.4f quanta/net=%.3f Xbar/net=%.3f (NR-count law %.3f); harmonic ctrl %.1e; gp_leftovers eps=0.1 check f=%.3f vs 8.09'
          % (ad['j_Jr_over_n'], ad['eps_NR'], ad['quanta_per_net_charge'], ad['Xbar_per_net_charge'], ad['NR_count_law_at_eps_q']['Xbar_per_net_charge'],
             ad['harmonic_control_max_rel_err'], ad['gp_leftovers_check_eps0.1']['f_adiabatic_here']))
    print('energy cap: quanta/net %.3f (no V0 %.3f); Xbar/net <= %.3f (%.3f)' % (cap['quanta_per_net_charge_with_V0'], cap['quanta_per_net_charge_without_V0'],
                                                                                 cap['Xbar_per_net_charge_max_with_V0'], cap['Xbar_per_net_charge_max_without_V0']))
    print('template: Xbar/Q0=%.5f s/Q0=%.5f Xesc/Q0=%.4f cores=%.4f; energy carried %.3f of release' % (tmpl['fractions']['Xbar_over_Q0'], tmpl['fractions']['s_over_Q0'],
                                                                                                        tmpl['fractions']['X_escaped_over_Q0'], tmpl['fractions']['cores_and_halo_over_Q0'], tmpl['scaled_to_point']['energy_fraction_carried']))
    print('bose series vs quad max rel %.1e; Gamma_s=%.5e GeV tri=%.6f' % (max(series_vs_quad), width, tri))
    for r in reservoir:
        print('  T/m=%5.1f n_s/n-=%.6f  scalar fraction=%.6f' % (r['T_over_m'], r['n_s_over_nminus'], r['scalar_fraction_of_gas_energy']))
    for k, st in states.items():
        print('  %-24s T/m=%.5f mu/m=%.5f cond=%.6f (%.4f of n_Q) Xbar/net=%.4f s/net=%.4f scalar frac=%.5f sigma/net=%.2f'
              % (k, st['T_over_m'], st['mu_over_m'], st['condensed_charge'], st['condensed_fraction'], st['Xbar_per_net_charge'], st['s_per_net_charge'], st['scalar_energy_fraction'], st['entropy_per_net_charge']))


if __name__ == '__main__':
    main(sys.argv[1])
