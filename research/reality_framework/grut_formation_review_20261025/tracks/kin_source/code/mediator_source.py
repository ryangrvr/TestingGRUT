#!/usr/bin/env python3
"""kin_source (Claude, own code): the mediator-debris sources (R = 1.14) as late sources with their own timing.

(1) Y_L -> X X X  or  Xbar Xbar Xbar (half each; Y_L is the real CP-even eigenmode) [D, sp_yledger]:
    m_L^2 = M_Y^2 - sqrt2 mu h, Gamma_L = 3 kappa^2 m_L/(256 pi^3), parent population N_L = N_YH + N_YL (the heavy
    quanta cascade to Y_L first).  Three-body contact decay with massless daughters: dN/dE = 8 E/m_L^2 on
    [0, m_L/2] per daughter (uniform Dalitz density), mean E = m_L/3.  Time profile dN_L/dln a = gamma_L a^2 N_L(a)
    (own ODE).  Daughter comoving momentum q = a E.  Spectrum per daughter (derived here, checked numerically):
        dN/dq = (4 q gamma_L / m^2) E_1(2 gamma_L q^2 / m^2),  <q> = m <a>/3 with <a> = a_unit sqrt(pi/2),
    with m = m_L/H_i and a_unit = gamma_L^(-1/2).  R_Y = 3 N_L <q> / E_X = N_L (m_L/H_i) <a>_L / E_X.
(2) Cascade vectors: one A per Y_H (Y_H -> Y_L + A_L, Goldstone-equivalence width estimate [E]); they are
    nonrelativistic long before their decay and follow the vector law of decay_source.py (case cascade_vectors_*).
(3) sigma -> A A (two vectors per Sigma radial quantum); negligible.
Usage: python3 -B mediator_source.py ../results/decay_source_own.json ../results/mediator_source.json
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp, quad
from scipy.special import exp1

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kincommon as K

X_GRID = np.geomspace(1e-4, 10.0, 241)


def yl_source():
    ym = K.y_masses()
    H = K.H_GEV
    m = ym['mYL_GeV'] / H
    gam = ym['gamma_L']
    a_unit = ym['a_unit_L']
    NL0 = K.N_YH + K.N_YL
    nd_sign = 1.5                              # daughters of one sign per decay (3 x 1/2)

    def rhs(l, y):
        a = math.exp(l)
        r = gam * a * a * y[0]
        # per sign: 1.5 daughters per decay, each with mean rest-frame energy m/3, conformal a m/3
        return [-r, nd_sign * r, nd_sign * r * a * m / 3]
    sol = solve_ivp(rhs, (0.0, math.log(10 * a_unit)), [NL0, 0.0, 0.0], method='DOP853', rtol=1e-12, atol=1e-30,
                    dense_output=True)
    NL_end, Ns_end, Sq = sol.y[:, -1]
    mean_q = Sq / Ns_end
    # analytic spectrum per daughter in y = 2 gamma q^2/m^2: pdf(y) = E1(y)
    a_g = X_GRID * a_unit
    y_sol = sol.sol(np.log(a_g))
    q_grid = np.geomspace(1e-4, 10.0, 721) * m * a_unit / 3  # q grid around <q>
    pdf_q = (4 * q_grid * gam / m**2) * exp1(2 * gam * q_grid**2 / m**2)   # per daughter, dN/dq
    dN_sign_dlnq = nd_sign * NL0 * pdf_q * q_grid
    # numerical convolution check: dN/dq = int da P(a) (1/a) 8 (q/a)/m^2 for a > 2q/m
    def conv(q):
        f = lambda a: gam * a * math.exp(-gam * a * a / 2) * 8 * q / (a * a * m * m)
        return quad(f, 2 * q / m, 10 * a_unit, limit=200, epsabs=0, epsrel=1e-10)[0]
    conv_check = max(abs(conv(q) / p - 1) for q, p in zip(q_grid[::40], pdf_q[::40]))
    norm_check = float(np.trapezoid(pdf_q, q_grid))
    mean_q_spec = float(np.trapezoid(pdf_q * q_grid, q_grid))
    R_Y = 3 * NL0 * mean_q / K.E_X
    cum = np.cumsum(np.r_[0, 0.5 * (pdf_q[1:] + pdf_q[:-1]) * np.diff(q_grid)])
    return dict(
        parents=dict(N_L=NL0, N_YH=K.N_YH, N_YL=K.N_YL, mYL_GeV=ym['mYL_GeV'], mYL_over_H=m, Gamma_L_GeV=ym['Gamma_L_GeV'],
                     gamma_L=gam, a_unit=a_unit, a_mean_decay=ym['a_mean_L'], t_unit_over_Hinv=0.5 * a_unit**2,
                     T_vis_at_a_mean_GeV=K.T_vis_of_a(ym['a_mean_L']),
                     branching='1/2 to XXX, 1/2 to XbarXbarXbar [D, real eigenmode]; Y_H cascades to Y_L first'),
        daughters=dict(per_sign_per_decay=nd_sign, N_sign=Ns_end, N_sign_per_net_charge=Ns_end / abs(K.N_X),
                       mean_q_over_H=mean_q, mean_q_analytic=m * ym['a_mean_L'] / 3, mean_q_from_spectrum=mean_q_spec,
                       mean_E_rest_frame_over_mL=1 / 3, mean_p_at_reference_over_m=mean_q * K.H_GEV / (K.A_REF * K.MX),
                       R_Y=R_Y, R_Y_closed_form=NL0 * m * ym['a_mean_L'] / K.E_X,
                       f_rest_if_unprocessed=2 * K.MX * K.YQ * (Ns_end / abs(K.N_X)) / K.DM_PER_S,
                       number_closure=(NL_end + Ns_end / nd_sign) / NL0 - 1,
                       spectrum_norm_check=norm_check, convolution_vs_analytic_max_rel=conv_check,
                       rest_frame_energy_pdf='dN/dE = 8 E / m_L^2 on [0, m_L/2] (massless daughters, contact |M|^2)'),
        time_profile=dict(x=X_GRID, a=a_g, N_L_remaining=y_sol[0], dN_sign_dlna=nd_sign * gam * a_g**2 * y_sol[0],
                          dN_sign_dt_over_H=nd_sign * gam * y_sol[0], cumulative_fraction=y_sol[1] / (nd_sign * NL0),
                          T_vis_GeV=[K.T_vis_of_a(a) for a in a_g[::30]], T_vis_a=a_g[::30]),
        spectrum=dict(q_over_H=q_grid, dN_sign_dlnq=dN_sign_dlnq, dE_conf_sign_dlnq=dN_sign_dlnq * q_grid,
                      dR_dlnq_both_signs=2 * dN_sign_dlnq * q_grid / K.E_X, cumulative_fraction=cum,
                      p_over_mX_at_a_ref=q_grid * K.H_GEV / (K.A_REF * K.MX)),
        two_d_rule='inject at time a: three daughters, energies E ~ 8E/m_L^2 dE on [0, m_L/2] (rest frame = lab, parent cold), q = a E; sign XXX or XbarXbarXbar with probability 1/2 each',
    )


def main(decay_path, out):
    dec = json.loads(Path(decay_path).read_text())
    cases = {c['label']: c for c in dec['cases']}
    ym = K.y_masses()
    ylsrc = yl_source()
    casc = {}
    for e in (1e-7, 1e-6):
        c = cases['cascade_vectors_e1e-7' if e == 1e-7 else 'cascade_vectors_e1e-6']
        vp = K.vector_params(e)
        # cascade vector: emitted at a_emit ~ few (Gamma_H->LA / H = const x a^2), relativistic at birth
        GH = ym['Gamma_H_to_LA_GeV_estimate']
        a_eq = math.sqrt(K.H_GEV / GH)
        EA = ym['E_A_from_YH_GeV']
        q_A = a_eq * EA / K.H_GEV              # comoving momentum of the emitted vector, H_i units (order of magnitude)
        a_NR = q_A / vp['mA_over_H']            # when it becomes nonrelativistic
        casc['e1e-7' if e == 1e-7 else 'e1e-6'] = dict(
            N_A=K.N_YH, R=c['moments']['R_from_energy_moment'], pairs_per_net_charge=c['moments']['pairs_per_net_charge'],
            Gamma_YH_to_YL_A_GeV_estimate=GH, a_where_Gamma_equals_H=a_eq, E_A_at_emission_GeV=EA,
            q_A_over_H_at_emission=q_A, a_nonrelativistic=a_NR, a_decay_mean=vp['a_mean'],
            decayed_fraction_while_relativistic=-math.expm1(-vp['gamma'] * a_NR**2 / 2),
            spectrum='identical Rayleigh law to the primordial vectors of this e (decay_source case cascade_vectors_*)',
            population_caveat='N_YH is the late undamped classical count (a = 30-40) [Q]; the population present when Gamma_H->LA = H (a ~ %.0f) is not re-derived here (bracket, hw_gates/sp_yledger)' % a_eq)
    sig = dict(N_sigma=K.N_SIG, vectors=2 * K.N_SIG, msigma_GeV=ym['msigma_GeV'], E_A_GeV=ym['msigma_GeV'] / 2,
               R_e1e_7=cases['sigma_vectors_e1e-7']['moments']['R_from_energy_moment'], note='negligible; included')
    totals = {}
    for e in (1e-7, 1e-6):
        lab = 'e1e-7' if e == 1e-7 else 'e1e-6'
        RY = ylsrc['daughters']['R_Y']; RC = casc[lab]['R']
        # sigma vectors are vectors of the same e: their deposit scales exactly like the cascade vectors'
        RS = sig['R_e1e_7'] * (casc['e1e-6']['R'] / casc['e1e-7']['R'] if e == 1e-6 else 1.0)
        totals[lab] = dict(R_Y=RY, R_cascade=RC, R_sigma=RS, R_mediator_total=RY + RC + RS,
                                   pairs_per_net_charge_total=ylsrc['daughters']['N_sign_per_net_charge'] + casc[lab]['pairs_per_net_charge'] + 2 * K.N_SIG / abs(K.N_X),
                                   codex_225028b=dict(R_Y=0.8055700748053104, R_cascade=0.338883874466159 if e == 1e-7 else 0.10716449056156375,
                                                      R_Sigma=1.240711650236604e-07 if e == 1e-7 else 3.923474734253848e-08))
    ordering = dict(
        a_YH_to_YL_A_prompt=casc['e1e-7']['a_where_Gamma_equals_H'],
        a_YL_decay_mean=ylsrc['parents']['a_mean_decay'],
        a_vector_decay_mean_e1e_7=cases['e1e-7_band']['decay']['a_mean_decay'],
        a_vector_decay_mean_e1e_6=cases['e1e-6_band']['decay']['a_mean_decay'],
        a_case_L_reference=K.A_REF,
        a_vector_daughters_nonrelativistic_e1e_7=cases['e1e-7_band']['moments']['mean_q_over_H'] / (K.MX / K.H_GEV),
        T_vis=dict(YL_decay=ylsrc['parents']['T_vis_at_a_mean_GeV'], vector_decay_e1e_7=cases['e1e-7_band']['decay']['T_vis_at_a_mean_GeV'],
                   vector_decay_e1e_6=cases['e1e-6_band']['decay']['T_vis_at_a_mean_GeV'], case_L_reference=K.T_vis_of_a(K.A_REF)))
    res = dict(method=__doc__, YL_three_body=ylsrc, cascade_vectors=casc, sigma_vectors=sig, totals=totals, time_ordering=ordering)
    K.save(out, res)
    d = ylsrc['daughters']; p = ylsrc['parents']
    print('Y_L: m_L=%.4e GeV Gamma_L=%.4e GeV a_unit=%.1f <a>=%.1f T_vis=%.2e GeV; N_sign=%.5f (%.3e per net) <q>=%.1f (analytic %.1f, spectrum %.1f) R_Y=%.5f closure=%.1e conv=%.1e norm=%.6f'
          % (p['mYL_GeV'], p['Gamma_L_GeV'], p['a_unit'], p['a_mean_decay'], p['T_vis_at_a_mean_GeV'], d['N_sign'], d['N_sign_per_net_charge'],
             d['mean_q_over_H'], d['mean_q_analytic'], d['mean_q_from_spectrum'], d['R_Y'], d['number_closure'], d['convolution_vs_analytic_max_rel'], d['spectrum_norm_check']))
    for e, c in casc.items():
        print('cascade %s: Gamma_H->LA~%.2e GeV (=H at a=%.1f), E_A=%.2e GeV, q_A~%.0f, a_NR~%.1e, decayed-while-relativistic %.1e, R=%.5f' %
              (e, c['Gamma_YH_to_YL_A_GeV_estimate'], c['a_where_Gamma_equals_H'], c['E_A_at_emission_GeV'], c['q_A_over_H_at_emission'], c['a_nonrelativistic'], c['decayed_fraction_while_relativistic'], c['R']))
    for e, t in totals.items():
        print('totals %s: R_Y=%.5f R_casc=%.5f R_sig=%.2e total=%.5f pairs/net=%.4e' % (e, t['R_Y'], t['R_cascade'], t['R_sigma'], t['R_mediator_total'], t['pairs_per_net_charge_total']))
    print('ordering', {k: v for k, v in ordering.items() if k != 'T_vis'}, ordering['T_vis'])


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
