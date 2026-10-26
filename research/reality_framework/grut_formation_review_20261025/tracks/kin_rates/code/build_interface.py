#!/usr/bin/env python3
"""kin_rates: assemble results/rates_interface.json for the kinetic solver (own code).
usage: python3 -B build_interface.py <kernels_n2.json> <xsec_2n.json|none> [supplement.json|none] [balance.json|none]
Reads results/thermal_masses.json, results/kernels_22.json, results/mc_checks.json and the MC tables; recomputes
the MB thermal averages of the n<->2 tables on a common T grid (tables.py).  Units: GeV throughout; n = int d^3p
f/(2pi)^3 per degree; rate densities gamma in GeV^4 (reactions per volume per time, MB, mu = 0, S_in and S_out
included); <sigma v> = gamma S_in/(n1 n2) in GeV^-2; sigma in GeV^-2; Gt(s) dimensionless."""
import json
import sys
import numpy as np
import model as M
import tables as TB

RES = '../results'
TGRID = sorted(set(np.round(np.concatenate([np.geomspace(0.03, 30, 40), [0.3, 1.0, 3.0, 7.344, 8.118, 12.0, 13.06]]), 6)))
CHAIN = {  # the sequential chain each subtracted line represents, in the network's own channels
    'A_XXXbXb_XXb': 'XXb->s (inverse decay) then s XXb->XXb (E); two pairs: XXb->s twice then ss->XXb; nested: XXb->s, s->X Xb*, Xb* X->s ...; XXb->X Xb* (elastic) then Xb* X->s',
    'B_XXXXb_XX': 'XXb->s then s XX->XX (E2); X X->X X* (elastic) etc.',
    'C_XXXbXb_ss': 'XXb->s then s XXb->ss (G); XXb->s twice then ss->ss (elastic)',
    'H_XXbss_XXb': 'XXb->s then sss->XXb (L) / s ss->ss ...; s s->s s elastic chains',
    'J_XXXbs_Xs': 'XXb->s then s X s->X s (G2); Xs->Xs elastic chains',
    'K_XXbss_ss': 'XXb->s then sss->ss (F)',
    'D_XXXb_Xs': 'XXb->s then s X->s X (elastic); X Xb->X Xb* then Xb* X->s',
    'E_XXbs_XXb': 'XXb->s then s s->XXb; s->XXb then ...',
    'E2_XXs_XX': 'none below 5.7 m_X (no neutral pair in the initial state); X s->X* ... lines at higher energy',
    'G_XXbs_ss': 'XXb->s then ss->ss (elastic)',
    'G2_Xss_Xs': 'Xs->Xs elastic with the u-channel X pole (s->X Xb, Xb X->s)',
    'F_sss_ss': 'none (no neutral X Xb pair; s s->s* off shell)', 'L_sss_XXb': 'none at tree level below the X-line window', 'N_ssss_ss': 'none'}


def cols(rows, keys):
    return {k: [r.get(k) for r in rows] for k in keys}


def thermal_block(tab, tab04=None):
    out = dict(T_over_mX=list(TGRID), gamma_sub=[], gamma_sub_err=[], gamma_full_reg=[], gamma_sub_reg0_04=[])
    for t in TGRID:
        T = t * M.MX
        g = tab.thermal(T, 'sub'); out['gamma_sub'].append(g)
        out['gamma_sub_err'].append(abs(tab.thermal(T, 'sub', shift=1.0) - g))
        out['gamma_full_reg'].append(tab.thermal(T, 'full'))
        out['gamma_sub_reg0_04'].append(tab04.thermal(T, 'sub') if tab04 else None)
    # resolved flag per T: |gamma_sub| > 2 err and regulator spread < 50 %; bound otherwise
    out['remainder_resolved'] = [bool(abs(g) > 2 * e and (g4 is None or abs(g4 - g) < 0.5 * (abs(g) + e))) for g, e, g4 in zip(out['gamma_sub'], out['gamma_sub_err'], out['gamma_sub_reg0_04'])]
    out['remainder_bound'] = [max(abs(g) + e, abs(g4) if g4 is not None else 0.0) for g, e, g4 in zip(out['gamma_sub'], out['gamma_sub_err'], out['gamma_sub_reg0_04'])]
    return out


def merge_supplement(chan, supp):
    if not supp:
        return chan
    chan = dict(chan)
    chan['rows'] = sorted(chan['rows'] + [r for r in supp['rows'] if 'freg' in r], key=lambda r: (r['freg'], r['sqrt_s_over_mX']))
    return chan


def main(n2path, twonpath, supppath='none', balpath='none'):
    thermal = json.load(open(f'{RES}/thermal_masses.json'))
    k22 = json.load(open(f'{RES}/kernels_22.json'))
    mc = json.load(open(f'{RES}/mc_checks.json'))
    n2 = json.load(open(n2path))
    twon = json.load(open(twonpath)) if twonpath != 'none' else dict(channels={})
    supp = json.load(open(supppath)) if supppath != 'none' else {}
    bal = json.load(open(balpath)) if balpath != 'none' else {}
    P = k22['params']
    out = {}
    out['meta'] = dict(
        track='kin_rates (Claude, own code), 26 October 2026; BRIEF_KINETIC.md item: collision kernels of the light sector',
        model='U = lP/4 (phi^2-v^2)^2 + g^2 phi^2 |X|^2 + lX |X|^4, phi = v + s; tree amplitudes (Berends-Giele recursion, exact vertices of U)',
        parameters=dict(gX=M.GX, lambdaX=M.LX, lambdaPhi=M.LP, v_GeV=M.V, mX_GeV=M.MX, ms_GeV=M.MS, lam3_sXXb_GeV=M.LAM3,
                        lam_s3_GeV=M.LAM_S3, lam_s4=M.LAM_S4, lam_ssXXb=M.LAM_SSXX, Gamma_s_vac_GeV=M.GAMMA_S_VAC),
        units='GeV; n = int d^3p f/(2pi)^3 per bosonic degree (X, Xb, s: one each); gamma = reactions per volume per time (GeV^4) '
              'with S_in and S_out included; <sigma v> = gamma S_in/(n1 n2) (GeV^-2); sigma (GeV^-2) includes 1/S_out; Gt(s) = '
              'int dPhi_n dPhi_2 |M|^2 (dimensionless), gamma_{n->2}(T) = (1/(S_in S_out)) (T/8pi^3) int ds sqrt(s) K1(sqrt(s)/T) Gt(s); '
              'gamma_{2->n}(T) = (1/S_in) (T/32pi^4) int ds lambda(s,m1^2,m2^2) sigma(s) K1(sqrt(s)/T)/sqrt(s)',
        statistics='all thermal tables are Maxwell-Boltzmann at mu = 0 (tabulated per channel); Bose/MB correction factors at '
                   'T/mX = 0.3, 1, 3, 7.344, 12 from the lab-frame MC with Bose f and 1+f (mc_checks); a momentum-grid solver should '
                   'use sigma(s)/Gt(s) with its own f; a moments closure should multiply the MB table by the interpolated Bose factor '
                   'and by exp(sum mu_in/T) for chemical potentials (exact for MB; Bose at mu != 0 needs the full integral)',
        symmetry_factors='S_in = product of n! over identical incoming species, S_out likewise; included in gamma and sigma; the '
                         'collision term for species i then carries the multiplicity of i in the channel',
        resonances='RIS (narrow-width) subtraction: for every internal line that can be on shell in the physical region (s from a '
                   'neutral X Xb pair, X from a charge +-1 subset) |D|^2 -> Re(D^2) with D = 1/(p^2 - m^2 + i m Gamma_reg), '
                   'including nested and disjoint multi-line terms by inclusion-exclusion; the removed on-shell pieces are the '
                   'sequential chains that the network carries explicitly (1<->2, 2<->2, 3<->2 with real s/X).  The 2->2 '
                   'quadrature uses the physical width (XXb_XXb s-pole) and Gamma_X,reg = 1e-7 m_X (Xs_Xs u-pole); the MC '
                   'uses Gamma_reg = 0.02 m_s and 0.05 m_X (primary) and 0.04/0.10 (secondary): the spread is the regulator '
                   'systematic of the finite part; the regulated full column scales as 1/Gamma_reg where on-shell lines dominate',
        in_medium='kernels at vacuum masses; masses(T) and the phase in masses_T; validity in masses_T.validity',
        provenance='Codex packets read only via git show/archive (9877518 KINETIC_SPEC.md, MATHEMATICS.md; 225028b INPUTS.json); no Codex code run; '
                   'own earlier tracks reused: Berends-Giele method of gp_leftovers/amp42.py (re-implemented), threshold constants [Q own] for comparison')
    # masses(T)
    mt = {}
    for scheme in ('oneloop', 'daisy'):
        rows = thermal[scheme]['rows']
        mt[scheme] = cols(rows, ['T_over_mX', 'T_GeV', 'vT_over_v', 'mX_GeV', 'ms_GeV', 'lam3_GeV', 'lam_s3_GeV', 'broken', 'ms_over_2mX',
                                 'Gamma_s_rest_medium_GeV', 'bose_enh_rest', 'decay_open', 'Gamma_s_vacformula_medium_masses_GeV'])
        mt[scheme].update(T_c_over_mX=thermal[scheme]['T_c_over_mX'], T_close_s_decay_over_mX=thermal[scheme]['T_close_s_decay_over_mX'],
                          vacuum_mass_deviation_crossings=thermal[scheme]['vacuum_mass_deviation_crossings'])
    mt['validity'] = dict(vacuum_masses_within_1pct='T <= %.2f m_X (m_s), %.2f m_X (m_X)' % (
        thermal['oneloop']['vacuum_mass_deviation_crossings']['ms_dev_1pct_T_over_mX'], thermal['oneloop']['vacuum_mass_deviation_crossings']['mX_dev_1pct_T_over_mX']),
        vacuum_masses_within_5pct='T <= %.2f m_X (m_s), %.2f m_X (m_X)' % (
        thermal['oneloop']['vacuum_mass_deviation_crossings']['ms_dev_5pct_T_over_mX'], thermal['oneloop']['vacuum_mass_deviation_crossings']['mX_dev_5pct_T_over_mX']),
        broken_phase='T < T_c = %.2f m_X (one loop), %.2f m_X (daisy); no objects and no s <-> X Xb above' % (thermal['oneloop']['T_c_over_mX'], thermal['daisy']['T_c_over_mX']),
        s_decay_open='T < %.2f m_X (one loop), %.2f m_X (daisy): m_s(T) > 2 m_X(T)' % (thermal['oneloop']['T_close_s_decay_over_mX'], thermal['daisy']['T_close_s_decay_over_mX']),
        scheme='one-loop Dolan-Jackiw, no daisy (primary) / Parwani daisy (bracket); mu = 0 in the loops; CW piece dropped (6e-6)')
    mt['width_vs_energy_oneloop'] = thermal['width_vs_energy_oneloop']
    mt['formulas'] = dict(hiT=thermal['hiT_estimate']['formula'], width_medium='Gamma(E,T) = (m_s/E) Gamma_vac(m_s(T), m_X(T), lam3(T)) <1 + f(E1) + f(E2)>_angle; '
                          'inverse decay by detailed balance (thermal.width_s_medium, thermal_mc.inverse_decay_bose)')
    out['masses_T'] = mt
    # decay
    dec = k22['decay']
    bose_dec = {}
    for row in mc['thermal_rates_and_balance']:
        k = 'db_s_XXb_BE_mup0.0_mum0.0'
        if k in row:
            bose_dec[str(row['T_over_mX'])] = dict(decay_BE_over_MB=row[k]['decay_mc'] / row[k]['decay_MB_analytic'], inverse_BE_quadrature_over_MB=row[k]['inverse_quad'] / row[k]['decay_MB_analytic'])
    out['decay_s_XXb'] = dict(Gamma_vac_GeV=dec['Gamma_s_vac_GeV'], formula=dec['formula'], MB=cols(dec['rows'], ['T_over_mX', 'n_s_eq', 'gamma_decay_MB', 'Gamma_eff_over_Gamma_vac']),
                              bose_over_MB_mu0=bose_dec, amplitude='M = -lam3 (constant); Gamma_vac = lam3^2 beta/(16 pi m_s)')
    # 2->2
    c22 = {}
    for ch, c in k22['sigma'].items():
        keys = [k for k in c['rows'][0] if k != 's_GeV2']
        th = k22['thermal'][ch]['rows']
        tkeys = [k for k in th[0]]
        bose = {}
        for row in mc['thermal_rates_and_balance']:
            key = ch
            bose[str(row['T_over_mX'])] = dict(bose_over_MB=row[key]['bose_over_MB'], mode=row[key]['mode'], mc_over_quad_MB=row[key]['mc_over_quad'])
            if ch + '_nores' in row:
                bose[str(row['T_over_mX'])]['bose_over_MB_nores'] = row[ch + '_nores']['bose_over_MB']
        c22[ch] = dict(Sin=c['Sin'], Sout=c['Sout'], masses_GeV=c['masses'], sqrt_s_th_over_mX=c['sqrt_s_th_over_mX'], sigma_table=cols(c['rows'], keys),
                       thermal_MB=cols(th, tkeys), bose_over_MB_mu0=bose, high_energy=c['high_energy'])
    c22['XXb_XXb']['prescription'] = ('sigma = sigma0 + sigma1 Re D + sigma2 R, D = 1/(s - ms^2 + i ms Gamma_s): R = |D|^2 (full, includes the on-shell '
                                      'sequential piece = inverse decay x decay) or R = Re(D^2) (subtracted; thermal average via finite-part quadrature, '
                                      'gamma_sub, may be negative and is <= 2e-3 of gamma_ris); gamma_ris == gamma_decay_MB (checked)')
    c22['Xs_Xs']['prescription'] = ('u-channel X pole inside the angular range for 3.606 < sqrt s/m_X < 5.0 (s -> X Xb then Xb X -> s): sigma_sub is the '
                                    'finite part (sign-indefinite there); for transport use sigma_nores (+ the chain through decay/inverse decay)')
    out['channels_22'] = c22
    # n->2
    cn2 = {}
    for ch, c in n2['channels'].items():
        c = merge_supplement(c, supp.get(ch))
        tab = TB.Table(c, 'n2'); tab04 = TB.Table(c, 'n2', freg=0.04)
        rows02 = sorted([r for r in c['rows'] if abs(r['freg'] - 0.02) < 1e-12], key=lambda r: r['sqrt_s_over_mX'])
        rows04 = sorted([r for r in c['rows'] if abs(r['freg'] - 0.04) < 1e-12], key=lambda r: r['sqrt_s_over_mX'])
        cn2[ch] = dict(legs=c['legs'], Sin=c['Sin'], Sout=c['Sout'], m_in=c['m_in'], m_out=c['m_out'], sqrt_s_th_over_mX=c['sqrt_s_th'] / M.MX,
                       Gt_table_reg0_02=dict(sqrt_s_over_mX=[r['sqrt_s_over_mX'] for r in rows02], G_sub=[r['G'] for r in rows02], G_sub_err=[r['err'] for r in rows02],
                                             G_full_reg=[r['G_full'] for r in rows02], G_full_err=[r.get('err_full') for r in rows02], lines_flagged=[r['lines_flagged'] for r in rows02], N=[r['N'] for r in rows02]),
                       Gt_table_reg0_04=dict(sqrt_s_over_mX=[r['sqrt_s_over_mX'] for r in rows04], G_sub=[r['G'] for r in rows04], G_sub_err=[r['err'] for r in rows04], G_full_reg=[r['G_full'] for r in rows04]),
                       chain_threshold_sqrt_s_over_mX=tab.chain_threshold(c), chain=CHAIN.get(ch, ''), threshold_power_law_exponent_full=tab.k_th,
                       thermal_MB=thermal_block(tab, tab04), extrapolation='below first point: power law of G_full between the first two points; above 400 m_X: constant')
    out['channels_n2'] = cn2
    # 2->n
    c2n = {}
    for ch, c in twon['channels'].items():
        tab = TB.Table(c, '2n'); tab04 = TB.Table(c, '2n', freg=0.04) if any(abs(r['freg'] - 0.04) < 1e-12 for r in c['rows']) else None
        rows02 = sorted([r for r in c['rows'] if abs(r['freg'] - 0.02) < 1e-12], key=lambda r: r['sqrt_s_over_mX'])
        c2n[ch] = dict(legs=c['legs'], Sin=c['Sin'], Sout=c['Sout'], m_in=c['m_in'], m_out=c['m_out'], sqrt_s_th_over_mX=c['sqrt_s_th'] / M.MX,
                       sigma_table_reg0_02=dict(sqrt_s_over_mX=[r['sqrt_s_over_mX'] for r in rows02], sigma_sub_GeV2=[r['sigma'] for r in rows02], sigma_sub_err=[r['err'] for r in rows02],
                                                sigma_full_reg_GeV2=[r['sigma_full'] for r in rows02], lines_flagged=[r['lines_flagged'] for r in rows02]),
                       thermal_MB=thermal_block(tab, tab04), high_energy=dict(sigma_times_s_at_top=rows02[-1]['sigma_full'] * (rows02[-1]['sqrt_s_over_mX'] * M.MX)**2,
                                                                              top_sqrt_s_over_mX=rows02[-1]['sqrt_s_over_mX']))
    out['channels_2n'] = c2n
    # MC checks, hard particle, energy exchange
    out['hard_particle_on_bath'] = dict(note='test X of lab energy E on a Bose bath at T (mu = 0); rate per bath density = <sigma v_rel> (GeV^-2); dE/dt per density (GeV^-1); '
                                             'channels X (XX->XX), Xb (XXb->XXb non-resonant), s (Xs->Xs subtracted; above 5 m_X no pole)', rows=mc['hard_particle'])
    out['energy_exchange_Xs'] = dict(note='MB species at T_X, T_s; Q_s = energy gained by s per volume per time (GeV^5) per n_X n_s (non-resonant kernel)', rows=mc['energy_exchange_Xs_MB'])
    out['checks'] = dict(point_kinematics=mc['point_kinematics'], elastic_energy_moment=mc['elastic_energy_moment'],
                         bose_detailed_balance=[{k: v for k, v in row.items() if k.startswith('db_') or k == 'T_over_mX'} for row in mc['thermal_rates_and_balance']],
                         table_detailed_balance=bal.get('detailed_balance'), nr_limits=bal.get('nr_limits'), old_s_route=bal.get('old_s_route'))
    out['usage'] = dict(
        solver_recipe=['masses: vacuum for T <= 3 m_X; masses_T tables beyond; s <-> X Xb switched off above T_close',
                       '1<->2: Gamma(E,T,mu) from masses_T.formulas; inverse decay by detailed balance',
                       '2->2 elastic: sigma_nores (XXb, Xs) / sigma_sub (XX, ss) tables; transport columns for momentum relaxation',
                       '2<->2 inelastic ss<->XXb: sigma_sub (no pole); Bose factors from the solver',
                       'n<->2: the on-shell chains are the explicit 1<->2, 2<->2, 3<->2 channels; add the subtracted remainder G_sub where remainder_resolved, '
                       'otherwise treat it as zero with remainder_bound as the systematic; below the chain threshold G_sub is the whole kernel',
                       'late limit: thermal_MB tables -> A_old n^4/S_in as T -> 0 (checks.nr_limits)'])
    json.dump(out, open(f'{RES}/rates_interface.json', 'w'), indent=1, default=float)
    import os
    print('written', os.path.getsize(f'{RES}/rates_interface.json') / 1e6, 'MB')


if __name__ == '__main__':
    main(*sys.argv[1:5])
