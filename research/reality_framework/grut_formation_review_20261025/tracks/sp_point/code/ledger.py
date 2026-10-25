"""sp_point: Y-excitation ledger at the leading point from our own integration (repro.json, longrun*.json).

Re-derives: Gamma(Y->3X), a at decay (H = Gamma and t = 1/Gamma conventions), the Y/X energy ratio at decay,
pairs per net charge (Codex J_Y, eigenmode J_Y, N_X or N_eff denominator), ratio to the quoted 1% allowance,
the energy-equivalent quartic fraction, and two decay-channel estimates Codex did not state:
  * condensate-assisted two-body Y -> X X via 3 kappa X_bg (rate / H, constant in radiation),
  * parametric (narrow) resonance parameter of the Y condensate acting on X fluctuations.
Usage: python3 ledger.py ../results/repro.json ../results/ensembles.json ../results/ledger.json
"""
import json
import sys
import numpy as np
import spmodel as m

ALLOW = 0.0049296413205678555   # quoted pair allowance per net charge at 1% (2793dab / b5174e5) [Q]


def main(rp, ep, out):
    R = json.load(open(rp))
    w = R['runs']['Y8_c0.00625']['windows']['600-1000']
    H = m.H_STAR
    M, h, mu = m.pars(H)
    MY = m.M_GEV
    G3 = 3 * m.KAP**2 * MY / (256 * np.pi**3)
    a_dec_H = np.sqrt(H / G3)            # H(a) = H_i/a^2 = Gamma
    a_dec_t = np.sqrt(2 * H / G3)        # t = 1/(2H) = 1/Gamma
    NX, EX, JY, JYe, Neff = w['NX'], w['EX'], w['JY'], w['JYeig'], w['Neff']
    ratio_E = lambda a, J: J * a * M / EX
    out_d = dict(
        Gamma_Y3_GeV=G3, a_decay_H_eq_Gamma=a_dec_H, a_decay_t_eq_1_over_Gamma=a_dec_t,
        EY_over_EX_at_decay_codexJ=ratio_E(a_dec_H, JY), EY_over_EX_at_decay_eigJ=ratio_E(a_dec_H, JYe),
        EY_over_EX_at_decay_t_convention=ratio_E(a_dec_t, JY),
        pairs_codex=1.5 * JY / abs(NX), pairs_eig=1.5 * JYe / abs(NX), pairs_eig_Neff=1.5 * JYe / abs(Neff),
        allowance=ALLOW)
    out_d['ratio_to_allowance_codex'] = out_d['pairs_codex'] / ALLOW
    out_d['ratio_to_allowance_eig'] = out_d['pairs_eig'] / ALLOW
    out_d['ratio_to_allowance_eig_Neff'] = out_d['pairs_eig_Neff'] / ALLOW
    eps = w['epsX']
    out_d['energy_equivalent_eps'] = eps * (1 / (1 + ratio_E(a_dec_H, JY)))**0.75
    # condensate-assisted two-body channel: L ~ 3 kappa X_bg Y* dX^2 + h.c.; amplitude 6 kappa X_bg, identical pair:
    # Gamma2 = 9 kappa^2 |X_bg|^2 / (8 pi M_Y); <|X_bg|^2> ~ A^2/2 for the quartic oscillation [E]
    chiX = (EX / m.LX)**0.25                # conformal quartic amplitude
    Xphys_over_a = chiX * H                 # |X| = chiX H / a  (GeV, times 1/a)
    G2_over_H = 9 * m.KAP**2 * Xphys_over_a**2 / (8 * np.pi * MY) / H   # both scale as a^-2
    out_d['two_body_condensate_Gamma_over_H_rms_amplitude'] = G2_over_H / 2  # <|X|^2> ~ A^2/2
    out_d['two_body_integrated_to_decay'] = G2_over_H / 2 * np.log(a_dec_H / 10.)  # from a = 10 (end of strong coupling) to decay
    # narrow parametric resonance from the Y condensate: q ~ 6 kappa |d||X| / M_Y^2 (physical), |d|^2 = J_Y/(a M)
    a = 1000.
    d_conf = np.sqrt(JY / (a * M))
    q = 12 * m.KAP * (d_conf / a * H) * (chiX / a * H) / MY**2   # Mathieu q for resonance at M_Y/2 [E, factor-2 level]
    out_d['parametric_q_at_a1000'] = q
    # visible-sector comparison: rho_X / rho_vis at release, and Y energy at decay as dark radiation
    rX = EX * H**2 / (3 * m.MRED**2)
    out_d['rhoX_over_rhovis'] = rX
    out_d['rhoY_over_rhovis_at_decay'] = rX * ratio_E(a_dec_H, JY)
    try:
        E = json.load(open(ep))['results']
        ens = {k: v['stats'].get('pairs_per_net_charge') for k, v in E.items() if 'stats' in v}
        out_d['ensemble_pairs_per_net_charge'] = ens
        out_d['ensemble_ratio_to_allowance'] = {k: (v / ALLOW if v else None) for k, v in ens.items()}
    except Exception as ex:   # ensembles optional
        out_d['ensemble_error'] = str(ex)
    for k, v in out_d.items():
        print(k, v)
    json.dump(out_d, open(out, 'w'), indent=1)


if __name__ == '__main__':
    main(*sys.argv[1:4])
