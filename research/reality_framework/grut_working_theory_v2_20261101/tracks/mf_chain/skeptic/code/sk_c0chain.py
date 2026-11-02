#!/usr/bin/env python3
"""mf_chain SKEPTIC: the full chain (mf_ideal's engine, vectors, ledger, p map; read-only, outputs in the skeptic's scratch)
on a bench history with the holder amplitude c0 != 1, i.e. mu_eff(slice) = mu + c0 (10i - mu) (c0 = 1 is M1r's 10i
calibration, which mf_hybrid copied into |g_h| s_* = 12.30 H_i; the holder strength itself is a free model choice).
The preparation is the held stationary point at that u (checked against the skeptic's own continuation solver:
amplitudes agree to 1e-6).  Usage: python3 -B sk_c0chain.py TABLE_TAG c0 [c0 ...]
  TABLE_TAG: tag of an existing bench_run in the skeptic's scratch (its MODE 5 table benchtab_<tag>.txt and release are reused).
Labels: [C] given the model choices.
"""
import json, math, sys
import sk_common as S


def run(base, c0):
    rel0 = json.load(open(S.ENG / 'release' / ('release_%s.json' % base)))
    sp = rel0['spec']
    orig = S.rm.spec_env

    def spec_env(*a, **kw):
        env, gam = orig(*a, **kw)
        env['MFI_SM2TAB'] = '1'
        return env, gam
    S.rm.spec_env = spec_env
    tag = '%s_C%.3f' % (base, c0)
    try:
        S.rm.release(tag, sp['K_target'], sp['m_s'], 0.0, gs=sp['Gamma_s'], sstar=sp['s_star_GeV'], c0=c0, bg=sp['bg'],
                     Tl=rel0['result']['T_long'], bcap=sp['step_cap_b'])
    finally:
        S.rm.spec_env = orig
    S.vm.main(tag)
    o = S.cm.analyse(tag)
    S.dump(S.ENG / 'chain' / ('chain_%s.json' % tag), o)
    a, b = o['e1e-6'], o['e3e-6']
    print('%s: K %.3f eps %.4f E/Ec %.3f N_YL %.4g N_YH %.4g | 1e-6: Rv %.3g RY %.3g Rc %.3g S/Q %.3f %s M %.3g H_req %.3e | 3e-6: S/Q %.3f %s | inf %.3f | e_crit %s / %s' % (
        tag, o['K'], o['eps'], o['E_over_Ec'], o['N_YL'], o['N_YH'], a['R_vec'], a['R_Y'], a['R_casc'], a['S_per_Q'], a['branch'], a['M'],
        a['H_i_req_GeV'], b['S_per_Q'], b['branch'], o['S_per_Q_e_inf'], o['e_crit_bulk'], o['e_crit_Q2e5']), flush=True)
    return o


if __name__ == '__main__':
    base = sys.argv[1]
    for c in sys.argv[2:]:
        run(base, float(c))
