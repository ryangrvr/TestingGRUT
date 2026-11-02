#!/usr/bin/env python3
"""mf_chain SKEPTIC, attack on mf_ideal's map: mf_ideal scanned the holder amplitude c0 only in 0.9-1.2 and at one corner
(K = 1, m_s = 3, t_on = 0.5), and calls Q-matter at 1e-6 a 'narrow pocket needing a holder 3-5% stronger'.  Here mf_ideal's
own idealised histories (built-in dust fluid tuned to K, smootherstep turn-on) are re-run with mf_ideal's engine, vectors and
ledger (read-only; outputs in the skeptic's scratch) at c0 = 1.3-1.5.
Usage: python3 -B sk_idealc0.py K MS TON c0 [c0 ...]
Labels: [C] given the model choices.
"""
import json, sys
import sk_common as S


def run(K, ms, ton, c0):
    tag = 'id_K%g_m%g_t%g_C%.3f' % (K, ms, ton, c0)
    S.rm.release(tag, K, ms, ton, gs=0.0, sstar=1e16, c0=c0, bg='fluid', Tl=1.0e4)
    S.vm.main(tag)
    o = S.cm.analyse(tag)
    S.dump(S.ENG / 'chain' / ('chain_%s.json' % tag), o)
    a, b = o['e1e-6'], o['e3e-6']
    print('%s: eps %.4f E/Ec %.3f N_YL %.4g N_YH %.4g | 1e-6: Rv %.3g RY %.3g Rc %.3g S/Q %.3f %s H_req %.3e | 3e-6 %.3f | inf %.3f | e_crit %s / %s' % (
        tag, o['eps'], o['E_over_Ec'], o['N_YL'], o['N_YH'], a['R_vec'], a['R_Y'], a['R_casc'], a['S_per_Q'], a['branch'][:8],
        a['H_i_req_GeV'], b['S_per_Q'], o['S_per_Q_e_inf'], o['e_crit_bulk'], o['e_crit_Q2e5']), flush=True)
    return o


if __name__ == '__main__':
    K, ms, ton = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
    for c in sys.argv[4:]:
        run(K, ms, ton, float(c))
