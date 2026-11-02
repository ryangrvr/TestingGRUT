#!/usr/bin/env python3
"""mf_higgs SKEPTIC: R < 0 windows after eps_H = 1 on the raw homogeneous R of three histories (bench export, mf_chain
unblocked and Pauli rebuilds), with the tachyonic exponent at xi(H_i) and at large xi, the Mathieu q of the 2 m_chi
line from the raw R (fit of the oscillating part), and the R-window exponent that would compound only if q >~ 1.
[C on the exported histories]  Output: results/sk_windows.json
"""
import sys
sys.dont_write_bytecode = True
import math

import numpy as np

import sk_common as C

HP = {'EXP': C.TRACKS / 'mf_hybrid' / 'results' / 'bench_HFa.json',
      'MC_TR': C.TRACKS / 'mf_chain' / 'results' / 'histories' / 'bench_itHFa_track_it1.json',
      'MC_P0': C.TRACKS / 'mf_chain' / 'results' / 'histories' / 'bench_itHFa_P0_it2.json'}


def main():
    out = {}
    for lab, p in HP.items():
        d, h = C.load_hist(p)
        te = d['waterfall']['t_eps1']
        tsw = float(d['reheating'].get('t_sw', 0.6))
        m = (h['t'] >= te) & (h['t'] <= tsw)
        t, R, H = h['t'][m], h['R'][m], h['H'][m]
        tf = np.arange(t[0], t[-1], 2e-5)
        Rf = np.interp(tf, t, R)
        Hf = np.interp(tf, t, H)
        neg = Rf < 0
        wins = []
        i = 0
        while i < len(tf):
            if neg[i]:
                j = i
                while j < len(tf) and neg[j]:
                    j += 1
                seg = slice(i, j)
                k = i + int(np.argmin(Rf[seg]))
                wins.append(dict(t0=float(tf[i]), dur=float(tf[min(j, len(tf) - 1)] - tf[i]), Rmin_over_H2=float(Rf[k] / Hf[k] ** 2),
                                 Rmin=float(Rf[k]), I1=float(np.trapezoid(np.sqrt(np.maximum(-Rf[seg], 0)), tf[seg]))))
                i = j
            else:
                i += 1
        xi = 0.08773
        I1 = np.array([w['I1'] for w in wins])
        # oscillating part of R about its running mean -> amplitude and dominant frequency
        Rm = np.convolve(Rf, np.ones(650) / 650, mode='same')
        osc = (Rf - Rm)[700:-700]
        spec = np.abs(np.fft.rfft(osc * np.hanning(len(osc))))
        fr = np.fft.rfftfreq(len(osc), 2e-5) * 2 * math.pi
        om = float(fr[np.argmax(spec[5:]) + 5])
        amp = float(np.sqrt(2) * np.std(osc))
        o = dict(t_eps1=te, t_sw=tsw, n_windows=len(wins), first_windows=wins[:4],
                 frac_time_R_negative=float(np.mean(neg)), Rmin_over_H2_min=float(min(w['Rmin_over_H2'] for w in wins)),
                 max_window_exponent_xiHi=float(math.sqrt(xi) * I1.max()), sum_exponent_xiHi=float(math.sqrt(xi) * I1.sum()),
                 max_window_exponent_xi1000=float(math.sqrt(1000) * I1.max()),
                 R_line_omega=om, R_line_amp_rms_based=amp,
                 mathieu_q_xiHi=float(2 * xi * amp / om ** 2), mathieu_q_xi100=float(2 * 100 * amp / om ** 2),
                 mathieu_q_xi1000=float(2 * 1000 * amp / om ** 2))
        out[lab] = o
        print(lab, {k: (round(v, 5) if isinstance(v, float) else v) for k, v in o.items() if k != 'first_windows'})
        print('    first', [(round(w['t0'], 4), round(w['dur'], 4), round(w['Rmin_over_H2'], 2)) for w in wins[:4]])
    C.jdump(out, C.RES / 'sk_windows.json')


if __name__ == '__main__':
    main()
