#!/usr/bin/env python3
"""inf_chain (Claude, own code): ESTIMATE [E] of inflaton -> A_L A_L through the induced phi_I - sigma mixing, with the
dq_heat long-lived-vector ledger (each A_L deposits m_A <a_A>, not its production energy).

Induced portal (dq_entropy DQ2 form): V ⊃ c (|u|^2 - mu^2)(|Sigma|^2 + 2|Y|^2), u = mu + G phihat, so
A phihat |Sigma|^2 with A = 2 mu Re(G) c;  Sigma = (h + sigma)/sqrt2 gives the mixing mass^2 A h, angle
theta = A h/(m_sigma^2 - m_I^2), and sigma (d pi)^2/h gives phi -> pi pi (pi = A_L, Goldstone equivalence) with
    Gamma_leak = A^2 m_I^3 / (32 pi (m_sigma^2 - m_I^2)^2)        (inf_reheat SKEPTIC form, m_sigma = 40 H_i).
c is the supplied renormalisation condition [model choice]: we bracket
    (a) c = 1/(16 pi^2) (inf_reheat skeptic's 'natural size', L = 1),  (b) c = 1/(4 pi^2) (dq_entropy's unit, L = 1:
    threshold size when the portals vanish at M_Y; below M_Y the Y loops decouple),  (c) c = 13.31/(4 pi^2) (portals vanish
    at M_P).
Number of inflaton quanta that decay: TB* (matter era) N_I = x0/(eta2 m_I); Gamma = 0.5 runs N_I ≈ 3 K^2/(eta2 m_I a_d)
with a_d = a(tau = 1/Gamma) from the dense run [E].  Leak heat R_leak = 2 (Gamma_leak/Gamma_I) N_I m_A <a_A>/E_X at
e = 1e-6 (deposit ∝ e^{-1/2}).  The effect on S/Q, p, M and H_i,req is given at fixed Gamma_I (no re-calibration).
Output ../results/leak_estimate.json.  Label: [E] throughout.
"""
import json, math
import numpy as np
import ch_common as ch
import chain as cn

MSIG = 40.0
CONV = {'a_skeptic_L1': 1 / (16 * math.pi ** 2), 'b_MY_L1': 1 / (4 * math.pi ** 2), 'c_MP': 13.31 / (4 * math.pi ** 2)}
CASES = {'HT6t': 'HT6t', 'HT6t_cal1': 'HT6t_cal1', 'T1': 'T1', 'T1_cal1': 'T1_cal1', 'TB1': 'TB1_g7.512e-09', 'TB2': 'TB2_g2.113e-09'}


def main():
    out = dict(method=__doc__, rows={})
    for key, vt in CASES.items():
        rel = json.load(open(ch.RES / ('release_%s.json' % key)))
        ch_ = json.load(open(ch.RES / ('chain_%s.json' % vt)))
        b = rel['bench']
        G = complex(*b['G']); mI = b['m_I_over_Hi']; eta2 = b['eta2']; gam = ch_['gam'] or b['Gamma_I_over_Hi']
        EX = rel['result']['E_X']; K = ch_['K']
        v6 = [x for x in json.load(open(ch.RES / ('vectors_%s.json' % vt)))['by_e'] if abs(x['e'] - 1e-6) < 1e-12][0]
        DA = v6['deposit_per_vector']
        if 'x0' in rel['result']:
            NI = rel['result']['x0'] / (eta2 * mI)
        else:
            S = np.load(ch.RUNS / (key + '_dense.npz'))['S']
            ad = float(np.interp(1.0 / gam, S[:, 0], S[:, 1]))
            NI = 3 * K ** 2 / (eta2 * mI * ad)
        r = dict(m_I=mI, Gamma_I=gam, N_I=NI, deposit_A_e1e6=DA, E_X=EX, R_total_e1e6=ch_['e1e6']['R_total'], S_per_Q_e1e6=ch_['e1e6']['S_per_Q'])
        for lab, c in CONV.items():
            A = 2 * ch.MU * G.real * c
            Gl = A * A * mI ** 3 / (32 * math.pi * (MSIG ** 2 - mI ** 2) ** 2)
            BR = Gl / gam
            Rl = 2 * BR * NI * DA / EX
            Rt = ch_['e1e6']['R_total'] + Rl
            sq = cn.SQ(ch_['E_over_Ec'] * (1 + Rt))
            p = cn.p_of_SQ(sq) * cn.sc_ratio(sq)
            Hreq = cn.H_of(1 + 2 * p, ch_['Y_over_YQ_at_HM1'])
            r[lab] = dict(c=c, A=A, Gamma_leak=Gl, BR=BR, R_leak_e1e6=Rl, S_per_Q_e1e6=sq, p=p, H_i_req_fixed_Gamma=Hreq)
        out['rows'][key] = r
        print('%-10s m_I=%.2f Gamma=%.2e N_I=%.2e | ' % (key, mI, gam, NI) + ' | '.join(
            '%s: BR=%.1e R_leak=%.3g S/Q=%.4g p=%.3g H_req=%.3e' % (l, r[l]['BR'], r[l]['R_leak_e1e6'], r[l]['S_per_Q_e1e6'], r[l]['p'], r[l]['H_i_req_fixed_Gamma'])
            for l in CONV))
    json.dump(out, open(ch.RES / 'leak_estimate.json', 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
