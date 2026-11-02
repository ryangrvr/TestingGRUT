#!/usr/bin/env python3
"""dq_heat SKEPTIC (Claude, own code, written fresh): semi-analytic longitudinal-vector abundance.

A method different from dq_heat/code/vector_modes.py (which integrates every mode from k/aH = 1000 through the recorded
background).  Here:
  1. Frozen Goldstone angle.  For k well outside the horizon at the slice the inflationary common-phase fluctuation is
     |theta_0|^2 = H^2 / (2 k^3 F_inf^2) (massless Goldstone in de Sitter, m_A << H), and a global U(1) rotation is a
     symmetry of the full background dynamics, so theta stays frozen through the release whatever F(eta) does.
  2. Universal radiation-era evolution with F = F_vac, a = K (eta - eta_0), M = a m (m = e F_vac).  With s = (K m)^{1/2},
     x = s (eta - eta_0), q = k/s, the longitudinal (Stueckelberg) equation is
         d/dx [ x^2/(q^2 + x^2) theta_x ] + x^2 theta = 0,   theta(0) = 1, theta_x(0) = 0,
     and the late number is n_k = (F_vac/F_inf)^2 K^2 H^2 C(q)^2 / (4 s^4 q), C^2 = lim (theta_x^2 + x^2 theta^2)/x
     (read with the canonical u = sqrt(P) theta adiabatic invariant; C^2 -> q^-3 for q >> 1).
  3. N_A = int d^3k/(2pi)^3 n_k = (F_vac/F_inf)^2 K^{3/2} H^2 m^{-1/2}/(8 pi^2) * int q C^2 dq   [D given 1-2].
This predicts N_A proportional to K^{3/2} (F_vac/F_inf)^2 e^{-1/2} exactly for the frozen part, and gives n_k mode by mode
for comparison with dq_heat's spectra (valid for k << 1, i.e. modes far outside the horizon at the slice).
Output: ../results/sk_frozen.json.   Labels: [C] computed, [D] derived.
"""
import json, math, time
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

HERE = Path(__file__).resolve().parent
RES = HERE.parent / 'results'
DQ = HERE.parents[1] / 'results'


def C2_of_q(q, rtol=1e-11):
    """late adiabatic invariant C^2(q) of the universal radiation-era longitudinal mode (theta(0) = 1)."""
    x0 = 1e-3 * min(q, 1.0, 1.0 / max(q, 1e-30))
    # series start: theta = 1 - q^2 x^2/6 (x << q) or 1 - x^4/12 (q << x << 1); both tiny at x0
    P0 = x0 * x0 / (q * q + x0 * x0)
    th0 = 1.0 - q * q * x0 * x0 / 6.0 if x0 < q else 1.0 - x0 ** 4 / 12.0
    thx0 = -q * q * x0 / 3.0 if x0 < q else -x0 ** 3 / 3.0
    y0 = [th0, P0 * thx0]

    def rhs(x, y):
        P = x * x / (q * q + x * x)
        return [y[1] / P, -x * x * y[0]]

    def readout(x, y):
        P = x * x / (q * q + x * x)
        th, thx = y[0], y[1] / P
        sP = math.sqrt(P); dsP = (q * q / (q * q + x * x) ** 1.5)          # d sqrt(P)/dx = q^2/(q^2+x^2)^{3/2}
        u, ux = sP * th, sP * thx + dsP * th
        w = math.sqrt(q * q + x * x)
        return (ux * ux + w * w * u * u) / (2 * w)                         # = (A/2)^{-1} n (A = K^2 F^2 q^2/s)

    # end: well into the adiabatic regime.  For q >> 1 the mode is relativistic inside the horizon at q x >> 1 and the
    # u-invariant is conserved through the later NR transition; for q <~ 1 wait until x >> 1 (H << m).
    if q > 3:
        xs = [200.0 / q, 400.0 / q]          # q x = 200, 400: deep inside the horizon; (sqrt P)''/sqrt P / w^2 <~ (q x)^-2
    else:
        xs = [max(20.0, 20 * q), max(40.0, 40 * q)]
    sol = solve_ivp(rhs, (x0, xs[1]), y0, method='DOP853', rtol=rtol, atol=1e-14, dense_output=True)
    vals = []
    for xx in xs:
        y = sol.sol(xx)
        # the u-number in x units equals (theta_x^2 + x^2 theta^2)/(2x) late; C^2 = 2 * that
        vals.append(2 * readout(xx, y))
    return vals[1], abs(vals[1] - vals[0]) / abs(vals[1])


def main():
    t0 = time.time()
    qg = np.geomspace(1e-4, 1e4, 321)
    C2 = []; chg = []
    for q in qg:
        c, d = C2_of_q(q)
        C2.append(c); chg.append(d)
    C2 = np.array(C2)
    lq = np.log(qg)
    integ = qg * qg * C2                     # q C^2 dq = q^2 C^2 dlnq
    # analytic tails: q << 1 : C^2 -> const (checked below), q >> 1 : C^2 -> q^-3
    I_num = float(np.trapezoid(integ, lq))
    I_tail_hi = 1.0 / qg[-1]                 # int_{qmax}^inf q q^-3 dq
    I = I_num + I_tail_hi
    out = dict(method=__doc__, q=qg.tolist(), C2=C2.tolist(), max_readout_change=float(max(chg)),
               C2_small_q=float(C2[0]), q3C2_large_q=float(qg[-1] ** 3 * C2[-1]), I=I, I_peak_q=float(qg[np.argmax(integ)]))
    print('I = int q C^2 dq = %.6f ; C2(q->0)=%.5f ; q^3 C2 (q=1e4) = %.5f ; max readout change %.1e' % (I, C2[0], qg[-1] ** 3 * C2[-1], max(chg)))

    # predictions in the dq_heat frame
    bg = json.load(open(DQ / 'background.json'))
    K = bg['late']['K_end']; Fi = bg['prep']['F_inf']; Fv = 240.0
    vm = json.load(open(DQ / 'vector_modes.json'))
    KG = np.array(vm['kgrid'])
    Fnox = math.sqrt(bg['prep']['F2_inf'] * (1 - bg['prep']['X_part_of_F2']))
    cmp = []
    for c in vm['cases']:
        e, var = c['e'], c['variant']
        KK = 1.0 if var == 'inst' else K
        FI = Fnox if var == 'noX' else Fi
        m = e * Fv; s = math.sqrt(KK * m)
        pref = (Fv / FI) ** 2 * KK ** 1.5 / (8 * math.pi ** 2 * math.sqrt(m))
        # frozen-part number with the same k <= 1 cut as dq_heat
        qcut = 1.0 / s
        sel = qg <= qcut
        I_cut = float(np.trapezoid(integ[sel], lq[sel])) + float(np.interp(math.log(qcut), lq, integ)) * (math.log(qcut) - lq[sel][-1])
        N_frozen_all = pref * I
        N_frozen_kle1 = pref * I_cut
        # mode-by-mode comparison n_k (theirs) / n_k (frozen) for k <= 0.3
        rows = []
        for r in c['spectrum']:
            k = r['k']; q = k / s
            C2k = math.exp(np.interp(math.log(q), lq, np.log(C2)))
            n_fr = (Fv / FI) ** 2 * KK ** 2 * C2k / (4 * s ** 4 * q)
            rows.append(dict(k=k, n_theirs=r['n_k'], n_frozen=n_fr, ratio=r['n_k'] / n_fr))
        lowk = [x['ratio'] for x in rows if 2e-3 <= x['k'] <= 0.05]
        midk = [x['ratio'] for x in rows if 0.05 < x['k'] <= 0.3]
        # their N_A over k <= 0.3 vs frozen over k <= 0.3
        ks = np.array([x['k'] for x in rows]); nt = np.array([x['n_theirs'] for x in rows]); nf = np.array([x['n_frozen'] for x in rows])
        s3 = ks <= 0.3
        N_t03 = float(np.trapezoid(ks[s3] ** 3 * nt[s3] / (2 * math.pi ** 2), np.log(ks[s3])))
        N_f03 = float(np.trapezoid(ks[s3] ** 3 * nf[s3] / (2 * math.pi ** 2), np.log(ks[s3])))
        cc = dict(e=e, variant=var, K=KK, F_inf=FI, N_theirs_total=c['N_A'], N_theirs_k_le_1=c['N_A_k_le_1'],
                  N_frozen_all_k=N_frozen_all, N_frozen_k_le_1=N_frozen_kle1, ratio_k_le_1=c['N_A_k_le_1'] / N_frozen_kle1,
                  N_theirs_k_le_0p3=N_t03, N_frozen_k_le_0p3=N_f03, ratio_k_le_0p3=N_t03 / N_f03,
                  ratio_n_k_2em3_to_0p05=[min(lowk), max(lowk)], ratio_n_k_0p05_to_0p3=[min(midk), max(midk)], modes=rows)
        cmp.append(cc)
        print('e=%.1e %-4s N(k<=1): theirs %.5f frozen %.5f ratio %.4f | k<=0.3: theirs %.5f frozen %.5f ratio %.4f | n_k ratio 2e-3..0.05 [%.4f, %.4f], 0.05..0.3 [%.3f, %.3f]' % (
            e, var, c['N_A_k_le_1'], N_frozen_kle1, cc['ratio_k_le_1'], N_t03, N_f03, cc['ratio_k_le_0p3'], min(lowk), max(lowk), min(midk), max(midk)))
    out['comparison'] = cmp
    # K-scaling check: full vs inst at 1e-6 on modes k <= 0.05 (frozen region)
    out['seconds'] = time.time() - t0
    json.dump(out, open(RES / 'sk_frozen.json', 'w'), indent=1)
    print('seconds %.1f' % out['seconds'])


if __name__ == '__main__':
    main()
