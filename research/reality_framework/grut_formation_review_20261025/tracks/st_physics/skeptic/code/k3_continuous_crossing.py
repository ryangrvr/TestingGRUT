"""st_physics skeptic, check k3 (own code on Codex's published arrays).

(1) Continuous (between-sample) maxima of the fractional amplitude response at k = 2.13 near the published unit
    crossings.  From every published sample we propagate the X background and the X part of Codex's saved 12-vector
    over one sample interval with the X-only conformal quartic equations (the reviewer and we both find these
    reproduce the full tangent to < 0.2% over 600 tau; over one 0.5-tau interval the error is checked against the
    next sample).  The fine-grained maximum per pericentre passage is compared with the sampled maximum; the first
    pericentre whose continuous maximum exceeds 1 is the true first crossing at that k.
(2) Size of the neglected nonlinear terms when s * delta chi_X is inserted (cos(kx) = 1 point), over the whole radial
    period containing the published 0.03 crossing:
        N = lam [2 chi |d|^2 + conj(chi) d^2 + |d|^2 d]   (complex notation; exact remainder of the cubic force)
    compared with k^2 |d| and with the full linear force |k^2 d + lam(2|chi|^2 d + chi^2 conj(d))|.
    The reviewer quotes "<= 7%"; we split the quadratic (2 chi|d|^2 + conj(chi) d^2) and cubic (|d|^2 d) parts.
(3) First sample where s |delta chi_X| >= r_max (k = 2.13).

usage: python3 -I k3_continuous_crossing.py <codex timing packet dir> <out.json>
"""
import sys, json
import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid

LAM = 1e-3
pkt, out = sys.argv[1:3]
rr = {r['tag']: r for r in json.load(open(f'{pkt}/RESPONSE_RESULTS.json'))['rows']}


def rhs(e, y, k):
    c = y[0] + 1j * y[1]; cp = y[2] + 1j * y[3]; d = y[4] + 1j * y[5]; dp = y[6] + 1j * y[7]
    cpp = -LAM * abs(c) ** 2 * c
    dpp = -k * k * d - LAM * (2 * abs(c) ** 2 * d + c * c * np.conj(d))
    return [cp.real, cp.imag, cpp.real, cpp.imag, dp.real, dp.imag, dpp.real, dpp.imag]


res = {}
for tag in ('quiet3', 'quiet5'):
    d = np.load(f'{pkt}/{tag}.npz')
    t, bg, sv, svk = d['times'], d['bg'], d['selected_vectors'], d['selected_vector_ks']
    j = int(np.argmin(abs(svk - 2.13))); k = float(svk[j])
    a = np.exp(bg[:, 15]); H = bg[:, 16]
    x = bg[:, 2] + 1j * bg[:, 3]; xd = bg[:, 8] + 1j * bg[:, 9]
    chi = np.sqrt(2) * a * x; chip = np.sqrt(2) * a * a * (xd + H * x)
    eta = np.r_[0, cumulative_trapezoid(1 / a, t)]
    dq = sv[:, j, 0] + 1j * sv[:, j, 1]; dp = sv[:, j, 6] + 1j * sv[:, j, 7]
    A_s = np.real(dq * np.conj(chi)) / abs(chi) ** 2           # sampled amplitude response per radian
    pub = d['amplitude_per_radian'][:, int(np.argmin(abs(d['ks'] - 2.13)))]
    out_tag = dict(k=k, sampled_vs_published_amp_maxreldiff=float(np.max(abs(A_s[t >= 80] - pub[t >= 80]) / np.max(abs(pub[t >= 80])))))
    cr = {c['phase_seed_radians']: c for c in rr[tag]['crossings']}
    rmax_late = float(np.max(abs(chi[t >= 200])))
    for s, lo, hi in ((0.03, 80, 160), (1e-6, 680, 800)):
        idx = np.flatnonzero((t >= lo) & (t < hi))
        fine_t, fine_A, endchk = [], [], []
        for i in idx:
            if i + 1 >= len(t):
                break
            y0 = [chi[i].real, chi[i].imag, chip[i].real, chip[i].imag, dq[i].real, dq[i].imag, dp[i].real, dp[i].imag]
            ee = np.linspace(eta[i], eta[i + 1], 41)
            so = solve_ivp(rhs, (eta[i], eta[i + 1]), y0, args=(k,), method='DOP853', rtol=1e-12, atol=1e-12, t_eval=ee)
            c = so.y[0] + 1j * so.y[1]; dd = so.y[4] + 1j * so.y[5]
            A = np.real(dd * np.conj(c)) / abs(c) ** 2
            fine_t.append(np.interp(ee, eta, t)); fine_A.append(A)
            endchk.append(abs(dd[-1] - dq[i + 1]) / abs(dq[i + 1]))
        fine_t = np.concatenate(fine_t); fine_A = np.concatenate(fine_A)
        # group by pericentre: local minima of |chi| in the sampled series
        r = abs(chi)
        per = [i for i in idx[1:-1] if r[i] < r[i - 1] and r[i] <= r[i + 1]]
        passes = []
        for i in per:
            w = (t >= t[i] - 3) & (t <= t[i] + 3)
            wf = (fine_t >= t[i] - 3) & (fine_t <= t[i] + 3)
            passes.append(dict(tau=float(t[i]), sampled_max=float(s * np.max(abs(A_s[w]))),
                               continuous_max=float(s * np.max(abs(fine_A[wf])))))
        first_sampled = next((p for p in passes if p['sampled_max'] >= 1), None)
        first_cont = next((p for p in passes if p['continuous_max'] >= 1), None)
        out_tag[f'seed_{s}'] = dict(window=[lo, hi], max_interval_end_mismatch=float(np.max(endchk)),
                                    pericentre_passes=passes, first_sampled=first_sampled, first_continuous=first_cont,
                                    published=cr[s]['amplitude_unity'])
        print(tag, s, 'end mismatch %.2e' % np.max(endchk), 'first sampled', first_sampled, 'first continuous', first_cont)
        for p in passes:
            print('   peri tau %.1f sampled %.3f continuous %.3f' % (p['tau'], p['sampled_max'], p['continuous_max']))
    # (2) nonlinear-term sizes over the radial period containing the 0.03 crossing (sampled + fine)
    s = 0.03
    tc = cr[0.03]['amplitude_unity']['tau_bracket'][1]
    per_tau = 18.0 if tag == 'quiet3' else 18.5
    w = (t >= tc - per_tau) & (t <= tc + 1)
    dd = s * dq[w]; cc = chi[w]
    lin = abs(k * k * dd + LAM * (2 * abs(cc) ** 2 * dd + cc * cc * np.conj(dd)))
    quad = LAM * abs(2 * cc * abs(dd) ** 2 + np.conj(cc) * dd ** 2)
    cub = LAM * abs(dd) ** 3
    i0 = int(np.argmin(abs(t[w] - tc)))
    out_tag['nonlinear_terms_0.03'] = dict(
        tau_window=[float(tc - per_tau), float(tc + 1)],
        at_crossing_sample=dict(tau=float(t[w][i0]), abs_d_over_r=float(abs(dd[i0]) / abs(cc[i0])),
                                cubic_over_k2d=float(cub[i0] / (k * k * abs(dd[i0]))),
                                quadratic_over_k2d=float(quad[i0] / (k * k * abs(dd[i0]))),
                                total_over_linear=float(abs(LAM * (2 * cc[i0] * abs(dd[i0]) ** 2 + np.conj(cc[i0]) * dd[i0] ** 2 + abs(dd[i0]) ** 2 * dd[i0])) / lin[i0])),
        max_over_period=dict(cubic_over_k2d=float(np.max(cub / (k * k * abs(dd)))),
                             quadratic_over_k2d=float(np.max(quad / (k * k * abs(dd)))),
                             total_over_linear=float(np.max(abs(LAM * (2 * cc * abs(dd) ** 2 + np.conj(cc) * dd ** 2 + abs(dd) ** 2 * dd)) / lin)),
                             abs_d_over_r=float(np.max(abs(dd) / abs(cc)))),
        rms_over_period_total_over_linear=float(np.sqrt(np.mean((abs(LAM * (2 * cc * abs(dd) ** 2 + np.conj(cc) * dd ** 2 + abs(dd) ** 2 * dd)) / lin) ** 2))))
    print(tag, 'nonlinear', json.dumps(out_tag['nonlinear_terms_0.03']))
    # (3) |delta chi_X| >= r_max
    for s in (0.03, 1e-6):
        wh = np.flatnonzero(s * abs(dq) >= rmax_late)
        out_tag[f'dchi_ge_rmax_{s}'] = None if not len(wh) else dict(tau=float(t[wh[0]]), a=float(a[wh[0]]))
    out_tag['rmax_late'] = rmax_late
    print(tag, 'dchi>=rmax', out_tag['dchi_ge_rmax_0.03'], out_tag['sampled_vs_published_amp_maxreldiff'])
    res[tag] = out_tag
json.dump(res, open(out, 'w'), indent=1)
