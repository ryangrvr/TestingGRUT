#!/usr/bin/env python3
"""st_clocks_skeptic (Claude, own code): adversarial checks of the time-ordering statement (st_clocks item T1/N1/Q1).

1. Condensate-assisted Y_L -> XX [E]: re-integrate Gamma_2 = c2 |X|^2 (c2 = 9 kappa^2/(8 pi m_L), re-derived in
   widths_skeptic.py by a real-field vertex sum) on TWO data sets: Codex's own background (quiet*.npz, physical complex
   canonical X in columns 2,3) and st_clocks' rp run (absX column), with several lower limits.
2. 1e-9 rad seed: Codex's saved per-radian amplitude response (quiet*.npz 'amplitude_per_radian', DATA) gives the
   envelope max_k |dX/X per rad|(tau); we fit the late growth of its running maximum against a (conformal time ~ a/C in
   radiation) and extrapolate where 1e-9 x envelope = 1 [E, extrapolation of a linear solution, not a computation].
3. Wavelength cube: Q = |N_X| (2 pi/k)^3 at k = 2.13 and at the k that actually carries each unit crossing.
Output: ../results/ordering_skeptic.json
"""
import json, math
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / 'results' / 'ordering_skeptic.json'
SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad')
CODEX = SP / 'b73/research/reality_framework/grut_source_timing_20261026'
RP = SP / 'review/st_clocks/rp/runs'
W = json.loads((HERE.parent / 'results' / 'widths_skeptic.json').read_text())['scalar']['theta_minus_pi_over_2']
BG = json.loads((HERE.parent / 'results' / 'background_skeptic.json').read_text())['histories']
C2 = W['cond_assisted_coeff']; G3 = W['gamma_L']


def trap(x, y):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))


def cond_integral(t, X2, lo, hi):
    s = (t >= lo) & (t <= hi)
    return C2 * trap(t[s], X2[s])


def main():
    cr = {r['tag']: r for r in json.loads((CODEX / 'RESPONSE_RESULTS.json').read_text())['rows']}
    res = dict(c2=C2, gamma3=G3, histories={})
    for tag, f in [('quiet3', 'q3_nominal_800.npy'), ('quiet5', 'q5_nominal_800.npy')]:
        z = np.load(CODEX / (tag + '.npz'))
        t = z['times']; bg = z['bg']; a = np.exp(bg[:, 15])
        X2c = bg[:, 2] ** 2 + bg[:, 3] ** 2
        S = np.load(RP / f); tr = S[:, 0]; X2r = S[:, 13] ** 2
        unit = cr[tag]['crossings'][0]['amplitude_unity']['tau_bracket'][1]
        h = dict(unit_sample=unit, condensate_assisted={})
        for lo in (0.0, 1.0, 5.0):
            for hi in (unit, 800.0):
                key = 'from_%g_to_%g' % (lo, hi)
                h['condensate_assisted'][key] = dict(codex_data=cond_integral(t, X2c, lo, hi),
                                                     rp_data=cond_integral(tr, X2r, lo, hi),
                                                     three_body=G3 * (hi - lo))
        # rate ratio near the unit sample (window mean)
        s = (tr > unit - 10) & (tr < unit + 10)
        h['Gamma2_over_gamma3_near_unit'] = C2 * float(np.mean(X2r[s])) / G3
        # where does the [E] channel accumulate?  fraction of the 0->unit integral in tau < 5
        h['share_of_integral_before_tau5'] = cond_integral(tr, X2r, 0, 5) / cond_integral(tr, X2r, 0, unit)
        # ---- 1e-9 seed extrapolation from Codex's per-radian envelope
        amp = np.abs(z['amplitude_per_radian']); env = amp.max(axis=1)
        run = np.maximum.accumulate(env)
        late = t >= 500
        # fit ln(running max) vs a on tau 500-800
        p = np.polyfit(a[late], np.log(run[late]), 1)
        # also fit on the last 150 tau only
        last = t >= 650
        p2 = np.polyfit(a[last], np.log(run[last]), 1)
        need = math.log(1e9)
        ext = {}
        for name, pp in [('fit_500_800', p), ('fit_650_800', p2)]:
            a_cross = (need - pp[1]) / pp[0]
            C = BG[tag]['C_own']; ts = BG[tag]['tau_star_own']
            ext[name] = dict(slope_per_unit_a=float(pp[0]), a_cross=float(a_cross),
                             tau_cross=float(a_cross ** 2 / (2 * C) + ts))
        # growth rate per unit a in earlier windows, for the trend
        rates = {}
        for lo, hi in [(60, 140), (140, 300), (300, 500), (500, 800)]:
            s = (t >= lo) & (t <= hi)
            rates['%d_%d' % (lo, hi)] = float(np.polyfit(a[s], np.log(run[s]), 1)[0])
        h['seed_1e9'] = dict(envelope_at_800=float(run[-1]), response_1e9_at_800=float(1e-9 * run[-1]),
                             response_1e6_at_800=float(1e-6 * run[-1]), extrapolation=ext,
                             growth_per_unit_a=rates,
                             YL_1percent_tau=-math.log1p(-0.01) / G3)
        # decaying-rate extrapolation: fit ln(rate) ~ ln(a) and integrate (more pessimistic)
        mids, rr = [], []
        for lo, hi in [(140, 300), (300, 500), (500, 800)]:
            s = (t >= lo) & (t <= hi)
            mids.append(float(np.mean(a[s]))); rr.append(float(np.polyfit(a[s], np.log(run[s]), 1)[0]))
        q = np.polyfit(np.log(mids), np.log(rr), 1)      # rate = exp(q1) a^q0
        a0 = a[-1]; lnR0 = math.log(run[-1]); aa = a0
        while lnR0 < need and aa < 1e6:
            aa_new = aa * 1.001
            lnR0 += math.exp(q[1]) * (0.5 * (aa + aa_new)) ** q[0] * (aa_new - aa)
            aa = aa_new
        h['seed_1e9']['power_law_rate_extrapolation'] = dict(rate_exponent=float(q[0]), a_cross=aa,
                                                             tau_cross=aa ** 2 / (2 * BG[tag]['C_own']) + BG[tag]['tau_star_own'])
        # ---- wavelength cube
        NX = cr[tag]['late_NX_normalization_through80']
        kk = [cr[tag]['crossings'][i]['amplitude_unity']['k'] for i in range(2)]
        h['Q_cube'] = {('k=%g' % k): abs(NX) * (2 * math.pi / k) ** 3 for k in sorted(set([2.13] + kk))}
        h['unit_crossing_k'] = kk
        h['tenpercent_crossing_k_0.03'] = cr[tag]['crossings'][0]['amplitude_10percent']['k']
        res['histories'][tag] = h
        print(tag, json.dumps(h, indent=1), flush=True)
    OUT.write_text(json.dumps(res, indent=2) + '\n')


if __name__ == '__main__':
    main()
