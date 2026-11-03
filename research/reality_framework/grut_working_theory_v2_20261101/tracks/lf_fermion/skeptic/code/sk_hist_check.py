"""Skeptic S5: energy-conservation and continuity checks of results/hist_HFa_LF.json (read-only).
Claude, lf_fermion skeptic (issue #13 context). Checks: fractions sum to 1; f_rad = sum f_N + f_SM; integrated
Friedmann/continuity consistency ln H(t)/H(t0) = -int eps_H H dt and ln a = int H dt; window means of the numerical
-dH/dt/H^2 against the stored eps_H and 1.5(1+w); largest sample-to-sample jump of f_rad; lattice -> fluid splice.
Usage: python3 -B sk_hist_check.py -> ../results/sk_hist_check.json"""
import json, sys
from pathlib import Path
import numpy as np
TR = Path(__file__).resolve().parents[2]
d = json.load(open(TR / 'results' / 'hist_HFa_LF.json'))
trap = getattr(np, 'trapezoid', None) or np.trapz
out = {}


def chk(name, L, fl=None):
    t = np.array(L['t']); a = np.array(L['a']); H = np.array(L['H']); w = np.array(L['w']); eps = np.array(L['eps_H'])
    fs = sum(np.array(L[k]) for k in ['f_sig', 'f_chi', 'f_grad', 'f_s', 'f_rad'])
    fN = np.array(L['f_N']).sum(0) + np.array(L['f_SM'])
    epsn = -np.gradient(H, t) / H ** 2
    def wavg(y, lo, hi):
        s = (t >= lo) & (t < hi); return float(trap(y[s], t[s]) / (t[s][-1] - t[s][0]))
    I = np.concatenate([[0], np.cumsum(0.5 * (eps[1:] * H[1:] + eps[:-1] * H[:-1]) * np.diff(t))])
    Ia = np.concatenate([[0], np.cumsum(0.5 * (H[1:] + H[:-1]) * np.diff(t))])
    r = dict(sum_f_minus_1=float(np.abs(fs - 1).max()), frad_minus_sumN=float(np.abs(np.array(L['f_rad']) - fN).max()),
             lnH_dev=float(np.abs(np.log(H / H[0]) + I).max()), lna_dev=float(np.abs(np.log(a / a[0]) - Ia).max()),
             windows=[dict(lo=lo, hi=hi, eps_num=wavg(epsn, lo, hi), eps_H=wavg(eps, lo, hi), eps_w=wavg(1.5 * (1 + w), lo, hi))
                      for lo, hi in [(0.1, 0.3), (0.3, 0.5), (0.5, 1.0), (1.0, 2.0), (2.0, 2.99)]])
    fr = np.array(L['f_rad']); j = int(np.argmax(np.abs(np.diff(fr))))
    r['max_dfrad'] = float(np.abs(np.diff(fr)).max()); r['t_max_dfrad'] = float(t[j])
    if fl is not None:
        tf = np.array(fl['t']); Hf = np.array(fl['H']); ef = np.array(fl['eps_H'])
        If = np.concatenate([[0], np.cumsum(0.5 * (ef[1:] * Hf[1:] + ef[:-1] * Hf[:-1]) * np.diff(tf))])
        r['splice'] = dict(lat_end=dict(t=float(t[-1]), a=float(a[-1]), H=float(H[-1]), f_sig=L['f_sig'][-1],
                                        f_chi_tot=L['f_chi'][-1] + L['f_grad'][-1] + L['f_s'][-1], f_rad=L['f_rad'][-1]),
                           fluid_start=dict(t=fl['t'][0], a=fl['a'][0], H=fl['H'][0], f_sig=fl['f_sig'][0], f_chi=fl['f_chi'][0], f_rad=fl['f_rad'][0]),
                           fluid_lnH_dev=float(np.abs(np.log(Hf / Hf[0]) + If).max()),
                           fluid_sum_f_minus_1=float(np.abs(np.array(fl['f_sig']) + np.array(fl['f_chi']) + np.array(fl['f_rad']) - 1).max()))
    out[name] = r
    print(name, json.dumps(r, indent=None)[:1500])


chk('KFmid', d['lattice'], d['fluid'])
for k, v in d['alternates'].items():
    chk(k, v['lattice'])
(TR / 'skeptic' / 'results' / 'sk_hist_check.json').write_text(json.dumps(dict(label='C', note=__doc__, runs=out), indent=1))
