"""pt_relax skeptic (Claude, own code): robustness of pt_relax's N32 late-time fit A_-(s) -> A_inf.

Reads pt_relax results/long_run.json (windows of the N32 run to s = 1004, which this skeptic reproduced bit-for-bit
to s = 254 with scan.py) and refits A(s) = A_inf + dA exp(-G (s - s0)) from several start points, and with A_inf
fixed at the self-consistent Rayleigh-Jeans lattice end point (equilibria.py).  Also tabulates the IR occupations
(k <= 8 shells) that set whether the classical number-changing rate is trustworthy.  Writes results/refit_n32.json.
"""
import json, os
import numpy as np
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
PTR = os.path.join(HERE, '..', '..', 'results', 'long_run.json')
EQ = os.path.join(HERE, '..', 'results', 'equilibria.json')
OUT = os.path.join(HERE, '..', 'results', 'refit_n32.json')


def model(s, Ainf, dA, G):
    return Ainf + dA * np.exp(-G * (s - 299.0))


def main():
    d = json.load(open(PTR))
    res = {}
    W = d['IR12_s1_w1']['windows']
    s = np.array([np.mean(w['s']) for w in W]); A = np.array([w['total_all_modes'] for w in W])
    fits = []
    for s0 in [29, 199, 299, 399, 499, 599, 699]:
        m = s >= s0
        try:
            p, C = curve_fit(model, s[m], A[m], p0=[0.5, 0.9, 1.5e-3], maxfev=20000)
            r = A[m] - model(s[m], *p)
            fits.append(dict(s_start=s0, n=int(m.sum()), A_inf=p[0], A_inf_err=float(np.sqrt(C[0, 0])), Gamma=p[2],
                             Gamma_err=float(np.sqrt(C[2, 2])), rms=float(np.sqrt(np.mean(r * r)))))
        except Exception as e:
            fits.append(dict(s_start=s0, error=str(e)))
    res['free_fits'] = fits
    rj = {x['N']: x['Xbar_per_net'] for x in json.load(open(EQ))['rj_lattice']}
    fixed = []
    for Ainf in [rj[32], 0.4078, 0.51]:
        for s0 in [299, 499]:
            m = s >= s0
            f = lambda ss, dA, G: model(ss, Ainf, dA, G)
            p, C = curve_fit(f, s[m], A[m], p0=[0.9, 1.5e-3], maxfev=20000)
            r = A[m] - f(s[m], *p)
            fixed.append(dict(A_inf_fixed=Ainf, s_start=s0, Gamma=p[1], rms=float(np.sqrt(np.mean(r * r))),
                              max_abs_resid=float(np.max(np.abs(r)))))
    res['fixed_Ainf_fits'] = fixed
    # local log-slope of (A - A_RJ): if the approach to RJ were a single exponential this would be constant
    loc = []
    for a, b in [(299, 399), (399, 499), (499, 599), (599, 699), (699, 799), (799, 899), (899, 999)]:
        m = (s >= a) & (s <= b)
        loc.append(dict(range=[a, b], dlnA_ds=float(np.polyfit(s[m], np.log(A[m]), 1)[0]),
                        dln_A_minus_RJ_ds=float(np.polyfit(s[m], np.log(A[m] - rj[32]), 1)[0])))
    res['local_slopes_N32'] = loc
    # IR occupations, both pt_relax runs
    occ = {}
    for tag in ['IR12_s1_w1', 'IR12_fine_s1_w1']:
        O = d[tag]['occupations']; ks = d[tag].get('k_shells')
        for key, o in O.items():
            ox = np.array(o['occX_by_shell']); op = np.array(o['occphi_by_shell'])
            dk = 2 * np.pi / 12.0
            kk = np.arange(len(ox)) * dk
            band = (kk > 0) & (kk <= 8 + 1e-9)
            occ[f'{tag}:{key}'] = dict(X_band_min=float(ox[band].min()), X_band_median=float(np.median(ox[band])),
                                       X_band_max=float(ox[band].max()), phi_band_min=float(op[band].min()),
                                       phi_band_median=float(np.median(op[band])),
                                       X_all_min=float(ox[1:].min()), phi_all_min=float(op[1:].min()))
    res['occupations_k_le_8'] = occ
    json.dump(res, open(OUT, 'w'), indent=1, default=float)
    for f in fits: print('FREE', {k: (round(v, 5) if isinstance(v, float) else v) for k, v in f.items()})
    for f in fixed: print('FIXED', {k: (round(v, 6) if isinstance(v, float) else v) for k, v in f.items()})
    for l in loc: print('LOCAL', l)
    for k, v in occ.items(): print('OCC', k, {a: round(b, 3) for a, b in v.items()})


if __name__ == '__main__':
    main()
