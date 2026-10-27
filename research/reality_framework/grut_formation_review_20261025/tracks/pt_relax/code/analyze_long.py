"""Track pt_relax (Claude, own code): analysis of the long IR12 (N32) and IR12_fine (N48) continuations.

Reads the part*.npz files written by long_run.py from $RUNDIR/<tag>/ and writes results/long_run.json:
window series, drift fits, cascade diagnostics, occupations and the classical-reliability statement.
"""
import os, sys, json, glob
import numpy as np
from scipy.optimize import curve_fit

RUN = os.environ['RUNDIR']
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
L = 12.0
dk = 2 * np.pi / L


def load(tag):
    rows, specs, wins = [], [], []
    for f in sorted(glob.glob(os.path.join(RUN, tag, 'part*.npz'))):
        d = np.load(f)
        r = json.loads(str(d['rows'])); s = json.loads(str(d['specs'])); w = json.loads(str(d['wins']))
        # each restart repeats the boundary sample: drop duplicates by s
        rows += [x for x in r if not rows or x['s'] > rows[-1]['s'] + 1e-9]
        specs += [x for x in s if not specs or x['s'] > specs[-1]['s'] + 1e-9]
        wins += w
    return rows, specs, wins


def fit_exp(sc, A):
    f = lambda s, Ainf, dA, G: Ainf + dA * np.exp(-G * (s - sc[0]))
    p, cov = curve_fit(f, sc, A, p0=[0.5, A[0] - 0.5, 1e-3], maxfev=20000)
    return dict(A_inf=p[0], dA=p[1], Gamma=p[2], Gamma_err=float(np.sqrt(cov[2, 2])),
                rms_resid=float(np.sqrt(np.mean((f(sc, *p) - A) ** 2))))


def analyse(tag):
    rows, specs, wins = load(tag)
    W = [w for w in wins if w['window'][0] >= 24 - 1e-9]
    sc = np.array([np.mean(w['window']) for w in W])
    tot_all = np.array([w['total_minus_all'] for w in W]); tot8 = np.array([w['total_minus_band'] for w in W])
    coh = np.array([w['coherent_minus'] for w in W]); band = np.array([w['band_minus'] for w in W])
    uv = np.array([w['uv_minus'] for w in W]); chk = np.array([w['charge_check_all'] for w in W])
    res = dict(tag=tag, s_end=float(rows[-1]['s']),
               windows=[dict(s=w['window'], coherent=w['coherent_minus'], band_k_le_8=w['band_minus'],
                             k_gt_8=w['uv_minus'], total_codex_band=w['total_minus_band'], total_all_modes=w['total_minus_all'],
                             spectral_charge=w['charge_check_all'], mean_abs_w_minus_band=w['mean_abs_w_minus_band'])
                        for w in wins],
               max_energy_rel_drift=float(max(abs(r['rho'] / rows[0]['rho'] - 1) for r in rows)),
               max_charge_rel_drift=float(max(abs(r['Q'] / rows[0]['Q'] - 1) for r in rows)),
               spectral_charge_range=[float(chk.min()), float(chk.max())])
    # drift diagnostics
    def slope(lo, hi, y):
        m = (sc >= lo) & (sc <= hi)
        if m.sum() < 3: return None
        p = np.polyfit(sc[m], y[m], 1); return float(p[0])
    segs = [(29, 99), (99, 199), (199, 299), (299, 499), (499, 699), (699, 1000)]
    res['slopes_total_all_per_unit_s'] = {f'{a}-{b}': slope(a, b, tot_all) for a, b in segs}
    res['slopes_total_band_per_unit_s'] = {f'{a}-{b}': slope(a, b, tot8) for a, b in segs}
    res['fractional_change_total_all'] = dict(first=float(tot_all[0]), last=float(tot_all[-1]),
                                              first_s=W[0]['window'], last_s=W[-1]['window'],
                                              change=float(tot_all[-1] / tot_all[0] - 1))
    if sc[-1] > 300:
        res['exp_fit_total_all_s_ge_29'] = fit_exp(sc, tot_all)
        m = sc >= 299
        res['exp_fit_total_all_s_ge_299'] = fit_exp(sc[m], tot_all[m])
    # equal-time series (decimated)
    keys = ['rho', 'Q', 'r2', 'mH2', 'mphi2', 'outside_mean', 'phi_total_frac', 'phi_quad_frac', 'E_uv_frac',
            'ref_minus_all', 'ref_minus_band', 'ref_minus_uv', 'Q0_over_Q', 'cond_r']
    step = max(1, len(rows) // 400)
    res['equal_time'] = {k: [float(r[k]) for r in rows[::step]] for k in keys}
    res['equal_time']['s'] = [float(r['s']) for r in rows[::step]]
    # window means of equal-time quantities used elsewhere
    def wmean(k, lo, hi):
        v = [r[k] for r in rows if lo - 1e-9 <= r['s'] <= hi + 1e-9]
        return float(np.mean(v)) if v else None
    res['window_means'] = {f'{lo}-{hi}': {k: wmean(k, lo, hi) for k in ['r2', 'mH2', 'mphi2', 'phi_total_frac', 'E_uv_frac', 'outside_mean']}
                           for lo, hi in [(14, 24), (24, 34), (194, 204), (494, 504), (994, 1004), (244, 254)]}
    # spectra: cascade and occupations
    S = np.array([x['s'] for x in specs])
    EX = np.array([x['EX'] for x in specs]); EP = np.array([x['Ephi'] for x in specs])
    occX = np.array([x['occX'] for x in specs]); occP = np.array([x['occphi'] for x in specs])
    nsh = EX.shape[1]; kk = np.arange(nsh) * dk
    Etot = EX + EP
    kmean = (Etot * kk).sum(1) / Etot.sum(1)
    k2 = np.sqrt((Etot * kk ** 2).sum(1) / Etot.sum(1))
    # front: largest shell whose X energy exceeds 1e-3 of the X spectrum maximum
    front = np.array([kk[np.nonzero(e > 1e-3 * e.max())[0].max()] for e in EX])
    res['cascade'] = dict(s=S[::5].tolist(), k_mean_energy=kmean[::5].tolist(), k_rms_energy=k2[::5].tolist(),
                          k_front_1em3=front[::5].tolist())
    for lo, hi in [(10, 100), (34, 200), (100, 300)]:
        m = (S >= lo) & (S <= hi)
        if m.sum() > 5 and S[m].max() >= hi - 1:
            p = np.polyfit(np.log(S[m]), np.log(k2[m]), 1)
            res['cascade'][f'power_k_rms_vs_s_{lo}-{hi}'] = float(p[0])
    # occupations at selected times
    occ = {}
    for t in [0, 34, 200, 500, 1000, 250]:
        i = int(np.argmin(abs(S - t)))
        if abs(S[i] - t) > 1: continue
        ox, op = occX[i], occP[i]
        good = np.arange(nsh) > 0
        frac_phi_low = float(EP[i][good & (op < 1)].sum() / (EX[i].sum() + EP[i].sum()))
        frac_X_low = float(EX[i][good & (ox < 1)].sum() / (EX[i].sum() + EP[i].sum()))
        occ[f's{S[i]:.0f}'] = dict(occX_min=float(ox[good].min()), occX_max=float(ox[good].max()),
                                   occX_at_k8=float(ox[int(round(8 / dk))]), occX_last_shell=float(ox[-1]),
                                   occphi_min=float(op[good].min()), occphi_max=float(op[good].max()),
                                   quad_energy_fraction_phi_in_shells_occ_lt_1=frac_phi_low,
                                   quad_energy_fraction_X_in_shells_occ_lt_1=frac_X_low,
                                   occX_by_shell=ox.tolist(), occphi_by_shell=op.tolist(),
                                   EX_by_shell=EX[i].tolist(), Ephi_by_shell=EP[i].tolist())
    res['occupations'] = occ
    res['k_shells'] = kk.tolist()
    return res


def main():
    out = {}
    for tag in sys.argv[1:] or ['IR12_s1_w1', 'IR12_fine_s1_w1']:
        if not os.path.isdir(os.path.join(RUN, tag)): continue
        r = analyse(tag); out[tag] = r
        print(tag, 's_end', r['s_end'], 'drift E', r['max_energy_rel_drift'], 'Q', r['max_charge_rel_drift'], 'charge', r['spectral_charge_range'])
        print('  slopes all', r['slopes_total_all_per_unit_s'])
        print('  change', r['fractional_change_total_all'])
        for k in ['exp_fit_total_all_s_ge_29', 'exp_fit_total_all_s_ge_299']:
            if k in r: print(' ', k, r[k])
        print('  cascade powers', {k: v for k, v in r['cascade'].items() if k.startswith('power')})
        for k, v in r['occupations'].items():
            print('  occ', k, {a: (round(b, 4) if isinstance(b, float) else None) for a, b in v.items() if not isinstance(b, list)})
        print('  window means', r['window_means'])
    json.dump(out, open(os.path.join(OUT, 'long_run.json'), 'w'), indent=None, separators=(',', ':'))


if __name__ == '__main__':
    main()
