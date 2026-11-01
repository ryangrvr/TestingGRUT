#!/usr/bin/env python3
"""fx_holder track (Claude, own code) -- items 3 and 4 of the repair brief: evolution through reheating with the
inflaton-trilinear holder (rptri.c), compared with the unheld M1 preparation ("original") and Option 1A, for the (3, 0.5)
and (5, 1) inflaton histories, at lambda_Y = lambda_Sigma = 1 and 0.5; window statistics of the late X charge (the
X <-> Sigma exchange drift), yields, the H_i recalibration with Y_Q ∝ H_i^{3/2} at fixed dimensionless history, and the
light-Y debris through our decay ledger.  Output: ../results/tri_runs.json (+ scratch .npy series).
"""
import json, sys, time
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import fx_common as fc

HIST = {'b3_g0.5': (3.0, 0.5), 'b5_g1': (5.0, 1.0)}


def summarize(rr, p):
    end = rr['end']
    hi = min(end, 1000.0)
    nx = fc.window_stats(rr, 'NX', 400.0, hi, 80.0)
    out = dict(tag=rr['tag'], lam=rr['lam'], b=rr['b'], gam=rr['gam'], end=end, rtol=rr['rtol'], G=rr['G'], lock=rr['lock'], I0=rr['I0'],
               NX_400_1000=nx['mean'], NX_run80_min=nx['run_min'], NX_run80_max=nx['run_max'], NX_run80_std=nx['run_std'],
               NX_640_800=fc.window_mean(rr, 'NX', 640.0, 800.0), NX_800_1000=fc.window_mean(rr, 'NX', 800.0, 1000.0),
               NX_W80=rr['W']['W80']['NX'], checks=rr['END'], seconds=rr['seconds'])
    if end >= 2000:
        out['NX_1000_2000'] = fc.window_mean(rr, 'NX', 1000.0, 2000.0)
        out['NX_1500_2000'] = fc.window_mean(rr, 'NX', 1500.0, 2000.0)
        s2 = fc.window_stats(rr, 'NX', 1000.0, 2000.0, 80.0)
        out['NX_run80_1000_2000'] = [s2['run_min'], s2['run_max'], s2['run_std']]
    out['running80'] = nx['running']
    # exchange partner: N_X + 6 N_Sigma and 3 N_Y over the same windows (Gauss: N_X + 3N_Y + 6N_S = 0)
    B = rr['bins']; t0 = B[:, 0]
    sel = (t0 >= 400 - 1e-9) & (B[:, 1] <= hi + 1e-9)
    nxb, nsb, nyb = fc.bin_col(rr, 'NX')[sel], fc.bin_col(rr, 'NS')[sel], fc.bin_col(rr, 'NY')[sel]
    out['exchange'] = dict(corr_NX_6NS=float(np.corrcoef(nxb, 6 * nsb)[0, 1]), corr_NX_3NY=float(np.corrcoef(nxb, 3 * nyb)[0, 1]),
                           std_NX=float(nxb.std()), std_NX_plus_6NS=float((nxb + 6 * nsb).std()), std_3NY=float((3 * nyb).std()))
    for k in ('NYL', 'NYH', 'Nsig', 'EX', 'eps'):
        out[k] = fc.window_mean(rr, k, 400.0, hi)
    Y, ratio, Rr = fc.yield_over_target(rr, out['NX_400_1000'])
    out['yield_400_1000'] = Y; out['yield_over_YQ_400_1000'] = ratio; out['Rr_end'] = Rr
    out['yield_over_YQ_640_800'] = fc.yield_over_target(rr, out['NX_640_800'])[1]
    out['yield_over_YQ_run80'] = [fc.yield_over_target(rr, out['NX_run80_min'])[1], fc.yield_over_target(rr, out['NX_run80_max'])[1]]
    # recalibrated H_i at fixed dimensionless history: Y ∝ H^{3/2}
    out['H_new_GeV_400_1000'] = fc.H_GEV * ratio ** (-2 / 3)
    out['H_new_GeV_640_800'] = fc.H_GEV * out['yield_over_YQ_640_800'] ** (-2 / 3)
    out['H_new_GeV_run80_range'] = sorted([fc.H_GEV * v ** (-2 / 3) for v in out['yield_over_YQ_run80']])
    out['ledger_e1e-7'] = fc.ledger(out['NYH'], out['NYL'], out['Nsig'], out['EX'], Rr, p, e=1e-7)
    out['ledger_e1e-6'] = fc.ledger(out['NYH'], out['NYL'], out['Nsig'], out['EX'], Rr, p, e=1e-6)
    # early build-up (fraction of the late charge) and the u(t) excursion
    S = rr['S']; C = fc.COLS
    out['charge_fraction_at_tau'] = {str(t): float(np.interp(t, S[:, 0], S[:, C.index('NX')]) / out['NX_400_1000']) for t in (0.4, 0.8, 1.2, 2, 3, 5, 10, 20, 50)}
    early = S[S[:, 0] <= 6.0]
    out['u_excursion_tau_le_6'] = dict(min_Im_u=float(early[:, C.index('uim')].min()), max_abs_u=float(np.hypot(early[:, C.index('ure')], early[:, C.index('uim')]).max()),
                                       min_mL2_forced=float(early[:, C.index('mL2')].min()))
    return out


def main():
    t0 = time.time()
    res = dict(runs={}, comparisons={})
    for lam in (1.0, 0.5):
        p = fc.params(lam)
        st = fc.stationary(fc.U_PREP, p)
        zh = st['z']
        zo, info = fc.prepare_held(p, 3.0)
        for hk, (b, gam) in HIST.items():
            G, I0, ph0 = fc.trilinear_G(zh, p, b)
            end = 2000.0 if hk == 'b3_g0.5' else 1000.0
            ns = int(end / 0.2)
            key = f'lam{lam:g}_{hk}'
            runs = {}
            runs['original'] = fc.run(key + '_original', p, b, gam, z0=zo, end=end, nsamp=ns)
            runs['option1A'] = fc.run(key + '_option1A', p, b, gam, z0=zh, lock=fc.LOCK_1A, end=end, nsamp=ns)
            runs['trilinear'] = fc.run(key + '_trilinear', p, b, gam, z0=zh, G=G, end=end, nsamp=ns)
            if lam == 1.0 and hk == 'b3_g0.5':
                runs['trilinear_rtol12'] = fc.run(key + '_trilinear_rtol12', p, b, gam, z0=zh, G=G, end=1000.0, nsamp=5000, rtol=1e-12, atol=1e-15)
                # CP conjugate: u(0) = -10 i (G -> conj(10i - mu)... i.e. G' with mu + G' phi0 = -10 i), conjugate preparation
                Gc, _, _ = fc.trilinear_G(np.conj(zh), p, b, u_target=-10j)
                runs['trilinear_conjugate'] = fc.run(key + '_trilinear_conj', p, b, gam, z0=np.conj(zh), G=Gc, end=1000.0, nsamp=5000)
                # opposite inflaton branch: same Lagrangian, phi_I(0) < 0  <=>  G -> -G with phi_I(0) > 0 (inflaton sector is Z2 even)
                uo = 2 * p['mu'] - 10j
                sto = fc.stationary(uo, p)
                runs['trilinear_opposite_branch'] = fc.run(key + '_trilinear_opp', p, b, gam, z0=sto['z'], G=-G, end=1000.0, nsamp=5000)
                # clock sensitivity: preparation slice at phi_I = (1 +- 1%) phi_I(0) (re-prepared conditional minimum)
                for sgn in (+1, -1):
                    s = 1 + 0.01 * sgn
                    us = p['mu'] + G * ph0 * s
                    sts = fc.stationary(us, p, guess=zh)
                    runs[f'trilinear_slice_{s:g}'] = fc.run(key + f'_trilinear_s{s:g}', p, b, gam, z0=sts['z'], G=G * s, end=1000.0, nsamp=5000)
                # same clock test for Option 1A (R slice 12 (1 +- 1%))
                for sgn in (+1, -1):
                    s = 1 + 0.01 * sgn
                    us = p['mu'] + fc.LOCK_1A * 12 * s
                    sts = fc.stationary(us, p, guess=zh)
                    runs[f'option1A_slice_{s:g}'] = fc.run(key + f'_option1A_s{s:g}', p, b, gam, z0=sts['z'], lock=fc.LOCK_1A * s, end=1000.0, nsamp=5000)
            summ = {k: summarize(v, p) for k, v in runs.items()}
            res['runs'][key] = summ
            o, a, t = summ['original'], summ['option1A'], summ['trilinear']
            comp = {}
            for w in ('NX_400_1000', 'NX_640_800', 'NX_800_1000', 'NX_W80'):
                comp[w] = dict(option1A_over_original=a[w] / o[w], trilinear_over_original=t[w] / o[w], trilinear_over_option1A=t[w] / a[w])
            if 'NX_1000_2000' in o:
                comp['NX_1000_2000'] = dict(option1A_over_original=a['NX_1000_2000'] / o['NX_1000_2000'], trilinear_over_original=t['NX_1000_2000'] / o['NX_1000_2000'],
                                            trilinear_over_option1A=t['NX_1000_2000'] / a['NX_1000_2000'])
            # running-window ratio spread (tri / original, 1A / original) over 80-wide windows
            ro = np.array([v for _, v in o['running80']]); ra = np.array([v for _, v in a['running80']]); rt = np.array([v for _, v in t['running80']])
            comp['running80_trilinear_over_original'] = [float((rt / ro).min()), float((rt / ro).max())]
            comp['running80_option1A_over_original'] = [float((ra / ro).min()), float((ra / ro).max())]
            if 'trilinear_slice_1.01' in summ:
                comp['dlnNX_dln_phiI_slice_trilinear'] = (np.log(abs(summ['trilinear_slice_1.01']['NX_400_1000'])) - np.log(abs(summ['trilinear_slice_0.99']['NX_400_1000']))) / (np.log(1.01) - np.log(0.99))
                comp['dlnNX_dlnR_slice_option1A'] = (np.log(abs(summ['option1A_slice_1.01']['NX_400_1000'])) - np.log(abs(summ['option1A_slice_0.99']['NX_400_1000']))) / (np.log(1.01) - np.log(0.99))
                comp['conjugate_over_trilinear'] = summ['trilinear_conjugate']['NX_400_1000'] / t['NX_400_1000']
                comp['opposite_branch_over_trilinear'] = summ['trilinear_opposite_branch']['NX_400_1000'] / t['NX_400_1000']
                comp['rtol12_rel_change'] = summ['trilinear_rtol12']['NX_400_1000'] / t['NX_400_1000'] - 1
            res['comparisons'][key] = comp
            print(key, json.dumps(comp, default=float))
            for k, v in summ.items():
                print('   %-28s NX(400-1000)=%11.3f run80=[%.3f, %.3f] NX(640-800)=%.3f Y/YQ=%.6f H_new=%.6e NYL=%.4f R_Y=%.4f checks=%s' % (
                    k, v['NX_400_1000'], v['NX_run80_min'], v['NX_run80_max'], v['NX_640_800'], v['yield_over_YQ_400_1000'], v['H_new_GeV_400_1000'],
                    v['NYL'], v['ledger_e1e-7']['R_Y'], v['checks']))
    # drop the long running-window lists from the saved json except for the main three at lambda = 1, (3, 0.5)
    for key, summ in res['runs'].items():
        for k, v in summ.items():
            if not (key == 'lam1_b3_g0.5' and k in ('original', 'option1A', 'trilinear')):
                v.pop('running80', None)
    res['seconds'] = time.time() - t0
    json.dump(res, open(fc.RES / 'tri_runs.json', 'w'), indent=1, default=float)
    print('total seconds', res['seconds'])


if __name__ == '__main__':
    main()
