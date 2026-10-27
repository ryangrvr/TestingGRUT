"""en_particles: collect headline numbers, dispersion fits and verdicts into results/summary.json (own code).
Usage: python3 -I summarize_particles.py <scratch_dir> <track_results_dir>
"""
import sys, json, shutil
from pathlib import Path
import numpy as np


def disp_fit(shells, key, kmin=3.5, kmax=10.):
    """Fit w^2 = k_fd4^2 + m^2 on the branch peak (F) for kmin<=k<=kmax; returns m^2 mean and spread."""
    v = []
    for r in shells:
        k = r['k_fd4_rms']; w = r[key].get('omega')
        if w and kmin <= k <= kmax: v.append(w*w-k*k)
    v = np.array(v); return dict(m2_mean=float(v.mean()), m2_min=float(v.min()), m2_max=float(v.max()), shells=len(v))


def main(scr, res):
    scr = Path(scr); res = Path(res)
    files = {'spectral_W_early.json': scr/'W_early/spectral.json', 'spectral_W_late.json': scr/'W_late/spectral.json',
             'spectral_W_late_h2.json': scr/'W_late/spectral_h2.json', 'spectral_W_ext.json': scr/'W_ext/spectral.json',
             'response_eps005.json': scr/'resp_005.json', 'response_eps02.json': scr/'resp_02.json',
             'subwindow_seed2.json': scr/'ctrl/sub_seed2.json', 'subwindow_box16_N64.json': scr/'ctrl/sub_box16.json'}
    for k, v in files.items(): shutil.copy(v, res/k)
    S = {k: json.load(open(res/k)) for k in files}
    sub = json.load(open(res/'subwindow_totals.json'))['rows']+json.load(open(res/'subwindow_totals_first.json'))['rows']
    sub += S['subwindow_seed2.json']['rows']+S['subwindow_box16_N64.json']['rows']
    tot = [r['n_minus_all'] for r in sub]
    out = dict(source='own code (analysis) on reproduction-type regenerated histories (Codex engine unmodified)')
    out['windows'] = {}
    for name in ['spectral_W_early.json', 'spectral_W_late.json', 'spectral_W_late_h2.json', 'spectral_W_ext.json']:
        d = S[name]; w = d['window']; t = d['totals']; e = d['equal_time']; z = d['zero_mode']
        out['windows'][name] = dict(
            s=[w['s_start'], w['s_end']],
            tau=[20+8.869218738806*w['s_start']+1.9857962750634526*w['s_start']**2/2,
                 20+8.869218738806*w['s_end']+1.9857962750634526*w['s_end']**2/2],
            resolution=w['frequency_resolution'],
            spectral_nminus_nonzero=t['spectral_n_minus_per_net'], spectral_nplus_nonzero=t['spectral_n_plus_per_net'],
            branchfreq_nminus_nonzero=t['branchfreq_n_minus_per_net'],
            charge_check=[t['spectral_charge_nonzero_per_net'], t['equal_time_charge_nonzero_per_net']],
            zero_mode=dict(omega0=z['omega0_lab'], nminus=z['n_minus'], nplus=z['n_plus'], r_min=z['radius_min'], r_max=z['radius_max']),
            nminus_all_modes=t['spectral_n_minus_per_net']+z['n_minus'],
            equal_time_massless=e['massless'], equal_time_hartree=e['hartree'], mH2=e['mH2'],
            phi_quanta_spectral=t['phi_quanta_spectral_per_net'],
            dispersion_pos=disp_fit(d['shells'], 'pos'), dispersion_neg=disp_fit(d['shells'], 'neg'),
            lab_symmetry_low_k=[(r['k_fd4_rms'], r['pos'].get('omega'), r['neg'].get('omega')) for r in d['shells'][:6]])
    rp = S['response_eps005.json']
    out['response'] = dict(s0=rp['s0'], T=rp['T'], eps=[rp['eps'], S['response_eps02.json']['eps']],
                           Z_plus_range=[min(r['Z_plus'] for r in rp['shells'][:22]), max(r['Z_plus'] for r in rp['shells'][:22])],
                           Z_minus_range=[min(r['Z_minus'] for r in rp['shells'][:22]), max(r['Z_minus'] for r in rp['shells'][:22])],
                           Z_sum_range=[min(r['Z_sum'] for r in rp['shells'][:22]), max(r['Z_sum'] for r in rp['shells'][:22])],
                           n_count_minus=rp['n_count_minus_total'], nbar_minus_FDT=rp['nbar_minus_total'],
                           nbar_minus_FDT_eps02=S['response_eps02.json']['nbar_minus_total'],
                           rho_m2_k_ge_4=[r['rho_peak_pos']**2-(r['k_cont'])**2 for r in rp['shells'] if 4 <= r['k_cont'] <= 10],
                           lab_symmetry_rho=[(r['k_cont'], r['rho_peak_pos'], r['rho_peak_neg']) for r in rp['shells'][:8]])
    out['total_nminus_all_modes'] = dict(min=min(tot), max=max(tot), subwindows=len(tot),
                                         tau_range=[min(r['tau_start'] for r in sub), max(r['tau_end'] for r in sub)])
    out['controls_tau_550_800'] = {r['window']: dict(nminus_nonzero=r['n_minus_nonzero'], nminus_zero=r['n_minus_zero'],
                                                     eq_hartree=r['eq_hartree'], eq_massless=r['eq_massless'])
                                   for r in sub if r['window'] in ('W_late_seed2', 'W_late_box16') and r['tau_end'] > 790}
    json.dump(out, open(res/'summary.json', 'w'), indent=1)
    print(json.dumps(out, indent=1)[:6000])


if __name__ == '__main__':
    main(*sys.argv[1:])
