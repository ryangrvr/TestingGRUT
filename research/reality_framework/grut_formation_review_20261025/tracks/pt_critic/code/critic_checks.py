#!/usr/bin/env python3
"""pt_critic: quick completeness checks on Codex 5f5513b and on the four review tracks.

INDEPENDENT code (own arithmetic). It imports no Codex module and no track module. Inputs are
Codex's published JSON/npz files (git archive of 5f5513b) and the tracks' stored result JSONs.

Usage: python3 -I critic_checks.py <codex packet dir> <early packet dir> <tracks dir> <out.json>
"""
import json, sys, math
from pathlib import Path
import numpy as np

PK = Path(sys.argv[1]); EARLY = Path(sys.argv[2]); TR = Path(sys.argv[3]); OUT = Path(sys.argv[4])
out = {}

# ---------------------------------------------------------------- A. packet arithmetic
R = json.load(open(PK / 'RESULTS.json'))
rows = {r['case']['tag']: r for r in R['rows']}
tags = list(rows)
A = {}
for tag, r in rows.items():
    A[tag] = {}
    for w in r['windows']:
        key = f"{int(w['s'][0])}_{int(w['s'][1])}"
        tail = w['omitted_Hartree_reference_above8']
        rm = w['reference_mass_means']
        pos = [rm[k] for k in ('7', '12', '20', '26')]
        A[tag][key] = dict(coh=w['coherent_negative_action'], wave=w['wave_negative_action'],
                           tot=w['total_negative_action'], tail=tail,
                           tot_incl=w['total_negative_action'] + tail,
                           charge_identity=abs(w['spectral_charge_ratio'] - 1.0),
                           refmass_spread_7_26=(max(pos) - min(pos)) / min(pos),
                           refmass_massless_over_H=rm['0'] / rm['H'])
def rng(key, field):
    v = [A[t][key][field] for t in tags]
    return [min(v), max(v)]
packet = {}
for key in ('14_24', '24_34'):
    packet[key] = dict(total_k8_range=rng(key, 'tot'), total_incl_tail_range=rng(key, 'tot_incl'),
                       tail_range=rng(key, 'tail'), charge_identity_max=max(A[t][key]['charge_identity'] for t in tags),
                       refmass_spread_7_26_range=rng(key, 'refmass_spread_7_26'),
                       massless_over_Hartree_range=rng(key, 'refmass_massless_over_H'),
                       per_case={t: A[t][key] for t in tags})
# declines window to window, per case
packet['decline_percent_k8'] = {t: 100 * (A[t]['24_34']['tot'] / A[t]['14_24']['tot'] - 1) for t in tags}
packet['decline_percent_incl_tail'] = {t: 100 * (A[t]['24_34']['tot_incl'] / A[t]['14_24']['tot_incl'] - 1) for t in tags}
# spatial refinement IR12_s1_w1 (N32) -> IR12_fine_s1_w1 (N48), both k<=8 and incl tail, and wave part
c, f = A['IR12_s1_w1'], A['IR12_fine_s1_w1']
packet['refinement_fine_minus_coarse_percent'] = {
    key: dict(total_k8=100 * (f[key]['tot'] / c[key]['tot'] - 1), total_incl_tail=100 * (f[key]['tot_incl'] / c[key]['tot_incl'] - 1),
              wave_k8=100 * (f[key]['wave'] / c[key]['wave'] - 1), tail_ratio_fine_over_coarse=f[key]['tail'] / c[key]['tail'])
    for key in ('14_24', '24_34')}
packet['codex_stated_refinement'] = R['IR_spatial_refinement']
# base48 -> IR12_fine isolates the IR addition (same N48, same seed band)
b = A['base48']
packet['base48_to_IR12_fine_percent'] = {key: dict(k8=100 * (f[key]['tot'] / b[key]['tot'] - 1), incl_tail=100 * (f[key]['tot_incl'] / b[key]['tot_incl'] - 1)) for key in ('14_24', '24_34')}
# time refinement
h = A['base48_half']
packet['time_refinement_abs'] = {key: abs(h[key]['tot'] - b[key]['tot']) for key in ('14_24', '24_34')}
packet['codex_time_refinement'] = R['time_refinement_maximum_total_action_absolute']
# drifts from RUN files
drift = {}
for tag in tags:
    run = json.load(open(PK / f'RUN_{tag}.json'))
    drift[tag] = dict(energy=run['energy_max_relative_drift'], charge=run['charge_max_relative_drift'],
                      export_max_charge=max(e['charge_relative'] for e in run['export_checks']),
                      export_max_energy=max(e['energy_relative'] for e in run['export_checks']), seconds=run.get('seconds'))
packet['drift'] = drift
packet['drift_max'] = dict(energy=max(d['energy'] for d in drift.values()), charge=max(d['charge'] for d in drift.values()),
                           export_charge=max(d['export_max_charge'] for d in drift.values()), export_energy=max(d['export_max_energy'] for d in drift.values()))
packet['codex_conservation'] = R['conservation']
packet['codex_action_ranges'] = R['negative_action_ranges_by_window']
out['A_packet_arithmetic'] = packet

# ---------------------------------------------------------------- B. legacy orbit averages (own trapezoid + Simpson + cycle mean)
z = np.load(EARLY / 'FIELD_CURVES.npz')
L = json.load(open(PK / 'LEGACY_ORBIT_AVERAGES.json'))
legacy = {}
def trap(y, x):  # own trapezoid
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * (x[1:] - x[:-1])))
def simpson(y, x):
    from scipy.integrate import simpson as sp
    return float(sp(y, x=x))
for row in L['rows']:
    tag = row['tag']
    t = z[tag + '_tau']; s = z[tag + '_eta']
    sel = (t >= 700 - 1e-9) & (t <= 800 + 1e-9)
    tt, ss = t[sel], s[sel]
    rec = dict(samples=int(sel.sum()), stats={})
    for name in ('energy_outside_spatial_mean_fraction', 'phi_plus_interaction_energy_fraction', 'spatial_mean_charge_fraction'):
        v = z[tag + '_' + name][sel]
        cm = trap(v, tt) / (tt[-1] - tt[0]); cf = trap(v, ss) / (ss[-1] - ss[0])
        cms = simpson(v, tt) / (tt[-1] - tt[0])
        cx = row['statistics'][name]
        rec['stats'][name] = dict(cosmic_mean=cm, cosmic_mean_simpson=cms, conformal_mean=cf, min=float(v.min()), max=float(v.max()), last=float(v[-1]),
                                  codex_cosmic=cx['cosmic_time_mean'], codex_conformal=cx['conformal_time_mean'],
                                  diff_cosmic=cm - cx['cosmic_time_mean'], diff_conformal=cf - cx['conformal_time_mean'],
                                  diff_min=float(v.min()) - cx['minimum'], diff_max=float(v.max()) - cx['maximum'], diff_last=float(v[-1]) - cx['last'])
    # integer-cycle mean of the energy-outside fraction: between successive maxima of the curve inside the window
    v = z[tag + '_energy_outside_spatial_mean_fraction'][sel]
    pk = [i for i in range(1, len(v) - 1) if v[i] >= v[i - 1] and v[i] > v[i + 1] and v[i] > 0.6]
    # keep peaks separated by at least a third of the typical period (period ~0.89 in eta ~ 44-48 in tau)
    peaks = []
    for i in pk:
        if not peaks or tt[i] - tt[peaks[-1]] > 15: peaks.append(i)
    if len(peaks) >= 2:
        i0, i1 = peaks[0], peaks[-1]
        rec['cycle_mean_energy_outside'] = dict(n_cycles=len(peaks) - 1, tau_from=float(tt[i0]), tau_to=float(tt[i1]),
                                                cosmic_mean=trap(v[i0:i1 + 1], tt[i0:i1 + 1]) / (tt[i1] - tt[i0]),
                                                period_tau=float((tt[i1] - tt[i0]) / (len(peaks) - 1)))
    legacy[tag] = rec
out['B_legacy_orbit_averages'] = legacy

# ---------------------------------------------------------------- C. scale chain (own arithmetic from the brief's constants)
lam_X, g2, lam_phi = 1e-3, 1e-2, 3e-2
a_m, C = 8.869218738806, 1.9857962750634526
kappa = 1.44067458242127e-24          # lambda_phi v^2 / H_i^2  [Q: pt_relax scales.json; brief's clock constants]
v_GeV = 0.924                           # [Q: gauge-point packet, via pt_matching_skeptic]
H_i_GeV = math.sqrt(lam_phi * v_GeV ** 2 / kappa)
scale = dict(kappa=kappa, implied_H_i_GeV=H_i_GeV, rows=[])
for r2, label in ((3600., 'brief |chi|=60'), (3927.19, 'equilibrium r_c^2+<|dpsi|^2> (pt_relax)'), (10104., 'measured medium s24-34'), (12100., 'brief |chi|=110')):
    a_on = math.sqrt(g2 * r2 / kappa)            # lambda_phi (v a/H_i)^2 = g^2 <r^2>
    s_on = (a_on - a_m) / C                      # a = a_m + C s (radiation era, dtau/ds = a)
    M2 = g2 * r2 / lam_phi                       # conformal phi vev^2 at onset
    mX = math.sqrt(g2 * M2); mphi = math.sqrt(2 * lam_phi * M2)
    H_conf = C / a_on                            # a'/a in conformal units
    # phi -> X Xbar at tree level: H = g^2/2 phi^2 r^2 with r^2 = chi1^2+chi2^2 (two real fields, coupling (g^2 M) delta_phi chi_i^2 / 2 each)
    # Gamma = sum_i mu^2/(32 pi m_phi) * beta with mu = g^2 M, identical-particle factor 1/2 per real field
    beta = math.sqrt(max(0.0, 1 - 4 * mX ** 2 / mphi ** 2))
    Gam_conf_realfields = 2 * (g2 * math.sqrt(M2)) ** 2 / (32 * math.pi * mphi) * beta
    # equivalent complex-field form: Phi = psi/sqrt2, vertex 2 g^2 M dphi |Phi|^2 -> Gamma = (2g^2 M)^2 beta/(16 pi m_phi)
    Gam_conf_complex = (2 * g2 * math.sqrt(M2)) ** 2 * beta / (16 * math.pi * mphi)
    scale['rows'].append(dict(label=label, r2=r2, a_onset=a_on, s_onset=s_on, Gamma_crit=1 / s_on, m_X_conf=mX, m_phi_conf=mphi, m_phi_over_m_X=mphi / mX,
                              H_conf=H_conf, Gamma_phi_XX_conf=Gam_conf_realfields, Gamma_phi_XX_over_mX=Gam_conf_realfields / mX,
                              Gamma_complex_form=Gam_conf_complex, Gamma_over_H=Gam_conf_realfields / H_conf,
                              eps_at_onset_omega2=(1 / s_on) / 2.0, eps_at_onset_omega8=(1 / s_on) / 7.71))
scale['phi_to_XX_open_condition'] = dict(requires='lambda_phi > 2 g^2', lambda_phi=lam_phi, two_g2=2 * g2, open=lam_phi > 2 * g2,
                                        candidate1_lambda_phi_0p01_ratio=math.sqrt(2 * 0.01 / g2))
# e-fold table with the skeptic's c_eff (own recomputation of the same V^4 T formula)
T = 17.556; s_on = 2.629e12
c_lo, c_hi = 3.023e-8, 4.111e-7
def ef(V4, c): return c * V4 * T * s_on
scale['efolds_recomputed'] = dict(
    X_only_V4=(4 * lam_X) ** 4, X_only=[ef((4 * lam_X) ** 4, c_hi)],
    g2_only=[ef((2 * g2) ** 4, c_lo), ef((2 * g2) ** 4, c_hi)],
    mixed=[ef((2 * g2) ** 2 * (6 * lam_phi) ** 2, c_lo), ef((2 * g2) ** 2 * (6 * lam_phi) ** 2, c_hi)],
    phi_self=[ef((6 * lam_phi) ** 4, c_lo), ef((6 * lam_phi) ** 4, c_hi)],
    phi_self_at_lambda_phi_0p01=[ef((6 * 0.01) ** 4, c_lo), ef((6 * 0.01) ** 4, c_hi)],
    mixed_at_lambda_phi_0p01=[ef((2 * g2) ** 2 * (6 * 0.01) ** 2, c_lo), ef((2 * g2) ** 2 * (6 * 0.01) ** 2, c_hi)],
    note='c_eff from pt_relax_skeptic mc24 (single real scalar, Bose, mu=0, m/T 0.36 and 0.16); V^4 T formula as in pt_relax; combinatorics O(1-10) uncomputed')
out['C_scale_chain'] = scale

# ---------------------------------------------------------------- D. cross-reads of stored track results (not prose)
X = {}
cs = json.load(open(TR / 'pt_relax/skeptic/results/cutoff_scan.json'))
def wins(N):
    return {f"{int(w['s'][0])}_{int(w['s'][1])}": w for w in cs['runs'][str(N)]['windows']}
W = {N: wins(N) for N in (32, 40, 48, 64)}
X['cutoff_all_mode_total'] = {N: {k: W[N][k]['hann'] for k in ('24_34', '94_104', '194_204', '244_254', '294_304', '394_404', '494_504') if k in W[N]} for N in (32, 40, 48, 64)}
drop = {N: W[N]['24_34']['hann'] - W[N]['244_254']['hann'] for N in (32, 40, 48, 64)}
X['drop_24_34_to_244_254'] = drop
X['early_rate_per_unit_s'] = {N: drop[N] / 220 for N in (32, 40, 48, 64)}
X['aitken_drop_29_249'] = cs['continuum_extrapolation']['drop_29_to_249']
X['continuum_early_rate_per_unit_s'] = cs['continuum_extrapolation']['drop_29_to_249']['aitken_continuum'] / 220
X['continuum_early_rate_over_Gamma_crit'] = X['continuum_early_rate_per_unit_s'] / 3.8e-13
X['late_N48_vs_N64'] = cs['late_N48_vs_N64']
X['split_494_504'] = cs['split_494_504']
mt = json.load(open(TR / 'pt_matching/skeptic/results/skC_transition.json'))
X['transition_pooled'] = {k: dict(nminus=mt['pooled'][k]['nminus'], Nphi=mt['pooled'][k]['Nphi'], nminus_over_J0=mt['pooled'][k]['nminus_over_J0'], Nphi_over_J0=mt['pooled'][k]['Nphi_over_J0']) for k in mt['pooled']}
X['transition_late_both_codes'] = dict(nminus=mt['late_eps1e-3_rand_both_codes'], Nphi=mt['late_eps1e-3_rand_both_codes_Nphi'])
X['transition_tilt_beta0p3'] = {k: mt['bias'][k]['nminus_over_J0'] for k in mt['bias'] if k.startswith('beta0.3')}
X['transition_eps_trend'] = mt['eps_trend']['tau20']; X['transition_seed_trend'] = mt['seed_trend']['tau20']
X['transition_uniform_test'] = mt['uniform_test']['tau20']
ms = json.load(open(TR / 'pt_matching/results/summary.json'))['key_numbers']
X['pt_matching_key_numbers'] = dict(J_r=ms['J_r_per_q_tau20'], A_minus=ms['A_minus_per_q_tau20'], late_J_r=ms['late_orbit_A_match']['J_r_per_q'], slow_limit_tau20=ms['slow_limit_tau20'], slow_limit_late=ms['slow_limit_late'])
pp = json.load(open(TR / 'pt_response/skeptic/results/phase_probes_all.json')); pp2 = json.load(open(TR / 'pt_response/skeptic/results/phase_probes_2092.json'))
ratios = {}
for src in (pp, pp2):
    for name, rec in src.items():
        for shell, rr in rec['rows'].items():
            # find a ratio-like key
            keyr = [k for k in rr if 'ratio' in k.lower()]
            ratios[f"s0={rec['s0']}|{shell}"] = {k: rr[k] for k in keyr} if keyr else list(rr.keys())
X['response_ratio_rows'] = ratios
sv = json.load(open(TR / 'pt_rerun/skeptic/results/sk_verdicts.json'))['tail_arithmetic']
X['pt_rerun_skeptic_tail'] = sv
eq = json.load(open(TR / 'pt_relax/skeptic/results/equilibria.json'))
X['equilibria'] = dict(quantum_hartree_phi_eq=eq['quantum_hartree']['phi_equilibrium']['Xbar_per_net'], quantum_hartree_phi_absent=eq['quantum_hartree']['phi_absent']['Xbar_per_net'],
                       T_phi_eq=eq['quantum_hartree']['phi_equilibrium']['T'], rj_lattice_Xbar_per_net={r['N']: r['Xbar_per_net'] for r in eq['rj_lattice']})
out['D_cross_reads'] = X

OUT.write_text(json.dumps(out, indent=1, default=float) + '\n')
print(json.dumps({k: (v if k != 'D_cross_reads' else {kk: vv for kk, vv in v.items() if kk != 'response_ratio_rows'}) for k, v in out.items() if k in ('A_packet_arithmetic', 'C_scale_chain', 'D_cross_reads')}, indent=1, default=float)[:200])
print('written', OUT)
