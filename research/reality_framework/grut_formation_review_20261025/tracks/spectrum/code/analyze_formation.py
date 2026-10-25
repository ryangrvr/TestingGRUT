#!/usr/bin/env python3
"""Spectrum-track analysis of one independent formation run (kg_radial.c output).
Usage: analyze_formation.py RUNDIR NAME OUTJSON
All numbers in Qbar units (physical quanta = value / g_X^2)."""
import json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spec_lib as S

run, name, outp = sys.argv[1], sys.argv[2], sys.argv[3]
meta = json.load(open(f'{run}/{name}_meta.json'))
Rext, Rin = meta['rext'], meta['rin']
t, d = S.load_ext(f'{run}/{name}', len(Rext))
H = S.load_hist(f'{run}/{name}', len(Rin))
th, E, Q = H[:, 0], H[:, 1], H[:, 2]
Q0, E0 = Q[0], E[0]
Qin = {R: H[:, 3 + 2 * j] for j, R in enumerate(Rin)}
Ein = {R: H[:, 4 + 2 * j] for j, R in enumerate(Rin)}
T = t[-1]
out = dict(run=name, h=meta['h'], dt=meta['dt'], T=T, L=meta['L'], Q0=Q0, E0=E0,
           max_rel_energy_error=float(np.max(np.abs(E / E0 - 1))), max_rel_charge_error=float(np.max(np.abs(Q / Q0 - 1))))

# ---- 1. core retention and core frequency
def late_mean(R, a, b):
    m = (th >= a) & (th <= b); return float(np.mean(Qin[R][m]) / Q0)
out['retention'] = {f'Q{int(R)}_mean_{a}_{b}': late_mean(R, a, b) for R in [30, 60, 90, 120] for (a, b) in [(9600, 12000), (21600, 24000), (45600, 48000)] if b <= T + 1}
out['retention_Q60_slope_per_1000_9600_12000'] = float(np.polyfit(th[(th >= 9600) & (th <= 12000)], Qin[60][(th >= 9600) & (th <= 12000)] / Q0, 1)[0] * 1000)
for (a, b) in [(21600, 24000), (45600, 48000)]:
    if T >= b:
        m = (th >= a) & (th <= b)
        out[f'retention_Q60_slope_per_1000_{a}_{b}'] = float(np.polyfit(th[m], Qin[60][m] / Q0, 1)[0] * 1000)
U5 = d[:, Rext.index(5.0), 0] + 1j * d[:, Rext.index(5.0), 1]
ph = np.unwrap(np.angle(U5))
core_w = {}
for a in range(2000, int(T), 2000):
    m = (t >= a) & (t < a + 2000)
    core_w[str(a)] = float(-np.polyfit(t[m], ph[m], 1)[0])
out['core_frequency_from_phase_at_r5'] = core_w

# ---- 2. per-radius time-domain fluxes, fixed-radius spectra, nonlinearity criteria
kb = [0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 12.0]
def binned(w, dN, dw, sign):
    k = np.sqrt(np.maximum(w * w - 1, 0)); res = []
    for lo, hi in zip(kb[:-1], kb[1:]):
        m = (np.sign(w) == sign) & (k >= lo) & (k < hi) & (np.abs(w) > 1)
        res.append(float(np.sum(dN[m]) * dw))
    return res
out['k_bins'] = kb
fr = {}
for j, R in enumerate(Rext):
    if R < 120: continue
    U, Ut, Ur, Z, Zt, Zr = S.series(d, j)
    fl = S.time_fluxes(t, U, Ut, Ur, Z, Zt, Zr)
    psi2 = np.abs(U)**2 / R**2; Sm = 1 + Z / R
    dm2 = Sm**2 - 1 + 0.1 * psi2
    # slow (dressing) part of the X mass shift: 40-unit running mean removes the neutral-radiation
    # oscillations (w >= sqrt 6), which are non-resonant for X; the oscillating remainder is reported separately
    nav = int(round(40.0 / (t[1] - t[0])))
    dm2_slow = np.convolve(dm2, np.ones(nav) / nav, mode='same')
    dm2_osc = dm2 - dm2_slow
    rec = dict(time_domain=dict(Qpass_T=float(fl['Qpass'][-1]), Epass_T=float(fl['Epass'][-1]), Es_pass_T=float(fl['Es_pass'][-1])))
    if R in Qin:
        rec['time_domain']['Q0_minus_Qinside_T'] = float(Q0 - Qin[R][-1]); rec['time_domain']['E0_minus_Einside_T'] = float(E0 - Ein[R][-1])
    arrived = t > R / 1.0
    late = arrived & (t > R / 1.0 + 100)
    rec['nonlinearity'] = dict(
        max_lambda_psi2=float(np.max(0.1 * psi2[arrived])), max_abs_S2m1=float(np.max(np.abs(Sm**2 - 1)[arrived])),
        max_abs_dm2=float(np.max(np.abs(dm2)[arrived])), max_abs_dm2_slow=float(np.max(np.abs(dm2_slow)[late])),
        max_abs_dm2_osc=float(np.max(np.abs(dm2_osc)[late])),
        max_abs_dm2_slow_after_t3000=float(np.max(np.abs(dm2_slow)[t > 3000])))
    windows = {}
    for Tw in [6000, 12000, 18000, 24000, 36000, 48000]:
        if Tw > T + 1: continue
        m = t <= Tw; tt = t[m]
        wt = S.taper(tt, -1e9, Tw, 500.)
        cs = S.charged_spectrum(tt, U[m], Ur[m], wt)
        ns = S.neutral_spectrum(tt, Z[m], Zr[m], wt)
        flw = S.time_fluxes(tt, U[m] * wt, Ut[m] * wt, Ur[m] * wt, Z[m] * wt, Zt[m] * wt, Zr[m] * wt)
        kk = np.sqrt(np.maximum(cs['w']**2 - 1, 0)); posm = cs['w'] > 1
        cumX = np.cumsum((cs['dNout'] * posm)[np.argsort(kk)]) * cs['dw']; ks = np.sort(kk)
        k10 = float(ks[np.searchsorted(cumX, 0.1 * cumX[-1])]); k50 = float(ks[np.searchsorted(cumX, 0.5 * cumX[-1])])
        dms = rec['nonlinearity']['max_abs_dm2_slow']
        fracX_k2_below_10dm2 = float(np.sum((cs['dNout'] * posm)[kk**2 < 10 * dms]) * cs['dw'] / cs['NX'])
        windows[str(Tw)] = dict(k10_X=k10, k50_X=k50, fracX_k2_below_10dm2slow=fracX_k2_below_10dm2,
            phase_beyond_R_at_k10=float(dms * R / (2 * k10)), NX=cs['NX'], NXbar=cs['NXbar'], EX=cs['EX'], EXbar=cs['EXbar'],
            NX_incoming=cs['NX_in'], NXbar_incoming=cs['NXbar_in'], flux_subthreshold=cs['flux_subthreshold'],
            tapered_time_domain_Qpass=float(flw['Qpass'][-1]), untapered_Qpass=float(S.time_fluxes(tt, U[m], Ut[m], Ur[m], Z[m], Zt[m], Zr[m])['Qpass'][-1]),
            Ns=ns['Ns'], Es=ns['Es'], Ns_incoming=ns['Ns_in'], neutral_Eflux_subthreshold=ns['Eflux_subthreshold'],
            X_bins=binned(cs['w'], cs['dNout'], cs['dw'], +1), Xbar_bins=binned(cs['w'], cs['dNout'], cs['dw'], -1))
    rec['windows'] = windows
    fr[str(R)] = rec
out['fixed_radius'] = fr

# ---- 3. spatial projections of snapshots
sp = {}
for Ts in [6000, 12000, 18000, 24000, 36000, 48000]:
    f = f'{run}/{name}_snap{Ts}.bin'
    if not os.path.exists(f): continue
    snap = S.load_snap(f)
    for Rw in [90., 120., 180., 240., 360.]:
        c = S.spatial_projection(snap, Rw, 30.)
        n = S.spatial_projection(snap, Rw, 30., m2=6., field='neutral')
        kbins_X = [float(c['NXn'][(c['k'] >= lo) & (c['k'] < hi)].sum()) for lo, hi in zip(kb[:-1], kb[1:])]
        kbins_Xb = [float(c['NXbn'][(c['k'] >= lo) & (c['k'] < hi)].sum()) for lo, hi in zip(kb[:-1], kb[1:])]
        Rn = int(round(Rw)); qin = Qin.get(Rw)
        sp[f'T{Ts}_Rw{Rn}'] = dict(NX=c['NX'], NXbar=c['NXbar'], EX=c['EX'], EXbar=c['EXbar'], Q_window_direct=c['Q_window_direct'],
                                   Ns=n['Ns'], Es=n['Es'], X_bins=kbins_X, Xbar_bins=kbins_Xb,
                                   Q0_minus_Qinside=float(Q0 - qin[np.argmin(np.abs(th - Ts))]) if qin is not None else None)
out['spatial'] = sp

# ---- 4. late-emission spectrogram at R=240 and 480: charge passing per k-band in 2000-unit blocks
sg = {}
for R in [240., 480.]:
    j = Rext.index(R); U, Ut, Ur, Z, Zt, Zr = S.series(d, j)
    blocks = {}
    for a in range(2000, int(T) - 2499, 2000):
        m = (t >= a - 500) & (t < a + 2500); tt = t[m]
        wt = S.taper(tt, a - 500, a + 2500, 500.)
        cs = S.charged_spectrum(tt, U[m], Ur[m], wt)
        blocks[str(a)] = dict(X_bins=binned(cs['w'], cs['dNout'], cs['dw'], +1), Xbar_bins=binned(cs['w'], cs['dNout'], cs['dw'], -1),
                              NX=cs['NX'], NXbar=cs['NXbar'])
    sg[str(int(R))] = blocks
out['spectrogram'] = sg
out['Qinside_final'] = {str(int(R)): float(Qin[R][-1]) for R in Rin}
json.dump(out, open(outp, 'w'), indent=1)
print('wrote', outp)
