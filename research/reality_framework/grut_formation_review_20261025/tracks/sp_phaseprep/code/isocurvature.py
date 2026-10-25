#!/usr/bin/env python3
"""sp_phaseprep: observable-universe coherence, CMB isocurvature and order-of-magnitude new-source estimates (own code).

Reads results/spectrum.json, results/transfer.json, results/ppi.json.

1. Coherence of the source phase across the observable universe (baseline action, light mode m = 0.61 H):
   long-wavelength rms on the observable scale x_* = e^{-N_*}:  sigma_L = C_nu/(2 pi f) x_*^{(3-2nu)/2}/sqrt(3-2nu)
   (linear, Gaussian), memory factor e^{-(3/2-nu) N_*} of a value set when the observable scale left the horizon,
   and the Gaussian exponent for the observable-scale mean to sit inside the band around pi/2.
2. CDM isocurvature (uncorrelated) S = T dpsi_CMB, dpsi_CMB = C_nu/(2 pi f) e^{-(3/2-nu) N_*} per ln k (real nu;
   m >= 3H/2 is suppressed as e^{-3N/2}), against S_max = 9.1e-6 [Q: Planck 2018 beta_iso < 0.038, as used by our
   affleck_dine and gp_matching tracks].  T from transfer.json (local, and smoothed over sigma_S = 0.03 / 0.29).
3. Order-of-magnitude [E] new sources from a PPI: sudden-quench Bogoliubov production when the PPI masses are removed,
   and a Planck-suppressed inflaton decay channel; tensor ratio implied by H_CMB >= H_end (retrospective).

Usage: python3 isocurvature.py ../results/isocurvature.json
"""
import json
import sys
from math import gamma
from pathlib import Path
import numpy as np

R = Path(__file__).resolve().parents[1] / 'results'
SP = json.load(open(R / 'spectrum.json'))
TRF = json.load(open(R / 'transfer.json'))
PPI = json.load(open(R / 'ppi.json'))
SMAX = 9.1e-6
H_GEV = 1.1211167304108865e13
MRED = 2.435e18
AS = 2.1e-9


def Cnu(nu):
    return 2**(nu - 1.5) * gamma(nu) / gamma(1.5)


def dpsi_cmb(m, f, N):
    nu2 = 2.25 - m * m
    if nu2 <= 0:
        return np.exp(-1.5 * N) / (2 * np.pi * f)   # heavy: amplitude ~ e^{-3N/2} (order of magnitude)
    nu = np.sqrt(nu2)
    return Cnu(nu) / (2 * np.pi * f) * np.exp(-(1.5 - nu) * N)


def main(out):
    res = {}
    L = SP['H_star']['light']
    m0, f0 = L['mass_over_H'], L['f_eff_over_H']
    nu0 = np.sqrt(2.25 - m0 * m0)
    band_gap = np.pi / 2 - 0.09
    coh = {}
    for N in (50, 55, 60):
        x = np.exp(-N)
        sL = Cnu(nu0) / (2 * np.pi * f0) * x**((3 - 2 * nu0) / 2) / np.sqrt(3 - 2 * nu0)
        coh[f'N{N}'] = dict(sigma_L_rad=float(sL), memory_factor=float(np.exp(-(1.5 - nu0) * N)),
                            gaussian_exponent_for_band=float(band_gap**2 / (2 * sL**2)),
                            release_patches_in_observable_universe_log10=float(3 * (N - np.log(40)) / np.log(10)))
    res['coherence'] = dict(m_over_H=m0, f_eff=f0, nu=float(nu0), rows=coh,
                            correlation_power=float(3 - 2 * nu0),
                            note='linear Gaussian extrapolation; the exponent only shows that the probability is nil')
    # transfers
    T = dict(local=abs(TRF['local']['T_local_lnNX']))
    for k in ('smoothed_sigma0.03', 'smoothed_sigma0.29'):
        if k in TRF:
            T[k] = abs(TRF[k]['T_smoothed']); T[k + '_sem'] = TRF[k]['T_smoothed_sem']
    res['transfer'] = T
    # isocurvature versus phase mass
    cases = [('baseline light phase (displacement imposed, not an attractor)', m0, f0)]
    for key in ('tuned_cS0.5', 'tuned_cS1.0', 'tuned_cS2.0'):
        d = PPI['designs'][key]
        cases.append((f'PPI {key}', d['spectrum']['light_phase']['mass_over_H'], d['spectrum']['light_phase']['f_eff']))
    rows = []
    for name, m, f in cases:
        r = dict(case=name, m_over_H=m, f_eff=f)
        for N in (50, 60):
            dp = dpsi_cmb(m, f, N)
            r[f'dpsi_cmb_N{N}'] = float(dp)
            for tk in ('local', 'smoothed_sigma0.03', 'smoothed_sigma0.29'):
                if tk in T:
                    r[f'S_over_Smax_N{N}_{tk}'] = float(T[tk] * dp / SMAX)
        rows.append(r)
    res['isocurvature_cases'] = rows
    # minimal mass for S < S_max
    ms = np.linspace(0.3, 1.5, 1201)
    mins = {}
    for N in (50, 60):
        for tk in ('local', 'smoothed_sigma0.03', 'smoothed_sigma0.29'):
            if tk in T:
                ok = [m for m in ms if T[tk] * dpsi_cmb(m, f0, N) < SMAX]
                mins[f'N{N}_{tk}'] = float(min(ok)) if ok else None
    res['minimal_mass_for_isocurvature'] = mins
    # 3a. sudden-quench production when PPI masses vanish at release [E]
    #     single mode: |beta_k|^2 = (w1 - w2)^2 / (4 w1 w2), w_i = sqrt(k^2 + m_i^2); n = int d^3k/(2pi)^3 |beta|^2
    #     (flat-space estimate at a = 1; modes k < H are not in a flat vacuum, so this is an order of magnitude).
    q = []
    for key in ('tuned_cS0.5', 'tuned_cS1.0', 'tuned_cS2.0'):
        mb = np.array(PPI['designs'][key]['mass_over_H'])
        ma = np.array(PPI['designs']['tuned_cS0.0']['mass_over_H'])
        n = 0.0
        k = np.linspace(1e-4, 30, 300001)
        for m1, m2 in zip(np.sort(mb)[1:], np.sort(ma)[1:]):   # skip the gauge zero mode
            w1 = np.sqrt(k * k + m1 * m1); w2 = np.sqrt(k * k + m2 * m2)
            b2 = (w1 - w2)**2 / (4 * w1 * w2)
            n += np.trapezoid(k * k * b2, k) / (2 * np.pi**2)
        q.append(dict(design=key, n_quanta_over_H3=float(n), per_net_charge=float(n / 238.65)))
    res['sudden_quench_estimate'] = dict(rows=q, net_charge_NX=238.65, allowance_pairs_per_charge=0.00493,
                                         note='flat-space sudden approximation, masses before (PPI) vs after (no PPI) '
                                              'paired by order; order of magnitude only [E]')
    # 3b. Planck-suppressed inflaton decay through the A-term coupling [E]: g = |muP| H / M_red (dimensionless
    #     coupling of I Sigma* Y^2 if |muP| H = g <I> with <I> ~ M_red), Gamma3 ~ g^2 m_I / (256 pi^3),
    #     branching B = Gamma3 / Gamma_I with Gamma_I ~ H (instant reheating), dark quanta per entropy B * 3 T_R/(4 m_I).
    TR = (90 * MRED**2 * H_GEV**2 / (np.pi**2 * 106.75))**0.25
    yx = 65.98 * 4.8294e-9
    dec = []
    for key in ('tuned_cS1.0', 'tuned_cS2.0'):
        g = PPI['designs'][key]['abs_muP_over_H'] * H_GEV / MRED
        for mI_over_H in (1.0, 10.0):
            mI = mI_over_H * H_GEV
            G3 = g * g * mI / (256 * np.pi**3)
            B = G3 / H_GEV
            nds = B * 3 * TR / (4 * mI)
            dec.append(dict(design=key, g=float(g), m_I_over_H=mI_over_H, branching=float(B), n_over_s=float(nds),
                            pairs_per_net_charge=float(3 * nds / (2 * yx))))
    res['inflaton_decay_estimate'] = dict(T_R_GeV=float(TR), rows=dec, note='model-dependent order of magnitude [E]')
    # 3c. tensor-to-scalar ratio if H_CMB >= H_end (retrospective) [D]
    res['tensor_ratio_lower'] = float(2 * H_GEV**2 / (np.pi**2 * MRED**2 * AS))
    Path(out).write_text(json.dumps(res, indent=1) + '\n')
    print(json.dumps(res, indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
