#!/usr/bin/env python3
"""kin_source (Claude, own code): the broader-grid parent spectrum (our hw_gates gate 1) as the second vector source.

Reads our own hw_gates results (tracks/hw_gates/results/gate1_vector.json: 49 modes k = 1e-5..1 H_i, n_k per mode,
profiles P0 constant vacuum mass and P1 evolving Higgs norm; e = 1e-7 and 1e-6) and Codex's 225028b
VECTOR_RESULTS.json modes (read through git show) for the in-band comparison only.
  N_A = int dln k  k^3 n_k / (2 pi^2)   (one longitudinal polarisation, H_i units).
Also: the in-band (k/k* = 0.2-5) part of our spectrum against Codex's band number; the number-weighted parent
momentum <k^2>^(1/2); and the bound on the daughter comoving-momentum smearing from the parent motion,
delta q = (k/2) cos(theta*) for a nonrelativistic parent, relative to q0.
Usage: python3 -B broad_spectrum.py ../results/broad_source.json
"""
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import kincommon as K

HWFILE = K.TRACKS / 'hw_gates' / 'results' / 'gate1_vector.json'


def integrand(k, n):
    return k**3 * n / (2 * math.pi**2)


def band_integral(k, dNdlnk, kmin, kmax):
    """Trapezoid in ln k with log-log linear interpolation at the band edges."""
    lk, li = np.log(k), np.log(dNdlnk)
    edges = np.log([kmin, kmax])
    inside = (lk > edges[0]) & (lk < edges[1])
    xs = np.concatenate([[edges[0]], lk[inside], [edges[1]]])
    ys = np.concatenate([[np.interp(edges[0], lk, li)], li[inside], [np.interp(edges[1], lk, li)]])
    return float(np.trapezoid(np.exp(ys), xs))


def main(out):
    hw = json.loads(HWFILE.read_text())
    vr = K.git_show_json(K.COMMIT_PARENT, 'research/reality_framework/grut_heavy_vector_budget_20261025/VECTOR_RESULTS.json')
    res = dict(method=__doc__, hw_point=hw['point'], spectra={})
    for sp in hw['spectra']:
        e, prof = sp['e'], sp['profile']
        if e not in (1e-7, 1e-6):
            continue
        k = np.array([m['k'] for m in sp['modes']]); n = np.array([m['n_k'] for m in sp['modes']])
        dn = integrand(k, n)
        lk = np.log(k)
        NA = float(np.trapezoid(dn, lk))
        kstar = math.sqrt(6 * e * 40.0)                 # k* = sqrt(m_A,vac/H_i), h/H_i = 40 exactly
        band = band_integral(k, dn, 0.2 * kstar, 5 * kstar)
        below = band_integral(k, dn, k[0], 0.2 * kstar)
        above = band_integral(k, dn, 5 * kstar, k[-1])
        k2 = float(np.trapezoid(dn * k * k, lk) / NA)
        vp = K.vector_params(e)
        # Codex modes at this e (profile 'actual' for 1e-7; e_sensitivity rows carry e=1e-6)
        cm = [m for m in vr['modes'] if abs(m['e'] / e - 1) < 1e-9 and m['profile'] in ('actual',)]
        comp = None
        if cm:
            ck = np.array([m['k_over_H'] for m in cm]); cd = np.array([m['dN_dlnk'] for m in cm])
            o = np.argsort(ck); ck, cd = ck[o], cd[o]
            ours_at = np.exp(np.interp(np.log(ck), lk, np.log(dn)))
            ratio = ours_at / cd
            codex_band = float(np.trapezoid(cd, np.log(ck)))
            comp = dict(codex_nodes=len(ck), codex_profile='actual', codex_band_trapezoid=codex_band,
                        ratio_ours_over_codex_min=float(ratio.min()), ratio_max=float(ratio.max()),
                        ratio_number_weighted=float(np.trapezoid(ours_at, np.log(ck)) / codex_band),
                        samples=[dict(k_over_kstar=float(kk / kstar), codex=float(c), ours_interp=float(u))
                                 for kk, c, u in zip(ck[::12], cd[::12], ours_at[::12])])
        key = f'e{e:.0e}_{prof}'
        res['spectra'][key] = dict(
            e=e, profile=prof, N_A_full_grid=NA, k_min=float(k[0]), k_max=float(k[-1]), k_star=kstar,
            N_A_in_codex_band=band, N_A_below_band=below, N_A_above_band=above, band_fraction=band / NA,
            codex_band_N_A=K.CODEX_BAND[e], in_band_ratio_ours_over_codex=band / K.CODEX_BAND[e],
            full_over_codex_band=NA / K.CODEX_BAND[e],
            parent_k_rms_over_H=math.sqrt(k2), parent_p_over_mA_at_1pct_decay=math.sqrt(k2) / (vp['mA_over_H'] * 0.1419 * vp['a_unit']),
            daughter_q_smearing_rms_over_q0=0.5 * math.sqrt(k2 / 3) / vp['q0_over_H'],
            daughter_q_smearing_max_over_q0=0.5 * float(k[-1]) / vp['q0_over_H'],
            relativistic_parent_decay_fraction_bound=-math.expm1(-vp['gamma'] * (10 * float(k[-1]) / vp['mA_over_H'])**2 / 2),
            dN_dlnk=[dict(k_over_H=float(kk), k_over_kstar=float(kk / kstar), n_k=float(nn), dN_A_dlnk=float(dd))
                     for kk, nn, dd in zip(k, n, dn)],
            codex_comparison=comp)
    # Codex's own constant-vacuum control vs our P0, in band
    res['codex_constant_vacuum_band'] = vr['spectra']['constant_vacuum']['N_A_band']
    K.save(out, res)
    for key, s in res['spectra'].items():
        print('%-14s N_A=%.5f band=%.5f (%.3f of full; ours/Codex in band %.4f; full/Codex band %.4f) below=%.4f above=%.5f k_rms=%.4f dq/q0=%.1e'
              % (key, s['N_A_full_grid'], s['N_A_in_codex_band'], s['band_fraction'], s['in_band_ratio_ours_over_codex'],
                 s['full_over_codex_band'], s['N_A_below_band'], s['N_A_above_band'], s['parent_k_rms_over_H'],
                 s['daughter_q_smearing_rms_over_q0']))
        if s['codex_comparison']:
            c = s['codex_comparison']
            print('   vs Codex modes: ratio %.4f..%.4f, number-weighted %.4f' % (c['ratio_ours_over_codex_min'], c['ratio_max'], c['ratio_number_weighted']))


if __name__ == '__main__':
    main(sys.argv[1])
