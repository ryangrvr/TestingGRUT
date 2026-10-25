#!/usr/bin/env python3
"""Track sp_critic (Claude, own code): consolidation checks for the review of Codex 373d37e.

Parts
  A  Y-decay energy conventions (Codex H = Gamma, sp_point t = 1/Gamma, sp_yledger decay average) [D, C].
  B  Combined leftover ledger assembled from the track JSONs, with the pair -> f conversion re-derived [C].
  C  Independent check of the frozen superhorizon gauge angle (sp_vector piece A): longitudinal Proca mode
     (z^2 B')' + a^2 m^2 B = 0, z^2 = a^2 m^2/(k^2 + a^2 m^2) (Graham-Mardon-Rajendran action), BD start at
     k/aH = 1000, de Sitter -> instant radiation, constant m_vac (P0) versus m_inf during inflation with
     m^2 = m_vac^2 + (m_inf^2 - m_vac^2)/a^2 after release (P1) and an instant step (P2) [C].
  E  Codex robustness ensembles: ensemble-averaged Y energy against the case-L lump threshold [Q inputs, C].
  D  Scan-wide Y-debris energy: for every published Codex grid release (SCAN, REFINEMENT, LOCAL_WINDOW csv,
     read from a git-archive extraction), the decay-averaged Y energy R_Y = <a> J_Y (M_Y/H_i)/E_X and the
     Y pairs per net charge, for releases with |eps_X| >= 0.19 [Q inputs, C].
Usage: python3 -B critic_checks.py <packet_dir> <tracks_dir> <out.json> [rmax_relic*.json ...]
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma as Gfun

HI = 1.1211167304108865e13
MY, H_S, MU, KAPPA, E_G = 1e12, 1e11, 5e11, 1e-3, 1e-7
YQ, BDM, MX = 4.8294009076988175e-9, 4.4e-10, 0.1 * 0.924089637375362


def part_a(L):
    G = 3 * KAPPA**2 * MY / (256 * math.pi**3)
    mL = math.sqrt(MY**2 - 2 * MU * H_S / math.sqrt(2))
    GL = G * mL / MY
    EX = 741.64988
    JY, NYe = 13.26903, L['numbers_late']['NY_modes']
    a_HG = math.sqrt(HI / G)
    a_t = math.sqrt(2 * HI / G)
    a_avg = Gfun(1.5) * math.sqrt(2 * HI / GL)
    R = lambda a, N, m: a * N * (m / HI) / EX
    out = dict(Gamma=G, Gamma_L=GL, m_L=mL,
               a_H_eq_Gamma=a_HG, a_t_eq_1_over_Gamma=a_t, a_decay_average_L=a_avg,
               R_codex_H_eq_Gamma=R(a_HG, JY, MY),
               R_sp_point_t_eq_1_over_Gamma=R(a_t, JY, MY),
               R_decay_average_JY_MY=R(Gfun(1.5) * a_t, JY, MY),
               R_decay_average_NYeig_mL=R(a_avg, NYe, mL),
               note='Energy deposited into a radiation-like hidden sector by an exponential decay in radiation '
                    'domination is N m <a>, <a> = Gamma(3/2) sqrt(2 H_i/Gamma); H = Gamma is low by 1/1.2533 and '
                    't = 1/Gamma high by 1/0.8862. Correct value: 349 E_X (sp_yledger).')
    return out


def part_b(L, V):
    allow = 0.01 * BDM / (2 * MX * YQ)               # pairs per net charge at 1% of calibrated DM
    f_per_pair = 2 * MX * YQ / BDM
    eps = abs(L['windows']['2000-3000']['eps'])
    NX = abs(L['numbers_late']['NX'])
    p = L['pairs_per_net_charge']
    pY, pS = p['Y_daughters'], p['sigma_vectors']
    pYH = p['YH_vectors_range']
    pA = [V['piece_A']['P1']['pairs_actual_superhorizon'], V['piece_A']['P1']['pairs_actual']]
    pB2_sup = V['piece_B2']['codex_0.03']['range_pairs_per_net']
    pB2_nat = V['piece_B2']['natural_0.489']['range_pairs_per_net']
    f_rem = (MX / (BDM / YQ)) * (1 / eps - 1)
    deb = [pY + pS + pYH[0] + pA[0], pY + pS + pYH[1] + pA[1]]
    led = dict(allowance_pairs_per_net=allow, f_per_pair=f_per_pair, eps_final=eps, NX_final=NX,
               pairs=dict(Y_daughters=pY, sigma_to_AA=pS, YH_to_YL_A=pYH, inflationary_vector_pieceA=pA,
                          release_goldstone_B2_supplied_sigma003=pB2_sup, release_goldstone_B2_natural=pB2_nat,
                          release_goldstone_B1='unresolved'),
               debris_pairs_computed=deb, debris_f=[x * f_per_pair for x in deb],
               debris_over_allowance=[x / allow for x in deb],
               X_remainder_frozen_f=f_rem,
               frozen_total_f=[f_rem + x * f_per_pair for x in deb],
               frozen_total_with_B2_supplied_f=[f_rem + (deb[0] + pB2_sup[0]) * f_per_pair,
                                                f_rem + (deb[1] + pB2_sup[1]) * f_per_pair],
               natural_spread_remainder_f={str(e): (MX / (BDM / YQ)) * (1 / e - 1) for e in (0.048, 0.058, 0.072)},
               no_double_count=('Pairs are disjoint particle populations: Y_L daughters (1.5 per Y quantum, heavy quanta '
                                'counted once in N_Y and their vector separately), two vectors per sigma quantum, one '
                                'pair per vector. The X remainder is the cold condensate symmetric part (gp_leftovers '
                                'number count). In equilibrium only energies enter, from disjoint Hamiltonian pieces.'))
    return led


def proca_ratio(k, mvac, minf, profile, a_end):
    """Late adiabatic number of a longitudinal mode (complex B, Pi = z^2 B')."""
    def m2(a, infl):
        if profile == 'P0':
            return mvac**2
        if infl:
            return minf**2
        if profile == 'P1':
            return mvac**2 + (minf**2 - mvac**2) / a**2
        return mvac**2                                            # P2 instant step

    def rhs(t, y, infl):
        a = -1 / t if infl else t + 2
        mm = m2(a, infl)
        z2 = a * a * mm / (k * k + a * a * mm)
        B = y[0] + 1j * y[1]; P = y[2] + 1j * y[3]
        dB = P / z2; dP = -a * a * mm * B
        return [dB.real, dB.imag, dP.real, dP.imag]

    t0 = -1000 / k
    a0 = -1 / t0
    mm = m2(a0, True)
    w2 = k * k + a0 * a0 * mm
    z = a0 * math.sqrt(mm / w2)
    v = np.exp(-1j * k * t0) / math.sqrt(2 * k)
    dv = -1j * k * v
    zp_over_z = (-1 / t0) * (k * k / w2)          # z'/z = (a'/a) k^2/w^2 at constant m, a'/a = -1/t
    B = v / z
    P = z * (dv - zp_over_z * v)
    y0 = [B.real, B.imag, P.real, P.imag]
    s1 = solve_ivp(rhs, (t0, -1.0), y0, args=(True,), method='DOP853', rtol=1e-11, atol=1e-300)
    s2 = solve_ivp(rhs, (-1.0, a_end - 2), s1.y[:, -1], args=(False,), method='DOP853', rtol=1e-11, atol=1e-300)
    a = a_end
    mm = m2(a, False)
    w = math.sqrt(k * k + a * a * mm)
    z2 = a * a * mm / w**2
    B = s2.y[0, -1] + 1j * s2.y[1, -1]; P = s2.y[2, -1] + 1j * s2.y[3, -1]
    n = (abs(P)**2 / (2 * z2) + a * a * mm * abs(B)**2 / 2) / w - 0.5
    return n, s1.nfev + s2.nfev


def part_c():
    mvac = 6 * E_G * H_S / HI
    Fg = 34.8322
    minf = E_G * Fg
    a_star = 1 / math.sqrt(mvac)
    out = dict(m_vac=mvac, m_inf=minf, mass_ratio_sq=(mvac / minf)**2, a_star=a_star, rows=[])
    for k in (3e-6, 1e-5, 3e-5):
        rows = {}
        for prof in ('P0', 'P1', 'P2'):
            for aend in (20 * a_star, 30 * a_star):
                n, nf = proca_ratio(k, mvac, minf, prof, aend)
                rows[f'{prof}_aend{aend / a_star:.0f}'] = dict(n=n, nfev=nf)
        r = dict(k=k, **rows,
                 ratio_P1_P0=rows['P1_aend30']['n'] / rows['P0_aend30']['n'],
                 ratio_P2_P0=rows['P2_aend30']['n'] / rows['P0_aend30']['n'],
                 conv_P0_20_30=rows['P0_aend20']['n'] / rows['P0_aend30']['n'] - 1)
        out['rows'].append(r)
        print('C', k, r['ratio_P1_P0'], r['ratio_P2_P0'], r['conv_P0_20_30'], flush=True)
    return out


def regularity(pkt, tol=1e-3):
    """Codex's own 1e-9 rad test: coarse rows come in (delta, delta + 1e-9) pairs; refinement rows carry
    phase_shift_abs_delta_epsX. Returns {(H rounded, cH, delta rounded): bool}."""
    reg = {}
    with open(Path(pkt) / 'SCAN_RESULTS.csv') as fh:
        rows = list(csv.DictReader(fh))
    groups = {}
    for r in rows:
        groups.setdefault((round(float(r['H_GeV']), -3), float(r['cH']), round(float(r['delta']), 6)), []).append(r)
    for k, v in groups.items():
        if len(v) == 2:
            v.sort(key=lambda r: float(r['delta']))
            ok = abs(float(v[1]['epsX']) - float(v[0]['epsX'])) < tol
            for r in v:
                reg[(k[0], k[1], round(float(r['delta']), 5))] = ok
    with open(Path(pkt) / 'REFINEMENT_RESULTS.csv') as fh:
        for r in csv.DictReader(fh):
            if r.get('phase_shift_abs_delta_epsX'):
                reg[(round(float(r['H_GeV']), -3), float(r['cH']), round(float(r['delta']), 5))] = \
                    float(r['phase_shift_abs_delta_epsX']) < tol
    return reg


def part_d(pkt, rmax):
    G = 3 * KAPPA**2 * MY / (256 * math.pi**3)
    rows = []
    for name in ('SCAN_RESULTS.csv', 'REFINEMENT_RESULTS.csv', 'LOCAL_WINDOW_RESULTS.csv'):
        with open(Path(pkt) / name) as fh:
            for r in csv.DictReader(fh):
                try:
                    H, eps, JY, EX = float(r['H_GeV']), float(r['epsX']), float(r['JY_approx']), float(r['EX'])
                    NX, Y = float(r['NX']), float(r['YX_over_YQ'])
                except (KeyError, ValueError):
                    continue
                a_avg = Gfun(1.5) * math.sqrt(2 * H / G)
                rows.append(dict(src=name, H=H, cH=float(r['cH']), delta=float(r['delta']), eps=eps, YX=Y,
                                 JY=JY, EX=EX, NX=NX, R_Y=a_avg * JY * (MY / H) / EX,
                                 pairs_Y=1.5 * JY / abs(NX) if NX else float('inf')))
    big = [r for r in rows if abs(r['eps']) >= 0.19]
    both = [r for r in big if abs(r['YX']) >= 1]
    out = dict(n_releases=len(rows), n_eps_ge_019=len(big), n_eps_ge_019_and_Y_ge_YQ=len(both))
    for lab, sel in (('eps_ge_019', big), ('eps_ge_019_Y_ge_YQ', both)):
        if not sel:
            continue
        rmin = min(sel, key=lambda r: r['R_Y'])
        pmin = min(sel, key=lambda r: r['pairs_Y'])
        out[lab] = dict(min_R_Y=rmin, min_pairs_Y=pmin,
                        R_Y_quantiles=np.quantile([r['R_Y'] for r in sel], [0, 0.1, 0.5, 0.9, 1]).tolist(),
                        pairs_Y_quantiles=np.quantile([r['pairs_Y'] for r in sel], [0, 0.1, 0.5, 0.9, 1]).tolist())
    allr = np.array([r['R_Y'] for r in rows])
    out['all_R_Y_quantiles'] = np.quantile(allr, [0, 0.1, 0.5, 0.9, 1]).tolist()
    out['min_R_Y_any_release'] = min(rows, key=lambda r: r['R_Y'])
    if rmax:
        out['R_max_benchmark'] = rmax
        # case H: f = 0.01 threshold, log-interpolated in eps; case L: bracket of the condensate (lump) threshold
        hpts = sorted((float(k.split('_')[0]), v['R_max']) for k, v in rmax['H'].items() if v.get('R_max'))
        e_h = np.array([e for e, _ in hpts]); lr_h = np.log([r for _, r in hpts])
        def rmax_h(e):
            return float(np.exp(np.interp(e, e_h, lr_h)))
        for r in big:
            r['R_max_H'] = rmax_h(abs(r['eps'])); r['R_over_RmaxH'] = r['R_Y'] / r['R_max_H']
            r['eps_beyond_Hgrid'] = abs(r['eps']) > e_h[-1]
        out['case_H'] = dict(thresholds=hpts,
                             min_R_over_RmaxH=min(r['R_over_RmaxH'] for r in big),
                             argmin=min(big, key=lambda r: r['R_over_RmaxH']),
                             n_beyond_grid=sum(r['eps_beyond_Hgrid'] for r in big),
                             n_pass=sum(r['R_over_RmaxH'] < 1 for r in big))
        lb = rmax.get('L_brackets', {})
        out['case_L'] = dict(brackets=lb)
        if lb:
            # per-release lump threshold, linearly interpolated in eps between the bracket ends; beyond the largest
            # scanned eps the last lower bound is used (conservative for "pass", the threshold rises with eps)
            es = sorted(float(e) for e in lb)
            lo_v = [lb[str(e) if str(e) in lb else repr(e)][0] for e in es]
            hi_v = [lb[str(e) if str(e) in lb else repr(e)][1] for e in es]
            hi_v = [h if h is not None else float('inf') for h in hi_v]
            def thr(e, v):
                return float(np.interp(e, es, v))
            reg = regularity(pkt)
            for r in rows:
                e = abs(r['eps'])
                # below the smallest scanned eps the threshold is lower still: lower bound 0, upper bound = first upper
                r['RL_lo'] = thr(e, lo_v) if e >= es[0] else 0.0
                r['RL_hi'] = (thr(e, hi_v) if all(np.isfinite(hi_v)) or e <= es[-2] else float('inf')) \
                    if e >= es[0] else hi_v[0]
                r['regular_1e9'] = reg.get((round(r['H'], -3), r['cH'], round(r['delta'], 5)))
            prod = [r for r in rows if abs(r['eps']) >= 0.01 and abs(r['YX']) >= 1 and r['regular_1e9'] is True]
            out['case_L_regular_productive'] = dict(
                selection='|eps| >= 0.01, |Y| >= Y_Q, passes Codex 1e-9 test (coarse pairs or refinement shift)',
                n=len(prod), min_R_Y=min(r['R_Y'] for r in prod),
                n_RY_below_lower_threshold=sum(r['R_Y'] < r['RL_lo'] for r in prod),
                n_RY_below_upper_threshold=sum(r['R_Y'] < r['RL_hi'] for r in prod),
                closest=[dict(H=r['H'], cH=r['cH'], delta=r['delta'], eps=r['eps'], YX=r['YX'], R_Y=r['R_Y'],
                              RL_lo=r['RL_lo'], RL_hi=r['RL_hi'], pairs_Y=r['pairs_Y'], src=r['src'])
                         for r in sorted(prod, key=lambda x: x['R_Y'] / max(x['RL_lo'], 1e-9))[:8]])
            below = [r for r in big if r['R_Y'] < r['RL_lo']]
            out['case_L'].update(
                n_RY_below_lower_threshold=len(below),
                n_RY_below_upper_threshold=sum(r['R_Y'] < r['RL_hi'] for r in big),
                below_lower=[dict(H=r['H'], cH=r['cH'], delta=r['delta'], eps=r['eps'], R_Y=r['R_Y'],
                                  RL_lo=r['RL_lo'], pairs_Y=r['pairs_Y'], regular_1e9=r['regular_1e9'], src=r['src'])
                             for r in sorted(below, key=lambda x: x['R_Y'])],
                n_below_lower_and_regular=sum(1 for r in below if r['regular_1e9'] is True),
                n_below_lower_regularity_unknown=sum(1 for r in below if r['regular_1e9'] is None))
    out['note'] = ('Codex single-mass J_Y and E_X at a = 600-1000 [Q]; decay-averaged <a> with Gamma = 377.95 GeV. '
                   'Eigenmode counting (+3%) and m_L (-3.6%) nearly cancel. Sigma-radial and vector debris are not in '
                   'the published columns, so R_Y is a lower bound on the injected energy.')
    return out


def part_e(pkt, rmax):
    """Codex's published robustness ensembles (ENSEMBLE_MEMBERS.csv; sigma = 0.001 for 23 candidates, 0.03 for the
    leading point): ensemble-averaged decay-averaged Y energy R_ens = <a> <J_Y> (M_Y/H)/<E_X> against the case-L
    lump threshold at |<eps>| [Q inputs, C]."""
    G = 3 * KAPPA**2 * MY / (256 * math.pi**3)
    lb = rmax['L_brackets']
    es = sorted(float(e) for e in lb)
    lo = [lb[str(e)][0] for e in es]
    hi = [lb[str(e)][1] if lb[str(e)][1] is not None else float('inf') for e in es]
    S = json.load(open(Path(pkt) / 'ENSEMBLE_SUMMARY.json'))
    meta = {c['ensemble_id']: c for key in ('initial_ensembles', 'candidate_ensembles', 'high_fraction_controls')
            for c in (S.get(key) or []) if isinstance(c, dict) and 'ensemble_id' in c}
    mem = {}
    with open(Path(pkt) / 'ENSEMBLE_MEMBERS.csv') as fh:
        for r in csv.DictReader(fh):
            mem.setdefault(r['ensemble_id'], []).append(r)
    out = []
    for eid, rows in mem.items():
        H = float(rows[0]['H_GeV'])
        JY = np.array([float(r['JY_approx']) for r in rows]); EX = np.array([float(r['EX']) for r in rows])
        NX = np.array([float(r['NX']) for r in rows]); eps = np.array([float(r['epsX']) for r in rows])
        Y = np.array([float(r['YX_over_YQ']) for r in rows])
        R = Gfun(1.5) * math.sqrt(2 * H / G) * JY.mean() * (MY / H) / EX.mean()
        e = abs(eps.mean())
        RLlo = float(np.interp(e, es, lo)) if e >= es[0] else 0.0
        RLhi = (float(np.interp(e, es, hi)) if e <= es[-2] else float('inf')) if e >= es[0] else hi[0]
        out.append(dict(ensemble_id=eid, sigma=meta.get(eid, {}).get('seed_sigma_rad'), H=H, cH=float(rows[0]['cH']),
                        delta_nominal=meta.get(eid, {}).get('delta'), N=len(rows), eps_mean=float(eps.mean()),
                        Y_mean=float(Y.mean()), R_ens=R, RL_lo=RLlo, RL_hi=RLhi,
                        pairs_Y_ens=float(1.5 * JY.mean() / abs(NX.mean())), below_lo=R < RLlo,
                        below_hi=(R < RLhi) if math.isfinite(RLhi) else None))   # None: upper bound not scanned
    out.sort(key=lambda x: x['R_ens'] / max(x['RL_lo'], 1e-9))
    for x in out:
        x['below_lo'] = bool(x['below_lo'])
        x['below_hi'] = None if x['below_hi'] is None else bool(x['below_hi'])
    return dict(n_ensembles=len(out), n_below_lower=int(sum(x['below_lo'] for x in out)),
                n_below_upper=int(sum(x['below_hi'] is True for x in out)),
                n_upper_unknown=int(sum(x['below_hi'] is None for x in out)), rows=out,
                note='Y-only energy (Sigma radial and vector debris not in the published columns): a lower bound.')


def main():
    pkt, tracks, out = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
    rmax = None
    if len(sys.argv) > 4:                        # one or more rmax_relic*.json files
        rmax = {'H': {}, 'L_brackets': {}}
        rowsL = {}
        for fn in sys.argv[4:]:
            d = json.load(open(fn))
            for k, v in d.get('R_max', {}).items():
                if k.endswith('_H'):
                    rmax['H'][k] = v
            for r in d['rows']:
                if r.get('case') == 'L' and 'f_extra' in r and r.get('ymult', 1.0) == 1.0:
                    rowsL.setdefault(round(r['eps'], 4), []).append((r['R'], r['f_extra']))
        for e, pts in rowsL.items():
            pts.sort()
            ok = [R for R, f in pts if f <= 0.01]; bad = [R for R, f in pts if f > 0.01]
            rmax['L_brackets'][str(e)] = [max(ok) if ok else None, min(bad) if bad else None]
    L = json.load(open(tracks / 'sp_yledger/results/point_ledger.json'))
    V = json.load(open(tracks / 'sp_vector/results/vector_ledger.json'))
    res = dict(method=__doc__)
    res['A_Y_energy'] = part_a(L)
    print('A', json.dumps(res['A_Y_energy'], indent=1))
    res['B_ledger'] = part_b(L, V)
    print('B', json.dumps(res['B_ledger'], indent=1, default=float))
    res['D_scan_Y_energy'] = part_d(pkt, rmax)
    print('D', json.dumps({k: v for k, v in res['D_scan_Y_energy'].items()}, indent=1, default=float)[:3000])
    if rmax:
        res['E_ensembles_vs_lump_threshold'] = part_e(pkt, rmax)
        print('E', json.dumps({k: v for k, v in res['E_ensembles_vs_lump_threshold'].items() if k != 'rows'}))
    res['C_frozen_angle'] = part_c()
    Path(out).write_text(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
