#!/usr/bin/env python3
"""ls_ledger (Claude, own code): injection ledger for Codex b448f25 candidates 10 and 17 from OUR releases.

Inputs: results/ls_release.json (ls_release.py), results/ls_damped.json (ls_damped.py), and, for comparison only,
Codex's published JSON read from the git-archived packet (data, not code).
Ledger (disjoint; every population counted once):
  Y quanta (N_YH + N_YL): each ends as one Y_L -> XXX / XbXbXb; deposit m_L <a>_L per quantum (the heavy quantum's
      extra mass m_H - m_L leaves in its cascade vector, counted separately or dropped in the conservative ledger);
  Sigma radial quanta: sigma -> A_L A_L; each A -> X Xbar, deposit m_A <a>_A per vector (birth momentum redshifts);
  heavy cascade (central): N_YH extra vectors;  primordial vectors: Codex's value quoted [Q] (sp_critic confirmed).
Decay average: radiation era, <a> = Gamma(3/2) sqrt(2 H_i / Gamma) (sp_critic); checked here by direct
quadrature with the actual birth epoch and, for vectors, the two-stage sigma -> A -> X Xbar chain with time
dilation of the relativistic A.
Widths [D]: tree level, vacuum masses; Y_H -> Y_L A_L via Goldstone equivalence (coupling sqrt2 mu G Y1 Y2).
Usage: python3 -B ls_ledger.py <ls_release.json> <ls_damped.json> <codex_packet_dir> <out.json>
"""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad

MY, HV, MU, KAP, E = 1e12, 1e11, 5e11, 1e-3, 1e-7
MX = 0.092408963737536
YQ, DM = 4.8294009076988175e-9, 4.4e-10
MRED, GSTAR = 2.435e18, 106.75


def masses_widths():
    sv = HV / math.sqrt(2)
    mH = math.sqrt(MY**2 + 2 * MU * sv); mL = math.sqrt(MY**2 - 2 * MU * sv)
    ms = math.sqrt(2) * HV; mA = 6 * E * HV
    pHL = math.sqrt((mH**2 - (mL + mA)**2) * (mH**2 - (mL - mA)**2)) / (2 * mH)
    r = mA**2 / ms**2
    return dict(mYH=mH, mYL=mL, msig=ms, mA=mA, pA_YH=pHL, pA_sig=math.sqrt(ms**2 / 4 - mA**2),
                G_YL3=3 * KAP**2 * mL / (256 * math.pi**3), G_YH3=3 * KAP**2 * mH / (256 * math.pi**3),
                G_YHLA=(math.sqrt(2) * MU)**2 * pHL / (8 * math.pi * mH**2),
                G_sig=ms**3 / (32 * math.pi * HV**2) * math.sqrt(1 - 4 * r) * (1 - 4 * r + 12 * r * r),
                G_A=E**2 * mA / (48 * math.pi) * (1 - 4 * MX**2 / mA**2)**1.5)


def a_mean_closed(G, H):
    return math.gamma(1.5) * math.sqrt(2 * H / G)


def a_mean_from_birth(G, H, a_b):
    """<a> at decay for a NR particle born at a_b (t = a^2/2 in 1/H_i), clock starting at birth."""
    g = G / H; tb = a_b**2 / 2
    f = lambda t: g * math.exp(-g * (t - tb)) * math.sqrt(2 * t)
    sc = 1 / g
    return quad(f, tb, tb + 60 * sc, epsabs=0, epsrel=1e-11, limit=400, points=[tb + sc])[0]


def vector_deposit(GA, mA, H, a_b, p_phys):
    """Mean conformal energy (GeV, a E) deposited by A born at a_b with physical momentum p_phys (direct quadrature
    over a of the decay density with the time-dilated optical depth; independent of Codex's Laguerre inversion)."""
    q = a_b * p_phys; b = q / mA; g = GA / H
    def tau(a):   # (Gamma/H) int_{a_b}^a x^2/sqrt(x^2+b^2) dx, closed form
        P = lambda x: 0.5 * (x * math.sqrt(x * x + b * b) - b * b * math.asinh(x / b))
        return g * (P(a) - P(a_b))
    def dens(la):
        a = math.exp(la)
        return a * g * a * a / math.sqrt(a * a + b * b) * math.exp(-tau(a)) * math.sqrt(q * q + a * a * mA * mA)
    amid = math.sqrt(2 / g)
    pts = sorted([math.log(max(b, a_b * 1.0001)), math.log(amid)])
    lo, hi = math.log(a_b), math.log(amid * 12)
    s = 0.0; edges = [lo] + [p for p in pts if lo < p < hi] + [hi]
    for x0, x1 in zip(edges, edges[1:]):
        s += quad(dens, x0, x1, epsabs=0, epsrel=1e-10, limit=500)[0]
    return s


def vector_chain(mw, H):
    """sigma -> A A -> X Xbar: average the vector deposit over the sigma decay epoch (exponential in t from release)."""
    g = mw['G_sig'] / H
    def f(t):
        ab = math.sqrt(2 * t)
        return g * math.exp(-g * (t - 0.5)) * vector_deposit(mw['G_A'], mw['mA'], H, ab, mw['pA_sig'])
    return quad(f, 0.5, 0.5 + 40 / g, epsrel=1e-7, limit=200)[0]


def stats(x):
    x = np.asarray(x, float)
    return dict(mean=float(x.mean()), sem=float(x.std(ddof=1) / math.sqrt(len(x))) if len(x) > 1 else None, n=len(x))


def main(rel_path, damp_path, codex_dir, out):
    rel = json.loads(Path(rel_path).read_text()); H = rel['H_GeV']
    damp = json.loads(Path(damp_path).read_text()) if Path(damp_path).exists() else None
    cd = Path(codex_dir)
    cx_dec = json.loads((cd / 'DECAY_RESULTS.json').read_text())
    cx_mode = json.loads((cd / 'MODE_RESULTS.json').read_text())
    cx_vec = json.loads((cd / 'VECTOR_RESULTS.json').read_text())
    mw = masses_widths()
    T = (90 * MRED**2 * H**2 / (math.pi**2 * GSTAR))**0.25
    s_i = 2 * math.pi**2 / 45 * GSTAR * T**3
    yfac = H**3 / s_i
    aL = a_mean_closed(mw['G_YL3'], H); aA = a_mean_closed(mw['G_A'], H)
    checks = dict(
        a_mean_closed=dict(YL=aL, A=aA, sigma=a_mean_closed(mw['G_sig'], H), YH_cascade=a_mean_closed(mw['G_YHLA'], H)),
        YL_born_at_release_over_closed=a_mean_from_birth(mw['G_YL3'], H, 1.0) / aL,
        YL_born_at_a23_over_closed=a_mean_from_birth(mw['G_YL3'], H, 23.0) / aL,
        A_deposit_born_a10_sigma_over_mA_aA=vector_deposit(mw['G_A'], mw['mA'], H, 10.0, mw['pA_sig']) / (mw['mA'] * aA),
        A_deposit_born_a1000_YH_over_mA_aA=vector_deposit(mw['G_A'], mw['mA'], H, 1000.0, mw['pA_YH']) / (mw['mA'] * aA),
        A_chain_sigma_decay_averaged_over_mA_aA=vector_chain(mw, H) / (mw['mA'] * aA),
        Codex_vector_control_a10_sigma_over_mA_aA=cx_dec['vector_time_dilation_controls'][0]['mean_conformal_energy_GeV'] / (mw['mA'] * aA))
    width_compare = {k: (mw[k2], cx_dec['vacuum'][kc]) for k, k2, kc in [
        ('YL3', 'G_YL3', 'GammaYL3'), ('YHLA', 'G_YHLA', 'GammaYHLA'), ('sigmaAA', 'G_sig', 'GammaSigmaAA'),
        ('AXX', 'G_A', 'GammaAXX')]}
    pair_allow = 0.01 * DM / (2 * MX * YQ)
    cx_rows = {(r['candidate'], r['member']): r for r in cx_mode['endpoint_rows']}
    result = dict(H_GeV=H, masses_widths_GeV=mw, widths_ours_vs_codex=width_compare, decay_average_checks=checks,
                  yield_factor=yfac, pair_allowance=pair_allow, candidates={})
    for c in (10, 17):
        runs = [r for r in rel['runs'] if r['candidate'] == c and r['a_end'] == 1000.0 and r['rtol'] == 1e-11]
        mem = [r for r in runs if r['member'] >= 0]
        nom = [r for r in runs if r['member'] < 0][0]
        def late(r, src):
            if src == 'window':
                w = r['windows']['600-1000']
                return dict(NX=w['NX'], Ncirc=w['Ncirc'], EX=w['EX'], NYH=w['NYH'], NYL=w['NYL'], Nsig=w['Nsig'])
            s = r['snapshot']; w = r['windows']['600-1000']
            return dict(NX=w['NX'], Ncirc=w['Ncirc'], EX=w['EX'], NYH=s['NYH'], NYL=s['NYL'], Nsig=s['Nsig'])
        cres = {}
        for src in ('window', 'snapshot'):
            L = [late(r, src) for r in mem]
            mean = {k: float(np.mean([x[k] for x in L])) for k in L[0]}
            unitY = mw['mYL'] / H * aL; unitA = mw['mA'] / H * aA
            ny = mean['NYH'] + mean['NYL']
            RY = ny * unitY / mean['EX']; RS = 2 * mean['Nsig'] * unitA / mean['EX']; RH = mean['NYH'] * unitA / mean['EX']
            Rprim = float(np.mean(next(x for x in cx_vec['estimates'] if x['candidate'] == c)['corrected_R_range']))
            per = [dict(R_Y=(x['NYH'] + x['NYL']) * unitY / x['EX'], R_Sigma=2 * x['Nsig'] * unitA / x['EX'],
                        R_casc=x['NYH'] * unitA / x['EX']) for x in L]
            pairsY = 1.5 * ny / abs(mean['NX']); pairsS = 2 * mean['Nsig'] / abs(mean['NX']); pairsH = mean['NYH'] / abs(mean['NX'])
            cres[src] = dict(
                members=[r['member'] for r in mem], n=len(mem),
                eps_net=mean['NX'] / mean['Ncirc'], eps_mean=stats([x['NX'] / x['Ncirc'] for x in L]),
                yield_over_YQ=mean['NX'] * yfac / YQ, NX=mean['NX'], EX=mean['EX'],
                NYH=stats([x['NYH'] for x in L]), NYL=stats([x['NYL'] for x in L]), Nsig=stats([x['Nsig'] for x in L]),
                R_Y=RY, R_Sigma=RS, R_cascade_central=RH, R_cascade_max=ny * unitA / mean['EX'], R_primordial_quoted=Rprim,
                R_conservative=RY + RS, R_central=RY + RS + RH + Rprim,
                R_conservative_per_member=stats([p['R_Y'] + p['R_Sigma'] for p in per]),
                R_central_per_member=stats([p['R_Y'] + p['R_Sigma'] + p['R_casc'] + Rprim for p in per]),
                pairs=dict(Y=pairsY, Sigma=pairsS, cascade=pairsH, total_central=pairsY + pairsS + pairsH,
                           over_allowance=(pairsY + pairsS + pairsH) / pair_allow))
        # member-by-member comparison with Codex's projection of the same published member
        cmp = []
        for r in mem:
            cx = cx_rows[(c, r['member'])]
            w = r['windows']['600-1000']; s = r['snapshot']
            cmp.append(dict(member=r['member'], eps=(w['NX'] / w['Ncirc'], cx['eps_window']), NX=(w['NX'], cx['NX_window']),
                            EX=(w['EX'], cx['EX_window']), NYH=(s['NYH'], w['NYH'], cx['NYH']),
                            NYL=(s['NYL'], w['NYL'], cx['NYL']), Nsig=(s['Nsig'], w['Nsig'], cx['NSIG'])))
        rel_dev = {k: float(np.max([abs(x[k][0] / x[k][-1] - 1) for x in cmp])) for k in ('eps', 'NX', 'EX', 'NYH', 'NYL', 'Nsig')}
        # Codex's 48-member value restricted to OUR member subset (sampling vs method)
        sub = [cx_rows[(c, m)] for m in cres['window']['members']]
        cx_sub = dict(NYH=float(np.mean([x['NYH'] for x in sub])), NYL=float(np.mean([x['NYL'] for x in sub])),
                      Nsig=float(np.mean([x['NSIG'] for x in sub])),
                      eps_net=float(np.sum([x['NX_window'] for x in sub]) / np.sum([x['NX_window'] / x['eps_window'] for x in sub])))
        # Sigma radial decomposition (free vs static vs forced lines)
        dec = [r['sigma_decomposition'] for r in mem]
        sigdec = dict(N_formula=stats([d['Nsig_formula_mean'] for d in dec]), N_free=stats([d['N_free'] for d in dec]),
                      N_static=stats([d['N_static'] for d in dec]), N_forced_lines=stats([d['N_forced_lines_total'] for d in dec]),
                      dr_static=stats([d['dr_static_fit'] for d in dec]),
                      dr_centrifugal_pred=stats([d['dr_centrifugal_pred'] for d in dec]),
                      dr_Yforced_pred=stats([d['dr_Yforced_pred'] for d in dec]),
                      fit_residual=stats([d['fit_residual_rms_over_dr_rms'] for d in dec]))
        # stability and tolerance controls (nominal)
        ctrl = {}
        for r in rel['runs']:
            if r['candidate'] == c and (r['a_end'] != 1000.0 or r['rtol'] != 1e-11):
                ctrl[r['label']] = {k: {q: r['windows'][k][q] for q in ('eps', 'NX', 'EX', 'NYH', 'NYL', 'Nsig')}
                                    for k in r['windows'] if k in ('300-600', '600-1000', '1000-2000', '2000-3000')}
        nomw = nom['windows']
        ctrl[f'{c}_nominal'] = {k: {q: nomw[k][q] for q in ('eps', 'NX', 'EX', 'NYH', 'NYL', 'Nsig')} for k in nomw}
        cx = next(r for r in cx_dec['rows'] if r['candidate'] == c)
        cres.update(nominal_window=nomw['600-1000'], member_comparison=cmp, max_relative_deviation_vs_codex_same_member=rel_dev,
                    codex_on_our_member_subset=cx_sub, sigma_decomposition=sigdec, controls=ctrl,
                    codex=dict(eps_net=cx['epsilon_net'], yield_over_YQ=cx['ymult_raw'], NYH=cx['NYH'], NYL=cx['NYL'],
                               Nsig=cx['NSIG'], R_Y=cx['R_Y'], R_Sigma=cx['R_Sigma'], R_conservative=cx['R_lower'],
                               R_central=cx['R_central'], pairs_total_central=cx['pairs_total_central']))
        if damp:
            dr = [r for r in damp['runs'] if r['label'].startswith(f'{c}_')]
            pairs = []
            for r in dr:
                if r['label'].endswith('_G1') and '_m' in r['label']:
                    base = next(x for x in dr if x['label'] == r['label'][:-3] + '_G0')
                    pairs.append(dict(member=r['label'], N_decayed=r['N_sigma_decayed_final'],
                                      N_decayed_omegaE=r.get('N_sigma_decayed_omegaE_final'),
                                      ratio_omegaE=(r.get('N_sigma_decayed_omegaE_final') or float('nan')) / base['window_100_300']['Nsig'],
                                      N_late_undamped=base['window_100_300']['Nsig'],
                                      ratio=r['N_sigma_decayed_final'] / base['window_100_300']['Nsig'],
                                      eps_damped=r['window_100_300']['eps'], eps_undamped=base['window_100_300']['eps'],
                                      NY_damped=r['window_100_300']['NYH'] + r['window_100_300']['NYL'],
                                      NY_undamped=base['window_100_300']['NYH'] + base['window_100_300']['NYL'],
                                      frac_by_a=r['decayed_fraction_by_a']))
            nomd = {r['label']: dict(N_decayed=r['N_sigma_decayed_final'], N_decayed_omegaE=r.get('N_sigma_decayed_omegaE_final'),
                                     eps=r['window_100_300']['eps'], NX=r['window_100_300']['NX'],
                                     NY=r['window_100_300']['NYH'] + r['window_100_300']['NYL'],
                                     Nsig_left=r['window_100_300']['Nsig'], frac_by_a=r['decayed_fraction_by_a'])
                    for r in dr if 'nominal' in r['label']}
            ratio = stats([p['ratio'] for p in pairs]); ratio_w = stats([p['ratio_omegaE'] for p in pairs])
            epsr = stats([p['eps_damped'] / p['eps_undamped'] for p in pairs]); nyr = stats([p['NY_damped'] / p['NY_undamped'] for p in pairs])
            RS = cres['window']['R_Sigma']; RY = cres['window']['R_Y']
            cres['damped_sigma_model'] = dict(members=pairs, nominal=nomd, ratio_decayed_over_late=ratio,
                                              ratio_decayed_omegaE_over_late=ratio_w, eps_ratio=epsr, NY_ratio=nyr,
                                              R_Sigma_model_range=[RS * ratio_w['mean'], RS * ratio['mean']],
                                              R_Y_with_backreaction=RY * nyr['mean'],
                                              R_conservative_model_range=[RY * nyr['mean'] + RS * ratio_w['mean'],
                                                                          RY * nyr['mean'] + RS * ratio['mean']])
        result['candidates'][str(c)] = cres
    Path(out).write_text(json.dumps(result, indent=1, default=float))
    for c, v in result['candidates'].items():
        w = v['window']
        print(c, 'eps_net %.7f (Codex %.7f) yield %.5f (%.5f) NYH %.3f NYL %.3f Nsig %.3f | RY %.2f RS %.2f Rcons %.2f Rcen %.2f'
              % (w['eps_net'], v['codex']['eps_net'], w['yield_over_YQ'], v['codex']['yield_over_YQ'], w['NYH']['mean'],
                 w['NYL']['mean'], w['Nsig']['mean'], w['R_Y'], w['R_Sigma'], w['R_conservative'], w['R_central']))
    print(json.dumps(checks, indent=1)); print(json.dumps(width_compare, indent=1))


if __name__ == '__main__':
    main(*sys.argv[1:5])
