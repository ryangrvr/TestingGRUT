#!/usr/bin/env python3
"""hw_gates gate 1, piece B (Claude, own code): release-induced Goldstone angle at the spot point.

Mechanism (sp_vector, [D]): in each homogeneous patch N_X + 3N_Y + 6N_Sigma = 0, so the Sigma condensate carries
-N_X/6 and its phase winds, d theta_Sigma/da = -N_Sigma/(a^2 hhat^2) (vacuum |Sigma|); patch-to-patch differences of the
late Sigma phase freeze into the longitudinal mode while m_A < H (until a_* = 1/sqrt(m_vac)); theta_G -> theta_Sigma(inf)/6.
Conversion rule (sp_vector): N_A = int dlnk [dN_A^{P0}/dlnk] P_thetaG(k)/alpha_BD^2, alpha_BD = H_i/(2 pi F_vac).
Here: separate-universe linear responses of theta_Sigma(inf), theta_G, N_Sigma, eps to the source phase delta and to
the prepared moduli |X|, |Sigma|, |Y| (own hwgrel releases, relaxed preparation, a_end = 60, tail added);
natural amplitudes at k_* for the heavy inflationary modes (m > 3H/2: (H/2pi) e^{-3 N_*/2} per ln k, order of
magnitude), B1 = response x amplitude against the criterion set by the P0 spectrum and the remaining allowance;
B2 = Hubble-scale patch-noise model of sp_vector with the spot-point responses.
Usage: python3 -B gate1b_pieceB.py <gate1_vector.json> <out.json>
"""
import json, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import spotmodel as S

P = S.params(); D0 = np.pi / 2; AEND = 60.0


def job(args):
    tag, z0 = args
    r = S.run(np.array(z0, complex), P, a_end=AEND, rtol=1e-10, atol=1e-13, nsamp=30, win0=30.0)
    w = r['window']; s = r['samples'][-1]
    tail = -s['NS'] / (P['h']**2 * AEND)
    return dict(tag=tag, eps=w['eps'], NX=w['NX'], NS=w['NS'], thS_end=s['thS'], thS_inf=s['thS'] + tail, thG_end=s['thG'],
                thG_inf=(s['thS'] + tail) / 6, thX_end=s['thX'], hist=[(q['a'], q['thS'], q['NS']) for q in r['samples']])


def main(g1file, out):
    t0 = time.time()
    g1 = json.load(open(g1file))
    zb, info = S.prep('relaxed', D0, P)
    jobs = [('base', zb)]
    for dd in (1e-3, 1e-2):
        for sg in (+1, -1):
            z, _ = S.prep('relaxed', D0 + sg * dd, P); jobs.append((f'delta{sg * dd:+.0e}', z))
    for comp, name in ((0, 'X'), (2, 'S')):
        for sg in (+1, -1):
            z = zb.copy(); z[comp] = z[comp] * (1 + sg * 1e-3); jobs.append((f'ln{name}{sg * 1e-3:+.0e}', z))
    for sg in (+1, -1):
        z = zb.copy(); z[1] = z[1] * (1 + sg * 1e-3); jobs.append((f'lnY{sg * 1e-3:+.0e}', z))
    with Pool(4) as pool:
        rows = pool.map(job, jobs)
    R = {r['tag']: r for r in rows}
    base = R['base']
    resp = {}
    for key, pl, mi, dd in (('delta_1e-3', 'delta+1e-03', 'delta-1e-03', 1e-3), ('delta_1e-2', 'delta+1e-02', 'delta-1e-02', 1e-2),
                            ('lnX', 'lnX+1e-03', 'lnX-1e-03', 1e-3), ('lnS', 'lnS+1e-03', 'lnS-1e-03', 1e-3), ('lnY', 'lnY+1e-03', 'lnY-1e-03', 1e-3)):
        a, b = R[pl], R[mi]
        resp[key] = dict(dthS_inf=(a['thS_inf'] - b['thS_inf']) / (2 * dd), dthG_inf=(a['thG_inf'] - b['thG_inf']) / (2 * dd),
                         dNS=(a['NS'] - b['NS']) / (2 * dd), deps=(a['eps'] - b['eps']) / (2 * dd), dlnNX=(np.log(abs(a['NX'])) - np.log(abs(b['NX']))) / (2 * dd))
    # natural inflationary amplitudes at k_* (heavy modes, order of magnitude) and at horizon exit
    e = P['e']; m_vac = 6 * e * P['h']; k_star = np.sqrt(m_vac); N_star = np.log(1 / k_star)
    xa, ya, sa = S.aligned(P); za = np.array([xa, -ya, -sa])
    sp = S.spectrum(za, P)
    light = sp['light_phase']; f = light['f_eff']
    heavy_decay = np.exp(-1.5 * N_star)
    amp = dict(delta_exit=1 / (2 * np.pi * f), lnS_exit=1 / (2 * np.pi * np.sqrt(2) * sa), lnX_exit=1 / (2 * np.pi * np.sqrt(2) * xa),
               lnY_exit=1 / (2 * np.pi * np.sqrt(2) * ya), decay_to_kstar=heavy_decay, N_star=N_star, k_star=k_star,
               note='all inflationary modes at this point have m > 3H/2 (lightest 2.53 H): superhorizon amplitudes fall as e^{-3N/2} (times oscillation); e^{-3N_*/2} used as the order of magnitude')
    # criterion from the P0 spectrum and remaining allowance
    P0 = next(s for s in g1['spectra'] if s['profile'] == 'P0' and s['e'] == 1e-7)
    P1 = next(s for s in g1['spectra'] if s['profile'] == 'P1' and s['e'] == 1e-7)
    row = next(r for r in g1['ledger_rows'] if r['profile'] == 'P1' and r['e'] == 1e-7)
    dep = row['deposit_per_A_in_EX']; RL = g1['point']['R_L']
    R_scalar = row['R_Y'] + row['R_Sigma'] + row['R_cascade']
    N_allow_total = (RL - R_scalar) / dep
    N_allow_B = N_allow_total - P1['N_A']
    alpha_BD = 1 / (2 * np.pi * 6 * P['h'])
    crit = alpha_BD * np.sqrt(max(N_allow_B, 0) / P0['N_A'])     # rms Goldstone angle per ln k at k ~ k_* allowed for B1
    B1 = {}
    for key, a_exit in (('delta_1e-3', amp['delta_exit']), ('lnX', amp['lnX_exit']), ('lnS', amp['lnS_exit']), ('lnY', amp['lnY_exit'])):
        dthG = abs(resp[key]['dthG_inf']) * a_exit * heavy_decay
        B1[key] = dict(dthetaG_kstar_rad_per_lnk=dthG, ratio_to_criterion=dthG / crit,
                       N_A_B1=P0['N_A'] * (dthG / alpha_BD)**2, pairs_per_net=P0['N_A'] * (dthG / alpha_BD)**2 / abs(base['NX']))
    # B2: Hubble-scale patch noise (sp_vector model): modes k in (1/40, 1) with amplitude (exit) k^{3/2} (heavy decay),
    # delta N_Sigma = (dN_S/dq) dq, accumulated Sigma phase until horizon entry a_k = 1/k from the nominal history,
    # counted as (hhat^2/2) a_k (delta theta_Sigma)^2 quanta per ln k
    a_h = np.array([h[0] for h in base['hist']]); th_h = np.array([h[1] for h in base['hist']]); NS_nom = base['NS']
    ks = np.geomspace(1 / 40, 1.0, 100)
    B2 = {}
    for key, a_exit in (('delta_1e-3', amp['delta_exit']), ('lnX', amp['lnX_exit']), ('lnS', amp['lnS_exit']), ('lnY', amp['lnY_exit'])):
        dN = abs(resp[key]['dNS'])
        integ = []
        for k in ks:
            a_k = 1 / k
            dq = a_exit * k**1.5
            th_s = th_h[0]; th_k = np.interp(np.log(a_k), np.log(a_h), th_h)
            dth = dN * dq * (th_k - th_s) / NS_nom
            integ.append(P['h']**2 / 2 * a_k * dth**2)
        NA = float(np.trapezoid(integ, np.log(ks)))
        B2[key] = dict(N_A_B2=NA, pairs_per_net=NA / abs(base['NX']), ratio_to_remaining_allowance=NA / max(N_allow_B, 1e-300))
    res = dict(method=__doc__, base=dict(eps=base['eps'], NX=base['NX'], NS=base['NS'], thS_release=base['hist'][0][1], thS_inf=base['thS_inf'],
                                        total_winding_rad=base['thS_inf'] - base['hist'][0][1], thG_inf=base['thG_inf'], thG_end=base['thG_end'],
                                        thX_end=base['thX_end']),
               responses=resp, amplitudes=amp, allowance=dict(deposit_per_A_in_EX=dep, R_L=RL, R_scalar=R_scalar, N_A_allowed_total=N_allow_total,
                                                               N_A_pieceA_P1=P1['N_A'], N_A_remaining_for_B=N_allow_B, alpha_BD=alpha_BD,
                                                               criterion_dthetaG_kstar_rad=crit, P0_N_A=P0['N_A']),
               B1=B1, B2=B2, runs=rows)
    Path(out).write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps({k: v for k, v in res.items() if k not in ('method', 'runs')}, indent=1, default=float))
    print('done', time.time() - t0)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
