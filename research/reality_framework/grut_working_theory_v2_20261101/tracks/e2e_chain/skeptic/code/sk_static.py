#!/usr/bin/env python3
"""e2e_chain_skeptic (Claude): static re-computations from e2e_chain / e2e_relic outputs (read-only) and our engine runs.

  A. mu_eff at the slice from the c0 definition mu_eff = mu + c0 e^{i delta}(10i - mu), mu = 7.159456 H_i, vs the text.
  B. The final histories (e2e_chain scratch, read-only): K = a^2 H at the end, convergence, the m_s^2 and R jumps at t_sw.
  C. N_* and r recomputed with our own entropy bookkeeping (not mf_hybrid's Nstar_from_K); s_* floor and beta recomputed.
  D. Branch-aware self-consistency.  e2e_chain's fixpoint.py interpolates f = ln(H_req/H_i) between the two nodes with
     each node's OWN branch (M ~ 1 for free X, ~0.80 for Q-matter) and then classifies the root by the interpolated S/Q.
     When the nodes sit on different branches, f mixes two relic functions.  Here: for each branch b in {free X,
     Q-matter} f_b(node) = ln(H_req,node (M_b/M_node)^(-2/3) / H_node), with M_b evaluated at the node's S/Q
     (free X: e2e_relic p_selfconsistent; Q-matter: e2e_relic relic(), calibrated, extrapolated to S/Q <= 3.4716);
     the root of f_b (linear in ln H, as e2e_chain) is self-consistent iff the interpolated S/Q there is on branch b.
     Also: the same with interpolation linear in H (interpolation-rule sensitivity).
     GMg03 (fragment spread) has only the 1.3e11 node: its root uses the plain-GM 8.4e10/1.3e11 slopes of f and S/Q at
     the same (c0, delta) [E].
  E. Relic arithmetic: object mass/radius/flux/impulse, the SHM mean speed (analytic), a flat-detector geometry factor,
     free-X cross sections (independent amplitude bookkeeping).
  F. e2e_chain's Starobinsky-Yokoyama numbers re-derived (for the record; their premise is tested in sk_natural.py).
Labels: [C] given the model choices; [E] where flagged.
"""
import json, math, glob
import numpy as np
from scipy.special import erf
import sk_common as S

E = S.E
import pmap as PM                                            # e2e_relic, read-only
import relic as RL                                           # e2e_relic, read-only

SQB = 3.4716191256166837
MU = E.mc.MU
CL = 'cL0_vanish_at_MY'


def key(c0, d):
    return '%.3f:%+.1f' % (round(float(c0), 3), round(float(d), 1))


def nodes():
    """(hist, key) -> node dict from e2e_chain scans and our sk_engine runs."""
    out = {}
    for f in glob.glob(str(S.E2C_RES / 'scan_*.json')):
        tag = f.split('scan_')[1][:-5]
        for k, o in json.load(open(f))['runs'].items():
            b = o['by_e'][CL]
            out[(tag, key(o['c0'], o['delta_deg']))] = dict(src='e2e_chain', H=o['H_i_GeV'], K=o['K'], SQ=b['e1e-6']['S_per_Q'],
                                                           SQ3=b['e3e-6']['S_per_Q'], M=b['e1e-6']['M'], R=b['e1e-6']['H_req_over_H_i'],
                                                           ecrit=b['e_crit_bulk'], NYL=o['N_YL'])
    mine = S.load(S.RES / 'sk_engine.json', {})
    alias = {'rHFa_GM': 'hHFa_GM', 'rH84_GM': 'hH84_GM', 'rhHFa_P0': 'hHFa_P0', 'rhH84_P0': 'hH84_P0', 'rHFa_GMg03': 'hHFa_GMg03',
             'rhHFa_GM_anh': 'hHFa_GM_anh', 'rhH84_GM_anh': 'hH84_GM_anh'}
    for base, blk in mine.items():
        tag = alias.get(base, base)
        for k, o in blk['runs'].items():
            a = o['e1e-6']
            kk = key(o['c0'], o['delta_deg'])
            if (tag, kk) in out and tag in ('hHFa_GM', 'hH84_GM', 'hHFa_P0', 'hH84_P0', 'hHFa_GMg03'):
                out[(tag, kk)]['regression_SQ_skeptic'] = a['S_per_Q']
                continue
            out[(tag, kk)] = dict(src='skeptic', H=o['H_i_GeV'], K=o['K'], SQ=a['S_per_Q'], SQ3=o['e3e-6']['S_per_Q'], M=a['M'],
                                  R=a['H_req_over_H_i'], ecrit=o['e_crit_bulk'], NYL=o['N_YL'])
    return out


def M_branch(b, sq):
    if b == 'free X':
        p, _ = PM.p_selfconsistent(max(sq, SQB + 1e-6))
        return 1 + 2 * p
    return RL.relic(min(sq, SQB - 1e-6))['M_central']


def solve_branch(lo, hi, b, rule='lnH'):
    x = (lambda H: math.log(H)) if rule == 'lnH' else (lambda H: H)
    f = lambda n: math.log(n['R'] * (M_branch(b, n['SQ']) / n['M']) ** (-2.0 / 3.0))
    x0, x1, f0, f1 = x(lo['H']), x(hi['H']), f(lo), f(hi)
    w = -f0 / (f1 - f0)
    sq = lo['SQ'] + w * (hi['SQ'] - lo['SQ'])
    Hs = math.exp(x0 + w * (x1 - x0)) if rule == 'lnH' else x0 + w * (x1 - x0)
    on = (sq > SQB) == (b == 'free X')
    ec = None
    if lo['ecrit'] and hi['ecrit']:
        ec = math.exp(math.log(lo['ecrit']) + w * (math.log(hi['ecrit']) - math.log(lo['ecrit'])))
    return dict(branch=b, H_star=Hs, w=w, extrapolated=not (0 <= w <= 1), SQ=sq, self_consistent=bool(on), e_crit=ec,
                K=lo['K'] + w * (hi['K'] - lo['K']), M=M_branch(b, sq))


def fix_all(N):
    res = {}
    for clo, (tlo, thi) in {'GM': ('hH84_GM', 'hHFa_GM'), 'P0': ('hH84_P0', 'hHFa_P0')}.items():
        keys = sorted({k for (t, k) in N if t == tlo} & {k for (t, k) in N if t == thi})
        for k in keys:
            if k.startswith('1.000'):
                continue
            lo, hi = N[(tlo, k)], N[(thi, k)]
            e2c_w = None
            fe = lambda n: math.log(n['R'])
            w = -fe(lo) / (fe(hi) - fe(lo))
            sq_e2c = lo['SQ'] + w * (hi['SQ'] - lo['SQ'])
            M_e2c = lo['M'] + w * (hi['M'] - lo['M'])
            r = dict(nodes=dict(lo=lo, hi=hi), e2e_chain_rule=dict(w=w, SQ=sq_e2c, M_interp=M_e2c,
                                                                  H_star=math.exp(math.log(lo['H']) + w * (math.log(hi['H']) - math.log(lo['H'])))))
            for b in ('free X', 'Q-matter'):
                r[b] = solve_branch(lo, hi, b)
                r[b + ' (linear in H)'] = solve_branch(lo, hi, b, rule='H')
            r['n_selfconsistent_roots'] = sum(r[b]['self_consistent'] for b in ('free X', 'Q-matter'))
            res['%s %s' % (clo, k)] = r
    # GMg03: single node, plain-GM slopes at the same key [E]
    for (t, k), n in N.items():
        if t != 'hHFa_GMg03' or k.startswith('1.000'):
            continue
        lo, hi = N.get(('hH84_GM', k)), N.get(('hHFa_GM', k))
        if not (lo and hi):
            continue
        dl = math.log(hi['H'] / lo['H'])
        out = dict(node=n, slopes_from='plain GM at the same (c0, delta)')
        for b in ('free X', 'Q-matter'):
            fb = lambda m: math.log(m['R'] * (M_branch(b, m['SQ']) / m['M']) ** (-2.0 / 3.0))
            sf = (fb(hi) - fb(lo)) / dl
            ssq = (hi['SQ'] - lo['SQ']) / dl
            x = -fb(n) / sf                                   # ln(H*/H_node)
            sq = n['SQ'] + ssq * x
            out[b] = dict(H_star=n['H'] * math.exp(x), dlnH=x, SQ=sq, self_consistent=bool((sq > SQB) == (b == 'free X')),
                          slope_f=sf, slope_SQ=ssq, M=M_branch(b, sq), e_crit_node=n['ecrit'])
        res['GMg03 %s' % k] = out
    return res


def hist_checks():
    out = {}
    for f in sorted(glob.glob(str(S.E2C_SP / 'hist' / 'bench_h*.json'))):
        d = json.load(open(f)); h = d['history']
        g = lambda k: np.array([np.nan if v is None else float(v) for v in h[k]])
        t, a, H, R, ms, fr = g('t'), g('a'), g('H'), g('R'), g('m_s'), g('f_rad')
        K = a * a * H
        j = int(np.searchsorted(t, 3.0 + 1e-9))
        lat = (t > 2.0) & (t <= 3.0)
        out[d['name']] = dict(K_end=float(K[-1]), K_stored=d['reheating']['K'], t_end=float(t[-1]),
                              K_rel_change_last_decade=float(K[-1] / np.interp(t[-1] / 10, t, K) - 1), f_rad_end=float(fr[-1]),
                              ms2_last_lattice=float(ms[j - 1] ** 2), ms2_first_kinetic=float(ms[j] ** 2),
                              ms2_jump=float(ms[j] ** 2 / ms[j - 1] ** 2 - 1),
                              ms2_lattice_mean_t2_3=float(np.mean(ms[lat] ** 2)), ms2_lattice_rel_std_t2_3=float(np.std(ms[lat] ** 2) / np.mean(ms[lat] ** 2)),
                              R_jump=float(R[j] / R[j - 1] - 1), H_jump=float(H[j] / H[j - 1] - 1))
    return out


def nstar(H, K, g=106.75, gs0=3.931, T0=2.7255 * 8.617333e-14, kstar=0.05):
    """N_* = ln(a_slice H_i / k_*), a_slice from entropy conservation; rho_r a^4 = 3 K^2 M_P^2 H_i^2 (a = 1 at slice)."""
    MP = E.MP_GEV
    Ta = (90 * K ** 2 * MP ** 2 * H ** 2 / (math.pi ** 2 * g)) ** 0.25
    a_sl = T0 * gs0 ** (1 / 3) / (Ta * g ** (1 / 3))
    k = kstar * 0.1973269804e-15 / 3.0856775814913673e22       # GeV
    return math.log(a_sl * H / k)


def shm_mean_speed_analytic(v0=238.0, vE=250.5):
    """<|v|> of an untruncated Maxwellian (dispersion v0/sqrt2 per axis) boosted by vE [D]."""
    eta = vE / v0
    return v0 * (math.exp(-eta ** 2) / math.sqrt(math.pi) + (eta + 1 / (2 * eta)) * erf(eta))


def relic_checks(m2p):
    out = {}
    q = m2p['M2Q']['relic']; o = q['observables']
    Q = q['Q_mean']; MX = RL.MX
    m_GeV = 0.775866 * Q * MX
    m_kg = m_GeV * 1.78266192e-27
    R_fm = (3 * Q / (4 * math.pi * 449.31336788454803)) ** (1 / 3) * 0.1973269804 / MX
    fobj = 1 - q['f_v']['perfect'] / (0.775866 + 0.224135 * q['f_v']['perfect'])
    n_m3 = 0.3 * fobj / m_GeV * 1e6
    v_an = shm_mean_speed_analytic()
    flux = n_m3 * v_an * 1e3 * 3.15576e7
    I = 2 * 6.67430e-11 * m_kg * 1e-6 / (1e-3 * v_an * 1e3)
    out['M2Q_objects'] = dict(m_obj_ug=m_kg * 1e9, track_ug=o['object_mass_ug'], R_nm=R_fm * 1e-6, track_R_nm=o['radius_nm'],
                              n_m3=n_m3, track_n_m3=o['n_local_per_m3'], mean_speed_analytic_untruncated=v_an,
                              track_mean_speed_MC_truncated=o['mean_speed_km_s'], flux_sphere_or_normal_plane_m2_yr=flux,
                              flux_isotropic_flat_plane_one_side_m2_yr=flux / 2,
                              track_flux=o['flux_per_m2_per_yr'], impulse=I, track_impulse=o['impulse_on_sensor']['impulse_kg_m_s'],
                              M_from_fv=0.775866 + 0.224135 * 0.5 * (q['f_v']['perfect'] + q['f_v']['beliaev']), track_M=q['M_over_mX'],
                              mass_over_MPl=m_kg / 2.176434e-8)
    # free-X tree amplitudes, own bookkeeping: lam (X*X)^2 contact -4 lam; phi exchange with vertex 2 g^2 v = 2 g m_X
    lam, g = 1e-3, 0.1
    y2 = (2 * g) ** 2                               # (2 g^2 v)^2 / m_X^2 with g v = m_X
    mphi2 = 6.0
    A_XX = -4 * lam - y2 * (1 / (0 - mphi2) + 1 / (0 - mphi2))          # t = u = 0 at threshold
    A_XXb = -4 * lam - y2 * (1 / (4 - mphi2) + 1 / (0 - mphi2))         # s = 4 m^2, t = 0
    conv = 0.3893794e-27
    sXX = A_XX ** 2 / (32 * math.pi * 4) / MX ** 2 * conv
    sXXb = A_XXb ** 2 / (16 * math.pi * 4) / MX ** 2 * conv
    mg = MX * 1.78266192e-24
    out['freeX_sigma'] = dict(sigma_over_m_XX=sXX / mg, sigma_over_m_XXbar=sXXb / mg,
                              track=m2p['M2X']['relic']['observables']['self_interaction'])
    return out


def main():
    m2p = S.load(S.E2C_RES / 'm2_points.json')
    res = {}
    # A
    A = {}
    for nm, c0 in (('M2X', 1.55), ('M2Q', 1.45)):
        u = MU + c0 * (10j - MU)
        A[nm] = dict(mu_eff=[u.real, u.imag], abs=abs(u), arg_deg=math.degrees(math.atan2(u.imag, u.real)),
                     json=m2p[nm]['holder']['mu_eff_slice'])
    A['NOTES_text'] = dict(M2X='-2.38 + 15.81i (|mu| 15.99, 104.3 deg)', M2Q='-1.66 + 14.85i (102.5 deg)',
                           text_components_abs_arg=dict(M2X=[abs(-2.38 + 15.81j), math.degrees(math.atan2(15.81, -2.38))],
                                                        M2Q=[abs(-1.66 + 14.85j), math.degrees(math.atan2(14.85, -1.66))]))
    res['A_mu_eff'] = A
    res['B_histories'] = hist_checks()
    # C
    C = {}
    for nm in ('M2X', 'M2Q'):
        sc = m2p[nm]['self_consistent']; H, K = sc['H_i_GeV'], sc['K']
        D = m2p[nm]['isocurvature']['D']
        smin = abs(D) * H / (2 * math.pi * math.sqrt(0.038 * 2.1e-9 / 0.962))
        C[nm] = dict(Nstar_own=nstar(H, K), Nstar_track=m2p[nm]['cmb']['Nstar'], r_own=2 * H ** 2 / (math.pi ** 2 * 2.1e-9 * E.MP_GEV ** 2),
                     r_track=m2p[nm]['cmb']['r'], s_floor_own=smin, s_floor_track=m2p[nm]['isocurvature']['s_star_min_GeV'],
                     dNstar_dlnK=0.5)
    # D cross-check by finite differences over the scanned c0 grid at 1.3e11 GM (free-X, Omega ~ H_req^-3/2 at fixed H)
    sc = json.load(open(S.E2C_RES / 'scan_hHFa_GM.json'))['runs']
    Rr = {k: sc[k]['by_e'][CL]['e1e-6']['H_req_over_H_i'] for k in ('1.500:+0', '1.550:+0', '1.600:+0', '1.534:+0', '1.565:+0')}
    C['D_fd_GM'] = {'1.50-1.55': -1.5 * math.log(Rr['1.550:+0'] / Rr['1.500:+0']) / math.log(1.55 / 1.50),
                    '1.55-1.60': -1.5 * math.log(Rr['1.600:+0'] / Rr['1.550:+0']) / math.log(1.60 / 1.55),
                    '1.534-1.565': -1.5 * math.log(Rr['1.565:+0'] / Rr['1.534:+0']) / math.log(1.565 / 1.534)}
    res['C_cmb_iso'] = C
    # D
    N = nodes()
    res['D_fixpoints'] = fix_all(N)
    res['regressions'] = {'%s %s' % k: dict(track=v['SQ'], skeptic=v['regression_SQ_skeptic']) for k, v in N.items() if 'regression_SQ_skeptic' in v}
    # E
    res['E_relic'] = relic_checks(m2p)
    # F
    F = {}
    for m2 in (3e-5, 4.5e-7, 2.5e-7, 4.5e-8):
        srms = math.sqrt(3 / (8 * math.pi ** 2 * m2))
        F['%.2g' % m2] = dict(s_rms_over_H=srms, N_relax=3 / m2, M2X_ratio=m2p['M2X']['isocurvature']['s_star_min_GeV'] / m2p['M2X']['self_consistent']['H_i_GeV'] / srms)
    res['F_SY'] = F
    S.dump('sk_static.json', res)
    # print summary
    print(json.dumps(A, indent=0)[:900])
    for k, v in res['B_histories'].items():
        print(k, {q: round(x, 5) if isinstance(x, float) else x for q, x in v.items()})
    print(json.dumps(C, indent=0, default=float))
    for k, r in res['D_fixpoints'].items():
        if k.startswith('GMg03'):
            print('%-16s node SQ %.3f R %.3f | fx: H* %.3e SQ %.3f sc %s | Q: H* %.3e SQ %.3f sc %s' % (
                k, r['node']['SQ'], r['node']['R'], r['free X']['H_star'], r['free X']['SQ'], r['free X']['self_consistent'],
                r['Q-matter']['H_star'], r['Q-matter']['SQ'], r['Q-matter']['self_consistent']))
            continue
        e = r['e2e_chain_rule']
        print('%-14s e2c: w %.2f SQ %.3f M %.3f | fx: H* %.3e w %.2f SQ %.3f sc %-5s | Q: H* %.3e w %.2f SQ %.3f sc %-5s | n %d | linH fx %.3e Q %.3e' % (
            k, e['w'], e['SQ'], e['M_interp'], r['free X']['H_star'], r['free X']['w'], r['free X']['SQ'], r['free X']['self_consistent'],
            r['Q-matter']['H_star'], r['Q-matter']['w'], r['Q-matter']['SQ'], r['Q-matter']['self_consistent'], r['n_selfconsistent_roots'],
            r['free X (linear in H)']['H_star'], r['Q-matter (linear in H)']['H_star']))
    print(json.dumps(res['E_relic'], indent=0, default=float))
    print(json.dumps(res['regressions']))


if __name__ == '__main__':
    main()
