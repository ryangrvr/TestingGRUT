#!/usr/bin/env python3
"""lf_chain step 4 (Claude): windows, representative points, D / s_* floor / r, the KFlo/KFhi brackets, Q3 (B, D, E, F),
H1 and Q4 (a, b, c) on the LF closure.  Reads results/fixpoints_LF.json, results/scan_*.json, results/pr_edge.json,
results/pr_phase_soft.json (pr_principle's edge.py / phase_soft.py, run verbatim from this track) and the release
timelines in scratch.  Writes results/windows.json.
Definitions (fixed before reading this output; see NOTES §2): M2-Q = self-consistent Q-matter root at e = 1e-6;
M2-X proper = self-consistent free-X root with 1e-6 < e_crit <= 4e-6 (BRIEF_E2E "e_crit ~ 2e-6"; e2e_chain skeptic's
window range 1.3-3.8e-6); M2-X loose = self-consistent free-X root with any e_crit <= 1e-3 (e2e_chain's operational
window).  D = -(3/2) dln H_req/dln c0 at a fixed node (the e2e_chain skeptic's finite-difference rule, Omega ~ H_req^-3/2),
same branch on both sides; s_* floor = D H_i / (2 pi sqrt(beta A_s/(1 - beta))), beta = 0.038, A_s = 2.1e-9.
Labels: [C] at nodes, [E] between nodes, [A] where an [A] node (nearest-H lattice block) is used.
"""
import json, math, glob
import numpy as np
import lfc_common as L
import fp_lf as F

E = L.E
cm = E.cm
SQB = F.SQB
BETA, AS = 0.038, 2.1e-9
SCR = L.SCR / 'lf_chain'


def runs(tag):
    return L.load(L.RES / ('scan_%s.json' % tag), {}).get('runs', {})


def keyf(k):
    c, d = k.split(':'); return float(c), float(d)


def D_at(tag, c0, d, branch):
    R = runs(tag)
    pts = []
    for k, o in R.items():
        c, dd = keyf(k)
        if dd != d or k.startswith('1.000'):
            continue
        b = o['by_e'][F.CL]['e1e-6']
        if (b['S_per_Q'] > SQB) == (branch == 'free X'):
            pts.append((c, b['H_i_req_GeV']))
    pts.sort()
    cs = [p[0] for p in pts]
    if c0 not in cs:
        return None
    i = cs.index(c0)
    out = {}
    for j, lab in ((i - 1, 'left'), (i + 1, 'right')):
        if 0 <= j < len(pts) and abs(pts[j][0] - c0) <= 0.051:
            out[lab] = -1.5 * math.log(pts[j][1] / pts[i][1]) / math.log(pts[j][0] / c0)
            out[lab + '_c0'] = pts[j][0]
    vals = [v for kk, v in out.items() if not kk.endswith('_c0')]
    out['D'] = float(np.mean(vals)) if vals else None
    return out


def floor(D, H):
    return D * H / (2 * math.pi * math.sqrt(BETA * AS / (1 - BETA)))


def rep(fx, k, branch, tag_for_D):
    r = [a for a in fx[k]['all_nodes'] if a['self_consistent'] and a['branch'] == branch]
    if not r:
        return None
    r = r[0]
    c0, d = keyf(k)
    Dd = D_at(tag_for_D, c0, d, branch)
    H = r['H_star']
    u = E.mc.MU + c0 * complex(math.cos(math.radians(d)), math.sin(math.radians(d))) * (10j - E.mc.MU)
    out = dict(c0=c0, delta_deg=d, mu_eff=[u.real, u.imag], arg_mu_eff_deg=math.degrees(math.atan2(u.imag, u.real)), H_i_GeV=H,
               w=r['w'], pair=r['pair'], uses_A_node=r['uses_A_node'], K=r['K'], SQ_1e6=r['SQ'], SQ_3e6=r['SQ3'], e_crit=r['e_crit'],
               M=r['M'], eps=r['eps'], N_YL=r['NYL'], r=cm.r_of(H), D=Dd, D_node=tag_for_D,
               s_star_floor_GeV=floor(Dd['D'], H) if Dd and Dd['D'] else None,
               primary_rule=fx[k].get('primary', {}).get(branch), primary_linearH=fx[k].get('primary_linearH', {}).get(branch),
               e2e_chain_rule=fx[k].get('e2e_chain_rule'))
    return out


def brackets():
    """KFlo / KFhi nodes at 1.3e11: root with the KFmid (8.4e10, 1.3e11) slopes at the same key [E], as the e2e_chain skeptic
    did for GMg03."""
    out = {}
    lo, mid = runs('hLF84_mid'), runs('hLF_mid')
    for tag in ('hLF_lo', 'hLF_hi'):
        for k, o in runs(tag).items():
            if k.startswith('1.000') or k not in lo or k not in mid:
                continue
            n, a, b_ = F.nd(o), F.nd(lo[k]), F.nd(mid[k])
            dl = math.log(b_['H'] / a['H'])
            rec = dict(node=dict(SQ=n['SQ'], R=n['R'], ecrit=n['ecrit'], K=n['K']))
            for br in ('free X', 'Q-matter'):
                sf = (F.fb(b_, br) - F.fb(a, br)) / dl
                ssq = (b_['SQ'] - a['SQ']) / dl
                x = -F.fb(n, br) / sf
                sq = n['SQ'] + ssq * x
                rec[br] = dict(H_star=n['H'] * math.exp(x), SQ=sq, self_consistent=bool((sq > SQB) == (br == 'free X')))
            out['%s %s' % (tag, k)] = rec
    return out


def t_release(tag, frac=0.9):
    """time at which |N_X| first reaches frac of its final value (release timeline of the engine run) [C]."""
    fs = glob.glob(str(SCR / 'sp*' / 'engine' / 'release' / ('release_%s.json' % tag)))
    if not fs:
        return None
    tl = json.load(open(fs[0]))['dense']['timeline']
    t = np.array([r['tau'] for r in tl]); nx = np.abs(np.array([r['NX'] for r in tl], float))
    fin = nx[-1]
    j = int(np.argmax(nx >= frac * fin))
    if j == 0:
        return None
    return float(t[j - 1] + (frac * fin - nx[j - 1]) * (t[j] - t[j - 1]) / (nx[j] - nx[j - 1]))


def cross(fx, k1, k2, branch, thr=SQB, field='SQ'):
    """c0 at which the interpolated root S/Q (or ln e_crit) of `branch` crosses thr between keys k1 < k2 (same delta) [E]."""
    g = lambda k: [a for a in fx[k]['all_nodes'] if a['branch'] == branch][0]
    a, b = g(k1), g(k2)
    va, vb = a[field], b[field]
    if field == 'e_crit':
        va, vb, thr = math.log(va), math.log(vb), math.log(thr)
    c1, c2 = keyf(k1)[0], keyf(k2)[0]
    return c1 + (c2 - c1) * (va - thr) / (va - vb)


def summary(fx, W):
    S = {}
    S['bistable_band_delta0'] = [cross(fx, '1.374:+0', '1.378:+0', 'Q-matter'), cross(fx, '1.374:+0', '1.378:+0', 'free X')]
    S['M2X_proper_pocket_delta-4'] = [cross(fx, '1.080:-4', '1.090:-4', 'free X', 4e-6, 'e_crit'), cross(fx, '1.100:-4', '1.120:-4', 'free X', 4e-6, 'e_crit')]
    soft = {'-4': '1.080:-4', '+0': '1.130:+0', '+4': '1.170:+4'}
    h1 = {}
    for d, k in soft.items():
        r = [a for a in fx[k]['all_nodes'] if a['self_consistent']]
        ec = r[0]['e_crit'] if r else None
        h1[d] = dict(key=k, root=r[0] if r else None, class_all=fx[k]['class_all'], class_primary=fx[k]['class_primary'],
                     strict_M2X=bool(isinstance(ec, float) and 1e-6 < ec <= 4e-6), loose_M2X=bool(isinstance(ec, float) and ec <= 1e-3),
                     M2Q=any(a['branch'] == 'Q-matter' for a in r), root_on_computed_lattice_range=bool(r and 8.4e10 <= r[0]['H_star'] <= 1.3e11))
    S['H1'] = h1
    S['H1_verdict_strict'] = 'FAILS' if not any(v['strict_M2X'] or v['M2Q'] for v in h1.values()) else 'SUCCEEDS'
    S['H1_verdict_loose'] = '%d of 3 points in a free-X-with-e_crit region, none on the computed lattice H range' % sum(v['loose_M2X'] for v in h1.values())
    return S


def main():
    fx = L.load(L.RES / 'fixpoints_LF.json')
    W = dict(definitions=__doc__, closure='LF (lf_fermion KFmid calibrated local blocking with capacity kappa_j(t))')
    W['grid'] = {k: dict(class_all=v['class_all'], class_primary=v['class_primary'], n_roots_all=v['n_roots_all'],
                         roots=[dict(branch=a['branch'], H_star=a['H_star'], w=a['w'], SQ=a['SQ'], e_crit=a['e_crit'], K=a['K'],
                                     uses_A_node=a['uses_A_node'], extrapolated=a['extrapolated'])
                                for a in v['all_nodes'] if a['self_consistent']]) for k, v in fx.items()}
    W['M2Q_rep_c1.45'] = rep(fx, '1.450:+0', 'Q-matter', 'hLF84_mid')
    W['M2Q_rep_c1.40'] = rep(fx, '1.400:+0', 'Q-matter', 'hLF84_mid')
    W['M2X_proper_rep_c1.37'] = rep(fx, '1.370:+0', 'free X', 'hLF84_mid')
    W['M2X_proper_c1.40_d+4'] = rep(fx, '1.400:+4', 'free X', 'hLF84_mid')
    W['old_points_on_LF'] = {'GM M2-X (1.55, 0)': fx.get('1.550:+0', {}).get('class_all'), 'P0 M2-Q (1.45, 0)': fx.get('1.450:+0', {}).get('class_all')}
    W['brackets_KFlo_KFhi'] = brackets()
    # Q3 -----------------------------------------------------------------------------------------------------------
    ps = L.load(L.RES / 'phase_soft.json', {}); ed = L.load(L.RES / 'edge.json', {})
    q3 = {}
    q3['B'] = dict(delta_max_torque=ps.get('torque'), evaluated={k: fx[k]['class_all'] for k in fx if keyf(k)[1] <= -10})
    epsl = sorted((keyf(k)[0], o['eps']) for k, o in runs('hLF_mid').items() if keyf(k)[1] == 0)
    q3['D'] = dict(eps_vs_c0_delta0_1p3e11=epsl, monotonic=bool(all(b[1] > a[1] for a, b in zip(epsl, epsl[1:]))))
    edge = {}
    for dd, rows in ed.items():
        lam = [r['eig'][1] for r in rows]
        edge[dd] = dict(min_lightest_nongauge=min(lam), c0_at_min=rows[int(np.argmin(lam))]['c0'], any_negative=bool(min(lam) < 0),
                        c0_range=[rows[0]['c0'], rows[-1]['c0']])
    q3['E'] = edge
    sq_all = [(F.nd(o)['SQ'], k) for k, o in runs('hLF_mid').items()]
    q3['F'] = dict(max_SQ_1p3e11=max(sq_all), min_SQ_1p3e11=min(sq_all), note='branch-defining by construction')
    W['Q3'] = q3
    # H1 ------------------------------------------------------------------------------------------------------------
    soft = ps.get('soft', {})
    W['H1'] = dict(c0_soft=soft, points={k: dict(class_all=fx[k]['class_all'], class_primary=fx[k]['class_primary'],
                                                  roots=W['grid'][k]['roots']) for k in fx if keyf(k)[0] < 1.2})
    # Q4 ------------------------------------------------------------------------------------------------------------
    trel = {}
    for k in runs('hLF_mid'):
        c, d = keyf(k)
        tag = 'hLF_mid' if k.startswith('1.000') else 'hLF_mid_C%.3f_d%+.0f' % (c, d)
        trel[k] = t_release(tag)
    head = json.load(open(L.LF_HIST))['head']
    W['Q4'] = dict(c_release_synchrony=dict(t_release90_1p3e11=trel, t_coh_end_KFmid=head.get('t_fcoh_lt_0p1'),
                                            note='t_coh,end = first t with f_sig + f_chi < 0.1 (lf_fermion head); never for KFhi and the 8.4e10 run'))
    W['summary'] = summary(fx, W)
    L.dump(L.RES / 'windows.json', W)
    print('summary', json.dumps(W['summary']['bistable_band_delta0']), W['summary']['M2X_proper_pocket_delta-4'], W['summary']['H1_verdict_strict'], '|', W['summary']['H1_verdict_loose'])
    for k in ('M2Q_rep_c1.45', 'M2Q_rep_c1.40', 'M2X_proper_rep_c1.37', 'M2X_proper_c1.40_d+4'):
        r = W[k]
        if r:
            print(k, 'H %.4g w %.2f K %.2f SQ %.3f/%.3f ec %s M %.4f r %.3g D %s floor %.3g' % (r['H_i_GeV'], r['w'], r['K'], r['SQ_1e6'], r['SQ_3e6'],
                  r['e_crit'], r['M'], r['r'], r['D'], r['s_star_floor_GeV'] or 0), 'linH', r['primary_linearH'] and '%.4g' % r['primary_linearH']['H_star'])
    print('brackets', {k: (round(v['node']['SQ'], 3), {b: ('%.3g' % v[b]['H_star'], round(v[b]['SQ'], 3), v[b]['self_consistent']) for b in ('free X', 'Q-matter')}) for k, v in W['brackets_KFlo_KFhi'].items()})
    print('Q3 D', q3['D']['monotonic'], 'E', edge, 'F', q3['F'])
    print('t_release', {k: (round(v, 3) if v else v) for k, v in sorted(trel.items())}, 't_coh', head.get('t_fcoh_lt_0p1'))


if __name__ == '__main__':
    main()
