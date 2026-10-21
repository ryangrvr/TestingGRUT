"""Aggregate the verification runs into results/v_charged_exit_summary.json."""
import json, os, glob
import numpy as np
from vlat import NBAR_C

HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, 'results')


def load(name):
    p = os.path.join(R, name + '.json')
    if not os.path.exists(p):
        return None
    return json.load(open(p))


def arr(run, key):
    return np.array([s.get(key, np.nan) for s in run['series']], float)


def fit_rate(tau, xi, lo, hi):
    m = (tau >= lo) & (tau <= hi) & np.isfinite(xi)
    if m.sum() < 3:
        return None
    c = np.polyfit(tau[m], xi[m], 1)
    return float(c[0])


def form_index(run, thr=0.5):
    fp, fm = arr(run, 'f_plus'), arr(run, 'f_minus')
    idx = np.where(fp + fm >= thr)[0]
    return int(idx[0]) if len(idx) else None


def run_metrics(run, kappa=None):
    tau = arr(run, 'tau'); a = arr(run, 'a')
    fp, fm = arr(run, 'f_plus'), arr(run, 'f_minus')
    Lb = arr(run, 'L_basin'); Lr = arr(run, 'L_raw'); Lbn = arr(run, 'L_basin_noexcl')
    out = dict(name=run['name'], L_box=run['L_box'], max_Q_rel_dev=run.get('max_Q_rel_dev'),
               wall_seconds=run['wall_seconds'])
    i0 = form_index(run)
    if i0 is None:
        out['formed'] = False
        out['f_plus_min'] = float(fp.min()); out['f_minus_max'] = float(fm.max())
        return out
    out['formed'] = True
    out['a_form'] = float(a[i0]); out['tau_form'] = float(tau[i0])
    if kappa:
        out['n_form_over_nc'] = kappa / a[i0] ** 3
    out['p_plus_form'] = float(fp[i0] / (fp[i0] + fm[i0]))
    i8 = form_index(run, 0.8)
    if i8 is not None:
        out['p_plus_at_vac0.8'] = float(fp[i8] / (fp[i8] + fm[i8])); out['a_vac0.8'] = float(a[i8])
    out['p_plus_end'] = float(fp[-1] / (fp[-1] + fm[-1]))
    out['f_plus_end'] = float(fp[-1]); out['f_minus_end'] = float(fm[-1])
    post = slice(i0, None)
    sp = arr(run, 'f_seed_plus')[post]; sm = arr(run, 'f_seed_minus')[post]
    out['walls_every_output_after_form'] = bool(np.all((Lb[post] > 0) & (sp > 0) & (sm > 0)))
    out['min_frac_minority_after_form'] = float(np.min(np.minimum(fp[post], fm[post])))
    out['maxcl_plus_end'] = run['series'][-1].get('maxcl_plus'); out['maxcl_minus_end'] = run['series'][-1].get('maxcl_minus')
    out['ncl_plus_end'] = run['series'][-1].get('ncl_plus'); out['ncl_minus_end'] = run['series'][-1].get('ncl_minus')
    xi = 1.0 / Lb
    out['xi_c_form'] = float(xi[i0]); out['xi_c_end'] = float(xi[-1]); out['xi_phys_end'] = float(xi[-1] * a[-1])
    out['xi_c_end_over_L'] = float(xi[-1] / run['L_box'])
    out['dyn_range_xi_c'] = float(xi[-1] / xi[i0])
    out['dyn_range_xi_phys'] = float(xi[-1] * a[-1] / (xi[i0] * a[i0]))
    t_lo = tau[i0] + 5
    boxsafe = run['L_box'] / 8
    ok = np.where(xi <= boxsafe)[0]
    t_hi = tau[ok[-1]] if len(ok) else tau[-1]
    out['rate_boxsafe'] = fit_rate(tau, xi, t_lo, t_hi)
    out['rate_boxsafe_window'] = [float(t_lo), float(t_hi)]
    out['rate_all_post'] = fit_rate(tau, xi, t_lo, tau[-1])
    tm = 0.5 * (t_lo + tau[-1])
    out['rate_first_half'] = fit_rate(tau, xi, t_lo, tm); out['rate_second_half'] = fit_rate(tau, xi, tm, tau[-1])
    out['Lraw_over_Lbasin_end'] = float(Lr[-1] / Lb[-1]); out['Lnoexcl_over_Lbasin_end'] = float(Lbn[-1] / Lb[-1])
    out['A_area_end'] = float(Lb[-1] * tau[-1] / 2)
    s = run['series'][-1]
    for k in ('qfrac_wallcomp', 'qfrac_within_2', 'qfrac_within_4', 'qfrac_within_6', 'qfrac_Qmatter',
              'qfrac_Qmatter_on_wallcomp', 'qfrac_vac', 'qfrac_core', 'area_within_2', 'area_within_4', 'area_within_6',
              'f_wallcomp'):
        if k in s:
            out[k + '_end'] = s[k]
    if 'qfrac_wallcomp' in s:
        q = arr(run, 'qfrac_wallcomp')[post]
        out['qfrac_wallcomp_range_post'] = [float(np.nanmin(q)), float(np.nanmax(q))]
    return out


def cont_metrics(run, box_xi=None):
    tau = arr(run, 'tau'); Lb = arr(run, 'L_basin'); a = arr(run, 'a')
    xi = 1.0 / Lb
    o = dict(name=run['name'], tau0=float(tau[0]), tau_end=float(tau[-1]), xi_c_start=float(xi[0]), xi_c_end=float(xi[-1]),
             rate_all=fit_rate(tau, xi, tau[0], tau[-1]))
    if box_xi:
        ok = np.where(xi <= box_xi)[0]
        o['rate_boxsafe'] = fit_rate(tau, xi, tau[0], tau[ok[-1]]) if len(ok) else None
        o['boxsafe_tau_end'] = float(tau[ok[-1]]) if len(ok) else None
    n = len(tau); third = n // 3
    o['rate_first_third'] = fit_rate(tau, xi, tau[0], tau[third]); o['rate_last_third'] = fit_rate(tau, xi, tau[n - 1 - third], tau[-1])
    o['walls_every_output'] = bool(np.all(Lb > 0))
    o['p_plus_end'] = float(arr(run, 'f_plus')[-1] / (arr(run, 'f_plus')[-1] + arr(run, 'f_minus')[-1]))
    o['Lraw_over_Lbasin_end'] = float(arr(run, 'L_raw')[-1] / Lb[-1])
    s = run['series'][-1]
    for k in ('qfrac_wallcomp', 'qfrac_within_4', 'qfrac_Qmatter'):
        if k in s:
            o[k + '_start'] = run['series'][0][k]; o[k + '_end'] = s[k]
    o['series_tau_xi'] = [[float(t), float(x)] for t, x in zip(tau, xi)]
    return o


summ = dict(runs={}, groups={}, decharge={}, controls={})
for p in sorted(glob.glob(os.path.join(R, '*.json'))):
    nm = os.path.basename(p)[:-5]
    if nm in ('checks', 'v_charged_exit_summary'):
        continue
    run = load(nm)
    if nm.startswith('A_') and 'to1.' not in nm or nm.startswith('D3_'):
        summ['runs'][nm] = run_metrics(run, kappa=run['config'].get('kappa'))
    elif nm.startswith('A_') and 'to1.' in nm:
        summ['runs'][nm] = run_metrics(run, kappa=run['config'].get('kappa'))
    elif nm.startswith('B_'):
        summ['runs'][nm] = run_metrics(run)
    elif nm.startswith('C1_'):
        S_mean = arr(run, 'S_mean'); fp = arr(run, 'f_plus'); fm = arr(run, 'f_minus')
        summ['controls'][nm] = dict(meta=run['meta'], min_f_plus=float(fp.min()), max_f_minus=float(fm.max()),
                                    min_S_mean=float(S_mean.min()), any_negative_vacuum_cell=bool(fm.max() > 0),
                                    max_frac_Qmatter=float(np.nanmax(arr(run, 'f_Qmatter'))),
                                    max_Q_rel_dev=run['max_Q_rel_dev'], a_end=float(arr(run, 'a')[-1]))

# groups
for g, pref in (('A_k2_t40_N1024', 'A_k2_t40_N1024_s'), ('B_t40_N1024', 'B_t40_N1024_s')):
    rs = [v for k, v in summ['runs'].items() if k.startswith(pref)]
    if rs:
        G = {}
        for key in ('a_form', 'p_plus_form', 'p_plus_at_vac0.8', 'p_plus_end', 'rate_boxsafe', 'rate_all_post', 'xi_c_end',
                    'qfrac_wallcomp_end', 'qfrac_within_4_end', 'qfrac_Qmatter_end', 'dyn_range_xi_phys', 'dyn_range_xi_c'):
            v = np.array([r[key] for r in rs if r.get(key) is not None], float)
            if len(v):
                G[key] = dict(mean=float(v.mean()), std=float(v.std(ddof=1)) if len(v) > 1 else None, vals=v.tolist())
        G['walls_every_output_all'] = all(r['walls_every_output_after_form'] for r in rs)
        summ['groups'][g] = G

for tag, src, cc, dc, L in (('t150_from1.6', 'A_k2_t150_N2048_s301_to1.6', 'CC_k2_t150_N2048_s301_from1.6',
                             'DC_k2_t150_N2048_s301_from1.6', 512.0),
                            ('t40_from1.8', 'A_k2_t40_N2048_s302_to1.8', 'CC_k2_t40_N2048_s302_from1.8',
                             'DC_k2_t40_N2048_s302_from1.8', 512.0),
                            ('frozen_t150_at1.6', 'A_k2_t150_N2048_s301_to1.6', 'F_k2_t150_N2048_s301_freeze1.6',
                             'F_dc_k2_t150_N2048_s301_freeze1.6', 512.0),
                            ('frozen_t150_at1.6_s303', 'A_k2_t150_N2048_s303_to1.6', 'F_k2_t150_N2048_s303_freeze1.6',
                             'F_dc_k2_t150_N2048_s303_freeze1.6', 512.0)):
    rc, rd = load(cc), load(dc)
    if rc and rd:
        mc, md = cont_metrics(rc, L / 8), cont_metrics(rd, L / 8)
        d = dict(charged=mc, decharged=md, ratio_all=mc['rate_all'] / md['rate_all'])
        if mc.get('rate_boxsafe') and md.get('rate_boxsafe'):
            d['ratio_boxsafe'] = mc['rate_boxsafe'] / md['rate_boxsafe']
        # common box-safe window: up to the earlier box-safe end of the two
        te = min(mc.get('boxsafe_tau_end') or mc['tau_end'], md.get('boxsafe_tau_end') or md['tau_end'])
        tc = np.array(mc['series_tau_xi']); td = np.array(md['series_tau_xi'])
        rc_c = fit_rate(tc[:, 0], tc[:, 1], tc[0, 0], te); rd_c = fit_rate(td[:, 0], td[:, 1], td[0, 0], te)
        d['common_window_tau'] = [float(tc[0, 0]), float(te)]
        d['rate_common_charged'] = rc_c; d['rate_common_decharged'] = rd_c
        d['ratio_common'] = rc_c / rd_c if rc_c and rd_c else None
        # early-only window (first 40 conformal units) - least box-affected
        te2 = tc[0, 0] + 40
        a1 = fit_rate(tc[:, 0], tc[:, 1], tc[0, 0], te2); a2 = fit_rate(td[:, 0], td[:, 1], td[0, 0], te2)
        d['ratio_first40'] = a1 / a2 if a1 and a2 else None
        # half-run rates (late halves are box/stripe affected when xi_c > L/8)
        tm = 0.5 * (tc[0, 0] + tc[-1, 0])
        d['rate_halves_charged'] = [fit_rate(tc[:, 0], tc[:, 1], tc[0, 0], tm), fit_rate(tc[:, 0], tc[:, 1], tm, tc[-1, 0])]
        d['rate_halves_decharged'] = [fit_rate(td[:, 0], td[:, 1], td[0, 0], tm), fit_rate(td[:, 0], td[:, 1], tm, td[-1, 0])]
        if tag.startswith('frozen'):
            Lr = arr(rc, 'L_raw'); tt = arr(rc, 'tau')
            d['L_raw_charged_range_tau_ge_260'] = [float(Lr[tt >= 260].min()), float(Lr[tt >= 260].max())]
            d['L_basin_charged_fall_factor'] = float(arr(rc, 'L_basin')[0] / arr(rc, 'L_basin')[-1])
            d['L_raw_charged_fall_factor'] = float(Lr[0] / Lr[-1])
            d['L_raw_based_rate_charged_tau_ge_260'] = fit_rate(tt, 1.0 / Lr, 260, tt[-1])
        summ['decharge'][tag] = d

json.dump(summ, open(os.path.join(R, 'v_charged_exit_summary.json'), 'w'), indent=1)
print(json.dumps({k: v for k, v in summ.items() if k != 'runs'}, indent=1, default=str)[:6000])
for k, v in summ['runs'].items():
    print(k, {kk: (round(vv, 4) if isinstance(vv, float) else vv) for kk, vv in v.items() if kk in (
        'a_form', 'n_form_over_nc', 'p_plus_form', 'p_plus_at_vac0.8', 'p_plus_end', 'walls_every_output_after_form',
        'rate_boxsafe', 'rate_all_post', 'rate_first_half', 'rate_second_half', 'xi_c_end', 'xi_c_end_over_L', 'dyn_range_xi_phys',
        'qfrac_wallcomp_end', 'qfrac_within_2_end', 'qfrac_within_4_end', 'qfrac_Qmatter_end', 'qfrac_vac_end', 'Lraw_over_Lbasin_end',
        'max_Q_rel_dev', 'maxcl_plus_end', 'maxcl_minus_end')})
