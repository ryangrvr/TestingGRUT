"""e2e_lattice analysis: fluid continuation of each lattice end state (K), GW today, wall extrapolation with the bias
pressure, convergence table, and the chain exports results/hist_HFa.json, results/hist_HFb.json.
Usage: python3 -B analyze.py  ->  ../results/analysis.json, ../results/hist_HFa.json, ../results/hist_HFb.json

Closures (all [A], model choices):
  * continuation: HB.stage2 (mf_chain's copy of the mf_hybrid skeptic's per-species Pauli fluid) started from the box
    means at the lattice end; sigma energy -> rs, chi + gradient + spectator energy -> rc (matter-like fragments);
    N energies, numbers, blocking energies and SM radiation taken from the lattice fluids.  blocking = (closure is pauli).
  * GW: rho_GW a^4 is frozen once the source switches off (checked: constant to 1% after t = 1); Omega_GW in the
    late radiation era = (d rho_GW / d ln k) a^4 / (3 M_P^2 K^2); today with g_* = 106.75 and entropy conservation.
Every number is [C] on the lattice model unless labelled otherwise.
"""
import json, math, sys
import numpy as np
import latcommon as LC
import latpost as LP

TAGS = ['P0', 'P0GM', 'P10', 'P10GM', 'TR', 'HU', 'B_P0', 'B_P0GM', 'H84_P0', 'H84_P0GM', 'C_dt2', 'C_L63', 'C_N96', 'C_N96L63']
K_KINETIC = {'HFa': 17.698}           # e2e_kinetic/results/fluid_HFa.json summary K (read-only), for the GW bracket


def end_state(tag, D, v, i=-1):
    rho = D['rho'][i]
    return dict(t=float(D['t'][i]), a=float(D['a'][i]), rs=float(D['fsig'][i] * rho * v * v),
                rc=float((D['fchi'][i] + D['fgrad'][i] + D['fs'][i]) * rho * v * v),
                E=[float(D['rN%d' % j][i] * v * v) for j in (1, 2, 3, 4)],
                n=[float(D['nN%d' % j][i] * v * v) for j in (1, 2, 3, 4)],
                B=[float(D['B%d' % j][i] * rho * v * v) for j in (1, 2, 3, 4)], sm=float(D['rSM'][i] * v * v))


def runmean(t, y, w=0.05):
    c = np.concatenate([[0.0], np.cumsum(y)])
    lo = np.searchsorted(t, t - w / 2); hi = np.searchsorted(t, t + w / 2, side='right')
    return (c[hi] - c[lo]) / np.maximum(hi - lo, 1)


def walls(P):
    """Bias pressure vs wall tension (rho~ = rho/v^2 units, H_d = 1).  [D] given the potential; the network
    extrapolation is [E] (a scaling network with one wall per Hubble volume is assumed only to show what would happen
    if walls had formed)."""
    v = P['v']; lam = P['lam']; q = (v / P['mP']) ** 2
    dV = 2.0 * P['cb'] / v                                     # V(-v) - V(+v) at sigma = 0, leading order in the bias
    sig_w = (2.0 * math.sqrt(2.0) / 3.0) * math.sqrt(lam) * v   # wall tension / v^2
    delta = math.sqrt(2.0) / (math.sqrt(lam) * v)               # wall thickness (physical, H_d^-1)
    lb = sig_w / dV                                             # radius above which the pressure beats the tension
    H_dec = dV / sig_w                                          # network decay when H ~ dV / sigma_w (scaling, A ~ 1)
    H_dom = sig_w * q / 3.0                                     # rho_wall = sigma_w H equals 3 M_P^2 H^2
    return dict(dV_over_v2=dV, sigma_w_over_v2=sig_w, wall_thickness=delta, ell_bias=lb, H_dec_scaling=H_dec,
                H_dom_scaling=H_dom, eps_bias=dV / (0.25 * lam * v * v),
                verdict_if_network='dominates first' if H_dom > H_dec else 'decays first')


def main():
    out = dict(runs={}, convergence={}, notes=__doc__)
    hist_src = {}
    for tag in TAGS:
        R = json.load(open(LC.RES / 'runs' / ('%s.json' % tag)))
        info = json.load(open(LC.RUNS / tag / 'info.json'))
        d, P = LC.model(info['name'], info['H_GeV'])
        v, q = P['v'], (P['v'] / P['mP']) ** 2
        D = LC.load_diag(tag)
        t = D['t']
        row = dict(R['head'])
        row['config'] = R['config']
        # fluid continuation
        rec = None
        if info['closure'] != 'hubble':
            st = end_state(tag, D, v)
            rec = LC.HB.stage2(P, st, blocking=(info['closure'] == 'pauli'), Gth=info['Gth'])
            row['K_cont'] = float(rec['K'][-1]); row['K_cont_converged'] = bool(rec['converged'])
            j = int(np.argmax(rec['frad'] >= 0.5)) if (rec['frad'] >= 0.5).any() else None
            row['t_frad_half_cont'] = float(rec['t'][j]) if j is not None else None
        # walls: last time with any -v site, max -v fraction, largest -v cluster
        m = D['fminus'] > 0
        row['t_last_minus'] = float(t[m][-1]) if m.any() else None
        row['t_last_minus_1em4'] = float(t[D['fminus'] > 1e-4][-1]) if (D['fminus'] > 1e-4).any() else None
        # GW
        gw = []
        for s in R['snaps']:
            if 'gw' not in s:
                continue
            k = np.array(s['gw']['k']); sp = np.array(s['gw']['dlnk'])
            rho_t = float(np.interp(s['t'], t, D['rho']))
            e = dict(t=s['t'], a=s['a'], rhoGW_a4=s['gw']['rho_tt'] * s['a'] ** 4, OmGW_e=s['gw']['rho_tt'] / rho_t,
                     kpeak=s['gw']['kpeak'], k_nyq=math.pi * info['N'] / info['L'],
                     frac_above_nyq=float(sp[k > math.pi * info['N'] / info['L']].sum() / max(sp.sum(), 1e-300)))
            Ks = {}
            if rec is not None:
                Ks['cont'] = row['K_cont']
            if info['name'] in K_KINETIC and info['H_GeV'] == LC.H_REF:
                Ks['kinetic'] = K_KINETIC[info['name']]
            for nm, K in Ks.items():
                f, Om, aux = LP.today(k, sp * s['a'] ** 4 * q / (3.0 * K * K), info['H_GeV'], K)
                ip = int(np.argmax(Om))
                e['today_K_%s' % nm] = dict(K=K, f_peak_Hz=float(f[ip]), Omh2_peak=float(Om[ip]),
                                            Omh2_total=float(s['gw']['rho_tt'] * s['a'] ** 4 * q / (3 * K * K) * aux['factor']),
                                            f_chi_scale_Hz=float(np.interp(s['chi_kpeak'], k, f)))
            gw.append(e)
        row['gw'] = gw
        row['walls'] = walls(P)
        out['runs'][tag] = row
        hist_src[tag] = (D, rec, info, P, R)
        print(tag, 'K_cont=%s' % row.get('K_cont'), 'tlast-=%s' % row['t_last_minus'],
              'OmGW_e(last)=%s' % (gw[-1]['OmGW_e'] if gw else None), flush=True)
    # convergence table at t = 1.5 against P0
    keys = ['frad_t1.5', 'fgrad_t1.5', 'w_cyc_t1.5', 'R_cyc_t1.5', 'R_rms_t1.5', 'chi_std_t1.5', 'ms2a_std_over_mean_t1.5',
            'ms2b_std_t1.5', 'fminus_max', 't_fminus_max', 'fgrad_max', 't_fcoh_lt_0p1', 'R_rms_max', 'budget_absmax']
    base = out['runs']['P0']
    D0, _, info0, P0_, _ = hist_src['P0']
    i15 = int(np.searchsorted(D0['t'], 1.4999))
    K15 = float(LC.HB.stage2(P0_, end_state('P0', D0, P0_['v'], i15), blocking=True, Gth=0.0)['K'][-1])
    base['K_cont_switch1.5'] = K15
    for tag in ('C_dt2', 'C_L63', 'C_N96', 'C_N96L63'):
        r = out['runs'][tag]
        c = {k: dict(P0=base[k], run=r[k], rel=(r[k] - base[k]) / base[k] if base[k] not in (0, None) and r[k] is not None
                     else None) for k in keys if k in r and k in base}
        g0 = [g for g in base['gw'] if abs(g['t'] - 1.5) < 0.01]
        g1 = [g for g in r['gw'] if abs(g['t'] - 1.5) < 0.01]
        if g0 and g1:
            c['rhoGW_a4_t1.5'] = dict(P0=g0[0]['rhoGW_a4'], run=g1[0]['rhoGW_a4'], rel=g1[0]['rhoGW_a4'] / g0[0]['rhoGW_a4'] - 1)
            c['kpeak_GW_t1.5'] = dict(P0=g0[0]['kpeak'], run=g1[0]['kpeak'], k_nyq_run=g1[0]['k_nyq'])
        c['K_cont_switch1.5'] = dict(P0=K15, run=r['K_cont'], rel=r['K_cont'] / K15 - 1)
        out['convergence'][tag] = c
    LC.dump(LC.RES / 'analysis.json', out)
    # exports
    for name, prim, alts in (('HFa', 'P0', ['P0GM', 'P10', 'P10GM', 'H84_P0', 'H84_P0GM']), ('HFb', 'B_P0', ['B_P0GM'])):
        hist = dict(name=name, label='e2e_lattice coherent-stage history (3D lattice 64^3) + fluid continuation',
                    note=('t, a = 1 and H_d at sigma = sigma_c of the homogeneous attractor (H_i/H_d = 1 + 2e-5). Part "lattice": '
                          'box means of the 3D run (t <= 3 H^-1); part "fluid": HB.stage2 continuation of the end state '
                          '[A]. The chain should splice e2e_kinetic fluid_%s.json at t_sw in [1.5, 3] where the fractions '
                          'are matched by moments; the lattice part is the coherent stage, the fluid part only gives K for '
                          'this closure.' % name),
                    columns=dict(t='H_d^-1 from the slice', a='scale factor', H='H/H_d', eps_H='-Hdot/H^2 (raw)',
                                 w='p/rho box mean (raw)', w_cyc='w averaged over 0.05 H^-1', R='box-mean Ricci scalar / H_d^2 = (rho - 3p)/M_P^2 (raw)',
                                 R_cyc='R averaged over 0.05 H^-1', R_rms='spatial rms of R(x) over the box (raw)',
                                 f_sig='sigma (kinetic + g^2 chi^2 (s - s0)^2 / 2)', f_chi='chi (rest of potential + kinetic)',
                                 f_grad='gradient energy', f_rad='N + SM', f_SM='SM radiation', f_s='spectators',
                                 s_mean='HFa: s_a/s_* box mean; HFb: s_b/s_*', s_std='spatial std of s/s_*',
                                 m_s2_mean='spectator mass^2 / H_d^2 box mean (HFa: (g_s v chi)^2; HFb: 9 - 3R/4)',
                                 m_s2_std='spatial std of m_s^2', fminus='-v vacuum volume fraction (chi < 0)',
                                 area='wall area per comoving volume from sign changes of the 0.05-averaged chi (H_d)',
                                 K_running='a^2 H'),
                    primary=prim, alternates={})
        for tag in [prim] + alts:
            D, rec, info, P, R = hist_src[tag]
            t = D['t']
            sel = np.unique(np.concatenate([np.searchsorted(t, np.arange(t[0], 1.0, 0.004)),
                                            np.searchsorted(t, np.arange(1.0, t[-1] + 1e-9, 0.02)), [len(t) - 1]]))
            sel = sel[sel < len(t)]
            sk = 'sa' if name == 'HFa' else 'sb'
            mk = 'ms2a' if name == 'HFa' else 'ms2b'
            f7 = lambda y: [float('%.7g' % x) for x in y]
            frad = D['fN1'] + D['fN2'] + D['fN3'] + D['fN4'] + D['fSM']
            lat = dict(t=f7(t[sel]), a=f7(D['a'][sel]), H=f7(D['H'][sel]), eps_H=f7(D['eps'][sel]), w=f7(D['w'][sel]),
                       w_cyc=f7(runmean(t, D['w'])[sel]), R=f7(D['R_mean'][sel]), R_cyc=f7(runmean(t, D['R_mean'])[sel]),
                       R_rms=f7(D['R_rms'][sel]), f_sig=f7(D['fsig'][sel]), f_chi=f7(D['fchi'][sel]),
                       f_grad=f7(D['fgrad'][sel]), f_rad=f7(frad[sel]), f_SM=f7(D['fSM'][sel]), f_s=f7(D['fs'][sel]),
                       s_mean=f7(D[sk + '_mean'][sel]), s_std=f7(D[sk + '_std'][sel]), m_s2_mean=f7(D[mk + '_mean'][sel]),
                       m_s2_std=f7(D[mk + '_std'][sel]), fminus=f7(D['fminus'][sel]), area=f7(D['area_av'][sel]),
                       K_running=f7(D['a'][sel] ** 2 * D['H'][sel]))
            ent = dict(config=R['config'], head=out['runs'][tag], lattice=lat)
            if rec is not None:
                fl = dict(t=f7(rec['t']), a=f7(rec['a']), H=f7(rec['H']), eps_H=f7(rec['eps']),
                          R_cyc=f7(rec['R']), f_sig=f7(rec['fsig']), f_chi=f7(rec['fchi']),
                          f_rad=f7(rec['frad']), f_SM=f7(rec['fsm']), K_running=f7(rec['K']))
                ent['fluid'] = fl
                ent['fluid_note'] = 'R_cyc here is (rho_s + rho_c)/M_P^2 in H_d^2 units (matter part only), as in hist_build'
            if tag == prim:
                hist.update(ent)
            else:
                hist['alternates'][tag] = dict(config=R['config'], head=out['runs'][tag], lattice=lat,
                                               K_cont=out['runs'][tag].get('K_cont'))
        LC.dump(LC.RES / ('hist_%s.json' % name), hist)
    return out


if __name__ == '__main__':
    main()
