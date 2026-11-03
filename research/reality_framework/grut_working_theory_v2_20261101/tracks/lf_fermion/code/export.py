"""Exports of the calibrated-closure lattice runs: results/hist_HFa_LF.json (e2e_lattice hist format: lattice block of
box means to t = 3, plus the HB.stage2 fluid continuation for K as in e2e_lattice [A]) and results/fluid_moments_HFa_LF.json
(hand-over moments for e2e_chain / e2e_kinetic at t_sw = 1.5, 2, 2.5, 3).  The functions runmean and end_state are
re-typed from tracks/e2e_lattice/code/analyze.py (unchanged logic).  Claude, lf_fermion (issue #13 context)."""
import json, math, sys
import numpy as np
sys.dont_write_bytecode = True
import lfcommon as LC

PRIM = 'KFmid'
ALTS = ['KFlo_nl', 'KFhi_nl', 'KFmid_nl_s2', 'KFmid_H84_nl']


def runmean(t, y, w=0.05):
    c = np.concatenate([[0.0], np.cumsum(y)])
    lo = np.searchsorted(t, t - w / 2); hi = np.searchsorted(t, t + w / 2, side='right')
    return (c[hi] - c[lo]) / np.maximum(hi - lo, 1)


def end_state(D, v, i=-1):
    rho = D['rho'][i]
    return dict(t=float(D['t'][i]), a=float(D['a'][i]), rs=float(D['fsig'][i] * rho * v * v),
                rc=float((D['fchi'][i] + D['fgrad'][i] + D['fs'][i]) * rho * v * v),
                E=[float(D['rN%d' % j][i] * v * v) for j in (1, 2, 3, 4)],
                n=[float(D['nN%d' % j][i] * v * v) for j in (1, 2, 3, 4)],
                B=[float(D['B%d' % j][i] * rho * v * v) for j in (1, 2, 3, 4)], sm=float(D['rSM'][i] * v * v))


f7 = lambda y: [float('%.7g' % x) for x in y]


def lattice_block(D):
    t = D['t']
    sel = np.unique(np.concatenate([np.searchsorted(t, np.arange(t[0], 1.0, 0.004)),
                                    np.searchsorted(t, np.arange(1.0, t[-1] + 1e-9, 0.02)), [len(t) - 1]]))
    sel = sel[sel < len(t)]
    frad = D['fN1'] + D['fN2'] + D['fN3'] + D['fN4'] + D['fSM']
    return dict(t=f7(t[sel]), a=f7(D['a'][sel]), H=f7(D['H'][sel]), eps_H=f7(D['eps'][sel]), w=f7(D['w'][sel]),
                w_cyc=f7(runmean(t, D['w'])[sel]), R=f7(D['R_mean'][sel]), R_cyc=f7(runmean(t, D['R_mean'])[sel]),
                R_rms=f7(D['R_rms'][sel]), f_sig=f7(D['fsig'][sel]), f_chi=f7(D['fchi'][sel]),
                f_grad=f7(D['fgrad'][sel]), f_rad=f7(frad[sel]), f_SM=f7(D['fSM'][sel]), f_s=f7(D['fs'][sel]),
                s_mean=f7(D['sa_mean'][sel]), s_std=f7(D['sa_std'][sel]), m_s2_mean=f7(D['ms2a_mean'][sel]),
                m_s2_std=f7(D['ms2a_std'][sel]), fminus=f7(D['fminus'][sel]), area=f7(D['area_av'][sel]),
                K_running=f7(D['a'][sel] ** 2 * D['H'][sel]),
                f_N=[f7(D['fN%d' % j][sel]) for j in (1, 2, 3, 4)])


def head(D):
    t = D['t']
    frad = D['fN1'] + D['fN2'] + D['fN3'] + D['fN4'] + D['fSM']
    g = lambda arr, x: float(np.interp(x, t, arr))
    i_coh = np.argmax(D['fcoh'] < 0.1) if np.any(D['fcoh'] < 0.1) else None
    return dict(f_rad={str(x): g(frad, x) for x in (0.3, 0.5, 1.0, 1.5, 2.0, 3.0)},
                w_cyc_3=float(runmean(t, D['w'])[-1]), R_cyc_3=float(runmean(t, D['R_mean'])[-1]),
                f_grad_max=float(D['fgrad'].max()), t_fgrad_max=float(t[np.argmax(D['fgrad'])]),
                t_fcoh_lt_0p1=None if i_coh is None else float(t[i_coh]), fminus_max=float(D['fminus'].max()),
                budget_absmax=float(np.abs(D['budget']).max()), m_s2_mean_3=float(D['ms2a_mean'][-1]),
                m_s2_rel_spread_1p5=g(D['ms2a_std'] / D['ms2a_mean'], 1.5))


def moments(D, v, P, ts):
    out = []
    for x in ts:
        i = int(np.argmin(np.abs(D['t'] - x)))
        rho = D['rho'][i]
        frad = sum(D['fN%d' % j][i] for j in (1, 2, 3, 4)) + D['fSM'][i]
        st = end_state(D, v, i)
        out.append(dict(t_sw=float(D['t'][i]), a=float(D['a'][i]), H=float(D['H'][i]), f_sig=float(D['fsig'][i]),
                        f_chi=float(D['fchi'][i]), f_grad=float(D['fgrad'][i]), f_s=float(D['fs'][i]), f_rad=float(frad),
                        f_SM=float(D['fSM'][i]), f_N=[float(D['fN%d' % j][i]) for j in (1, 2, 3, 4)],
                        E_N_Hd4=st['E'], n_N_Hd3=st['n'], B_N_Hd4=st['B'], rho_SM_Hd4=st['sm'], rho_tot_Hd4=float(rho * v * v),
                        mean_E_over_M=[float(st['E'][j] / max(st['n'][j], 1e-300) / M) for j, M in enumerate((250., 400., 200., 120.))],
                        m_s2_mean=float(D['ms2a_mean'][i]), m_s2_std=float(D['ms2a_std'][i]), s_mean=float(D['sa_mean'][i])))
    return out


def main():
    hist = dict(name='HFa', label='lf_fermion coherent-stage history, 3D lattice 64^3 with the lattice-fermion calibrated '
                'closure (local blocking, Pauli capacity x kappa_j(t)) + HB.stage2 fluid continuation for K [A]',
                note=('Same format as tracks/e2e_lattice/results/hist_HFa.json (t, a = 1 and H_d at sigma = sigma_c). Primary '
                      'KFmid: kappa from the mean line-fermion production; alternates KFlo / KFhi (Q -+ 1 s.e.) and a second '
                      'seed. f_N = per-species N energy fractions (N1..N4). The fluid part uses the skeptic Pauli fluid '
                      'with kappa = 1 and no thermalisation physics; splice e2e_kinetic by moments instead (fluid_moments).'),
                primary=PRIM, alternates={})
    mom = dict(label='C', note=__doc__, runs={})
    for tag in [PRIM] + ALTS:
        info = json.loads((LC.RUNS / tag / 'info.json').read_text())
        d, P = LC.model(info['name'], info['H_GeV'])
        v = P['v']
        D = LC.load_diag(tag)
        cfg = dict(tag=tag, N=info['N'], L=info['L'], dt=info['dt'], seed=info['seed'], kcut_mchi=info['kcut_mchi'],
                   kappa=[l for l in (LC.RUNS / tag / 'run.par').read_text().splitlines() if l.startswith('kappa')])
        h = head(D)
        rec = LC.HB.stage2(P, end_state(D, v), blocking=True, Gth=0.0)
        h['K_cont_stage2'] = float(rec['K'][-1])
        ent = dict(config=cfg, head=h, lattice=lattice_block(D))
        if tag == PRIM:
            ent['fluid'] = dict(t=f7(rec['t']), a=f7(rec['a']), H=f7(rec['H']), eps_H=f7(rec['eps']), R_cyc=f7(rec['R']),
                                f_sig=f7(rec['fsig']), f_chi=f7(rec['fchi']), f_rad=f7(rec['frad']), f_SM=f7(rec['fsm']),
                                K_running=f7(rec['K']))
            ent['fluid_note'] = 'HB.stage2 continuation at kappa = 1 (no kinetic physics); for K use the e2e_chain hand-over'
            hist.update(ent)
        else:
            hist['alternates'][tag] = ent
        mom['runs'][tag] = moments(D, v, P, (1.5, 2.0, 2.5, 3.0))
        print(tag, json.dumps(h['f_rad']), 'K_cont %.3f' % h['K_cont_stage2'], 'fgrad_max %.3f' % h['f_grad_max'],
              't_coh %s' % h['t_fcoh_lt_0p1'], 'budget %.1e' % h['budget_absmax'], 'ms2 %.3f' % h['m_s2_mean_3'], flush=True)
    LC.dump(LC.TRACK / 'results' / 'hist_HFa_LF.json', hist)
    LC.dump(LC.TRACK / 'results' / 'fluid_moments_HFa_LF.json', mom)


if __name__ == '__main__':
    main()
