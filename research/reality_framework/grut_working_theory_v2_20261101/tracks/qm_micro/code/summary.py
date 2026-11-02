"""qm_micro summary (Claude, 2 Nov 2026): collects items 1-4 into results/summary.json, compact coarsening and lattice
results; computes e(S/Q), zeta, the vapour velocity and its redshift, a free-streaming estimate, M/m_X and H_i.
Own code.  Inputs: results/*.json of this track and the scratch run files (coarsen_*.json, lat_*.json)."""
import json, math, os, sys
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
SCR = os.environ.get('QM_SCR', '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/qm_micro')
DQ = json.load(open(os.path.join(HERE, '..', '..', 'dq_entropy', 'results', 'summary.json')))['for_critic']['SQ_of_W']
OMB = 0.7758647074595397
H_M1 = 1.342842735e11
YRATIO = 0.991023          # Y/Y_Q of the frozen M1r release (dq_entropy)


def e_of_SQ(SQ, incr=0.0055, broad=False):
    """Reading B of dq_entropy (transferred vector heat, nu = 0.960): R(e) = 0.3786 + 0.0184 x^1/2 + 0.960 (13.0412/15.6036)
    (4.0538 band | 4.7905 broad) x, x = 1e-6/e; S/Q = A (W - B)^p + C + incr, W = 15.6036 (1 + R).  [E] (dq_heat refines)."""
    A, B, p, C = DQ['A'], DQ['B'], DQ['p'], DQ['C']
    W = B + ((SQ - C - incr) / A) ** (1 / p)
    R = W / 15.6036 - 1
    Rv = 4.7905 if broad else 4.0538
    f = lambda x: 0.3786 + 0.0184 * math.sqrt(x) + 0.960 * 13.0412 / 15.6036 * Rv * x - R
    if f(0.0) > 0:
        return dict(R=R, e=None, note='below the e -> infinity floor (R_med = 0.3786, S/Q = %.3f)' % (A * (15.6036 * 1.3786 - B) ** p + C + incr))
    x = brentq(f, 0.0, 100.0)
    return dict(R=R, e=1e-6 / x)


# background for the free-streaming estimate (Planck-2018-like, as fx_observables)
h = 0.674; Om = (0.12 + 0.0224) / h ** 2; Orad = 2.4728e-5 * (1 + 0.2271 * 3.044) / h ** 2; OL = 1 - Om - Orad
H0_Mpc = h / 2997.92458
T0_GeV = 2.7255 * 8.617333e-14


def gs(Tv):
    tab = np.array([[1e-5, 3.91], [5e-4, 3.91], [1e-3, 10.75], [0.05, 10.75], [0.1, 14.25], [0.2, 55.0], [1.0, 69.0]])
    return float(np.interp(math.log(Tv), np.log(tab[:, 0]), tab[:, 1]))


def vapour_kinematics(T_f, Tvis_f):
    v_rms = math.sqrt(3 * T_f)                      # NR Maxwellian at T_f (m_X units), fugacity << 1
    af = (3.91 / gs(Tvis_f)) ** (1 / 3) * T0_GeV / Tvis_f
    v0 = v_rms * af
    # comoving free-streaming of the rms speed: lambda = int_{a_f}^{1} v(a)/a da/(a H), v = v_rms a_f/a (NR throughout)
    Hz = lambda a: H0_Mpc * math.sqrt(Orad / a ** 4 + Om / a ** 3 + OL)
    lam = quad(lambda la: v_rms * af / math.exp(la) / (math.exp(la) * Hz(math.exp(la))) , math.log(af), 0.0, limit=400)[0]
    return dict(v_rms_at_freeze_c=v_rms, a_f_over_a0=af, v_rms_today_c=v0, v_rms_today_km_s=v0 * 2.998e5,
                T_vapour_today_over_mX=T_f * af ** 2,                   # free NR gas after decoupling: T_v ~ a^-2
                lambda_fs_comoving_Mpc=lam, lambda_fs_comoving_kpc=lam * 1e3)


def main():
    out = dict(notes='qm_micro summary; labels in NOTES.md')
    # ---------------- item 1
    st = json.load(open(os.path.join(RES, 'stick_thermal.json')))
    sc = json.load(open(os.path.join(RES, 'scatter1d.json')))
    rad = json.load(open(os.path.join(RES, 'radial.json')))
    th = {('%g' % r['T_h']): r for r in st['thermal']}
    out['item1'] = dict(
        profile=sc['profile'],
        X_normal={('%g' % r['T']): r['s'] for r in sc['X_normal'] if r['T'] in (1e-4, 0.001, 0.0017988, 0.02, 0.04, 0.1, 0.3, 0.52, 1.0) or abs(r['T'] - 0.1) < 1e-9},
        thermal_flux_avg={k: dict(s_X_boltz=v['s_X_boltz'], s_X_bose=v['s_X_bose_coex'], Xbar_destroyed=v['d_Xbar_boltz'],
                                  Xbar_to_X=v['Xbar_to_X'], Xbar_to_phi=v['Xbar_to_phi'], Xbar_to_bulk=v['Xbar_to_bulk'])
                          for k, v in th.items()},
        threshold_law=st['threshold_law'], max_unitarity=st['max_unitarity'],
        resolution_control=sc['resolution_control'], domain_control=sc['domain_control'],
        finite_object=dict(undamped_Q2e5=rad['undamped_Q2e5_k0.0684'], sk_xcore=rad['sk_xcore_reference'],
                           kscan=rad['undamped_Q2e5_kscan'], absorber=rad['absorber']))
    lat = []
    for f in ('lat_lin1d.json', 'lat_nl1d.json', 'lat_q1d.json', 'lat_q1d_b.json'):
        p = os.path.join(SCR, f)
        if os.path.exists(p):
            for c in json.load(open(p))['cases']:
                lat.append({k: c[k] for k in ('label', 'inc', 'T', 'eps', 'Ny', 'theta_deg', 'kp_lattice', 'sigma', 'x0', 'dx', 'dt', 't_end', 'Nx',
                                              'frac_dir', 'slab_charge_over_Nin', 'vac_charge_over_Nin', 'pred_packet', 'pred_carrier',
                                              'energy_drift_pk', 'charge_drift_pk', 'Lb', 'Lv')})
    json.dump(dict(cases=lat), open(os.path.join(RES, 'lattice.json'), 'w'), indent=1, default=float)
    out['item1']['lattice_vs_bvp'] = [dict(inc=c['inc'], T=c['T'], eps=c['eps'], Ny=c['Ny'], theta=c['theta_deg'],
                                           lattice_reflected=c['frac_dir']['X_right'] if c['inc'] == 'X' else c['frac_dir']['Xbar_right'],
                                           bvp_reflected=c['pred_packet']['out_X'] if c['inc'] == 'X' else c['pred_packet']['out_Xbar'],
                                           lattice_phi=c['frac_dir']['phi_right'], lattice_other=c['frac_dir']['Xbar_right'] if c['inc'] == 'X' else c['frac_dir']['X_right'],
                                           slab_charge=c['slab_charge_over_Nin']) for c in lat]
    # ---------------- item 2
    gt = json.load(open(os.path.join(RES, 'vapour_gt.json')))
    out['item2'] = dict(kinetic=[{k: r[k] for k in ('T', 'mu_c', 'n_X', 's_flux_avg', 'L', 'K_W')} for r in gt['kinetic']],
                        wagner_test=dict(Rc2_over_KWt_late=float(np.median([r['Rc2_over_KWt'] for r in json.load(open(os.path.join(RES, 'wagner_test.json')))['rows'][-20:]])),
                                         Rmax_over_Rc_late=float(np.median([r['Rmax_over_Rc'] for r in json.load(open(os.path.join(RES, 'wagner_test.json')))['rows'][-20:]])),
                                         Rbar_over_Rc_late=float(np.median([r['Rbar_over_Rc'] for r in json.load(open(os.path.join(RES, 'wagner_test.json')))['rows'][-20:]])),
                                         wagner_max_cutoff=0.5))
    # ---------------- items 3-4
    item3 = {}
    compact = {}
    for SQ in ('1.5', '2.5', '3.0', '3.2'):
        p = os.path.join(SCR, 'coarsen_%s.json' % SQ)
        if not os.path.exists(p):
            continue
        d = json.load(open(p))
        runs = {}
        for r in d['runs']:
            f, z, b = r['freeze'], r['final'], r['variant_b']
            s_f = b['sticking_at_freeze']
            kin = vapour_kinematics(f['T'], f['T_vis_GeV'])
            EQ = z['EoverQ_charge_weighted']
            M_a = (1 - z['f_v']) * EQ + z['f_v']
            M_b = (1 - b['f_v_final']) * EQ + b['f_v_final']
            runs[r['case']] = dict(stick_mode=r['stick_mode'], Q0=r['Q0'], T_f=f['T'], T_vis_f_GeV=f['T_vis_GeV'],
                                   f_v_at_freeze=f['f_v'], f_v_final_eqT=z['f_v'], f_v_final_decoupled=b['f_v_final'],
                                   Q_charge_weighted_freeze=f['Q_charge_weighted'], Q_capture_eff_freeze=f['Q_capture_eff'],
                                   sigma_lnQ=f['sigma_lnQ_charge'], Q_charge_weighted_final=z['Q_charge_weighted'],
                                   s_thermal_at_freeze=s_f,
                                   zeta_charge_weighted=s_f * f['Q_charge_weighted'] ** (-1 / 3),
                                   zeta_capture_eff=s_f * f['Q_capture_eff'] ** (-1 / 3),
                                   EoverQ_objects=EQ, M_per_charge=[M_b, M_a],
                                   H_i_GeV=[H_M1 * (0.9374 / (M * YRATIO)) ** (2 / 3) for M in (M_a, M_b)],
                                   vapour=kin, steps=r['steps'], runtime_s=r['runtime_s'])
            thin = r['rows'][::4] + [r['rows'][-1]]
            compact.setdefault(SQ, {})[r['case']] = dict(freeze=f, final=z, variant_b=b,
                                                         rows=[{k: row[k] for k in ('lna', 'T', 'T_vis_GeV', 'eta', 'f_v_eq', 'Q_charge_weighted', 'Q_capture_eff',
                                                                                    'sigma_lnQ_charge', 'Gamma_capture_over_H', 'coag_rate_over_H', 'phi_vol')} for row in thin],
                                                         freeze_hist=r.get('freeze_charge_histogram_log10Q'))
        item3[SQ] = dict(e_band=e_of_SQ(float(SQ)), e_broad=e_of_SQ(float(SQ), broad=True), T_spinodal=d['T_spinodal'],
                         hydro_stage=d['hydro_stage'], runs=runs)
        b = runs['Q0_4e11']
        print('S/Q=%s e(band)=%s e(broad)=%s | T_f=%.4f (T_vis %.3g MeV) f_v=%.4f [%.4f..%.4f] <Q>=%.3e Qeff=%.3e s(T_f)=%.3f zeta=%.2e | v_rms=%.3f c -> today %.2e km/s, lambda_fs=%.3g kpc | M=%.5f..%.5f H_i=%.4e..%.4e' % (
            SQ, item3[SQ]['e_band'].get('e'), item3[SQ]['e_broad'].get('e'), b['T_f'], b['T_vis_f_GeV'] * 1e3, b['f_v_at_freeze'],
            b['f_v_final_decoupled'], b['f_v_final_eqT'], b['Q_charge_weighted_freeze'], b['Q_capture_eff_freeze'], b['s_thermal_at_freeze'],
            b['zeta_charge_weighted'], b['vapour']['v_rms_at_freeze_c'], b['vapour']['v_rms_today_km_s'], b['vapour']['lambda_fs_comoving_kpc'],
            b['M_per_charge'][0], b['M_per_charge'][1], b['H_i_GeV'][1], b['H_i_GeV'][0]))
        for k, v in runs.items():
            if k != 'Q0_4e11':
                print('      %-26s <Q>=%.3e f_v=%.4f..%.4f T_f=%.4f zeta=%.2e M=%.5f..%.5f' % (k, v['Q_charge_weighted_freeze'], v['f_v_final_decoupled'],
                                                                                    v['f_v_final_eqT'], v['T_f'], v['zeta_charge_weighted'], v['M_per_charge'][0], v['M_per_charge'][1]))
    out['item3_4'] = item3
    out['item3_controls'] = json.load(open(os.path.join(RES, 'coarsen_controls.json')))
    json.dump(compact, open(os.path.join(RES, 'coarsen_runs.json'), 'w'), indent=None, default=float, separators=(',', ':'))
    json.dump(out, open(os.path.join(RES, 'summary.json'), 'w'), indent=1, default=float)


if __name__ == '__main__':
    main()
