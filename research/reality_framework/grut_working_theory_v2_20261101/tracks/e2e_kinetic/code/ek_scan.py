#!/usr/bin/env python3
"""e2e_kinetic scans: numerics, H_i (mf_chain's rescaling rule), the neutrino completion (m~), the SM rates, the
perturbative shell width, the lattice interface (restart from moments).  Each case writes $SP/run_<tag>.json;
'collect' writes ../results/scan.json.  Every case is a model choice; labels in NOTES.md.
Usage: python3 -B ek_scan.py run TAG [TAG ...];  python3 -B ek_scan.py collect;  python3 -B ek_scan.py list
"""
import json, math, sys, time
import numpy as np
import ek_common as C
import ek_run as R

BASE = dict(name='HFa', H_GeV=C.H_REF, mnu=C.MNU_REF, t0=0.19, t_coh=3.0, t_end=1e8, lam_group='SU2')
CASES = {
    # reference and numerics
    'ref': {}, 'refb': dict(name='HFb'),
    'num_dk035': dict(dk_m=0.035), 'num_cdt02': dict(cdt=0.02), 'num_nq1400': dict(nqlog=1400), 'num_nq350': dict(nqlog=350),
    'tcoh1': dict(t_coh=1.0), 'tcoh10': dict(t_coh=10.0), 'tcoh30': dict(t_coh=30.0),
    # neutrino completion (all four singlets at a common m~); per-species derivatives
    'mnu0.15': dict(mnu=0.15e-9), 'mnu0.015': dict(mnu=0.015e-9), 'mnu0.005': dict(mnu=0.005e-9), 'mnu0.5': dict(mnu=0.5e-9),
    'mnuN1x2': dict(mnu=[0.1e-9, 0.05e-9, 0.05e-9, 0.05e-9]), 'mnuN2x2': dict(mnu=[0.05e-9, 0.1e-9, 0.05e-9, 0.05e-9]),
    'mnuN3x2': dict(mnu=[0.05e-9, 0.05e-9, 0.1e-9, 0.05e-9]), 'mnuN4x2': dict(mnu=[0.05e-9, 0.05e-9, 0.05e-9, 0.1e-9]),
    # H_i (rescaling rule of mf_chain: masses and f fixed in H_i units, amplitude in M_P units)
    'H2.4e11': dict(H_GeV=2.4e11), 'H8.4e10': dict(H_GeV=8.4e10), 'H5.0e10': dict(H_GeV=5.0e10),
    'bH8.4e10': dict(name='HFb', H_GeV=8.4e10), 'bH5.0e10': dict(name='HFb', H_GeV=5.0e10),
    # SM rates
    'lamSU3': dict(lam_group='SU3'), 'KLx10': dict(KLscale=10.0), 'KLx0.1': dict(KLscale=0.1),
    'cBBx2': dict(cBBscale=2.0), 'cBBx0.5': dict(cBBscale=0.5), 'cBBx10': dict(cBBscale=10.0), 'G12hi': dict(G12=5.0e-4),
    'noinv': dict(winv=0.0),
    # perturbative shell width (fragmented / moving parents broaden it; lattice input)
    'dshell0.03': dict(dshell=0.03), 'dshell0.1': dict(dshell=0.1), 'dshell0.3': dict(dshell=0.3), 'dshell1': dict(dshell=1.0),
    'H8.4e10_dshell0.3': dict(H_GeV=8.4e10, dshell=0.3), 'H5.0e10_dshell0.3': dict(H_GeV=5.0e10, dshell=0.3),
    # lattice interface: restart at t_sw from moments only (see restart())
    'restart0.45': dict(restart=0.45), 'restart1.0': dict(restart=1.0),
}


def occ_from_moments(E, Nn, q, w, a, M):
    """occupation (both helicities) uniform up to k_max with the given physical energy and number density: the minimal
    information a lattice drain term provides (energy per species and mean energy). k_max from E/N (relativistic)."""
    if E <= 0 or Nn <= 0:
        return np.zeros(2 * len(q))
    kmax = 4.0 * E / (3.0 * Nn)
    k = q / a
    phase = np.sum(np.where(k < kmax, 2 * w / a ** 3, 0.0))
    nbar = min(Nn / max(phase, 1e-300), 0.5)
    o = np.where(k < kmax, nbar, 0.0)
    return np.r_[o, o]


def restart(tag, t_sw, **kw):
    """run the reference to t_sw (coherent), keep only the moments (per-species N energy and number, condensate and SM
    energies), rebuild occupations from moments, continue; compare with the uninterrupted reference."""
    kw0 = dict(BASE); kw0.update(kw)
    r1 = R.run_history(kw0['name'], tag + '_a', kw0['H_GeV'], kw0['mnu'], kw0['t0'], t_sw, t_sw * 1.0000001, kw0['lam_group'])
    rec = r1['rec']
    spec, mpar = C.species(kw0['name'], kw0['H_GeV'], kw0['mnu'])
    st = r1['info']['state_tcoh']
    rho = st['rho_tot']
    a = st['a']
    occ0 = []
    for j, s in enumerate(spec):
        E = float(rec['fN%d' % (j + 1)][-1] * rec['rho'][-1])
        Nn = float(rec['nN%d' % (j + 1)][-1])
        q, w = C.grid_linear(a * 3e4 * (1 if s['par'] else 1), 6000)
        occ0.append((q, w, occ_from_moments(E, Nn, q, w, a, s['M'])))
    state = dict(t=st['t'], a=a, H=st['H'], rho_sig=st['rho_sig'], rho_chi=st['rho_chi'], rho_h=st['rho_h'], n_h=st['n_h'], rho_th=st['rho_th'])
    r2 = R.run_history(kw0['name'], tag, kw0['H_GeV'], kw0['mnu'], t_sw, kw0['t_coh'] if kw0['t_coh'] > t_sw else t_sw + 2.0,
                       kw0['t_end'], kw0['lam_group'], state=state, occ0=occ0)
    r2['info']['restart_from'] = dict(t_sw=t_sw, moments='per-species E_N, n_N; condensates; SM (rho_h, n_h, rho_th)',
                                      first_leg_summary=r1['summary'])
    return r2


def run(tag):
    kw = dict(BASE); kw.update(CASES[tag])
    if 'restart' in kw:
        t_sw = kw.pop('restart')
        res = restart(tag, t_sw)
    else:
        name = kw.pop('name')
        res = R.run_history(name, tag, kw.pop('H_GeV'), kw.pop('mnu'), kw.pop('t0'), kw.pop('t_coh'), kw.pop('t_end'),
                            kw.pop('lam_group'), **kw)
    C.dump(C.SP / ('run_%s.json' % tag), R.jsonable(res, keep_rec=True))
    return res


def collect():
    out = dict(label='C on the stated rates (model choices)', note=__doc__, cases={})
    for tag in CASES:
        fn = C.SP / ('run_%s.json' % tag)
        if not fn.exists():
            continue
        d = json.load(open(fn))
        inf = d['info']
        out['cases'][tag] = dict(settings={k: inf.get(k) for k in ('name', 'H_GeV', 'mnu', 't0', 't_coh', 'lam_group', 'lam_eff', 'KLscale',
                                                                      'cBB', 'G12', 'dshell', 'nqlog', 'dk_m', 'cdt', 'winv')},
                                 summary=d['summary'], seconds=inf.get('seconds'), restart_from=inf.get('restart_from'))
    C.dump(C.RES / 'scan.json', out)
    for tag, v in out['cases'].items():
        s = v['summary']
        print('%-18s K %6.2f  t_SM50 %9.4g  t_th %9.4g  T_RH %.3g  Tmax %.3g  TK %.3g' % (tag, s['K'], s['t_SM_half'] or -1, s['t_thermal'] or -1,
              s['T_RH_GeV'] or -1, s['T_max_GeV'], s['T_K_GeV']))


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        R.BUDGET[0] = 470.0
        for tag in sys.argv[2:]:
            if (C.SP / ('run_%s.json' % tag)).exists() and not (C.SP / ('ckpt_%s.npz' % tag)).exists():
                print(tag, 'done already', flush=True)
                continue
            if time.time() - R.T_GLOBAL[0] > R.BUDGET[0]:
                print('BUDGET: stop before', tag, flush=True)
                break
            try:
                run(tag)
            except R.Incomplete:
                print('INCOMPLETE', tag, '(resume with the same command)', flush=True)
                break
    elif sys.argv[1] == 'collect':
        collect()
    else:
        print('\n'.join(CASES))
