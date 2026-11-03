#!/usr/bin/env python3
"""e2e_kinetic, resumed run: the perturbative-shell width study and the natural-width reference.

The dephased stage's golden-rule channels (sigma -> N1 N1, chi -> N4 N4, 2M < m) put all decay products into ONE
log bin when dshell = 0.  With fast N removal the shell occupation sits at the steady state n = R/(2R + G), R ~ 1/(bin
width), so the Pauli-blocked transfer rate depended on the grid (K = 20.7 / 22.4 / 24.3 at 350 / 700 / 1400 bins).
The physical width of the shell is the line of the damped final states: each N has removal rate G(k_res), the
pair-energy Lorentzian has HWHM G, i.e. HWHM_lnk = G / (m beta^2) [D, kinematics; A, Breit-Wigner final states].
ek_modes.c now computes it per species and per step when dshell <= -100 (fac = -dshell/1000 multiplies it).
Cases below: fixed Lorentzian 0.0063 (= G_N4/(m beta^2) at the reference, the previous session's choice) at two
grids, the self-consistent width at two grids, LCUT (Lorentzian cut in HWHM units) 20 / 80, width x 0.5 / x 2,
HFb, m~ and H_i.  Each case writes $SP/run_<tag>.json.  Every case is a model choice; labels in NOTES.md.
Usage: python3 -B ek_extra.py run TAG [TAG ...];  python3 -B ek_extra.py collect
"""
import json, os, shutil, sys, time
import ek_common as C
import ek_run as R

NAT = -1000.0
NATW = -1.0e5
CASES = {
    'natL063_nq6400': dict(dshell=-0.0063, nqlog=6400), 'natL063_nq12800': dict(dshell=-0.0063, nqlog=12800),
    'nat_nq6400': dict(dshell=NAT, nqlog=6400), 'nat_nq12800': dict(dshell=NAT, nqlog=12800),
    'nat_nq3200': dict(dshell=NAT, nqlog=3200),
    'nat_x0.5': dict(dshell=0.5 * NAT, nqlog=12800), 'nat_x2': dict(dshell=2 * NAT, nqlog=6400),
    'nat_HFb': dict(name='HFb', dshell=NAT, nqlog=6400),
    'nat_mnu0.15': dict(mnu=0.15e-9, dshell=NAT, nqlog=6400), 'nat_mnu0.015': dict(mnu=0.015e-9, dshell=NAT, nqlog=25600),
    'nat_H2.4e11': dict(H_GeV=2.4e11, dshell=NAT, nqlog=6400), 'nat_H5.0e10': dict(H_GeV=5.0e10, dshell=NAT, nqlog=6400),
    'nat_tcoh10': dict(t_coh=10.0, dshell=NAT, nqlog=6400),
    # omega-Lorentzian with the kinematic range of chi -> N N*(-> L H) instead of a cut in ln k (ek_modes.c, dshell <= -1e5)
    'natw_nq6400': dict(dshell=NATW, nqlog=6400), 'natw_nq12800': dict(dshell=NATW, nqlog=12800),
    'natw_HFb': dict(name='HFb', dshell=NATW, nqlog=6400), 'natw_x0.5': dict(dshell=0.5 * NATW, nqlog=12800),
    'natw_x2': dict(dshell=2 * NATW, nqlog=6400),
    'natw_mnu0.15': dict(mnu=0.15e-9, dshell=NATW, nqlog=6400), 'natw_mnu0.015': dict(mnu=0.015e-9, dshell=NATW, nqlog=25600),
    'natw_H2.4e11': dict(H_GeV=2.4e11, dshell=NATW, nqlog=6400), 'natw_H5.0e10': dict(H_GeV=5.0e10, dshell=NATW, nqlog=6400),
    'natw_tcoh10': dict(t_coh=10.0, dshell=NATW, nqlog=6400), 'natw_tcoh1': dict(t_coh=1.0, dshell=NATW, nqlog=6400),
}
# fx_: after the factor-2 fix of the shell rate (ek_modes.c line ~434); the nat*/natL*/natw* runs above are pre-fix
for _k in ('natw_nq6400', 'natw_nq12800', 'natw_HFb', 'natw_mnu0.15', 'natw_mnu0.015', 'natw_H2.4e11', 'natw_H5.0e10',
           'natw_tcoh10', 'natw_x2', 'natw_x0.5'):
    CASES['fx_' + _k] = CASES[_k]
CASES['fx_g0.1'] = dict(dshell=0.1, nqlog=3200)          # Gaussian ln-k spread of the products (moving / fragmented parent)
CASES['fx_g0.3'] = dict(dshell=0.3, nqlog=3200)
CASES['fx_dshell0'] = dict(dshell=0.0, nqlog=700)        # the one-bin shell (grid-dependent), for the record
EXPORT = {'fx_natw_nq6400': 'HFa', 'fx_natw_HFb': 'HFb'}
REFPAR = dict(name='HFa', H_GeV=C.H_REF, mnu=C.MNU_REF, t0=0.19, t_coh=3.0)
SEED = 'seedref'        # copy of the coherent-stage checkpoint (reference parameters, t_coh = 3) left by the interrupted run


def run(tag):
    kw = dict(REFPAR, t_end=1e8, lam_group='SU2'); kw.update(CASES[tag])
    ck, seed = C.SP / ('ckpt_%s.npz' % tag), C.SP / ('ckpt_%s.npz' % SEED)
    same = all(kw[k] == REFPAR[k] for k in REFPAR)
    if same and not ck.exists() and seed.exists():
        shutil.copy(seed, ck)          # the coherent stage does not depend on the dephased settings (dshell, nqlog)
    res = R.run_history(kw.pop('name'), tag, kw.pop('H_GeV'), kw.pop('mnu'), kw.pop('t0'), kw.pop('t_coh'), kw.pop('t_end'),
                        kw.pop('lam_group'), **kw)
    res['info']['lcut'] = float(os.environ.get('EK_LCUT', '40'))
    C.dump(C.SP / ('run_%s.json' % tag), R.jsonable(res, keep_rec=True))
    if tag in EXPORT and os.environ.get('EK_LCUT', '40') == '40':
        nm = EXPORT[tag]
        note = ('kinetic history from t0 = 0.19 H_i^-1 (homogeneous bench state at the first chi turning point); coherent mode '
                'stage to t_coh = 3 H_i^-1, dephased rate equations after, perturbative shell with the self-consistent natural line (omega-Lorentzian over the kinematic range)'
                'width. Lattice splice: NOTES.md section 6 (switch at t_sw in [0.45, 3] H_i^-1 by moments with ek_run.py --state, '
                'or splice this file at t >= t_sw when the lattice fractions agree with the history rows).')
        R.export_fluid(res, nm, C.RES / ('fluid_%s.json' % nm), note)
    return res


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        for tag in sys.argv[2:]:
            t1 = time.time()
            run(tag)
            print(tag, 'finished %.0fs' % (time.time() - t1), flush=True)
