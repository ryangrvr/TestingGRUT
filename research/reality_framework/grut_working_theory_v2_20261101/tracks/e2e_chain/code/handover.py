#!/usr/bin/env python3
"""e2e_chain item 1b (Claude): moments hand-over from e2e_lattice's coherent stage to e2e_kinetic's network at t_sw
(the kinetic track's preferred interface (a), its NOTES section 4), run with e2e_kinetic's own driver READ-ONLY
(code/ek_run.py export ... --state), its scratch redirected to ours (EK_SCRATCH) and its p-tables copied there.

State at t_sw (physical, H_i units) [model choices, A]:
  rho = 3 mP^2 H_lat^2; rho_sig = f_sig rho; rho_chi = (f_chi + f_grad) rho (the lattice's chi fragments are handed to
  the kinetic code's homogeneous chi condensate: the momentum spread of the fragments is dropped here and bracketed by
  the kinetic dshell = 0.3 run); SM: f_SM rho split into hard / thermal, and n_h, in the kinetic run's own ratios at
  t_sw; N: (f_rad - f_SM) rho split over the four species in the kinetic run's own ratios at t_sw, each with the
  kinetic run's own mean energy per particle (the lattice drain term carries total energy only).
Usage: python3 -B handover.py NAME LATKEY TAG [--H 1.3e11] [--kref fx_natw_nq6400] [--tsw 3]
Output: scratch kin/<TAG>_state.json, kin/fluid_<TAG>.json (the kinetic export, bench format), results/handover_<TAG>.json.
"""
import argparse, json, math, os, shutil, subprocess, sys
import numpy as np
import e2c_common as E
import splice as S

KIN = E.SP / 'kin'


def state(name, latkey, H_GeV, kref, t_sw=3.0):
    L = S.lattice_block(name, latkey)
    a = np.loadtxt(str(E.KIN_SP / ('%s_dep_rec.txt' % kref)), ndmin=2)
    r = {c: a[:, i] for i, c in enumerate(S.REC_COLS)}
    j = int(np.argmin(np.abs(r['t'] - t_sw)))
    mP = E.MP_GEV / H_GeV
    g = lambda k: float(np.interp(t_sw, L['t'], L[k]))
    H = g('H'); rho = 3 * mP ** 2 * H * H
    fsm, frad = g('f_SM'), g('f_rad')
    hard_share = r['f_hard'][j] / (r['f_hard'][j] + r['f_th'][j])
    rho_h = fsm * rho * hard_share
    st = dict(t=t_sw, a=g('a'), H=H, rho_sig=g('f_sig') * rho, rho_chi=(g('f_chi') + g('f_grad')) * rho,
              rho_h=rho_h, n_h=float(r['n_h'][j] * rho_h / (r['f_hard'][j] * r['rho'][j])), rho_th=fsm * rho * (1 - hard_share))
    fN = np.array([r['fN%d' % k][j] for k in (1, 2, 3, 4)])
    EkN = fN * r['rho'][j]
    nkN = np.array([r['nN%d' % k][j] for k in (1, 2, 3, 4)])
    EN = (frad - fsm) * rho * fN / fN.sum()
    nN = EN / np.maximum(EkN / np.maximum(nkN, 1e-300), 1e-300)
    # energy conservation: ek_scan.occ_from_moments caps the uniform occupation at 1/2 and would drop the excess energy.
    # Where the kinetic mean energy would need nbar > 1/2, hand over a degenerate sphere (nbar = 1/2) carrying the
    # lattice energy instead: k_max = (8 pi^2 E)^(1/4), n = k_max^3/(6 pi^2) (relativistic, both helicities) [A].
    over = []
    for k in range(4):
        if EN[k] <= 0 or nN[k] <= 0:
            continue
        km = 4 * EN[k] / (3 * nN[k])
        nbar = nN[k] * 3 * math.pi ** 2 / km ** 3
        if nbar > 0.5:
            km = (8 * math.pi ** 2 * EN[k]) ** 0.25
            over.append(dict(species=k + 1, nbar_needed=float(nbar), k_max=float(km)))
            nN[k] = km ** 3 / (6 * math.pi ** 2)
    st['E_N'] = [float(x) for x in EN]
    st['n_N'] = [float(x) for x in nN]
    st['degenerate_handover'] = over
    st['_provenance'] = dict(lattice=latkey, kinetic_split_from=kref, t_kin_row=float(r['t'][j]), f_sig=g('f_sig'), f_chi=g('f_chi'),
                             f_grad=g('f_grad'), f_rad=frad, f_SM=fsm, rho_tot=rho, kin_row=dict(f_sig=float(r['f_sig'][j]), f_chi=float(r['f_chi'][j]),
                             f_N=float(r['f_N'][j]), f_SM=float(r['f_hard'][j] + r['f_th'][j])))
    return st


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('name'); ap.add_argument('latkey'); ap.add_argument('tag')
    ap.add_argument('--H', type=float, default=E.H_REF); ap.add_argument('--kref', default='fx_natw_nq6400')
    ap.add_argument('--tsw', type=float, default=3.0); ap.add_argument('--dshell', default=None)
    a = ap.parse_args()
    KIN.mkdir(parents=True, exist_ok=True)
    for f in E.KIN_SP.glob('ptable_*.npz'):
        if not (KIN / f.name).exists():
            shutil.copy(f, KIN / f.name)
    st = state(a.name, a.latkey, a.H, a.kref, a.tsw)
    sp = KIN / ('%s_state.json' % a.tag)
    E.dump(sp, st)
    out = KIN / ('fluid_%s.json' % a.tag)
    # in-process call of e2e_kinetic's driver (read-only), scratch redirected; the coherent mode stage is NOT re-entered
    # after the hand-over (t_coh = t_sw = 3, where the kinetic reference itself switches to the dephased network): its
    # coherent grid tops out at (1.3 A + 4 M/m + 3) m and would truncate the handed-over N energy (found: 19% / 75% of
    # the N energy lost for P0GM / P0 when it was re-entered) [A; the kinetic reference shows t_coh 3 -> 10 moves K 0.2%].
    os.environ['EK_SCRATCH'] = str(KIN)
    sys.path.insert(0, str(E.TRK / 'e2e_kinetic' / 'code'))
    import ek_common as KC
    import ek_run as KR
    import ek_scan as KS
    assert str(KC.SP) == str(KIN), KC.SP
    spec_, _ = KC.species(a.name, a.H, 0.05e-9)
    stt = {k: float(st[k]) for k in ('t', 'a', 'H', 'rho_sig', 'rho_chi', 'rho_h', 'n_h', 'rho_th')}
    occ0 = []
    for j_, s_ in enumerate(spec_):
        q, w = KC.grid_linear(stt['a'] * 3e4, 6000)
        occ0.append((q, w, KS.occ_from_moments(float(st['E_N'][j_]), float(st['n_N'][j_]), q, w, stt['a'], s_['M'])))
    ds = float(a.dshell) if a.dshell else -1.0e5
    r = KR.run_history(a.name, a.tag, a.H, 0.05e-9, stt['t'], stt['t'], 1e8, 'SU2', dshell=ds, nqlog=3200 if a.dshell else 6400,
                       state=stt, occ0=occ0)
    KR.export_fluid(r, a.name, out, 'e2e_chain moments hand-over from e2e_lattice %s at t_sw = %g; dephased network from t_sw' % (a.latkey, a.tsw))
    d = json.load(open(out))
    E.dump(E.RES / ('handover_%s.json' % a.tag), dict(state=st, summary=d['summary'], info={k: d['info'].get(k) for k in
                                                                                    ('t0', 't_coh', 'H_GeV', 'dshell', 'nqlog', 'seconds', 'rebin_check')}))
    print('K', d['summary']['K'], 'T_RH', d['summary']['T_RH_GeV'], 't_SM_half', d['summary']['t_SM_half'])
