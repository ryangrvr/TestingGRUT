#!/usr/bin/env python3
"""lf_chain step 1 (Claude): final histories on the LF closure.
  (a) e2e_chain's handover.py (read-only, run unchanged via runpy as __main__, after lfc_common's redirection): the
      lattice -> kinetic moments hand-over at t_sw = 3 with its energy-conserving fixes (degenerate sphere where nbar > 1/2;
      dephased network from t_sw), N split over species in the kinetic reference's own ratios (handover.py default).
  (b) e2e_chain's splice.build (read-only) on that kinetic export.
  (c) the e2e_chain skeptic's correction 8 (sk_engine.py make_hist mode 'anh', re-typed): for t > t_sw the HFa kinetic-block
      spectator mass is m_s^2 = 9 (1 - 2x), x = rho_chi/(m_chi^2 v^2) [D there; effect on S/Q <= 7e-4 there].
Usage: python3 -B build_hist.py LATKEY TAG [--H 1.3e11]
Labels: [C] given the model choices.
"""
import sys, json, math, runpy, argparse
import numpy as np
import lfc_common as L

E, S = L.E, L.S


def col(h, k):
    return np.array([np.nan if v is None else float(v) for v in h[k]])


def anh(path):
    d = json.load(open(path))
    h = d['history']; p = d['params']
    t, H, fchi, ms = col(h, 't'), col(h, 'H'), col(h, 'f_chi'), col(h, 'm_s')
    tsw = float(d['reheating'].get('t_sw', 3.0))
    mP = E.MP_GEV / float(p['H_i_GeV'])
    v = float(p['v_MP']) * mP
    x = fchi * 3 * mP ** 2 * H ** 2 / (float(p['m_chi_over_Hi']) ** 2 * v ** 2)
    m = t > tsw + 1e-9
    new = ms.copy(); new[m] = np.sqrt(9.0 * (1 - 2 * x[m]))
    j = int(np.argmax(m))
    h['m_s'] = [None if not np.isfinite(q) else float('%.10g' % q) for q in new]
    d['lf_chain_mod'] = dict(mode='anh (e2e_chain skeptic correction 8)', ms2_before=float(ms[j - 1] ** 2),
                             ms2_after_old=float(ms[j] ** 2), ms2_after_new=float(new[j] ** 2), x_at_tsw=float(x[j]))
    json.dump(d, open(path, 'w'), indent=1, default=float)
    return d['lf_chain_mod']


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('latkey'); ap.add_argument('tag'); ap.add_argument('--H', type=float, default=E.H_REF)
    a = ap.parse_args()
    ho = 'ho' + a.tag
    fluid = E.SP / 'kin' / ('fluid_%s.json' % ho)
    if not fluid.exists():
        sys.argv = ['handover.py', 'HFa', a.latkey, ho, '--H', repr(a.H)]
        runpy.run_path(str(E.TRK / 'e2e_chain' / 'code' / 'handover.py'), run_name='__main__')
    o = S.build('HFa', a.tag, a.latkey, str(fluid), a.H)
    mod = anh(o['path'])
    chk = L.load(L.RES / 'splice_checks.json', {})
    sp = o['splice']
    chk[a.tag] = dict(lattice=a.latkey, H_design=a.H, K=o['reheating']['K'], H_jump_at_tsw=sp['H_jump'], lattice_state_tsw=sp['lattice'],
                      kinetic_first_row=sp['kinetic'], lattice_energy_budget=sp['lattice_budget'],
                      friedmann_max_rel_err_lattice_block=sp['friedmann_max_rel_err_off_joint'], K_drift_last_decade=sp['K_drift_last_decade'],
                      anh=mod, bench=o['path'])
    L.dump(L.RES / 'splice_checks.json', chk)
    print(a.tag, 'K %.4f' % o['reheating']['K'], 'H_jump %.2e' % sp['H_jump'], 'ms2 step', mod)
