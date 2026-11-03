#!/usr/bin/env python3
"""e2e_chain_skeptic (Claude): own driver for engine re-runs on (modified copies of) e2e_chain's final histories.

Subcommands
  hist SRC_BENCH OUT_TAG MODE   write a modified copy of a final history into this track's scratch (hist/bench_<OUT_TAG>.json)
      MODE = anh : for t > t_sw (the kinetic block) replace e2e_chain/e2e_kinetic's HFa spectator mass
                   m_s^2 = g_s^2 (v^2 + rho_chi/m_chi^2)   [a harmonic oscillation about v, <dchi> = 0]
                   by the leading anharmonic value of the same oscillation in V = lam/4 (chi^2 - v^2)^2:
                   <chi - v> = -(3/2v) <dchi^2>  ->  <chi^2> = v^2 - 2 rho_chi/m_chi^2,  m_s^2 = 9 (1 - 2x),
                   x = rho_chi/(m_chi^2 v^2)  [D, below].  The lattice block (t <= t_sw, measured box mean) is untouched.
      MODE = vac : m_s^2 = 9 exactly for t > t_sw (bracket).
  run BENCH BASETAG c0[:delta] ...   the e2e_chain chain on BENCH: base (c0 = 1) through mf_ideal's bench_run (MODE 5
      table, --R raw, --tlong 1e4, exactly as e2e_chain scan.py), then chain.release_c0 + chain.full for each point
      (e2e_chain code, read-only).  Results merged into results/sk_engine.json under BASETAG.

Derivation of the anharmonic shift [D]: chi = v + y, V = lam v^2 y^2 + lam v y^3 + lam y^4/4, m^2 = 2 lam v^2;
the time average of the equation of motion gives m^2 <y> + 3 lam v <y^2> = 0, so <y> = -(3/2v)<y^2>, and
<chi^2> = v^2 + 2v<y> + <y^2> = v^2 - 2<y^2>, with <y^2> = rho_chi/m^2 for the (nearly harmonic) oscillation.
Labels: [C] engine numbers given the model choices; the modification is a model choice [A].
"""
import json, math, sys, time, copy
import numpy as np
import sk_common as S

E = S.E
import chain as CH                                          # e2e_chain, read-only
import bench_run as br                                      # mf_ideal, read-only


def col(h, k):
    return np.array([np.nan if v is None else float(v) for v in h[k]])


def make_hist(src, out_tag, mode):
    d = json.load(open(src))
    h = d['history']; p = d['params']
    t, H, fchi, ms = col(h, 't'), col(h, 'H'), col(h, 'f_chi'), col(h, 'm_s')
    tsw = float(d['reheating'].get('t_sw', 3.0))
    mP = E.MP_GEV / float(p['H_i_GeV'])
    v = float(p['v_MP']) * mP
    mchi = float(p['m_chi_over_Hi'])
    rho = 3 * mP ** 2 * H ** 2
    x = fchi * rho / (mchi ** 2 * v ** 2)
    m = t > tsw + 1e-9
    new = ms.copy()
    if mode == 'anh':
        new[m] = np.sqrt(9.0 * (1 - 2 * x[m]))
    elif mode == 'vac':
        new[m] = 3.0
    else:
        raise SystemExit('mode?')
    h['m_s'] = [None if not np.isfinite(q) else float('%.10g' % q) for q in new]
    j = int(np.argmax(m))
    d['name'] = out_tag
    d['label'] = d['label'] + ' | e2e_chain_skeptic: m_s after t_sw replaced (%s)' % mode
    d['skeptic_mod'] = dict(mode=mode, t_sw=tsw, ms2_before=float(ms[j - 1] ** 2), ms2_after_old=float(ms[j] ** 2),
                            ms2_after_new=float(new[j] ** 2), x_at_tsw=float(x[j]))
    out = S.SKSP / 'hist' / ('bench_%s.json' % out_tag)
    out.parent.mkdir(parents=True, exist_ok=True)
    json.dump(d, open(out, 'w'), indent=1, default=float)
    print(json.dumps(d['skeptic_mod']))
    return out


def run(bench, base, pts):
    E.mc.build()
    allr = S.load(S.RES / 'sk_engine.json', {})
    blk = allr.setdefault(base, dict(bench=str(bench), runs={}))
    d = json.load(open(bench))
    blk['meta'] = dict(H_i_GeV=d['params']['H_i_GeV'], K=d['reheating']['K'], skeptic_mod=d.get('skeptic_mod'))
    if not (E.ENG / 'chain' / ('chain_%s.json' % base)).exists():
        t0 = time.time()
        sys.argv = ['bench_run.py', str(bench), '--tag', base, '--R', 'raw', '--tlong', '1e4']
        br.main()
        print('base %.0f s' % (time.time() - t0), flush=True)
    o = CH.full(base, str(bench), 1.0, 0.0)
    EX1 = o['E_X']
    blk['runs']['1.000:+0'] = summ(o)
    print(CH.line(o), flush=True)
    S.dump('sk_engine.json', allr)
    for p in pts:
        c0, dd = (float(q) for q in (p.split(':') + ['0'])[:2])
        key = '%.3f:%+.1f' % (c0, dd)
        tag = '%s_C%.3f_d%+.1f' % (base, c0, dd)
        t0 = time.time()
        if not (E.ENG / 'chain' / ('chain_%s.json' % tag)).exists():
            CH.release_c0(base, tag, c0, dd)
        o = CH.full(tag, str(bench), c0, dd, E_X_c1=EX1)
        r = summ(o); r['seconds'] = time.time() - t0
        blk['runs'][key] = r
        print(CH.line(o), '%.0fs' % r['seconds'], flush=True)
        S.dump('sk_engine.json', allr)


def summ(o):
    b = o['by_e']['cL0_vanish_at_MY']
    keep = {k: o[k] for k in ('c0', 'delta_deg', 'mu_eff', 'abs_mu_eff', 'arg_mu_eff_deg', 'H_i_GeV', 'K', 'eps', 'E_over_Ec', 'N_X',
                              'E_X', 'N_YL', 'N_YH', 'N_sigma', 'S_per_Q_e_inf', 'Y_over_YQ_at_HM1')}
    keep['checks'] = o['checks']
    keep['e_crit_bulk'] = b['e_crit_bulk']; keep['e_crit_Q2e5'] = b['e_crit_Q2e5']
    for lab in ('e1e-6', 'e3e-6'):
        r = b[lab]
        keep[lab] = {k: r.get(k) for k in ('R_vec', 'R_Y', 'R_casc', 'dR_grav', 'R_total', 'S_per_Q', 'branch', 'M', 'H_i_req_GeV',
                                           'H_req_over_H_i')}
    keep['spectator'] = o['spectator']
    return keep


if __name__ == '__main__':
    cmd = sys.argv[1]
    if cmd == 'hist':
        make_hist(sys.argv[2], sys.argv[3], sys.argv[4])
    elif cmd == 'run':
        run(sys.argv[2], sys.argv[3], sys.argv[4:])
