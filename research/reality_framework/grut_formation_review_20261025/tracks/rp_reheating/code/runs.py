#!/usr/bin/env python3
"""rp_reheating track (Claude, own code): the run campaign.
  reproduce   : Codex's six production runs, central tight, long (tau 500), RD reference, xi=0 control, phase +-1e-8
  gamma_scan  : b = 2, 10 x Gamma/H_i = 0.03 ... 100 (does a faster reheating restore the yield without retuning?)
  b_scan      : Gamma/H_i = 0.1, m_I/H_i = 1 ... 25 (does the oscillating R excite Y resonantly near 2 m_I = omega_YL?)
  models      : fluid inflaton with xi<R> (mode 1) and with -cH rho_I/3 (mode 2, hw_verify's model) vs the oscillator
Usage: python3 -B runs.py ../results <group> [<group> ...]      (~1-3 min per group on 4 cores)
"""
import json, sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import rp


def job(spec):
    tag, kw = spec
    try:
        res, S = rp.run(tag, **kw)
    except Exception as exc:
        return dict(tag=tag, error=repr(exc), kw={k: (v if not isinstance(v, np.ndarray) else v.tolist()) for k, v in kw.items()})
    # compact early profile (tau <= 12) and sparse late profile for the notes
    keep = ['tau', 'a', 'H', 'R', 'NX', 'EX', 'eps', 'NYH', 'NYL', 'Nsig', 'absX', 'absY', 'radfrac', 'w_I', 'm2curv', 'rho_r']
    idx = [rp.COLS.index(c) for c in keep]
    early = S[S[:, 0] <= 12.0][::max(1, int(len(S[S[:, 0] <= 12.0]) / 120))][:, idx]
    late = S[np.linspace(0, len(S) - 1, 60).astype(int)][:, idx]
    res['profile_columns'] = keep
    res['profile_early'] = early.tolist(); res['profile_late'] = late.tolist()
    res['NYL_max_over_run'] = float(np.nanmax(S[:, rp.COLS.index('NYL')]))
    res['NYH_max_over_run'] = float(np.nanmax(S[:, rp.COLS.index('NYH')]))
    print('%-22s b %-5g g %-6g eps %.5f NX %.1f NYH %.3e NYL %.4e Y/YQ %.4f (ratio %.3f / D %.3f) a %.2f Rr %.2f chk %s %.1fs'
          % (tag, res['b'], res['gamma'], res['late']['eps'], res['late']['NX'], res['late']['NYH'], res['late']['NYL'],
             res['yield_over_YQ'], res['charge_ratio_to_instant'], res['dilution'], res['final_a'], res['Rr_end'],
             res['checks'], res['seconds']), flush=True)
    return res


def specs(group):
    if group == 'reproduce':
        s = [('b%g_g%g' % (b, g), dict(mode=0, b=b, gam=g)) for b in (2.0, 10.0) for g in (0.03, 0.1, 0.3)]
        s += [('central_tight', dict(mode=0, b=2.0, gam=0.1, rtol=1e-12, atol=1e-15)),
              ('central_loose', dict(mode=0, b=2.0, gam=0.1, rtol=1e-8, atol=1e-11)),
              ('long_500', dict(mode=0, b=2.0, gam=0.1, end=500.0, nsamp=4000)),
              ('long_1500', dict(mode=0, b=2.0, gam=0.1, end=1500.0, nsamp=6000)),
              ('RD_reference', dict(mode=3, b=2.0, gam=0.1, end=799.5, eta2=0.0, nsamp=4000)),
              ('phase_plus', dict(mode=0, b=2.0, gam=0.1, delta=1e-8)),
              ('phase_minus', dict(mode=0, b=2.0, gam=0.1, delta=-1e-8))]
        z0, _ = rp.prepare(2.0)      # Codex's xi_zero_control keeps the xi = -5/6 preparation and runs with xi = 0
        s += [('xi_zero_control', dict(mode=0, b=2.0, gam=0.1, xi=0.0, z0=z0))]
        return s
    if group == 'gamma_scan':
        return [('gs_b%g_g%g' % (b, g), dict(mode=0, b=b, gam=g, end=max(25.0 / g, 80.0)))
                for b in (2.0, 10.0) for g in (0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)]
    if group == 'b_scan':
        bs = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 10.5, 11.0, 11.5, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0,
              18.0, 19.0, 20.0, 22.0, 25.0]
        return [('bs_b%g' % b, dict(mode=0, b=b, gam=0.1, end=250.0)) for b in bs]
    if group == 'models':
        s = []
        for g in (0.03, 0.1, 0.3, 1.0, 3.0, 10.0):
            s += [('fluidxiR_g%g' % g, dict(mode=1, b=2.0, gam=g, end=max(25.0 / g, 80.0))),
                  ('fluidcH_g%g' % g, dict(mode=2, b=2.0, gam=g, end=max(25.0 / g, 80.0)))]
        return s
    if group == 'extra':
        z0, _ = rp.prepare(2.0)      # de Sitter (xi R = -10) preparation, then released into each background
        s = [('RD_reference_dS', dict(mode=3, b=2.0, gam=0.1, end=799.5, eta2=0.0, nsamp=4000, z0=z0))]
        s += [('fluidxiR_dS_g%g' % g, dict(mode=1, b=2.0, gam=g, end=max(25.0 / g, 80.0), z0=z0)) for g in (0.03, 0.1, 0.3, 1.0, 3.0, 10.0)]
        # yield restoration grid with perturbative widths Gamma <= 0.3 m_I
        grid = [(3, 0.5), (3, 0.9), (5, 0.5), (5, 1.0), (5, 1.5), (10, 2.0), (15, 1.5), (15, 3.0), (15, 4.5), (20, 2.0), (20, 3.0),
                (20, 6.0), (30, 3.0), (30, 6.0), (30, 9.0), (50, 5.0), (50, 10.0), (50, 15.0)]
        s += [('yr_b%g_g%g' % (b, g), dict(mode=0, b=float(b), gam=g, end=max(25.0 / g, 80.0))) for b, g in grid]
        return s
    raise SystemExit('unknown group ' + group)


if __name__ == '__main__':
    out = Path(sys.argv[1]); out.mkdir(parents=True, exist_ok=True)
    rp.build()
    for group in sys.argv[2:]:
        sp = specs(group)
        rows = []
        with Pool(4) as p:
            for r in p.imap_unordered(job, sp):
                rows.append(r)
                (out / (group + '.json')).write_text(json.dumps(dict(method=__doc__, group=group, rows=sorted(rows, key=lambda x: x['tag'])),
                                                                indent=1, default=float))
        print('group', group, 'rows', len(rows), flush=True)
