#!/usr/bin/env python3
"""Escape-momentum sources for a range of imposed source epochs (Claude review, track `conversion`).

Monte Carlo (transport_mc.monte_carlo) where Gamma/H is moderate; local-equilibrium
drift where Gamma/H > 1500 (escape after a universal dln a = 0.0961).  Rate
table: quoted bb211d1 table with my log-log interpolation (and my own dense table
if 'mine' is passed).  Saves a_esc samples and metadata for class_screen.py.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import transport_mc as T


def main():
    table = sys.argv[1] if len(sys.argv) > 1 else 'codex'
    rate = T.codex_table() if table == 'codex' else T.my_table()
    rows = []
    for z in [1e5, 3e5, 1e6, 3e6, 1e7, 3e7, 1e8, 1e9]:
        af = 1 / (1 + z)
        g0 = float(T.lam_minus(rate, af, T.U0))
        fe = T.first_event(rate, af)
        if g0 < 1500:
            r = T.monte_carlo(rate, af, n=40000, ds_fac=0.004, seed=11, max_dlna=12.0)
            ae = r['a_esc'][np.isfinite(r['a_esc'])]
            meta = dict(method='MC', escaped_fraction=float(ae.size / r['n']),
                        anti_remaining_at_stop=float(r['anti_mask'].mean()), a_stop=r['a_end'],
                        forward=r['forward'], inverse=r['inverse'])
        else:
            eq = T.equilibrium_drift(rate, af)
            ae = np.full(40000, eq['a_escape'])
            meta = dict(method='equilibrium', escaped_fraction=1.0, anti_remaining_at_stop=0.0,
                        forward=eq['forward_events_per_particle'])
        u = T.UESC * ae
        grid = np.geomspace(ae.min() * 0.999, ae.max() * 1.001, 120)
        ht = np.array([T.horizon(x, T.UESC * x, 1.0) for x in grid])
        meta.update(z_form=z, a_form=af, Gamma_over_H_initial=g0, first_event_total_depth=fe['optical_depth_to_today'],
                    never_converted_today=fe['never_converted_today'],
                    z_escape_median=float(1 / np.median(ae) - 1), q_over_m_mean=float(u.mean()),
                    q_over_m_rel_spread=float(u.std() / u.mean()),
                    mean_horizon_today_Mpc=float(np.mean(np.interp(ae, grid, ht))),
                    daughter_horizon_today_Mpc=T.horizon(af, af * T.U0, 1.0))
        tag = f'epoch_{table}_z{z:.0e}'.replace('+', '')
        np.save(T.SCR / f'mc_aesc_{tag}.npy', ae)
        (T.SCR / f'mc_aesc_{tag}.json').write_text(json.dumps(meta))
        rows.append(meta)
        print(meta, flush=True)
    (T.RES / f'epoch_scan_{table}.json').write_text(json.dumps(rows, indent=1) + '\n')


if __name__ == '__main__':
    main()
