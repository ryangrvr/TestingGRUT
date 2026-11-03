#!/usr/bin/env python3
"""e2e_chain item 2 driver (Claude): holder scan on one spliced history.

Usage: python3 -B scan.py BENCH.json BASETAG c0[:delta_deg] [c0[:delta_deg] ...]
  - BASETAG at c0 = 1, delta = 0 is run first with mf_ideal's bench_run (MODE 5 table, --R raw, which here is the lattice's
    cycle-mean R and the kinetic cycle-mean R; spectator evolved with the table's m_s^2), unless its chain file exists.
  - each further point: chain.release_c0 on the same table, then chain.full (corrected relic, leak with c_L, grav).
Output: results/scan_<BASETAG>.json (merged; one entry per point).  [C] given the model choices.
"""
import json, sys, time
import e2c_common as E
import chain as CH
import bench_run as br                                    # mf_ideal, read-only

if __name__ == '__main__':
    bench, base = sys.argv[1], sys.argv[2]
    pts = list(sys.argv[3:])
    E.mc.build()
    outp = E.RES / ('scan_%s.json' % base)
    allr = E.load(outp, {})
    allr.setdefault('bench', bench)
    d = json.load(open(bench))
    allr['meta'] = dict(H_i_GeV=d['params']['H_i_GeV'], reheating=d.get('reheating'), splice=d.get('splice'))
    runs = allr.setdefault('runs', {})
    if not (E.ENG / 'chain' / ('chain_%s.json' % base)).exists():
        t0 = time.time()
        sys.argv = ['bench_run.py', bench, '--tag', base, '--R', 'raw', '--tlong', '1e4']
        br.main()
        print('base %.0f s' % (time.time() - t0), flush=True)
    o = CH.full(base, bench, 1.0, 0.0)
    EX1 = o['E_X']
    o['E_X_c1'] = EX1
    runs['1.000:+0'] = o
    print(CH.line(o), flush=True)
    E.dump(outp, allr)
    for p in pts:
        c0, dd = (float(x) for x in (p.split(':') + ['0'])[:2])
        key = '%.3f:%+.0f' % (c0, dd)
        tag = '%s_C%.3f_d%+.0f' % (base, c0, dd)
        t0 = time.time()
        if not (E.ENG / 'chain' / ('chain_%s.json' % tag)).exists():
            CH.release_c0(base, tag, c0, dd)
        o = CH.full(tag, bench, c0, dd, E_X_c1=EX1)
        o['seconds'] = time.time() - t0
        runs[key] = o
        print(CH.line(o), '%.0fs' % o['seconds'], flush=True)
        E.dump(outp, allr)
