#!/usr/bin/env python3
"""lf_chain step 2 driver (Claude): e2e_chain's scan.py run UNCHANGED (runpy, as __main__) after lfc_common's redirection
(scratch -> LFC_SCRATCH, results -> this track's results/scan_<BASE>.json).  Points already in the scan file are skipped.
Usage: LFC_SCRATCH=<scratch> python3 -B scan_run.py BENCH_TAG BASE c0[:delta] ...
"""
import sys, runpy
import lfc_common as L

if __name__ == '__main__':
    bench = str(L.HIST / ('bench_%s.json' % sys.argv[1]))
    base = sys.argv[2]
    have = L.load(L.RES / ('scan_%s.json' % base), {}).get('runs', {})
    pts = []
    for p in sys.argv[3:]:
        c0, dd = (float(x) for x in (p.split(':') + ['0'])[:2])
        if '%.3f:%+.0f' % (c0, dd) not in have:
            pts.append(p)
    sys.argv = ['scan.py', bench, base] + pts
    runpy.run_path(str(L.E.TRK / 'e2e_chain' / 'code' / 'scan.py'), run_name='__main__')
