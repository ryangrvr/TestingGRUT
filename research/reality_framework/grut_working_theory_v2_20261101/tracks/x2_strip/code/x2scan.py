#!/usr/bin/env python3
"""x2_strip (Claude): run e2e_chain's scan.py UNCHANGED on lf_chain's LF histories, writing to THIS track's results/.
Usage: LFC_SCRATCH=<scratch> python3 -B x2scan.py BENCH_TAG BASE c0:delta ..."""
import sys, os, runpy
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'lf_chain' / 'code'))
import lfc_common as L  # noqa: E402  (read-only; sets E2C_SCRATCH from LFC_SCRATCH)
NEW = HERE.parent / 'results'
L.RES = NEW; L.E.RES = NEW
if __name__ == '__main__':
    bench = str(L.HIST / ('bench_%s.json' % sys.argv[1]))
    sys.argv = ['scan.py', bench, sys.argv[2]] + sys.argv[3:]
    runpy.run_path(str(L.E.TRK / 'e2e_chain' / 'code' / 'scan.py'), run_name='__main__')
