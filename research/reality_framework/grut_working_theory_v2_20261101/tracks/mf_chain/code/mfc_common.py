#!/usr/bin/env python3
"""mf_chain (Claude, multi-field inflaton pass; GitHub issue #13 context): common module.

REUSE, NOT REWRITE.  The M1r chain engine and its drivers are mf_ideal's (tracks/mf_ideal/code: rpmfi.c MODE 5,
vmodes.c, mfi_common.py, release_mfi.py, vectors_mfi.py, chain_mfi.py, bench_run.py).  They are imported READ-ONLY and
unchanged; this module only redirects their output paths, which are module globals there, to this track
(results/engine/...) and to this track's scratch area.  Nothing in tracks/mf_ideal is written: Python runs with -B,
the binaries are compiled from mf_ideal's sources into OUR scratch (MFI_SCRATCH is set before the import).

Every history, rescaling rule, ledger entry and criterion used by this track is a model choice; none follows from closure.
"""
import os, sys, json, math
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRK = HERE.parents[1]                                    # .../tracks
TRACK = HERE.parent                                      # .../tracks/mf_chain
RES = TRACK / 'results'
SP = Path(os.environ.get('MFC_SCRATCH',
                         '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_chain'))
os.environ['MFI_SCRATCH'] = str(SP / 'mfi')             # mf_ideal's scratch root -> ours (read at import of mfi_common)
os.environ.setdefault('OMP_NUM_THREADS', '1')
MFI = TRK / 'mf_ideal' / 'code'
sys.path.insert(0, str(MFI))

import mfi_common as mc                                  # noqa: E402  (read-only)
import release_mfi as rm                                 # noqa: E402
import vectors_mfi as vm                                 # noqa: E402
import chain_mfi as cm                                   # noqa: E402

ENG = RES / 'engine'                                     # mf_ideal-format outputs of our runs
mc.RES = ENG
rm.OUT = ENG / 'release'
vm.OUT = ENG / 'vectors'
cm.OUT = ENG / 'chain'
assert str(mc.SP).startswith(str(SP)), mc.SP
for d in (ENG / 'release', ENG / 'vectors', ENG / 'chain', ENG / 'bench', SP / 'mfi' / 'runs', SP / 'hist'):
    d.mkdir(parents=True, exist_ok=True)

MP_GEV = 2.435e18
H_REF = 1.30e11
BENCH_DIR = TRK / 'mf_hybrid' / 'results'
SKC = TRK / 'mf_hybrid' / 'skeptic' / 'code'
SKR = TRK / 'mf_hybrid' / 'skeptic' / 'results'


def dump(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    json.dump(obj, open(path, 'w'), indent=1, default=float)


def load_bench(name_or_path):
    p = Path(name_or_path)
    if not p.exists():
        p = BENCH_DIR / ('bench_%s.json' % name_or_path)
    return json.load(open(p))
