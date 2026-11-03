#!/usr/bin/env python3
"""e2e_chain (Claude, end-to-end pass "M2"; GitHub issue #13 context): common module.

REUSE, NOT REWRITE.  The M1r chain (release -> vectors -> dq_heat ledger -> closed-form S/Q -> e_crit -> H_i) is
mf_ideal's engine and drivers (tracks/mf_ideal/code: rpmfi.c MODE 5, vmodes.c, mfi_common.py, release_mfi.py,
vectors_mfi.py, chain_mfi.py, bench_run.py), imported READ-ONLY and unchanged, exactly as mf_chain / mf_critic did.
Their output paths (module globals) are redirected to this track's scratch; the binaries are compiled from mf_ideal's
sources into that scratch.  The relic step uses e2e_relic's corrected p-map and relic() (read-only import) in place of
chain_mfi's inherited map (see chain.py).  Nothing in another track is written.  Run Python with -B, one thread.

Every history, splice rule, holder amplitude and phase, ledger entry and criterion is a model choice; none follows
from closure.
"""
import os, sys, json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRACK = HERE.parent                                      # .../tracks/e2e_chain
TRK = HERE.parents[1]                                    # .../tracks
RES = TRACK / 'results'
SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2')
SP = Path(os.environ.get('E2C_SCRATCH', str(SCR / 'e2e_chain')))
os.environ['MFI_SCRATCH'] = str(SP / 'mfi')             # read by mfi_common at import
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.path.insert(0, str(TRK / 'mf_ideal' / 'code'))
sys.path.insert(0, str(TRK / 'e2e_relic' / 'code'))

import mfi_common as mc                                  # noqa: E402  (read-only)
import release_mfi as rm                                 # noqa: E402
import vectors_mfi as vm                                 # noqa: E402
import chain_mfi as cm                                   # noqa: E402

ENG = SP / 'engine'
mc.RES = ENG
rm.OUT = ENG / 'release'
vm.OUT = ENG / 'vectors'
cm.OUT = ENG / 'chain'
assert str(mc.SP).startswith(str(SP)), mc.SP
for d in (ENG / 'release', ENG / 'vectors', ENG / 'chain', ENG / 'bench', SP / 'mfi' / 'runs', SP / 'hist'):
    d.mkdir(parents=True, exist_ok=True)

KIN_SP = SCR / 'e2e_kinetic'                             # e2e_kinetic's scratch: dephased records, read-only
MP_GEV = 2.435e18
H_REF = 1.30e11
A_S = 2.10e-9


def dump(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    json.dump(obj, open(path, 'w'), indent=1, default=float)


def load(path, default=None):
    path = Path(path)
    return json.load(open(path)) if path.exists() else default
