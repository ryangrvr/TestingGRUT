#!/usr/bin/env python3
"""mf_chain SKEPTIC (key mf_chain_skeptic; Claude, multi-field inflaton pass, issue #13 context): common module.

Where the skeptic needs the FULL chain (vectors, ledger, p map) on a modified history, it calls mf_ideal's engine and
drivers READ-ONLY (as mf_chain does), with every output redirected to the skeptic's scratch area; binaries are compiled from
mf_ideal's sources into that scratch.  The independent checks (release, ledger, S/Q, calibration, leak, isocurvature
amplitude) are in sk_release.py and sk_ledger.py and use none of this.
Every history and parameter is a model choice; nothing follows from closure.
"""
import os, sys, json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SKD = HERE.parent                                        # tracks/mf_chain/skeptic
TRK = HERE.parents[2]                                    # .../tracks
RES = SKD / 'results'
SP = Path(os.environ.get('SKC_SCRATCH',
                         '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_chain_skeptic'))
os.environ['MFI_SCRATCH'] = str(SP / 'mfi')
os.environ.setdefault('OMP_NUM_THREADS', '1')
MFI = TRK / 'mf_ideal' / 'code'
sys.path.insert(0, str(MFI))
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
for d in (ENG / 'release', ENG / 'vectors', ENG / 'chain', ENG / 'bench', SP / 'mfi' / 'runs', SP / 'hist', RES):
    d.mkdir(parents=True, exist_ok=True)
MFC_SP = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_chain')
MFC_RES = TRK / 'mf_chain' / 'results'


def dump(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    json.dump(obj, open(path, 'w'), indent=1, default=float)
