#!/usr/bin/env python3
"""mf_critic (Claude, final verdict of the multi-field inflaton pass; GitHub issue #13 context): common module.

REUSE, NOT REWRITE.  Where this track needs the full M1r chain (release -> vectors -> ledger -> S/Q, e_crit, relic, H_i),
it imports mf_ideal's engine drivers READ-ONLY (tracks/mf_ideal/code: rpmfi.c MODE 5, vmodes.c, mfi_common.py,
release_mfi.py, vectors_mfi.py, chain_mfi.py), exactly as mf_chain and its skeptic did.  Every output of those drivers is
redirected (module globals) to this track's scratch area; the binaries are compiled from mf_ideal's sources into that
scratch.  Nothing in another track is written.  Run Python with -B and one thread.

Every history, holder amplitude, ledger entry and criterion used here is a model choice; none follows from closure.
"""
import os, sys, json
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRACK = HERE.parent                                      # .../tracks/mf_critic
TRK = HERE.parents[1]                                    # .../tracks
RES = TRACK / 'results'
SP = Path(os.environ.get('MFCR_SCRATCH',
                         '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2/mf_critic'))
os.environ['MFI_SCRATCH'] = str(SP / 'mfi')             # read by mfi_common at import
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.path.insert(0, str(TRK / 'mf_ideal' / 'code'))

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
for d in (ENG / 'release', ENG / 'vectors', ENG / 'chain', ENG / 'bench', SP / 'mfi' / 'runs', SP / 'hist', SP / 'tab'):
    d.mkdir(parents=True, exist_ok=True)

SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2')
M1R_TABLE = SCR / 'mf_ideal' / 'table_m1r_dt0.001.txt'   # mf_ideal's R2 bridge table (M1r's own (3, 0.5) background), read-only
MP_GEV = 2.435e18
A_S = 2.10e-9


def dump(path, obj):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    json.dump(obj, open(path, 'w'), indent=1, default=float)


def short(o):
    """compact chain summary."""
    keep = ('K', 'N_X', 'E_X', 'eps', 'E_over_Ec', 'N_YL', 'N_YH', 'N_sigma', 'S_per_Q_e_inf', 'R_Y', 'e_crit_bulk', 'e_crit_Q2e5')
    out = {k: o.get(k) for k in keep}
    for lab in ('e1e-6', 'e3e-6'):
        a = o[lab]
        out[lab] = {k: a.get(k) for k in ('R_vec', 'R_Y', 'R_casc', 'W', 'S_per_Q', 'branch', 'M', 'asym_share', 'H_i_req_GeV', 'r')}
    return out
