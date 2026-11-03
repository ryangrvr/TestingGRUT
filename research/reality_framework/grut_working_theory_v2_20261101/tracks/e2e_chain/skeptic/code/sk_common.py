#!/usr/bin/env python3
"""e2e_chain_skeptic (Claude, end-to-end pass "M2"; GitHub issue #13 context): common module.

Adversarial checks of tracks e2e_chain and e2e_relic.  REUSE: e2e_chain's chain.py / e2c_common.py (and through them
mf_ideal's engine and drivers, e2e_relic's pmap / relic) are imported READ-ONLY and unchanged.  The only change is the
scratch redirection: E2C_SCRATCH is set to THIS track's scratch before e2c_common is imported, so every engine binary,
table, release, vector and chain file goes to  .../scratchpad/wt2/e2e_chain_skeptic/ .  Nothing in another track is
written.  e2e_chain's scan.py is NOT used (it writes into e2e_chain/results); sk_engine.py is our own driver.
Run Python with -B, one thread.

Every history modification, closure, holder amplitude and phase, interpolation rule and criterion used here is a model
choice; none follows from closure.  Labels: [C] computed here, [D] derived, [E] estimate, [A] argued, [Q] quoted.
"""
import os, sys, json
from pathlib import Path

HERE = Path(__file__).resolve().parent
SKTRACK = HERE.parent                                       # .../tracks/e2e_chain/skeptic
TRK = HERE.parents[2]                                       # .../tracks
RES = SKTRACK / 'results'
SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2')
SKSP = SCR / 'e2e_chain_skeptic'
os.environ['E2C_SCRATCH'] = str(SKSP)                       # read by e2e_chain's e2c_common at import
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.path.insert(0, str(TRK / 'e2e_chain' / 'code'))

import e2c_common as E                                      # noqa: E402  (read-only; scratch -> SKSP)
assert str(E.SP).startswith(str(SKSP)), E.SP
assert str(E.mc.SP).startswith(str(SKSP)), E.mc.SP

E2C_SP = SCR / 'e2e_chain'                                  # e2e_chain's scratch: final histories, read-only
E2C_RES = TRK / 'e2e_chain' / 'results'
REL_RES = TRK / 'e2e_relic' / 'results'


def dump(name, obj):
    RES.mkdir(parents=True, exist_ok=True)
    json.dump(obj, open(RES / name, 'w'), indent=1, default=float)


def load(path, default=None):
    path = Path(path)
    return json.load(open(path)) if path.exists() else default
