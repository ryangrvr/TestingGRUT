#!/usr/bin/env python3
"""lf_chain (Claude, lattice-fermion closure pass, BRIEF_LF.md; GitHub issue #13 context): common module.

REUSE, NOT REWRITE.  e2e_chain's code (e2c_common, splice, handover, scan, chain; and through it mf_ideal's engine,
e2e_kinetic's network, e2e_relic's pmap/relic) is imported READ-ONLY and unchanged.  Only two redirections:
  * E2C_SCRATCH -> this track's scratch (LFC_SCRATCH, default .../wt2/lf_chain/sp0), set before e2c_common is imported;
    E.RES (where scan.py and handover.py write their summaries) -> this track's results/.
  * splice.lattice_block (also used by handover.state) reads the lf_fermion export results/hist_HFa_LF.json for the
    closure keys 'KF*' (KFmid primary, KFlo_nl, KFhi_nl, KFmid_nl_s2, KFmid_H84_nl); same column logic as the original
    (re-typed, unchanged), other keys fall through to the original function.
Every history, splice and hand-over rule, closure, holder amplitude and phase, relic rule and criterion is a model
choice; none follows from closure.  Python -B, one thread per process.
"""
import os, sys, json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
TRACK = HERE.parent
TRK = HERE.parents[1]
RES = TRACK / 'results'
SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/wt2')
LSP = Path(os.environ.get('LFC_SCRATCH', str(SCR / 'lf_chain' / 'sp0')))
os.environ['E2C_SCRATCH'] = str(LSP)
os.environ.setdefault('OMP_NUM_THREADS', '1')
sys.path.insert(0, str(TRK / 'e2e_chain' / 'code'))
import e2c_common as E                                       # noqa: E402  (read-only)
assert str(E.SP).startswith(str(LSP)), E.SP
E.RES = RES
RES.mkdir(parents=True, exist_ok=True)
import splice as S                                            # noqa: E402  (read-only)

LF_HIST = TRK / 'lf_fermion' / 'results' / 'hist_HFa_LF.json'
LF_MOM = TRK / 'lf_fermion' / 'results' / 'fluid_moments_HFa_LF.json'
HIST = SCR / 'lf_chain' / 'sp0' / 'hist'                      # final histories of this track (scratch)
_orig_lb = S.lattice_block


def lattice_block(name, key):
    if not str(key).startswith('KF'):
        return _orig_lb(name, key)
    d = json.load(open(LF_HIST))
    if key == d['primary']:
        L, head, cfg = d['lattice'], d['head'], d['config']
    else:
        A = d['alternates'][key]
        L, head, cfg = A['lattice'], A['head'], A['config']
    g = lambda k: np.array([np.nan if x is None else float(x) for x in L[k]])
    return dict(t=g('t'), a=g('a'), H=g('H'), R=g('R'), R_cyc=g('R_cyc'), ms2=g('m_s2_mean'), ms2_std=g('m_s2_std'),
                f_sig=g('f_sig'), f_chi=g('f_chi'), f_grad=g('f_grad'), f_rad=g('f_rad'), f_SM=g('f_SM'), s_std=g('s_std'),
                eps=g('eps_H'), K=g('K_running'), head=head, cfg=cfg)


S.lattice_block = lattice_block


def dump(path, obj):
    E.dump(path, obj)


def load(path, default=None):
    return E.load(path, default)
