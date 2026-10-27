# REPRODUCTION (Codex's response.c via timing_common.integrate, run in our own work dir; nothing written into $SP/b73).
import sys, json, numpy as np
sys.dont_write_bytecode = True
SP = '/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad'
sys.path.insert(0, SP + '/b73/research/reality_framework/grut_source_timing_20261026/code')
import timing_common as T
C = np.load(SP + '/b73/research/reality_framework/grut_source_timing_20261026/quiet3.npz')
times = C['times']
exe = T.build()
out = {}
for rtol in (2e-10, 2e-12):
    bg, vec, ctl = T.integrate(exe, 'sk_tight_%g' % rtol, 3., .5, np.array([2.13]), times, e=1e-7, rtol=rtol)
    amp, dq, nx = T.responses(3., bg, vec)
    np.savez(SP + '/review/st_response_skeptic/codex_k213_rtol%g.npz' % rtol, times=times, amp=amp[:, 0], dq=dq[:, 0], nx=nx, a=np.exp(bg[:, 15]))
    print(rtol, ctl['check'], ctl['seconds'], flush=True)
