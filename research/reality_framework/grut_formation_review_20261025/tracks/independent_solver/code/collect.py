#!/usr/bin/env python3
"""Collect this track's run outputs into results/*.json and compare with Codex's recorded numbers.

Codex numbers are read from the ff49e5d / bb211d1 packets extracted with `git archive` into the scratchpad
(quoted [Q]); this script never runs Codex code.
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from drive_c import history, COLS

SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/independent_solver')
RUNS = SCR / 'runs'
FF = SCR / 'pk_ff49/research/reality_framework/grut_causal_receiver_20261024'
BB = SCR / 'pk_bb21/research/reality_framework/grut_core_conversion_20261023'
OUT = Path(__file__).resolve().parents[1] / 'results'
c = {k: i for i, k in enumerate(COLS)}


def codex_hist(name):
    x = json.loads((FF / f'{name}_RESULTS.json').read_text())
    return np.array(x['history']), x


def at(H, t):
    i = np.argmin(abs(H[:, 0] - t))
    return i if abs(H[i, 0] - t) < 1e-6 else None


def window_mean(H, col, Q0, lo, hi, step=None):
    tau = H[:, 0]
    m = (tau > lo + 1e-9) & (tau <= hi + 1e-9)
    if step:
        m &= np.isclose(np.mod(tau, step), 0) | np.isclose(np.mod(tau, step), step)
    if not m.any() or tau[-1] < hi - 1e-6:
        return None
    return float(np.mean(H[m, col]) / Q0)


def main():
    OUT.mkdir(exist_ok=True)
    runs = {}
    for f in sorted(RUNS.glob('*_history.bin')):
        name = f.name[:-len('_history.bin')]
        H = history(str(RUNS / name)); E0, Q0 = H[0, 1], H[0, 2]
        meta = json.loads((RUNS / f'{name}_meta.json').read_text()) if (RUNS / f'{name}_meta.json').exists() else {}
        r = dict(tau_end=float(H[-1, 0]), E0=E0, Q0=Q0, args=meta.get('args'),
                 max_rel_energy_plus_receiver_error=float(np.max(np.abs((H[:, 1] + H[:, 3]) / E0 - 1))),
                 max_rel_charge_plus_receiver_error=float(np.max(np.abs((H[:, 2] + H[:, 4]) / Q0 - 1))))
        for lo, hi in [(9600, 12000), (19200, 24000), (28800, 36000), (38400, 48000)]:
            v = window_mean(H, c['Q60'], Q0, lo, hi)
            if v is not None:
                r[f'Q60_mean_{lo}_{hi}'] = v
                r[f'Q60_mean_{lo}_{hi}_codex_sampling_40'] = window_mean(H, c['Q60'], Q0, lo, hi, 40.0)
                r[f'Q60trap_mean_{lo}_{hi}'] = window_mean(H, c['Q60_trap'], Q0, lo, hi)
                m = (H[:, 0] > lo) & (H[:, 0] <= hi + 1e-9)
                r[f'Q60_slope_per_1000tau_{lo}_{hi}'] = float(np.polyfit(H[m, 0], H[m, c['Q60']] / Q0, 1)[0] * 1000)
        r['Q60_at'] = {str(t): float(H[at(H, t), c['Q60']] / Q0) for t in [400, 800, 1600, 2400, 4000, 8000, 12000, 16000,
                                                                         20000, 24000, 36000, 48000] if at(H, t) is not None}
        res = RUNS / f'{name}_RESULT.json'
        if res.exists():
            rr = json.loads(res.read_text())
            r['samples'] = rr.get('samples'); r['modes'] = rr.get('modes')
        runs[name] = r
    codex = {}
    for n in ['P025_time', 'P0125', 'P00625', 'P025_remote', 'P025_remote_time', 'I025']:
        Hc, x = codex_hist(n)
        q0 = x['Qbar_initial']
        codex[n] = dict(core60_late_mean_fraction=x['core60_late_mean_fraction'],
                        Q60_at={str(int(t)): float(Hc[at(Hc, t), 5] / q0) for t in [400, 800, 1600, 2400, 4000, 8000, 12000, 16000, 20000, 24000] if at(Hc, t) is not None},
                        parameters=x['parameters'])
        if n == 'I025':
            codex[n]['final_net_Eflux'] = x['final_net_Eflux']; codex[n]['final_net_Qflux'] = x['final_net_Qflux']
    bb = {n: json.loads((BB / f'{n}_RESULTS.json').read_text())['core60_late_mean_fraction'] for n in
          ['N025', 'N0125', 'N00625', 'H0125_late']}
    (OUT / 'formation_runs.json').write_text(json.dumps(dict(runs=runs, codex_ff49e5d_quoted=codex,
                                                             codex_bb211d1_quoted=bb), indent=1))
    print(json.dumps({k: {kk: vv for kk, vv in v.items() if kk not in ('samples', 'modes', 'args', 'Q60_at')} for k, v in runs.items()}, indent=1))


if __name__ == '__main__':
    main()
