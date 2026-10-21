"""Dense Gaussian preparation (Codex's generator) with a forced box-scale S bias, reduced box.

Usage:
  python3 -I forced_bias.py predict <scratch_dir> <out_json>     # writes noise files + pre-registered predictions
  python3 -I forced_bias.py run <scratch_dir> <gaussian_binary>  # runs Codex's unchanged gaussian wrapper
Preparation [P]: N=128, dx=0.4 (L=51.2, Codex's resolution, smaller box), seed 101, Codex's filter exp(-k^2/2),
unit-rms S/X noise, then S noise -> S noise + mu (a uniform k=0 shift = forced box-scale bias of 0.01*mu in S).
"""
from __future__ import annotations
import json
import os
import subprocess
import sys
import time
from pathlib import Path
import numpy as np
from scipy.special import erf, erfinv

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gaussian_diagnostics as gd  # noqa: E402

N, DX, SEED = 128, 0.4, 101
MUS = [0.0, 0.02, 0.04, 0.07, 0.1, 0.15, 0.25]


def predict(scr, outp):
    scr.mkdir(parents=True, exist_ok=True)
    base = gd.codex_noise(N, DX, SEED)
    res = {'written_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), 'N': N, 'dx': DX, 'L': N * DX, 'seed': SEED,
           'rule': ('m_lin = <sign S_lin> of the linear field (incl. the forced zero mode, which grows with the k=0 rate) '
                    'at total rms 1.0. Bias law (Claude bias_threshold, gamma=1.41, x10=1.54): network removed within the '
                    'window if erfinv(m_f) K^gamma >= 1.54 with K the basin-area coarsening from formation to a=2.2 '
                    'measured in the unbiased (mu=0) run; predicted m_end = erf(erfinv(m_f) K^gamma). Threshold m_f for '
                    'takeover = erf(1.54 K^-gamma).'), 'runs': {}}
    for mu in MUS:
        q = base.copy()
        q[2] += mu
        path = scr / f'fb_N{N}_mu{mu:g}.noise.f64' if N != 128 else scr / f'fb_mu{mu:g}.noise.f64'
        q.astype('<f8').tofile(path)
        r = gd.analyze(q, N, DX, label=f'mu{mu:g}', field_marks=(0.5, 1.0))
        res['runs'][f'N{N}_mu{mu:g}'] = {'mu': mu, 'S_shift': 0.01 * mu, 'noise_file': str(path), 'predictions': r['predictions'],
                                    'C_uS': r['C_uS'], 'C_wS': r['C_wS']}
        p = r['predictions']['rms_1.0']
        print('mu', mu, 'a=%.3f m_lin(rms1)=%.4f [no bare cubic %.4f, fluct only %.4f] m_lin(rms0.5)=%.4f S0/sig=%.3f sd_m=%.4f' % (
            p['a'], p['mean_sign_with_zero_mode'], p['mean_sign_zero_mode_without_bare_cubic'], p['mean_sign_fluct_only'],
            r['predictions']['rms_0.5']['mean_sign_with_zero_mode'], p['S0_over_sigma_f'], p['ensemble_sd_of_mean_sign']), flush=True)
    for K in (1.2, 1.4, 1.6):
        res[f'takeover_threshold_m_f_K{K}'] = float(erf(1.54 * K ** -1.41))
    outp.write_text(json.dumps(res, indent=1))


def run(scr, exe):
    env = dict(os.environ, OMP_NUM_THREADS='2')
    log = []
    for mu in MUS:
        name = f'fb_N{N}_mu{mu:g}' if N != 128 else f'fb_mu{mu:g}'
        cmd = [str(exe), '--N', str(N), '--dx', str(DX), '--dt', '0.06', '--T', '48', '--noise-file',
               str(scr / f'{name}.noise.f64'), '--mirror', '1', '--out', str(scr / name), '--threads', '2']
        t0 = time.monotonic()
        r = subprocess.run(cmd, capture_output=True, text=True, env=env)
        log.append({'name': name, 'cmd': cmd, 'rc': r.returncode, 'stdout': r.stdout.strip(), 'stderr': r.stderr.strip(),
                    'wall_s': time.monotonic() - t0})
        print(json.dumps(log[-1]), flush=True)
        (scr / 'forced_bias_log.json').write_text(json.dumps(log, indent=1))


if __name__ == '__main__':
    # optional 5th/6th args: N and comma-separated mu list (paired larger-box controls; same seed 101 as Codex's
    # gaussian_s101 (N=192) and gaussian_large_s101 (N=288), so mu=0 there IS Codex's recorded run)
    if len(sys.argv) > 4:
        N = int(sys.argv[4])
        MUS = [float(x) for x in sys.argv[5].split(',')]
    if sys.argv[1] == 'predict':
        predict(Path(sys.argv[2]), Path(sys.argv[3]))
    else:
        run(Path(sys.argv[2]), Path(sys.argv[3]))
