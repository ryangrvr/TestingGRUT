#!/usr/bin/env python3
"""ls_ledger (Claude, own code): does the late undamped N_sigma count the sigma quanta that actually decay?  [E]

Physical sigma -> A_L A_L happens at <a> ~ 16 (Gamma_sig/H_i = 6.5e-3), while the undamped release only has a
well-defined (frozen) N_sigma from a ~ 90.  Model check: add the perturbative vacuum width as a radial friction
Gamma_sig * rdot on Sigma (ls_damped.c = our sp_yledger yrelease.c plus that one term and a counter of decayed
quanta, energy removed / m_sig).  Compare the decayed number with the undamped late N_sigma, and record the
backreaction on eps, N_X and N_Y.  The friction model is the standard phenomenological treatment of perturbative
decay of an oscillating field; it is not valid in the early anharmonic/nonadiabatic stage, so it is a model
estimate, not a kinetic solution.  Y_H -> Y_L A (a ~ 23) is not damped here (it conserves the Y number).
Usage: python3 -B ls_damped.py <out.json>   (~5 min, 4 cores)
"""
import json, subprocess, sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / 'sp_yledger' / 'code'))
import ymodel as Y   # noqa: E402  (our sp_yledger minimum)

SCR = Path('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/ls_ledger/build')
BIN = SCR / 'ls_damped'
H = 434867798206.6179
COLS = ['t', 'a', 'H', 'NX', 'NY', 'NS', 'EX', 'Ncirc', 'eps', 'NYH', 'NYL', 'JYcodex', 'Nsig', 'gauss', 'ledger',
        'a4rho', 'BoverA', 'EYexc', 'Esig', 'phi_frac', 'EYfull', 'Neff', 'Nsig_decayed', 'Nsig_decayed_omegaE']
HV, MA = 1e11, 6e-7 * 1e11
MSIG = np.sqrt(2) * HV
RA = MA**2 / MSIG**2
GSIG = MSIG**3 / (32 * np.pi * HV**2) * np.sqrt(1 - 4 * RA) * (1 - 4 * RA + 12 * RA * RA)   # GeV


def build():
    SCR.mkdir(parents=True, exist_ok=True)
    subprocess.run(['gcc', '-O3', '-o', str(BIN), str(HERE / 'ls_damped.c'), '-lm'], check=True)


def run(spec):
    label, delta, gfac, aend = spec
    z0, m = Y.initial(H, 10.0, delta)
    p = Y.params(H)
    g = lambda v: '%.17g' % float(v)
    args = [str(BIN), '0', '0', '10', g(p['M']), g(p['h']), g(p['mu']), g(0.5 * aend**2), '1e-11', '1e-14', '160',
            g(gfac * GSIG / H)] + [g(v) for c in z0 for v in (c.real, c.imag)]
    out = subprocess.run(args, capture_output=True, text=True, check=True).stdout.split('\n')
    S, W, END = [], {}, None
    for line in out:
        if line.startswith('S '):
            S.append(dict(zip(COLS, map(float, line.split()[1:]))))
        elif line.startswith('W '):
            v = line.split()
            W[f'{float(v[1]):g}-{float(v[2]):g}'] = dict(zip(COLS, map(float, v[4:])))
        elif line.startswith('END'):
            v = line.split(); END = dict(steps=int(v[1]), max_gauss=float(v[3]), max_ledger=float(v[4]))
    # where do the decays happen?  fraction of the decayed count accumulated by a
    nd = np.array([s['Nsig_decayed'] for s in S]); aa = np.array([s['a'] for s in S])
    ndw = np.array([s['Nsig_decayed_omegaE'] for s in S])
    frac = {f'a<={x}': float(np.interp(x, aa, nd) / nd[-1]) if nd[-1] > 0 else None for x in (5, 10, 17, 30, 50, 90)}
    w = W.get('100-300', {})
    return dict(label=label, delta=delta, Gamma_factor=gfac, a_end=aend, window_100_300={k: w.get(k) for k in
                ['eps', 'NX', 'EX', 'NYH', 'NYL', 'Nsig', 'Nsig_decayed', 'Ncirc']},
                N_sigma_decayed_final=float(nd[-1]), N_sigma_decayed_omegaE_final=float(ndw[-1]),
                decayed_fraction_by_a=frac, controls=END,
                history=[{k: s[k] for k in ('a', 'eps', 'NYH', 'NYL', 'Nsig', 'Nsig_decayed', 'Nsig_decayed_omegaE', 'EX')}
                         for s in S[::4]])


def main(out):
    build()
    specs = []
    for c, d0, seed in [(10, np.pi / 6, 2026102510), (17, 5 * np.pi / 24, 2026102517)]:
        deltas = d0 + np.random.default_rng(seed).normal(0, 1e-3, 48)
        for gf in (0.0, 1.0, 1 / 3, 3.0):
            specs.append((f'{c}_nominal_G{gf:.3g}', float(d0), gf, 300.0))
        for i in (0, 12, 24, 36):
            specs.append((f'{c}_m{i}_G1', float(deltas[i]), 1.0, 300.0))
            specs.append((f'{c}_m{i}_G0', float(deltas[i]), 0.0, 300.0))
    with Pool(4) as p:
        res = p.map(run, specs)
    Path(out).write_text(json.dumps(dict(method=__doc__, Gamma_sigma_GeV=GSIG, Gamma_sigma_over_H=GSIG / H,
                                         runs=res), indent=1, default=float))
    for r in res:
        w = r['window_100_300']
        print(r['label'], 'eps %.6f NX %.2f NYH+NYL %.3f Nsig_left %.4g Ndec %.3f NdecOmegaE %.3f frac(a<=17) %s' % (
            w['eps'], w['NX'], w['NYH'] + w['NYL'], w['Nsig'], r['N_sigma_decayed_final'], r['N_sigma_decayed_omegaE_final'],
            r['decayed_fraction_by_a']['a<=17']), flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
