#!/usr/bin/env python3
"""Run the full-ledger, conditional H55 spectral kinetic histories."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
from scipy.integrate import solve_ivp

from kinetic_engine import compile_engine, source as engine_source
from trace import ARCHIVE, COLUMNS, TRACKS

X3 = TRACKS / 'x3_actual'
TAB_INPUT = X3 / 'results/kinetic/actual_native_6400.in.gz'
CASES = {
    'matched': ('matched', 6400, 1.0),
    'matched_fine_q': ('matched', 12800, 1.0),
    'matched_half_step': ('matched', 6400, 0.5),
    'matched_omitted_s': ('omitted', 6400, 1.0),
    'harmonic': ('harmonic', 6400, 1.0),
    'harmonic_fine_q': ('harmonic', 12800, 1.0),
    'harmonic_half_step': ('harmonic', 6400, 0.5),
    'harmonic_to_interface': ('harmonic_stop', 6400, 1.0),
}
LEDGER_COLS = 't a rho_s rho_dr rho_sm rho exchange_error budget_error f_min f_max'.split()


def input_tables():
    """Load only the serialized probability tables from the archived X3 input."""
    tables = []
    with gzip.open(TAB_INPUT, 'rt') as f:
        for _ in range(8):
            next(f)
        for _ in range(4):
            next(f)  # channel parameters
            nq = int(next(f))
            for _ in range(3 * nq):
                next(f)  # momenta/weights and two occupations
            if int(next(f)):
                for _ in range(6 * nq):
                    next(f)
            for _ in range(2 * nq):
                next(f)  # frozen flags
            na, nx = map(int, next(f).split())
            la = np.array([float(next(f)) for _ in range(na)])
            pp = np.array([float(next(f)) for _ in range(na * nx)]).reshape(na, nx)
            assert np.all(np.diff(la) > 0) and np.min(pp) >= -1e-14 and np.max(pp) <= 2
            tables.append((la, pp))
        assert f.read().strip() == ''
    return tables


def load_source(path):
    metadata = json.loads(ARCHIVE.read_text())
    d = dict(zip(COLUMNS, np.loadtxt(path).T, strict=True))
    t = d['t']
    last = int(np.argmin(abs(t - metadata['termination']['last_diag_t'])))
    turns = np.flatnonzero(
        (d['chi_field'][1:] >= 1) & (d['chi_dot'][:-1] > 0) & (d['chi_dot'][1:] <= 0)
    ) + 1
    first_turn = int(turns[0])
    assert 0 < t[first_turn] < 0.4 < t[last]
    return metadata['run'], d, last, first_turn


def initialize(meta, d, last, first_turn, kind):
    harmonic = kind.startswith('harmonic')
    i = first_turn if harmonic else last
    v2 = meta['P']['v']**2
    energy = np.array([d[f'rho_{j}'][i] * v2 for j in range(1, 5)])
    native_counts = np.array([d[f'n_{j}'][i] * v2 for j in range(1, 5)])
    sigma = d['sigma'][i] * v2
    chi = d['chi'][i] * v2
    if harmonic:
        sigma += energy[0]
        chi += energy[1:].sum()
    # The same supplied hard/thermal split and mean hard-particle energy as X3.
    proxy = json.loads((TRACKS / 'lf_chain/results/handover_hohLF55_mid.json').read_text())['state']
    sm = d['sm'][i] * v2
    rho_h = sm * proxy['rho_h'] / (proxy['rho_h'] + proxy['rho_th'])
    st = dict(mP=meta['P']['mP'], t=d['t'][i], t_end=d['t'][last] if kind == 'harmonic_stop' else 1e8,
              a=d['a'][i], rho_sig=sigma, rho_chi=chi,
              rho_h=rho_h, rho_th=sm - rho_h, n_h=rho_h * proxy['n_h'] / proxy['rho_h'])
    spectator = d['spectator'][i] * v2 if kind != 'omitted' else 0.0
    return st, energy, native_counts, spectator, d['rho'][i] * v2, float(d['H'][i])


def check_spectator(binary):
    entries = []
    s0, dr0, dt, a0 = 2.0, 0.4, 0.9, 1.3
    for hubble, gamma in [(0, 0), (0, 0.3), (0.3, 0), (0.3, 0.3), (0.1, 0.5), (0.8, 0.1)]:
        args = [s0, dr0, gamma, dt, a0, a0 * math.exp(hubble * dt)]
        got = np.fromstring(subprocess.check_output(
            [str(binary), 'spectator-test', *map(str, args)], text=True
        ), sep=' ')
        ref = solve_ivp(
            lambda t, y: [-(3 * hubble + gamma) * y[0], -4 * hubble * y[1] + gamma * y[0]],
            (0, dt), [s0, dr0], method='DOP853', rtol=3e-14, atol=1e-14,
        ).y[:, -1]
        error = float(max(abs(got / ref - 1)))
        assert error < 5e-13
        entries.append(dict(H=hubble, Gamma=gamma, relative_error=error))
    return entries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    parser.add_argument('tape', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--cases', nargs='+', choices=CASES, default=list(CASES))
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--collect', action='store_true',
                        help='validate and collect already completed runs from this work directory')
    args = parser.parse_args()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    os.environ['EK_SCRATCH'] = str(work)
    sys.path.insert(0, str(TRACKS / 'e2e_kinetic/code'))
    import ek_common as common
    sys.path.insert(0, str(X3 / 'code'))
    from kinetic import sphere

    meta, data, last, turn = load_source(args.tape)
    scale = meta['H_GeV']
    spec, mpar = common.species('HFa', scale, 0.05e-9)
    sm = common.sm_coeffs()
    tables = input_tables()
    binaries = {}
    for factor in sorted({CASES[c][2] for c in args.cases}):
        if args.collect:
            binary = work / f'engine_dt{factor:g}'
            if Path(f'{binary}.c').read_text() != engine_source(factor):
                raise ValueError(f'Unexpected saved engine source for {binary}')
            binaries[factor] = binary
        else:
            binaries[factor] = compile_engine(work, factor)
    spectator_control = check_spectator(binaries[next(iter(binaries))])
    gamma = 1.3e-6 * scale / common.H_REF

    def run(name):
        kind, nq, dt_factor = CASES[name]
        st, energy, numbers, spectator, source_rho, source_H = initialize(meta, data, last, turn, kind)
        q, w = common.grid_log(3.0, 3e7, nq)
        occ = []
        represented_count = []
        for value, species in zip(energy, spec, strict=True):
            if kind.startswith('harmonic'):
                occ.append(np.zeros(2 * nq))
                represented_count.append(0.0)
            else:
                modes, moments = sphere(value, q, w, st['a'], species['M'], 0.5)
                occ.append(modes)
                represented_count.append(moments['N'])
        prefix = work / name
        infile = prefix.with_suffix('.in')
        generated = prefix.with_suffix('.collect.in') if args.collect else infile
        common.write_input(generated, st, spec, mpar, sm, dict(coherent=0), [(q, w)] * 4, occ,
                           tables=tables, dshell=-1e5)
        with generated.open('a') as f:
            f.write(f'{spectator:.17g} 0 {gamma:.17g}\n')
        if args.collect:
            if generated.read_bytes() != infile.read_bytes() or not Path(f'{prefix}_state.txt').is_file():
                raise ValueError(f'Missing or mismatched completed run for {name}')
            generated.unlink()
        else:
            proc = subprocess.run([str(binaries[dt_factor]), 'run', str(infile), str(prefix)],
                                  text=True, capture_output=True, timeout=2400)
            if proc.returncode:
                raise RuntimeError(f'{name}: {proc.stderr}')
        rec = common.read_rec(prefix)
        ledger = dict(zip(LEDGER_COLS, np.loadtxt(f'{prefix}_ledger.txt').T, strict=True))
        assert len(rec['t']) == len(ledger['t'])
        assert np.max(abs(rec['t'] / ledger['t'] - 1)) < 6e-11
        assert np.max(abs(rec['rho'] / ledger['rho'] - 1)) < 6e-7
        full_gap = float(ledger['rho'][0] / source_rho - 1)
        if kind != 'omitted':
            assert abs(full_gap) < 1e-6
        assert np.all(ledger['rho_s'] >= 0) and np.all(ledger['rho_dr'] >= 0)
        assert np.min(ledger['f_min']) >= -1e-14 and np.max(ledger['f_max']) <= 1 + 1e-14
        a, t = ledger['a'], ledger['t']
        K = a**2 * np.sqrt(ledger['rho'] / (3 * st['mP']**2))
        K_sm = a**2 * np.sqrt(ledger['rho_sm'] / (3 * st['mP']**2))
        ratio = ledger['rho_dr'][-1] / ledger['rho_sm'][-1]
        shift = 43 / 7 * (10.75 / common.GSTAR)**(1 / 3) * ratio
        balance = np.abs(ledger['rho_s'] * a**3 / (spectator * a[0]**3)
                         - np.exp(-gamma * (t - t[0]))) if spectator else np.zeros_like(t)
        assert max(balance) < 5e-11
        idx = sorted({0, len(t)-1, *(int(np.argmin(abs(t - x)))
                     for x in (0.3, 0.5, 1, 2, 3, 10, 100, 1e3, 1e4, 1e5, 1e6, 1e7))})
        result = dict(
            case=name, kind=kind, momentum_bins=nq, step_ceiling_factor=dt_factor,
            Gamma_s_over_H_i=gamma, initial=st, rho_s_initial=spectator,
            native_source_N_energy=list(energy), native_source_N_count=list(numbers),
            initialized_N_count=represented_count, full_source_energy_relative_gap=full_gap,
            Hubble_splice_relative_gap=float(K[0] / a[0]**2 / source_H - 1),
            maximum_fixed_volume_exchange_error=float(max(ledger['exchange_error'])),
            maximum_FRW_budget_error=float(max(ledger['budget_error'])),
            maximum_spectator_comoving_decay_error=float(max(balance)),
            occupation_range=[float(min(ledger['f_min'])), float(max(ledger['f_max']))],
            t_end=float(t[-1]), K_total=float(K[-1]), K_visible=float(K_sm[-1]),
            rho_DR_over_SM=float(ratio), Delta_Neff_if_standard_entropy=float(shift) if kind != "harmonic_stop" else None,
            final_spectator_fraction=float(ledger['rho_s'][-1] / ledger['rho'][-1]),
            samples=[dict(t=float(t[i]), a=float(a[i]), K_total=float(K[i]), K_visible=float(K_sm[i]),
                          f_N=float(rec['f_N'][i]), f_SM=float(ledger['rho_sm'][i]/ledger['rho'][i]),
                          f_s=float(ledger['rho_s'][i]/ledger['rho'][i]),
                          f_DR=float(ledger['rho_dr'][i]/ledger['rho'][i])) for i in idx],
            engine_sha256=hashlib.sha256(Path(f'{binaries[dt_factor]}.c').read_bytes()).hexdigest(),
            input_sha256=hashlib.sha256(infile.read_bytes()).hexdigest(),
            accepted_cosmological_point=False,
        )
        print(name, f'K={K[-1]:.9g} visible={K_sm[-1]:.9g} gap={full_gap:.3g}'
              f' budget={result["maximum_FRW_budget_error"]:.3g}', flush=True)
        return result

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(run, args.cases))
    out = dict(scope='conditional kinetic evolution, with a separate spectator and dark-radiation ledger',
               first_chi_turn=float(data['t'][turn]), final_archived_time=float(data['t'][last]),
               source_ptable_input_sha256=hashlib.sha256(TAB_INPUT.read_bytes()).hexdigest(),
               seesaw_eV=0.05, spectator_controls=spectator_control, cases=results)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, allow_nan=False) + '\n')


if __name__ == '__main__':
    main()
