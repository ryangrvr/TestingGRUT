#!/usr/bin/env python3
"""Verify a spectral restart without changing any particle distribution."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

from kinetic_replay import LEDGER_COLS, input_tables
from kinetic_engine import source
from trace import TRACKS


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path, help='completed kinetic_replay work directory')
    parser.add_argument('kinetic_results', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    work = args.work.resolve()
    os.environ['EK_SCRATCH'] = str(work)
    sys.path.insert(0, str(TRACKS / 'e2e_kinetic/code'))
    import ek_common as common

    results = json.loads(args.kinetic_results.read_text())
    cases = {c['case']: c for c in results['cases']}
    comparator = cases['harmonic']
    interface = cases['harmonic_to_interface']
    prefix = work / 'harmonic_to_interface'
    nq = interface['momentum_bins']
    state = common.read_state(prefix, [nq] * 4)
    old = dict(zip(LEDGER_COLS, np.loadtxt(f'{prefix}_ledger.txt').T, strict=True))
    for col in ('t', 'a'):
        assert state[col] == old[col][-1]
    assert abs(state['rho_tot'] / old['rho'][-1] - 1) < 1e-14
    assert state['t'] == results['final_archived_time']

    rho_s, rho_dr = old['rho_s'][-1], old['rho_dr'][-1]
    st = dict(state, mP=interface['initial']['mP'], t_end=1e8)
    spec, mpar = common.species('HFa', 5.5e10, results['seesaw_eV'] * 1e-9)
    q, w = common.grid_log(3.0, 3e7, nq)
    new = work / 'harmonic_restarted'
    inp = new.with_suffix('.in')
    common.write_input(inp, st, spec, mpar, common.sm_coeffs(), dict(coherent=0), [(q, w)] * 4,
                       state['occ'], frozen=state['frozen'], tables=input_tables(), dshell=-1e5)
    with inp.open('a') as f:
        f.write(f'{rho_s:.17g} {rho_dr:.17g} {interface["Gamma_s_over_H_i"]:.17g}\n')
    binary = work / 'engine_dt1'
    assert Path(f'{binary}.c').read_text() == source(1.0)
    subprocess.run([str(binary), 'run', str(inp), str(new)], check=True, timeout=2400)
    d = dict(zip(LEDGER_COLS, np.loadtxt(f'{new}_ledger.txt').T, strict=True))
    assert d['t'][0] == state['t'] and d['a'][0] == state['a']
    assert d['rho_s'][0] == rho_s and d['rho_dr'][0] == rho_dr
    a, mP = d['a'][-1], st['mP']
    K = float(a*a * np.sqrt(d['rho'][-1] / (3 * mP*mP)))
    K_sm = float(a*a * np.sqrt(d['rho_sm'][-1] / (3 * mP*mP)))
    init_error = float(d['rho'][0] / state['rho_tot'] - 1)
    total_difference = K / comparator['K_total'] - 1
    visible_difference = K_sm / comparator['K_visible'] - 1
    assert abs(init_error) < 1e-9
    assert max(abs(total_difference), abs(visible_difference)) < 0.01
    output = dict(
        scope='dephased harmonic reduction, with direct occupation-preserving restart',
        interface_time=state['t'], interface_a=state['a'], interface_rho=state['rho_tot'],
        interface_rho_s=float(rho_s), interface_rho_DR=float(rho_dr), restart_energy_relative_gap=init_error,
        K_total=K, K_visible=K_sm, K_total_relative_difference=total_difference,
        K_visible_relative_difference=visible_difference, accepted_cosmological_point=False,
        checkpoint_state_sha256=hashlib.sha256(Path(f'{prefix}_state.txt').read_bytes()).hexdigest(),
    )
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + '\n')
    print(json.dumps(output), flush=True)


if __name__ == '__main__':
    main()
