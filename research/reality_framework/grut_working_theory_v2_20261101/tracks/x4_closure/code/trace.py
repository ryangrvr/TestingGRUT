#!/usr/bin/env python3
"""Replay the archived H55 scalar source with passive, same-time-level output."""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

TRACK = Path(__file__).resolve().parents[1]
TRACKS = TRACK.parent
ARCHIVE = TRACKS / 'x2_strip/results/lower_H_lattice.json'
ENGINE = TRACKS / 'lf_fermion/code/lf3d.c'
VARS = ('rho', 'n', 'de', 'dn', 'de_pert', 'dn_pert')
COLUMNS = (
    't a H p rho sigma chi spectator sm sigma_field chi_field sigma_dot chi_dot '
    'kin_sigma kin_chi'
).split() + [f'{name}_{j}' for j in range(1, 5) for name in VARS]


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f'Expected exactly one source anchor: {old}')
    return text.replace(old, new, 1)


def instrument(source):
    patches = [
        ('A_NONE = A_TE0 + NL * NCH };',
         'A_PERT_E = A_TE0 + NL * NCH, A_PERT_N = A_PERT_E + NCH };'),
        ('    double acc[NACC];', '''
    char tape_path[640];
    snprintf(tape_path, sizeof tape_path, "%s/moments.txt", rundir);
    FILE *tape = fopen(tape_path, resumed ? "a" : "w");
    if (!tape) { perror("moments.txt"); return 1; }
    if (!resumed) fprintf(tape, "# '''+' '.join(COLUMNS)+r'''\n");
    double acc[NACC];'''),
        ('                    double Es = 0, Ex = 0;',
         '                    int pert[NCH] = {0};\n                    double Es = 0, Ex = 0;'),
        ('                            Gf[j] = gnp > gp ? gnp : gp;',
         '                            pert[j] = gp >= gnp;\n                            Gf[j] = gnp > gp ? gnp : gp;'),
        ('                            acc[A_INJN0 + j] += dE / Ei[j];', '''
                            acc[A_INJN0 + j] += dE / Ei[j];
                            if (pert[j]) {
                                acc[A_PERT_E + j] += dE;
                                acc[A_PERT_N + j] += dE / Ei[j];
                            }'''),
        ('        /* fluids */', r'''
        fprintf(tape, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g"
                " %.17g %.17g %.17g %.17g %.17g %.17g", S.t, a, H, p, rho,
                acc[A_ES], acc[A_EX] + Grd, acc[A_SKIN] + acc[A_SGRAD] + acc[A_SPOT],
                S.rSM, acc[A_SG], acc[A_XX], acc[A_SDOT], acc[A_XDOT], acc[A_KS], acc[A_KX]);
        for (int j = 0; j < NCH; j++)
            fprintf(tape, " %.17g %.17g %.17g %.17g %.17g %.17g", S.rN[j], S.nN[j],
                    acc[A_INJ0 + j], acc[A_INJN0 + j], acc[A_PERT_E + j], acc[A_PERT_N + j]);
        fputc('\n', tape);
        /* fluids */'''),
        ('    fclose(fd);', '    fclose(tape);\n    fclose(fd);'),
    ]
    for old, new in patches:
        source = replace_once(source, old, new)
    return source


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work', type=Path)
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--budget', type=float, default=650)
    args = parser.parse_args()
    config = json.loads(ARCHIVE.read_text())
    cfg = config['run']
    sys.path.insert(0, str(ENGINE.parent))
    import lfcommon as lattice

    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    lattice.SP = work
    lattice.RUNS = work / 'runs'
    lattice.BIN = work / 'lf3d'
    generated = work / 'lf3d_trace.c'
    code = instrument(ENGINE.read_text())
    if not generated.exists() or generated.read_text() != code:
        generated.write_text(code)
    lattice.SRC = generated
    tag = 'h55_trace'
    run = lattice.RUNS / tag
    if not (run / 'run.par').exists():
        lattice.make_run(
            tag, H_GeV=cfg['H_GeV'], N=cfg['N'], L=cfg['L'], dt=cfg['dt'],
            seed=cfg['seed'], kcut_mchi=cfg['kcut_mchi'], gw=0, snaps=[],
            extra=['tracers 0', *config['held_kappa']],
        )
    if not lattice.run(tag, budget=args.budget, threads=args.threads):
        print(f'Trace checkpointed in {run}. Resume with the same command.', flush=True)
        return 75
    data = np.loadtxt(run / 'moments.txt')
    assert len(data) == config['termination']['steps']
    assert data.shape[1] == len(COLUMNS) and np.all(np.isfinite(data))
    assert np.max(np.abs(np.diff(data[:, 0]) / cfg['dt'] - 1)) < 1e-10
    record = {
        'engine_sha256': hashlib.sha256(generated.read_bytes()).hexdigest(),
        'rows': len(data),
        'threads': args.threads,
        'first_t': float(data[0, 0]),
        'last_t': float(data[-1, 0]),
    }
    (run / 'trace_info.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
