#!/usr/bin/env python3
"""Measure the local mass-phase operator on a passive replay of the old H55 lattice."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from run import ARCHIVE, TRACKS, archived_history, read_columns

LATTICE = TRACKS / 'lf_fermion/code'
T_END = 0.8


def instrument(source, phase_scales):
    replacements = [
        ('A_NONE = A_TE0 + NL * NCH };',
         'A_PHASE = A_TE0 + NL * NCH, A_NONE = A_PHASE + 4 };'),
        ('    double acc[NACC];', r'''
    char phase_path[640];
    snprintf(phase_path, sizeof phase_path, "%s/phase.txt", rundir);
    FILE *phase = fopen(phase_path, resumed ? "a" : "w");
    if (!phase || A_NONE > NACC) return 1;
    if (!resumed) fprintf(phase, "# t a H cos_sigma sin_sigma cos_chi sin_chi sig_mean chi_mean\n");
    double acc[NACC];'''),
        ('                        acc[A_SA] += sav;', f'''
                        double sigma_phase = vH*s/{phase_scales[0]:.17g};
                        double chi_phase = vH*c/{phase_scales[1]:.17g};
                        acc[A_PHASE] += cos(sigma_phase);
                        acc[A_PHASE+1] += sin(sigma_phase);
                        acc[A_PHASE+2] += cos(chi_phase);
                        acc[A_PHASE+3] += sin(chi_phase);
                        acc[A_SA] += sav;'''),
        ('            double xm = acc[A_XX], sm_ = acc[A_SG];', r'''
            double xm = acc[A_XX], sm_ = acc[A_SG];
            fprintf(phase, "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                    S.t, a, H, acc[A_PHASE], acc[A_PHASE+1], acc[A_PHASE+2], acc[A_PHASE+3], sm_, xm);
            fflush(phase);'''),
        ('    fclose(fd);','    fclose(phase);\n    fclose(fd);'),
    ]
    for old,new in replacements:
        if source.count(old)!=1:
            raise ValueError(f'Expected one source anchor: {old}')
        source=source.replace(old,new,1)
    return source


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('work',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--threads',type=int,default=2)
    parser.add_argument('--budget',type=float,default=650)
    args=parser.parse_args()
    work=args.work.resolve()
    work.mkdir(parents=True,exist_ok=True)
    meta=json.loads(ARCHIVE.read_text())
    cfg=meta['run']
    scales=[cfg['channels'][i]['f'] for i in (0,1)]
    sys.path.insert(0,str(LATTICE))
    import lfcommon as L

    L.SP=work
    L.RUNS=work/'runs'
    L.BIN=work/'lf3d'
    L.SRC=work/'phase_lattice.c'
    code=instrument((LATTICE/'lf3d.c').read_text(),scales)
    if not L.SRC.exists() or L.SRC.read_text()!=code:
        L.SRC.write_text(code)
    tag='h55_phase'
    run=L.RUNS/tag
    if not (run/'run.par').exists():
        L.make_run(tag,H_GeV=cfg['H_GeV'],N=cfg['N'],L=cfg['L'],dt=cfg['dt'],
                   seed=cfg['seed'],kcut_mchi=cfg['kcut_mchi'],gw=0,snaps=[],t_end=T_END,
                   extra=['tracers 0',*meta['held_kappa']])
    if not L.run(tag,budget=args.budget,threads=args.threads):
        print('Phase replay paused at its wall-clock checkpoint.',flush=True)
        return 75
    path=run/'phase.txt'
    d=read_columns(path)
    legacy=archived_history()
    t=d['t']
    reference=slice(None,len(t))
    assert max(abs(t-legacy['t'][reference]))<5e-8
    archive_error={key:float(max(abs(d[key]/legacy[key][reference]-1))) for key in ('a','H')}
    assert max(archive_error.values())<0.005
    metrics={}
    for parent,mean,f in zip(('sigma','chi'),('sig_mean','chi_mean'),scales,strict=True):
        phase_mean=d[mean]*cfg['P']['v']/f
        raw=d[f'cos_{parent}']+1j*d[f'sin_{parent}']
        centered=raw*np.exp(-1j*phase_mean)
        assert np.max(abs(raw))<1+1e-12
        difference2=2*(1-centered.real)
        assert min(difference2)>-1e-12
        operator_rms=np.sqrt(np.maximum(difference2,0))
        coherence=abs(raw)
        peak=int(np.argmax(operator_rms))
        metrics[parent]=dict(
            peak_operator_rms=float(operator_rms[peak]), peak_t=float(t[peak]),
            minimum_phase_coherence=float(min(coherence)),
            scalar_coherent_fraction_at_peak=float(legacy['fcoh'][peak]),
            gradient_fraction_at_peak=float(legacy['fgrad'][peak]),
            samples=[dict(t=float(t[i]),operator_rms=float(operator_rms[i]),
                          phase_coherence=float(coherence[i]),
                          scalar_coherent_fraction=float(legacy['fcoh'][i]),
                          gradient_fraction=float(legacy['fgrad'][i]))
                     for i in sorted({0,len(t)-1,peak,*(int(np.argmin(abs(t-x))) for x in
                         (0.1,0.18,0.20,0.23,0.24,0.25,0.30,0.40,0.50,0.60,0.75))})],
        )
    result=dict(scope='passive old-drain 64-cubed scalar replay with local mass phases',
                last_t=float(t[-1]), rows=len(t), model_end=T_END, archive_relative_errors=archive_error,
                phase_scales=scales, engine_sha256=hashlib.sha256(code.encode()).hexdigest(),
                phase_history_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), parents=metrics)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:{key:v for key,v in y.items() if key!='samples'} for k,y in metrics.items()}),
          flush=True)
    return 0


if __name__=='__main__':
    raise SystemExit(main())
