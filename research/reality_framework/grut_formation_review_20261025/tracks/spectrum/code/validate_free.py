#!/usr/bin/env python3
"""Validate the extraction pipeline on a linear run with exactly known particle/antiparticle content.
Usage: validate_free.py RUNDIR OUTJSON  (expects pk.bin.json, pk_ext.bin, pk_hist.bin, pk_snap3000/6000.bin)."""
import json, sys
import numpy as np
sys.path.insert(0, __import__('os').path.dirname(__file__))
import spec_lib as S

run, outp = sys.argv[1], sys.argv[2]
truth = json.load(open(run + '/pk.bin.json'))
Rext = [300., 600., 1200.]
t, d = S.load_ext(run + '/pk', len(Rext))
hist = S.load_hist(run + '/pk', 4)
res = dict(truth={k: truth[k] for k in ['NX', 'NXbar', 'Q', 'E_charged', 'E_neutral', 'N_neutral']}, fixed_radius={}, spatial={})
for kR, R in enumerate(Rext):
    U, Ut, Ur, Z, Zt, Zr = S.series(d, kR)
    fl = S.time_fluxes(t, U, Ut, Ur, Z, Zt, Zr)
    iR = 1 + [60., 300., 600., 1200.].index(R)
    Qin_T = hist[-1, 1 + 2 * iR]
    wt = S.taper(t, -1e9, t[-1], 400.)
    cs = S.charged_spectrum(t, U, Ur, wt)
    ns = S.neutral_spectrum(t, Z, Zr, wt)
    res['fixed_radius'][str(R)] = dict(
        time_domain_Qpass=float(fl['Qpass'][-1]), Q0_minus_Qinside=float(truth['Q'] - Qin_T),
        NX=cs['NX'], NXbar=cs['NXbar'], NX_incoming=cs['NX_in'], NXbar_incoming=cs['NXbar_in'],
        EX_plus_EXbar=cs['EX'] + cs['EXbar'], flux_subthreshold=cs['flux_subthreshold'],
        Ns=ns['Ns'], Es=ns['Es'])
for T in [3000, 6000]:
    snap = S.load_snap(run + f'/pk_snap{T}.bin')
    for Rw in [0., 100.]:
        sp = S.spatial_projection(snap, Rw, 50.) if Rw > 0 else S.spatial_projection(snap, -1, 1e-9)
        spn = S.spatial_projection(snap, Rw, 50., m2=6., field='neutral') if Rw > 0 else S.spatial_projection(snap, -1, 1e-9, m2=6., field='neutral')
        res['spatial'][f'T{T}_Rw{Rw:.0f}'] = dict(NX=sp['NX'], NXbar=sp['NXbar'], E=sp['EX'] + sp['EXbar'], Ns=spn['Ns'], Es=spn['Es'])
# linearity / propagation residual: propagate R=300 outgoing spectrum to R=600 and compare
U1, _, Ur1, *_ = S.series(d, 0); U2, *_ = S.series(d, 1)
wt = S.taper(t, -1e9, t[-1], 400.)
cs1 = S.charged_spectrum(t, U1, Ur1, np.ones_like(t))
pred = S.propagate(t, cs1['aout'], cs1['w'], 300., len(t))
m = t < 5000
res['propagation_residual_300_to_600_t_lt_5000'] = float(np.linalg.norm((pred - U2)[m]) / np.linalg.norm(U2[m]))
json.dump(res, open(outp, 'w'), indent=1)
print(json.dumps(res, indent=1))
