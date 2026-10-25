#!/usr/bin/env python3
"""REPRODUCTION-DATA analysis: apply this track's own spectral extraction (spec_lib.py) to Codex's
r=240 probe traces, regenerated with Codex's unmodified code by the `reproduce` track and checked
against Codex's recorded raw sha256 hashes.  The arrays come from Codex's code; the analysis is ours.
Usage: analyze_codex_trace.py DATADIR RESULTSDIR OUTJSON"""
import hashlib, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import spec_lib as S

data, resdir, outp = sys.argv[1:4]
kb = [0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 12.0]
def binned(w, dN, dw, sign):
    k = np.sqrt(np.maximum(w * w - 1, 0)); res = []
    for lo, hi in zip(kb[:-1], kb[1:]):
        m = (np.sign(w) == sign) & (k >= lo) & (k < hi) & (np.abs(w) > 1)
        res.append(float(np.sum(dN[m]) * dw))
    return res
out = dict(k_bins=kb, probe_radius=240.0, runs={})
for name in ['P025_time', 'P0125', 'P00625', 'P025_remote_time', 'P025_remote']:
    tr = np.load(f'{data}/{name}_trace.npz')['trace']
    rj = json.load(open(f'{resdir}/{name}_RESULTS.json'))
    ok = hashlib.sha256(tr.tobytes()).hexdigest() == rj['raw_numeric_sha256']['trace']
    t = tr[:, 0]; R = 240.0
    psi = tr[:, 1] + 1j * tr[:, 2]; dpsi = tr[:, 4] + 1j * tr[:, 5]; pt = tr[:, 7] + 1j * tr[:, 8]
    U = R * psi; Ur = psi + R * dpsi; Ut = R * pt
    Z = R * tr[:, 3]; Zr = tr[:, 3] + R * tr[:, 6]; Zt = R * tr[:, 9]
    Q0 = rj['Qbar_initial']
    fl = S.time_fluxes(t, U, Ut, Ur, Z, Zt, Zr)
    rec = dict(hash_matches_codex_record=ok, Q0=Q0, Qpass_T_continuum_flux=float(fl['Qpass'][-1]), Epass_T=float(fl['Epass'][-1]),
               core60_late_mean_fraction_codex=rj['core60_late_mean_fraction'], windows={})
    for Tw in [6000, 12000, 18000, 24000]:
        if Tw > t[-1] + 1: continue
        m = t <= Tw; tt = t[m]; wt = S.taper(tt, -1e9, Tw, 500.)
        cs = S.charged_spectrum(tt, U[m], Ur[m], wt); ns = S.neutral_spectrum(tt, Z[m], Zr[m], wt)
        rec['windows'][str(Tw)] = dict(NX=cs['NX'], NXbar=cs['NXbar'], EX=cs['EX'], EXbar=cs['EXbar'], NX_incoming=cs['NX_in'],
                                       NXbar_incoming=cs['NXbar_in'], flux_subthreshold=cs['flux_subthreshold'], Ns=ns['Ns'], Es=ns['Es'],
                                       X_bins=binned(cs['w'], cs['dNout'], cs['dw'], 1), Xbar_bins=binned(cs['w'], cs['dNout'], cs['dw'], -1))
    out['runs'][name] = rec
json.dump(out, open(outp, 'w'), indent=1)
print(json.dumps({n: {k: v for k, v in r.items() if k != 'windows'} | {'w12000': {k: r['windows']['12000'][k] for k in ['NX', 'NXbar', 'EX', 'EXbar', 'NX_incoming', 'NXbar_incoming', 'Ns', 'Es']}} for n, r in out['runs'].items()}, indent=1))
