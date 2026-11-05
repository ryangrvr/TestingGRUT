#!/usr/bin/env python3
"""Recreate a lost scratch history from the exact committed hand-over state.

Only paths are adapted. The dephased solver, splice and anh correction are
the pinned lf_chain methods. No scan rows are read here.
"""
import argparse
import concurrent.futures as cf
import json
import math
import os
import shutil
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRACKS = HERE.parents[1]
NODES = {
    'hLF84_mid': (8.4e10, 'KFmid_H84_nl'),
    'hLF_mid': (1.3e11, 'KFmid'),
    'hLF55_mid': (5.5e10, 'KFmid_H84_nl'),
    'hLF270_mid': (2.7e11, 'KFmid'),
}
# The committed table summary uses the rounded 0.41667 for N1. The other
# channels' rounded values coincide with the exact input ratios.
TABLE_MASSES = (0.41667, 1.6, 0.8, 0.48)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('node', choices=(*NODES,'tables'))
    ap.add_argument('work', type=Path)
    a = ap.parse_args()
    work = a.work.resolve()
    root = work / 'histories' / a.node
    root.mkdir(parents=True, exist_ok=True)
    os.environ['LFC_SCRATCH'] = str(root)
    tables = work / 'original_ptables'
    os.environ['EK_SCRATCH'] = str(tables if a.node=='tables' else root/'kin')
    sys.path.insert(0, str(TRACKS / 'lf_chain' / 'code'))
    sys.path.insert(0, str(TRACKS / 'e2e_kinetic' / 'code'))
    import lfc_common as L
    import build_hist as B
    import ek_common as C
    import ek_ptable as P
    import ek_run as R
    import ek_scan as S
    import numpy as np

    if a.node=='tables':
        saved=json.loads((C.RES/'ptable_summary.json').read_text())['tables']
        C.compile_engine()
        def replay(mh):
            table=C.SP/f'ptable_{mh:.5f}.npz';raw=C.SP/f'tab_{mh:.5f}.txt'
            if not (table.exists() and raw.exists()):return P.build(mh)
            # Recheck a completed table without re-running its C integration.
            T=P.read_table(raw);assert list(sorted(T))==list(P.AGRID)
            rows=[];spectra=[]
            for A in sorted(T):
                k,pm,pp=T[A].T;tot=pm+pp
                sm=P.smooth(k,tot,min(2,.25*A+.3));om=np.sqrt(k*k+mh*mh)
                spectra.append(np.interp(np.linspace(0,1,P.NX)*(1.25*A+3*mh+2),k,sm,left=sm[0],right=0.))
                e=np.trapezoid(k*k*tot*om,k)/(2*math.pi**2)
                rows.append(dict(A=A,eps_per_period=e,eps_smoothed=np.trapezoid(k*k*sm*om,k)/(2*math.pi**2),
                                 h_plus_share=np.trapezoid(k*k*pp*om,k)/max(np.trapezoid(k*k*tot*om,k),1e-300),pmax=tot.max()))
            with np.load(table) as z:
                assert np.allclose(z['pp'],spectra,rtol=1e-13,atol=1e-15)
                assert np.allclose(z['lA'],np.log(sorted(T)),rtol=1e-13,atol=1e-15)
            return rows
        with cf.ThreadPoolExecutor(max_workers=4) as pool:
            fresh=list(pool.map(replay,TABLE_MASSES))
        errors={}
        for mh,rows in zip(TABLE_MASSES,fresh):
            k=f'{mh:.5f}';refs=saved[k];assert len(refs)==len(rows)
            errors[k]=float(max(abs(r[field]-s[field])/max(abs(s[field]),1e-30)
                                for r,s in zip(rows,refs) for field in r))
        ok=max(errors.values())<2e-6
        (tables/'check.json').write_text(json.dumps({'passed':ok,'relative_errors':errors},indent=2)+'\n')
        print('Original table summary replay', 'PASS' if ok else 'FAIL',errors,flush=True)
        if not ok:raise SystemExit(1)
        return

    H, lattice = NODES[a.node]
    ref = json.loads((L.RES / f'handover_ho{a.node}.json').read_text())
    spref = json.loads((L.RES / 'splice_checks.json').read_text())[a.node]
    raw = ref['state']
    assert raw['_provenance']['lattice'] == lattice
    assert ref['info']['H_GeV'] == H
    assert ref['info']['t0'] == ref['info']['t_coh'] == raw['t'] == 3.
    spec, _ = C.species('HFa', H, 0.05e-9)
    assert json.loads((tables/'check.json').read_text())['passed']
    # The first reconstruction attempt would have used the nearby mass ratios
    # in the species list. The original committed table summary instead uses
    # the four ek_ptable defaults; do not allow a nearer stale file to win.
    for p in C.SP.glob('ptable_*.npz'):p.unlink()
    for mh in TABLE_MASSES:
        p=tables/f'ptable_{mh:.5f}.npz';shutil.copyfile(p,C.SP/p.name)
    state = {k: float(raw[k]) for k in ('t','a','H','rho_sig','rho_chi','rho_h','n_h','rho_th')}
    occ = []
    for j, s in enumerate(spec):
        q, w = C.grid_linear(state['a'] * 3e4, 6000)
        occ.append((q, w, S.occ_from_moments(float(raw['E_N'][j]), float(raw['n_N'][j]), q, w, state['a'], s['M'])))
    res = R.run_history('HFa', 'ho'+a.node, H, 0.05e-9, 3., 3., 1e8, 'SU2',
                        dshell=ref['info']['dshell'], nqlog=ref['info']['nqlog'], state=state, occ0=occ)
    fluid = C.SP / f'fluid_ho{a.node}.json'
    R.export_fluid(res, 'HFa', fluid, f'e2e_chain moments hand-over from e2e_lattice {lattice} at t_sw = 3; dephased network from t_sw')
    export = json.loads(fluid.read_text())
    splice = L.S.build('HFa', a.node, lattice, str(fluid), H)
    mod = B.anh(splice['path'])
    # Independent saved summaries from the original history construction.
    pairs = {
        'kinetic_K': (export['summary']['K'], ref['summary']['K']),
        'kinetic_T_RH': (export['summary']['T_RH_GeV'], ref['summary']['T_RH_GeV']),
        'kinetic_t_SM_half': (export['summary']['t_SM_half'], ref['summary']['t_SM_half']),
        'splice_K': (splice['reheating']['K'], spref['K']),
        'splice_H_jump': (splice['splice']['H_jump'], spref['H_jump_at_tsw']),
        'splice_m_s2_after': (mod['ms2_after_new'], spref['anh']['ms2_after_new']),
    }
    checks = {k: {'replayed': x, 'archived': y, 'relative_difference': abs(x-y)/max(abs(y),1e-30)} for k,(x,y) in pairs.items()}
    ok = all(math.isclose(x,y,rel_tol=2e-6,abs_tol=1e-10) for x,y in pairs.values())
    (root / 'rebuild_check.json').write_text(json.dumps({'node':a.node,'passed':ok,'checks':checks},indent=2)+'\n')
    print(a.node,'history checks', 'PASS' if ok else 'FAIL', 'max relative difference',max(x['relative_difference'] for x in checks.values()),flush=True)
    if not ok:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
