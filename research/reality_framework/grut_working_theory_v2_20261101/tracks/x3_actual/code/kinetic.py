#!/usr/bin/env python3
"""X3 energy-preserving conditional network runs; original code stays read-only.

Usage: python -B kinetic.py WORK [CASE ...]
WORK/original_ptables must contain the verified original four tables.
"""
import argparse
import concurrent.futures
import gzip
import hashlib
import json
import math
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from moments import TRACK, TRACKS, RES, source, dump

CASES = {
    'actual_ref_3200':dict(nq=3200),
    'actual_ref_6400':dict(nq=6400),
    'actual_ref_12800':dict(nq=12800),
    'actual_ref_dt_half':dict(nq=6400,dt_scale=.5),
    'proxy_ref_6400':dict(nq=6400,proxy=True),
    'actual_native_6400':dict(nq=6400,native=True),
    'actual_dilute_6400':dict(nq=6400,cap=.25),
}


def sphere(energy,q,w,a,mass,cap):
    """Exact discrete energy match, two helicities and one partial edge bin."""
    weights=2*w/a**3
    eps=weights*np.hypot(q/a,mass)
    cumulative=np.cumsum(cap*eps)
    j=int(np.searchsorted(cumulative,energy))
    if j>=len(q):raise ValueError('spectrum exceeds grid support')
    n=np.zeros(len(q));n[:j]=cap
    n[j]=(energy-(cumulative[j-1] if j else 0))/eps[j]
    assert -1e-12<=n[j]<=cap*(1+1e-12)
    represented=float(np.dot(eps,n))
    assert abs(represented/energy-1)<1e-12
    return np.r_[n,n],dict(E=represented,N=float(np.dot(weights,n)),
        cutoff_physical=float(q[j]/a),cap=cap,energy_relative_error=abs(represented/energy-1))


def initial(C,options):
    ref=json.loads((TRACKS/'lf_chain/results/handover_hohLF55_mid.json').read_text())['state']
    meta,cols,data=source();d=dict(zip(cols,data[-1]));v=meta['run']['P']['v']
    Hscale=meta['run']['H_GeV'];mp=C.MP_GEV/Hscale
    if options.get('proxy'):
        st={k:float(ref[k]) for k in ['t','a','H','rho_sig','rho_chi','rho_h','n_h','rho_th']}
        energy=np.array(ref['E_N']);number=np.array(ref['n_N'])
        total=ref['_provenance']['rho_tot']
        spectator=None; rounding=None
    else:
        total=d['rho']*v*v
        fsm=d['fSM']; hard=ref['rho_h']/(ref['rho_h']+ref['rho_th'])
        rhoh=fsm*total*hard
        st=dict(t=float(d['t']),a=float(d['a']),H=float(d['H']),
                rho_sig=d['fsig']*total,rho_chi=(d['fchi']+d['fgrad'])*total,
                rho_h=rhoh,n_h=rhoh*ref['n_h']/ref['rho_h'],rho_th=fsm*total*(1-hard))
        native_e=np.array([d[f'rN{j}']*v*v for j in range(1,5)])
        native_n=np.array([d[f'nN{j}']*v*v for j in range(1,5)])
        energy=native_e if options.get('native') else native_e.sum()*np.array(ref['E_N'])/sum(ref['E_N'])
        number=native_n
        spectator=d['fs']*total
        rounding=total-(st['rho_sig']+st['rho_chi']+st['rho_h']+st['rho_th']+native_e.sum()+spectator)
    active=st['rho_sig']+st['rho_chi']+st['rho_h']+st['rho_th']+energy.sum()
    ledger=dict(source_total_energy=total,represented_energy=active,
        unrepresented_energy=total-active,spectator_energy=spectator,
        printed_column_residual=rounding,full_source_energy_gap=(active-total)/total,
        Hubble_relative_jump=math.sqrt(active/(3*mp*mp))/st['H']-1,
        native_number_input=list(number),E_N=list(energy),
        note='spectator is not evolved by this reheating network; no full-splice pass is implied')
    st['mP']=mp
    return st,energy,ledger


def compile_control(C,work,factor):
    original=C.HERE/'ek_modes.c';text=original.read_text()
    # Only timestep ceilings are changed. The physical right-hand sides are identical.
    for old,new in [('double lim = 0.02 / H;',f'double lim = {factor:.17g} * 0.02 / H;'),
                    ('0.01 * rhoP[P] / Pc[P]',f'{factor:.17g} * 0.01 * rhoP[P] / Pc[P]'),
                    ('0.2 * tk)',f'{factor:.17g} * 0.2 * tk)')]:
        assert text.count(old)==1,(old,text.count(old))
        if factor!=1:text=text.replace(old,new)
    path=work/f'engine_dt{factor:g}.c';path.write_text(text)
    binary=path.with_suffix('')
    subprocess.run(['gcc','-O3','-march=native','-o',str(binary),str(path),'-lm'],check=True)
    return binary,hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('work',type=Path)
    ap.add_argument('cases',nargs='*');ap.add_argument('--workers',type=int,default=3)
    args=ap.parse_args();work=args.work.resolve();work.mkdir(parents=True,exist_ok=True)
    os.environ['EK_SCRATCH']=str(work/'original_ptables')
    sys.path.insert(0,str(TRACKS/'e2e_kinetic/code'))
    import ek_common as C
    import ek_run as R
    assert json.loads((C.SP/'check.json').read_text())['passed']
    spec,mpar=C.species('HFa',5.5e10,.05e-9);sm=C.sm_coeffs();tables=R.tables(spec)
    bins={f:compile_control(C,work,f) for f in [1.,.5]}
    cases=args.cases or list(CASES)
    outdir=RES/'kinetic';outdir.mkdir(parents=True,exist_ok=True)

    def run(name):
        opt=dict({'cap':.5,'dt_scale':1.},**CASES[name])
        st,energy,ledger=initial(C,opt)
        q,w=C.grid_log(3.,3e7,opt['nq']);occ=[];spectra=[]
        for E,s in zip(energy,spec):
            o,chk=sphere(E,q,w,st['a'],s['M'],opt['cap']);occ.append(o);spectra.append(chk)
        prefix=work/name;inp=prefix.with_suffix('.in')
        C.write_input(inp,dict(st,t_end=1e8),spec,mpar,sm,dict(coherent=0),[(q,w)]*4,occ,tables=tables,dshell=-1e5)
        binary,digest=bins[opt['dt_scale']];t0=time.time()
        p=subprocess.run([str(binary),'run',str(inp),str(prefix)],text=True,capture_output=True,timeout=2400)
        if p.returncode:raise RuntimeError(name+': '+p.stderr)
        rec=C.read_rec(prefix);summary=R.summarise(rec,st['mP'],5.5e10,spec)
        K=rec['a']**2*rec['H'];tail=float(np.interp(rec['t'][-1]/10,rec['t'],K))
        neterr=3*st['mP']**2*rec['H'][0]**2/ledger['represented_energy']-1
        assert abs(neterr)<1e-6
        initial_numbers=[x['N'] for x in spectra]
        checks=dict(represented_network_energy_error=neterr,
            last_decade_K_relative_drift=abs(K[-1]/tail-1),
            E_match_max=max(x['energy_relative_error'] for x in spectra),
            initial_numbers=initial_numbers,
            N_ratio_to_native_or_proxy=[a/b for a,b in zip(initial_numbers,ledger['native_number_input'])],
            full_source_energy_continuity_pass=bool(abs(ledger['full_source_energy_gap'])<=1e-6))
        out=dict(case=name,settings=opt,state=st,splice_ledger=ledger,spectra=spectra,summary=summary,checks=checks,
            runtime_seconds=time.time()-t0,engine_sha256=digest,input_sha256=hashlib.sha256(inp.read_bytes()).hexdigest(),
            accepted_cosmological_point=False)
        for suffix in ['_rec.txt','_state.txt']:
            raw=Path(str(prefix)+suffix).read_bytes()
            (outdir/(name+suffix+'.gz')).write_bytes(gzip.compress(raw,mtime=0))
        (outdir/(name+'.in.gz')).write_bytes(gzip.compress(inp.read_bytes(),mtime=0))
        dump(outdir/(name+'.json'),out)
        print(name,'K',summary['K'],'network_energy_error',neterr,'seconds',out['runtime_seconds'],flush=True)
        return out

    failures=[];done=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures={pool.submit(run,name):name for name in cases}
        for f in concurrent.futures.as_completed(futures):
            try:done.append(f.result())
            except Exception as e:
                failures.append({'case':futures[f],'error':repr(e)});print('FAILED',failures[-1],flush=True)
    dump(RES/'kinetic_execution.json',dict(completed=[x['case'] for x in done],failures=failures,
         table_check=json.loads((C.SP/'check.json').read_text())))
    if failures:raise SystemExit(1)


if __name__=='__main__': main()
