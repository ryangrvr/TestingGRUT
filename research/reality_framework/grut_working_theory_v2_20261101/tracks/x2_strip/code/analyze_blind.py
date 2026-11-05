#!/usr/bin/env python3
"""Apply 6d836ca without changing fp_lf's node/branch functions.

Refuse to open a science analysis until every scheduled key is present on
both computed nodes. Extended nodes, if requested, are kept out of strict
classification. Original checkpoint rows take precedence over replays.
"""
import argparse
import gzip
import hashlib
import json
import math
import os
import sys
from pathlib import Path

from resume_blind import COMPUTED, EXTENDED, TRACK, TRACKS, grid, original, read, result

BAND = (1.03e-6,4.12e-6)
NODE_ORDER = ('hLF55_mid','hLF84_mid','hLF_mid','hLF270_mid')
BRANCHES = ('free X','Q-matter')


def collect(work,node):
    prior=read(original(node))
    initial=prior.get('runs',{})
    base=read(result(work/'baselines'/node,node))
    assert base['runs'].keys()=={'1.000:+0'}
    missing=[key for key in grid() if key not in initial and key!='1.000:+0']
    rows=dict(initial)
    rows.setdefault('1.000:+0',base['runs']['1.000:+0'])
    dirs=sorted((work/'shards'/node).glob('*'))
    for d in dirs:
        dd=read(result(d,node))
        for k,v in dd['runs'].items():
            if k=='1.000:+0': continue
            if k not in missing or k in rows:
                raise RuntimeError(f'{node}: unexpected/duplicate key {k}')
            rows[k]=v
    if set(rows)!=set(grid()):raise RuntimeError(f'{node}: incomplete grid')
    assert all(rows[k]==v for k,v in initial.items())
    for k,r in rows.items():
        assert k==f"{r['c0']:.3f}:{r['delta_deg']:+.0f}"
        assert r['H_i_GeV']==base['meta']['H_i_GeV']
    return {'node':node,'meta':base['meta'],'runs':{k:rows[k] for k in grid()},'preserved_checkpoint_rows':len(initial)}


def in_band(r):
    ec=r['e_crit'] if r else None
    return isinstance(ec,(int,float)) and math.isfinite(ec) and BAND[0]<=ec<=BAND[1]


def primary(r):
    return bool(r and r['self_consistent'] and 0<=r['w']<=1)


def runs_by_delta(keys):
    accepted=set(keys)
    out={str(d):[] for d in (-4,-2,0,2,4)}
    for d,seq in out.items():
        run=[]
        for k in [x for x in grid() if int(x.split(':')[1])==int(d)]+['end']:
            if k in accepted:run.append(float(k.split(':')[0]))
            elif run:
                seq.append({'c0_lo':run[0],'c0_hi':run[-1],'count':len(run),'isolated':len(run)==1,'c0_values':run})
                run=[]
    return out


def dump(path,data):
    path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('work',type=Path);ap.add_argument('--extended',action='store_true')
    ap.add_argument('--seal-only',action='store_true',help='Archive complete grids before evaluating any roots')
    a=ap.parse_args();work=a.work.resolve();out=work/'analysis'
    nodes=NODE_ORDER if a.extended else COMPUTED
    # Check counts first; there is no partial verdict.
    for n in nodes:
        keys=set(read(original(n)).get('runs',{}))
        keys.update(read(result(work/'baselines'/n,n)).get('runs',{}))
        for d in (work/'shards'/n).glob('*'):keys.update(read(result(d,n)).get('runs',{}))
        if keys!=set(grid()):raise SystemExit(f'{n}: {len(keys)}/135; analysis refused')
    scans={n:collect(work,n) for n in nodes}
    out.mkdir(parents=True,exist_ok=True)
    for n,data in scans.items():
        raw=json.dumps(data,separators=(',',':'),allow_nan=False).encode()
        archive=out/f'completed_scan_{n}.json.gz'
        sealed=gzip.compress(raw,mtime=0)
        if archive.exists() and archive.read_bytes()!=sealed:
            raise RuntimeError(f'{n}: grid does not match the sealed archive')
        archive.write_bytes(sealed)
    if a.seal_only:
        print('Sealed complete grids:',', '.join(nodes),'(no root analysis)')
        return
    os.environ['LFC_SCRATCH']=str(work/'analysis_import_cache')
    sys.path.insert(0,str(TRACKS/'lf_chain'/'code'))
    import fp_lf as F

    strictX=[];strictQ=[];extband=[];extold=[];extQ=[];overflow=[];rows={}
    def root(lo,hi,branch,key,pair):
        try:return F.root(lo,hi,branch)
        except OverflowError:
            # The preregistered function can only overflow in extrapolation.
            f0,f1=F.fb(lo,branch),F.fb(hi,branch)
            w=-f0/(f1-f0)
            if 0<=w<=1:raise
            overflow.append({'key':key,'pair':pair,'branch':branch,'w':w})
            return None

    for k in grid():
        nn={n:F.nd(scans[n]['runs'][k]) for n in nodes}
        lo,hi=(nn[n] for n in COMPUTED)
        roots={b:root(lo,hi,b,k,COMPUTED) for b in BRANCHES}
        x,q=(roots[b] for b in BRANCHES)
        okx=primary(x) and in_band(x);okq=primary(q)
        if okx:strictX.append(k)
        if okq:strictQ.append(k)
        ex=[]
        if a.extended:
            for i,(left,right) in enumerate(zip(nodes,nodes[1:])):
                for b in BRANCHES:
                    r=root(nn[left],nn[right],b,k,(left,right))
                    if not r:continue
                    inside=0<=r['w']<=1
                    if not (inside or i==0 and r['w']<0 or i==len(nodes)-2 and r['w']>1):continue
                    uses_A=left in EXTENDED or right in EXTENDED
                    if not (uses_A or r['extrapolated']):continue  # already in strict reading
                    ex.append(dict(r,uses_A_node=uses_A,node_names=[left,right]))
                    if not r['self_consistent']:continue
                    rec={'key':k,'root':ex[-1]}
                    if b=='Q-matter':extQ.append(rec)
                    else:
                        if in_band(r):extband.append(rec)
                        ec=r['e_crit']
                        if isinstance(ec,(int,float)) and 1e-6<ec<=4e-6:extold.append(rec)
        rows[k]={'primary':roots,'strict_M2_X':okx,'strict_M2_Q':okq,'extended_only':ex}

    report={'preregistration':'6d836ca5051b7bf3e9d677e26c3b7b2db0ae0f6c',
            'pause_commit':'800c3bfd21757f9d9826b6276686f6c3362be611','ecrit_band':BAND,
            'nodes':{n:{'grid_points':len(data['runs']),'preserved_checkpoint_rows':data['preserved_checkpoint_rows'],
                        'raw_scan_sha256':hashlib.sha256((out/f'completed_scan_{n}.json.gz').read_bytes()).hexdigest()}
                     for n,data in scans.items()},
            'strict_M2_X_keys':strictX,'strict_M2_X_strips':runs_by_delta(strictX),
            'strict_M2_Q_keys':strictQ,'strict_M2_Q_strips':runs_by_delta(strictQ),
            'strict_both_branch_keys':sorted(set(strictX)&set(strictQ)),
            'strict_M2_Q_keys_also_in_ecrit_band':[k for k in strictQ if in_band(rows[k]['primary']['Q-matter'])],
            'computed_pair_extrapolation_diagnostic_only':[k for k in grid() if (r:=rows[k]['primary']['free X'])
                                                            and r['self_consistent'] and r['extrapolated'] and in_band(r)],
            'extended':{'band':extband,'old_cutoff':extold,'Q_matter':extQ} if a.extended else None,
            'unrepresentable_extrapolations':overflow,'roots':rows}
    dump(out/'blind_verdict.json',report)
    print('Strict M2-X:',len(strictX),'of 135',runs_by_delta(strictX))
    print('Strict M2-Q:',len(strictQ),'of 135',runs_by_delta(strictQ))
    if a.extended:print('Extended-only: M2-X band',len(extband),'earlier cutoff',len(extold),'Q roots',len(extQ))


if __name__=='__main__':
    main()
