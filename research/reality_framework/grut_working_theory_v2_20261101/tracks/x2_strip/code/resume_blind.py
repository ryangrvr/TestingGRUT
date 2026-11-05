#!/usr/bin/env python3
"""Portable, restartable adapter around the unchanged preregistered scan.py.

Each subprocess owns its engine directory and result file. The parent reads
only completed keys, never classifies roots. Fixed assignment uses grid order
and the original checkpoint keys, with no dependence on calculated values.
"""
import argparse
import concurrent.futures as cf
import json
import os
import runpy
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
TRACK = HERE.parent
TRACKS = HERE.parents[1]
COMPUTED = ('hLF84_mid', 'hLF_mid')
EXTENDED = ('hLF55_mid', 'hLF270_mid')


def grid():
    cs = [1.,1.04,1.08,1.12,1.16,1.20] + [round(1.30+.01*i,2) for i in range(16)] + [1.355,1.365,1.375,1.385,1.395]
    return [f'{c:.3f}:{d:+.0f}' for d in (-4,-2,0,2,4) for c in sorted(set(cs))]


def read(path):
    return json.loads(path.read_text()) if path.exists() else {}


def original(node):
    return TRACK / 'results' / f'scan_x2_{node}.json'


def result(directory,node):
    return directory / 'results' / f'scan_x2_{node}.json'


def scan(work,node,engine_dir,keys):
    history = work / 'histories' / node
    assert read(history/'rebuild_check.json')['passed'] is True
    bench = history / 'hist' / f'bench_{node}.json'
    assert bench.is_file()
    os.environ['LFC_SCRATCH'] = str(engine_dir/'engine_cache')
    sys.path.insert(0,str(TRACKS/'lf_chain'/'code'))
    import lfc_common as L
    L.HIST = bench.parent
    L.RES = L.E.RES = engine_dir/'results'
    src = TRACKS/'e2e_chain'/'code'/'scan.py'
    sys.argv = [str(src),str(bench),f'x2_{node}'] + keys
    runpy.run_path(str(src),run_name='__main__')


def call(work,node,directory,keys):
    directory.mkdir(parents=True,exist_ok=True)
    cmd = [sys.executable,'-B',__file__,'worker',str(work),'--node',node,'--directory',str(directory),*keys]
    with (directory/'engine.log').open('a') as log:
        subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True)


def base(work,node):
    d = work/'baselines'/node
    if '1.000:+0' not in read(result(d,node)).get('runs',{}):
        call(work,node,d,[])
    prior = read(original(node)).get('runs',{}).get('1.000:+0')
    replay = read(result(d,node))['runs']['1.000:+0']
    # Do not print any science values. A platform/reconstruction check is
    # allowed before unblinding, and the original checkpoint row is retained.
    if prior is not None:
        probes = [('K',prior['K'],replay['K']),('eps',prior['eps'],replay['eps']),('E_X',prior['E_X'],replay['E_X'])]
        for field in ('S_per_Q','H_i_req_GeV'):
            probes.append((field,prior['by_e']['cL0_vanish_at_MY']['e1e-6'][field],replay['by_e']['cL0_vanish_at_MY']['e1e-6'][field]))
        errs={k:abs(a-b)/max(abs(a),1e-30) for k,a,b in probes}
        ok=max(errs.values())<2e-5
        (d/'replay_check.json').write_text(json.dumps({'passed':ok,'relative_errors':errs},indent=2)+'\n')
        print(node,'baseline replay', 'PASS' if ok else 'FAIL', 'maximum relative difference',max(errs.values()),flush=True)
        if not ok: raise RuntimeError(f'{node}: rebuilt baseline differs from checkpoint')
    else:
        print(node,'baseline built (no committed checkpoint)',flush=True)


def shard(work,node,index,keys):
    d = work/'shards'/node/f'{index:02d}'
    if not (d/'engine_cache').exists():
        shutil.copytree(work/'baselines'/node/'engine_cache',d/'engine_cache')
    # An interrupted write can only affect this shard. The per-point engine
    # cache remains available to the unchanged driver if summaries are replayed.
    done=set(read(result(d,node)).get('runs',{}))
    todo=[k for k in keys if k not in done]
    if todo:
        call(work,node,d,todo)
    have=set(read(result(d,node)).get('runs',{}))
    if not set(keys)<=have:
        raise RuntimeError(f'{node}/{index}: incomplete shard; see {d}/engine.log')
    print(node,'shard',index,'complete:',len(keys),'assigned keys',flush=True)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('mode',choices=('computed','extended','worker','progress'))
    ap.add_argument('work',type=Path)
    ap.add_argument('--node',choices=COMPUTED+EXTENDED)
    ap.add_argument('--directory',type=Path)
    ap.add_argument('--jobs',type=int,default=8)
    ap.add_argument('--shards',type=int,default=4)
    a,keys=ap.parse_known_args();work=a.work.resolve()
    if a.mode=='worker':
        scan(work,a.node,a.directory,keys);return
    if keys:ap.error(f'unexpected arguments {keys}')
    if a.mode=='progress':
        for node in COMPUTED+EXTENDED:
            keys=set(read(original(node)).get('runs',{}))
            for p in (work/'baselines'/node/'results',*(d/'results' for d in (work/'shards'/node).glob('*'))):
                try:keys.update(read(p/f'scan_x2_{node}.json').get('runs',{}))
                except json.JSONDecodeError:pass  # original driver is writing
            print(node,len(set(grid())&keys),'/',len(grid()))
        return
    nodes=COMPUTED if a.mode=='computed' else EXTENDED
    with cf.ThreadPoolExecutor(max_workers=min(a.jobs,len(nodes))) as pool:
        list(pool.map(lambda n:base(work,n),nodes))
    tasks=[]
    for node in nodes:
        done=set(read(original(node)).get('runs',{}))|{'1.000:+0'}
        missing=[key for key in grid() if key not in done]
        print(node,'scheduled',len(missing),'previously missing keys',flush=True)
        tasks.extend((node,i,missing[i::a.shards]) for i in range(a.shards) if missing[i::a.shards])
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as pool:
        for f in cf.as_completed([pool.submit(shard,work,*t) for t in tasks]):f.result()


if __name__=='__main__':
    main()
