"""en_rerun: compare regenerated Codex outputs with the published 7d9c93b copies (own comparison code).
Numeric leaves: max relative difference, ignoring wall-clock/path/hash/threads/bookkeeping keys and values below
an absolute floor of 1e-12 (round-off identities such as 2e-16 vs 0).
usage: python3 -I reproduction_compare.py <published_dir> <regenerated_dir> <out_json>"""
import sys,json,glob,os
import numpy as np
pub,reg,outp=sys.argv[1:4]
SKIP=('seconds','raw_sha256','check','snapshots','threads','wall','reproduced_by','chunked','chunked_restarts')
def walk(a,b,path,acc):
    if isinstance(a,dict):
        for k in a:
            if any(s==k for s in SKIP):continue
            if k not in b:acc['missing'].append(path+'/'+k);continue
            walk(a[k],b[k],path+'/'+k,acc)
    elif isinstance(a,list):
        if len(a)!=len(b):acc['length'].append(path);return
        for i,(x,y) in enumerate(zip(a,b)):walk(x,y,f'{path}[{i}]',acc)
    elif isinstance(a,bool) or a is None or isinstance(a,str):
        if a!=b:acc['other'].append((path,a,b))
    else:
        x,y=float(a),float(b);d=abs(x-y)
        if d>1e-12:
            r=d/max(abs(x),abs(y));
            if r>acc['max'][0]:acc['max']=(r,path,x,y)
out={}
for f in sorted(glob.glob(pub+'/*.json')):
    name=os.path.basename(f);g=reg+'/'+name
    if not os.path.exists(g) or name in('MANIFEST.json','RAW_FILES.json','SOURCES.json','AUDIT_RESULTS.json','BENCHMARK.json','EXECUTION.json'):continue
    A=json.load(open(f));B=json.load(open(g))
    if name=='EXPORT_CHECKS.json':  # row order differs (RUNS.json completion order); key by tag
        A['rows']=sorted(A['rows'],key=lambda r:r['tag']);B['rows']=sorted(B['rows'],key=lambda r:r['tag'])
    acc=dict(max=(0.,None,None,None),missing=[],length=[],other=[]);walk(A,B,'',acc);out[name]=acc
for f in sorted(glob.glob(pub+'/*.npz')):
    name=os.path.basename(f);g=reg+'/'+name
    if not os.path.exists(g):continue
    A=np.load(f);B=np.load(g);worst=(0.,None);shape=[]
    for k in A.files:
        if k not in B.files:shape.append(k+':missing');continue
        x,y=A[k],B[k]
        if x.shape!=y.shape:shape.append(f'{k}:{x.shape}->{y.shape}');continue
        if x.dtype.kind in 'fc':
            d=np.max(abs(x-y)) if x.size else 0.;s=np.max(abs(x)) if x.size else 1.
            r=float(d/s) if s else float(d)
            if r>worst[0]:worst=(r,k)
    out[name]=dict(max=worst,shape_mismatch=shape)
json.dump(out,open(outp,'w'),indent=1,default=str)
for k,v in out.items():print(k,v['max'],{x:v[x] for x in v if x!='max' and v[x]})
