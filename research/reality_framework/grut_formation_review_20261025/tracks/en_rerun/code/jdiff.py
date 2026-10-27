"""Recursive numeric diff of two JSON files: report max relative diff and any mismatched leaves."""
import json,sys,math
def walk(a,b,path,out):
    if isinstance(a,dict) and isinstance(b,dict):
        for k in sorted(set(a)|set(b)):
            if k not in a or k not in b: out.append((path+'/'+str(k),'missing',None,None));continue
            walk(a[k],b[k],path+'/'+str(k),out)
    elif isinstance(a,list) and isinstance(b,list):
        if len(a)!=len(b): out.append((path,'len',len(a),len(b)));return
        for i,(x,y) in enumerate(zip(a,b)):walk(x,y,path+f'[{i}]',out)
    elif isinstance(a,bool) or isinstance(b,bool) or a is None or b is None or isinstance(a,str) or isinstance(b,str):
        if a!=b: out.append((path,'value',a,b))
    else:
        x,y=float(a),float(b)
        d=abs(x-y)/max(abs(x),abs(y),1e-300) if (x or y) else 0.
        out.append((path,'num',d,(x,y)))
a=json.load(open(sys.argv[1]));b=json.load(open(sys.argv[2]));out=[]
walk(a,b,'',out)
thr=float(sys.argv[3]) if len(sys.argv)>3 else 1e-9
nums=[o for o in out if o[1]=='num'];other=[o for o in out if o[1]!='num']
big=[o for o in nums if o[2]>thr]
print(f'leaves={len(out)} numeric={len(nums)} max_rel={max([o[2] for o in nums],default=0):.3e} n_above_{thr:g}={len(big)} nonnumeric_mismatch={len(other)}')
skipkeys=('seconds','raw_sha256','wall','check','snapshots','sha')
for o in other[:20]:
    if not any(s in o[0] for s in skipkeys):print('  NONNUM',o)
for o in sorted(big,key=lambda o:-o[2])[:int(sys.argv[4]) if len(sys.argv)>4 else 15]:print('  ',o[0],f'{o[2]:.3e}',o[3])
