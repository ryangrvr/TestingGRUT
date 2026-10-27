import json,sys,numpy as np
pub=sys.argv[1];new=sys.argv[2]
def flat(x,p=''):
    if isinstance(x,dict):
        for k,v in x.items(): yield from flat(v,p+'/'+k)
    elif isinstance(x,list) and x and not isinstance(x[0],(int,float)):
        for i,v in enumerate(x): yield from flat(v,p+f'[{i}]')
    else: yield p,x
A=dict(flat(json.load(open(pub+'/PREPARATION_RESULTS.json'))));B=dict(flat(json.load(open(new+'/PREPARATION_RESULTS.json'))))
worst=0;wk=None
for k in A:
    a,b=A[k],B.get(k)
    if isinstance(a,str) or a is None: 
        if a!=b: print('DIFF',k)
        continue
    a=np.array(a,float);b=np.array(b,float)
    d=np.max(abs(a-b)/(abs(a)+1e-300)) if a.size else 0
    if d>worst: worst=d;wk=k
print('max relative diff PREPARATION_RESULTS',worst,wk)
for L in [12,24]:
    a=np.load(f'{pub}/PREPARATION_L{L}.npz');b=np.load(f'{new}/PREPARATION_L{L}.npz')
    for f in a.files:
        x=a[f];y=b[f];print(L,f,x.shape,float(np.max(abs(x-y))) if x.size else 0, float(np.max(abs(x))) if x.size else 0)
