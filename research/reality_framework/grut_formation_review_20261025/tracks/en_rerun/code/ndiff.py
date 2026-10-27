"""Compare two npz files key by key: max abs and rel diffs."""
import numpy as np,sys
a=np.load(sys.argv[1]);b=np.load(sys.argv[2])
ka=set(a.files);kb=set(b.files)
if ka!=kb:print('key mismatch',sorted(ka^kb)[:10])
worst=0
for k in sorted(ka&kb):
    x=a[k];y=b[k]
    if x.shape!=y.shape:print(k,'shape',x.shape,y.shape);continue
    if x.dtype.kind in 'fc' or y.dtype.kind in 'fc':
        d=np.max(abs(x-y)) if x.size else 0.;s=np.max(abs(x)) if x.size else 0.
        r=d/s if s else d;worst=max(worst,r)
        if r>float(sys.argv[3] if len(sys.argv)>3 else 0):print(f'  {k} shape={x.shape} maxabs={d:.3e} rel={r:.3e}')
    else:
        if not np.array_equal(x,y):print('  ',k,'nonfloat differs')
print('keys',len(ka&kb),'worst_rel',f'{worst:.3e}')
