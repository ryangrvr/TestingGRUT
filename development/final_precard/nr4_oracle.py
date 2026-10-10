"""Full-K once, exact weighted finite-chart ablation oracle; abstract fixtures."""
from fractions import Fraction as F


def measure(weights,mask,hits):
    if not len(weights)==len(mask)==len(hits): raise ValueError('array length mismatch')
    ws=tuple(F(x) for x in weights)
    if any(w<0 for w in ws) or sum(ws)<=0: raise ValueError('invalid reference measure')
    if any(type(x) is not bool for x in mask+hits): raise ValueError('boolean arrays required')
    denominator=sum(w for w,m in zip(ws,mask) if not m)
    numerator=sum(w for w,m,h in zip(ws,mask,hits) if not m and h)
    excluded=sum(w for w,m in zip(ws,mask) if m)/sum(ws)
    return {'fraction':None if denominator==0 else numerator/denominator,
            'excluded_fraction':excluded,'retained_weight':denominator}


def full_mask(images,references,tolerances):
    if not images or len(images)!=len(references) or not tolerances: raise ValueError('finite nonempty chart/domain required')
    ts={k:F(v) for k,v in tolerances.items()}
    if any(v<0 for v in ts.values()): raise ValueError('negative tolerance')
    out=[]
    for current,reference in zip(images,references):
        if current is None or reference is None or set(current)!=set(ts) or set(reference)!=set(ts):
            raise ValueError('unresolved full-K image/reference')
        out.append(all(abs(F(current[k])-F(reference[k]))<=ts[k] for k in ts))
    return out


def runs(weights,mask,variant_hits,threshold,fibers):
    if len(fibers)!=len(weights): raise ValueError('fiber length mismatch')
    threshold=F(threshold)
    if not 0<threshold<=1: raise ValueError('appearance threshold')
    result={}
    for name,hits in variant_hits.items():
        measured=measure(weights,mask,hits)
        per_fiber={}
        for f in dict.fromkeys(fibers):
            indices=[i for i,x in enumerate(fibers) if x==f]
            total=sum(F(weights[i]) for i in indices)
            if total==0: continue  # no reference mass: no 0/0 fiber verdict
            m=measure([weights[i] for i in indices],[mask[i] for i in indices],[hits[i] for i in indices])
            if m['fraction'] is not None: per_fiber[f]=m['fraction']
        result[name]=dict(measured,per_fiber=per_fiber,
            appears=None if measured['fraction'] is None else measured['fraction']>=threshold)
    return result
