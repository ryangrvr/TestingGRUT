"""Exact channel algebra and color-primitive classification; no trajectories."""
from fractions import Fraction as F
from itertools import combinations,permutations
from pathlib import Path
from typing import NamedTuple
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
checks=[]


def check(name,ok): checks.append({'name':name,'pass':bool(ok)})
def cabs2(z): return z[0]*z[0]+z[1]*z[1]
def cminus(z,a): return z[0]-a,z[1]


def norm_squared_at(a,d,z2,x):
    trace=a*x+d*(1-x)
    gap2=(a*x-d*(1-x))**2+4*z2*x*(1-x)
    return max(trace*trace,gap2)


def exact_success_norm2(u,v,c,pf):
    a=u-pf;d=v-pf;z2=cabs2(cminus(c,pf))
    candidates=[(F(0),d*d),(F(1),a*a)]
    K=(a+d)**2-4*z2
    if K<0:
        x=(d*(a+d)-2*z2)/K
        if 0<x<1:
            candidates.append((x,(a*x-d*(1-x))**2+4*z2*x*(1-x)))
    return max(candidates,key=lambda item:item[1])


def conditional2(u,v,c,x):
    p=u*x+v*(1-x)
    return x*(1-x)*((1-x)*cabs2(cminus(c,v))+x*cabs2(cminus(c,u)))/(p*p)


def polyadd(a,b):
    out=[F(0)]*max(len(a),len(b))
    for i,v in enumerate(a):out[i]+=v
    for i,v in enumerate(b):out[i]+=v
    return out


def polymul(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return out


samples=[(F(1,4),F(1,4),(F(1,4),F(0)),F(1,4)),
         (F(1,4),F(1,4),(F(0),F(0)),F(1,4)),
         (F(1,4),F(1,4),(F(-1,4),F(0)),F(1,4)),
         (F(3,10),F(1,5),(F(1,5),F(0)),F(1,4)),
         (F(1,4),F(1,4),(F(3,20),F(1,5)),F(1,4))]
xs=[F(0),F(1,7),F(1,2),F(4,5),F(1)]
for i,(u,v,c,pf) in enumerate(samples):
    check('Gram_positive_'+str(i),0<=u<=1 and 0<=v<=1 and cabs2(c)<=u*v)
    a=u-pf;d=v-pf;z2=cabs2(cminus(c,pf));vertex,maximum=exact_success_norm2(u,v,c,pf)
    check('quadratic_norm_closed_max_'+str(i),
          norm_squared_at(a,d,z2,vertex)==maximum and all(
          norm_squared_at(a,d,z2,x)<=maximum for x in xs))
    equality=True;bound=True
    for x in xs:
        p=u*x+v*(1-x)
        direct=(x*(u/p-1))**2+x*(1-x)*cabs2(cminus((c[0]/p,c[1]/p),F(1)))
        e2=u+v-2*c[0];det=u*v-cabs2(c)
        equality &= direct==conditional2(u,v,c,x)==x*(1-x)*(e2*p-det)/(p*p)
        bound &= conditional2(u,v,c,x)<=e2/(4*min(u,v))
    check('conditional_matrix_and_Gram_identity_'+str(i),equality)
    check('restricted_exchange_error_bound_'+str(i),bound)
    e0=cabs2(cminus(c,v));e1=cabs2(cminus(c,u));dd=u-v;eps2=F(1,10000)
    expanded=polyadd([eps2*z for z in polymul([v,dd],[v,dd])],
                     [-z for z in polymul([F(0),F(1),F(-1)],[e0,e1-e0])])
    check('cubic_certificate_coefficients_'+str(i),expanded==[
          eps2*v*v,2*eps2*v*dd-e0,eps2*dd*dd-e1+2*e0,e1-e0])

check('equal_probabilities_dephasing_error',conditional2(F(1,4),F(1,4),(F(0),F(0)),F(1,2))==F(1,4))
check('equal_probabilities_SWAP_error',conditional2(F(1,4),F(1,4),(F(-1,4),F(0)),F(1,2))==1)
check('complex_phase_exchange_bound_sharp',conditional2(F(1,4),F(1,4),(F(3,20),F(1,5)),F(1,2))==F(1,5))
check('reference_filtered_failure_not_scalar',
      (1-F(3,10))*F(1,2)!= (1-F(1,4))*F(1,2))
check('reference_flags_retain_failure_difference',
      abs(F(1,20)*F(1,2))+abs(F(-1,20)*F(1,2))==F(1,20)
      and F(1,20)*F(1,2)+F(-1,20)*F(1,2)==0)
# A common environment rotation: vectors (.5,0) and (0,.5), identically
# for both parity sectors, have equal Gram matrices but orthogonal support.
check('common_environment_rotation_instrument_equal',
      exact_success_norm2(F(1,4),F(1,4),(F(1,4),F(0)),F(1,4))[1]==0
      and F(1,2)**2+F(1,2)**2==F(1,2))
check('far_reference_but_direct_perfect_channel_possible',
      exact_success_norm2(F(1,2),F(1,2),(F(1,2),F(0)),F(1,4))[1]==F(1,16)
      and conditional2(F(1,2),F(1,2),(F(1,2),F(0)),F(1,2))==0)
check('rare_perfect_channel_fails_probability',F(1,100)<F(1,10)
      and conditional2(F(1,100),F(1,100),(F(1,100),F(0)),F(1,2))==0)
p0=F(514089,4000000);q0=F(434281,1000000);eta=F(1,1000)
check('unnormalized_map_output_margin',p0-eta==F(510089,4000000)>F(1,10))
check('retention_probability_margin',q0-eta==F(433281,1000000)>F(1,4))
check('normalization_cost_paid',eta/(p0-eta)==F(4000,510089)<F(1,100))
check('success_operation_budget_ceiling',(p0/F(101))/(p0-p0/F(101))==F(1,100))


class State(NamedTuple):
    word:tuple
    plus:tuple
    minus:tuple
    fuels:tuple


N=4;root=0;b=2;sites=(root,b)


def fields(s):
    charge=0;out=[]
    for x in s.word:
        charge+=int(x in s.plus)-int(x in s.minus);out.append(charge)
    assert charge==0
    return tuple(out[:-1])


def valid(s):
    return (len(s.plus)==len(s.minus)==2-len(s.fuels)
            and max(map(abs,fields(s)),default=0)<=1)


def carrier_order(s):
    tokens=[('p',x) for x in s.plus]+[('f',x) for x in s.fuels]
    def key(token):
        kind,x=token
        tie=(0 if kind=='f' else 1) if x==root else (1 if kind=='f' else 0)
        return s.word.index(x),tie
    return sorted(tokens,key=key)


def swap_parity(s,t,moved=None):
    tokens=carrier_order(s)
    mapped=[moved.get(z,z) if moved else z for z in tokens]
    target=carrier_order(t)
    assert set(mapped)==set(target) and len(mapped)==2
    return int(mapped!=target)


states=[]
for length in range(1,N+1):
    for suffix in permutations(range(1,N),length-1):
        word=(root,)+suffix
        allowed=tuple(x for x in sites if x in word)
        for nf in range(3):
            for f in combinations(allowed,nf):
                k=2-nf
                for p in combinations(word,k):
                    for m in combinations(word,k):
                        s=State(word,tuple(sorted(p)),tuple(sorted(m)),tuple(sorted(f)))
                        if valid(s):states.append(s)
catalog=set(states);exchange={s:set() for s in states}
plus_edges=minus_edges=conversion_edges=geometry_edges=0
classification=True;reverse=True;gate_support=True
for s in states:
    w,p,m,f=s
    for slots,species in ((p,'p'),(m,'m')):
        for x,y in zip(w,w[1:]):
            for old,new in ((x,y),(y,x)):
                if old not in slots or new in slots:continue
                moved=tuple(sorted(new if z==old else z for z in slots))
                t=State(w,moved,m,f) if species=='p' else State(w,p,moved,f)
                if not valid(t):continue
                parity=swap_parity(s,t,{('p',old):('p',new)} if species=='p' else None)
                expected=int(species=='p' and b in f and x==b)
                classification &= parity==expected
                if species=='p':plus_edges+=1
                else:minus_edges+=1
                if parity:
                    exchange[s].add(t)
                    gate_support &= b in f and root not in f and len(p)==1
    for x in sites:
        if x not in w:continue
        if x in f:
            if x in p or x in m:continue
            t=State(w,tuple(sorted(p+(x,))),tuple(sorted(m+(x,))),tuple(y for y in f if y!=x))
            moved={('f',x):('p',x)}
        else:
            if x not in p or x not in m:continue
            t=State(w,tuple(y for y in p if y!=x),tuple(y for y in m if y!=x),tuple(sorted(f+(x,))))
            moved={('p',x):('f',x)}
        classification &= swap_parity(s,t,moved)==0;conversion_edges+=1
    for v in range(N):
        if v in w:continue
        t=s._replace(word=w+(v,))
        classification &= swap_parity(s,t)==0;geometry_edges+=1
    if len(w)>1 and w[-1] not in p+m+f:
        t=s._replace(word=w[:-1])
        classification &= swap_parity(s,t)==0;geometry_edges+=1
for s,targets in exchange.items():
    for t in targets:reverse &= s in exchange[t]
check('all_color_primitives_classified',classification)
check('exchange_gate_partial_matching',all(len(v)<=1 for v in exchange.values()) and reverse)
check('crossing_requires_other_resource_in_matter',gate_support)
check('exchange_annihilates_retention_and_success',all(s.fuels
      and all(t.fuels for t in ts) for s,ts in exchange.items() if ts))
old=State((0,1,2),(2,),(2,),(2,))
grown=old._replace(word=(0,1,2,3))
check('color_identity_growth_can_enable_exchange',not exchange[old]
      and bool(exchange[grown]) and swap_parity(old,grown)==0)
initial=State((0,1,2),(),(),(0,2))
check('actual_two_fuel_seed_has_no_exchange',initial in catalog and not exchange[initial])

for item in json.loads((HERE/'SOURCE_PINS.json').read_text())['files']:
    check('pinned_'+item['path'].split('/')[-2]+'_'+item['path'].split('/')[-1],
          hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()==item['sha256'])
result={'status':'all_frozen_controls_pass' if all(c['pass'] for c in checks) else 'frozen_control_failure',
        'checks':checks,'kinematic_states':len(states),
        'directed_primitive_counts':{'plus':plus_edges,'minus':minus_edges,
                                    'conversion':conversion_edges,'geometry':geometry_edges,
                                    'exchange':sum(len(x) for x in exchange.values())},
        'arithmetic':'Fraction only; no floating-point decisions',
        'role':'exact Gram algebra and kinematic primitive-color classification',
        'actual_moving_Gram_entries_computed':False,'physical_dynamics_runs':0,
        'operating_gate':'unresolved','independent_review':'pending',
        'illustrative_channels_are_actual_B2_counterexamples':False}
print(json.dumps(result,indent=2,sort_keys=True))
raise SystemExit(0 if result['status']=='all_frozen_controls_pass' else 1)
