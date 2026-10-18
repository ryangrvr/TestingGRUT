"""Frozen exact finite-order calculation, not a time trajectory or new bath."""
from collections import Counter
from fractions import Fraction as F
from math import factorial,isqrt
from pathlib import Path
from typing import NamedTuple
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
N=4;MODES=64*N*N;ROOT_LABEL=0;RIGHT=1
JP=F(1,64);JM=F(1,96);J=F(1,3840);A=F(13,46080);OMEGA=2-2*A
checks=[]
def check(name,ok):
    checks.append({'name':name,'pass':bool(ok)})
def ratio(x):
    x=F(x)
    return {'numerator':str(x.numerator),'denominator':str(x.denominator)}
def bit(x,i):return (x>>i)&1
def wordmask(w):return sum(1<<i for i in w)

class State(NamedTuple):
    word:tuple
    q:int
    plus:int
    minus:int
    fuel:int
    photons:tuple
    parity:int

def fields(s):
    c=0;out=[]
    for x in s.word:
        c+=bit(s.plus,x)-bit(s.minus,x);out.append(c)
    return out
def valid(s):
    attached=wordmask(s.word)
    return (len(s.word)==len(set(s.word)) and s.word[0]==0
            and 1<=len(s.word)<=N and all(0<=i<N for i in s.word)
            and bit(s.q,0)==1 and 0<=s.q<1<<N
            and (s.plus|s.minus|s.fuel)&~attached==0
            and s.fuel&~3==0 and s.plus.bit_count()==s.minus.bit_count()==2-s.fuel.bit_count()
            and len(s.photons)==s.q.bit_count()-1
            and tuple(sorted(s.photons))==s.photons
            and all(0<=k<MODES for k in s.photons)
            and fields(s)[-1]==0 and max(map(abs,fields(s)),default=0)<=1
            and s.parity in (0,1))
def metric(s):
    out=1
    for count in Counter(s.photons).values():out*=factorial(count)
    return out
def colors(s):
    result=[]
    for x in s.word:
        if x==0 and bit(s.fuel,x):result.append(('f',x))
        if bit(s.plus,x):result.append(('p',x))
        if x!=0 and bit(s.fuel,x):result.append(('f',x))
    assert len(result)==2
    return result
def target(s,t,moved=None):
    before=[moved.get(z,z) if moved else z for z in colors(s)]
    after=colors(t)
    assert set(before)==set(after)
    t=t._replace(parity=s.parity^int(before!=after))
    assert valid(t)
    return t
def diagonal(s):
    attached=wordmask(s.word);q=s.q.bit_count();m=len(s.word)
    d=F(16*(q-m)+64*(s.q^attached).bit_count())
    d+=sum(v*v for v in fields(s)[:-1])
    d+=s.plus.bit_count()+s.minus.bit_count()+OMEGA*s.fuel.bit_count()
    for x,y in zip(s.word,s.word[1:]):
        d+=F(bit(s.plus,x)*bit(s.minus,y)+bit(s.minus,x)*bit(s.plus,y),4)
    for i,x in enumerate(s.word):
        deg=int(i>0)+int(i+1<len(s.word))
        d-=A*(2-deg)*bit(s.plus,x)*bit(s.minus,x)
    return d

def neighbors(s):
    """Every non-self primitive transition of the actual full finite law.

    diagonal(s) is retained separately. Self-loops cannot enter a shortest
    non-self path, so excluding them here does not delete a Taylor term at
    the minimal distance. No physical receiver mode or charged state cut.
    """
    w,q,p,m,f,photons,parity=s
    occ=Counter(photons)
    for old,count in occ.items():
        for new in (old-1,old+1):
            if not 0<=new<MODES:continue
            v=list(photons);v.remove(old);v.append(new)
            yield s._replace(photons=tuple(sorted(v))),F(-count),'receiver'
    for species,slots,coupling in (('p',p,JP),('m',m,JM)):
        for x,y in zip(w,w[1:]):
            for old,new in ((x,y),(y,x)):
                if not bit(slots,old) or bit(slots,new):continue
                moved_slots=slots^(1<<old)^(1<<new)
                t=s._replace(plus=moved_slots) if species=='p' else s._replace(minus=moved_slots)
                if not valid(t):continue
                moved={('p',old):('p',new)} if species=='p' else None
                yield target(s,t,moved),-coupling,'plus' if species=='p' else 'minus'
    for x in (0,RIGHT):
        if x not in w:continue
        mask=1<<x
        if bit(f,x):
            if bit(p,x) or bit(m,x):continue
            t=s._replace(plus=p|mask,minus=m|mask,fuel=f^mask)
            moved={('f',x):('p',x)}
        else:
            if not bit(p,x) or not bit(m,x):continue
            t=s._replace(plus=p^mask,minus=m^mask,fuel=f|mask)
            moved={('p',x):('f',x)}
        yield target(s,t,moved),-J,'conversion'
    head=w[-1]
    nh=bit(p,head)+bit(m,head)
    for v in range(N):
        if v in w:continue
        assert not bit(p|m|f,v)
        if not bit(q,v):
            t=s._replace(word=w+(v,),q=q|(1<<v),photons=tuple(sorted(photons+(0,))))
            yield target(s,t),F(-1,4),'assisted_append'
        if nh:
            t=s._replace(word=w+(v,))
            yield target(s,t),F(-nh,8),'bare_append'
    if len(w)>1 and not bit(p|m|f,head):
        if bit(q,head) and occ[0]:
            v=list(photons);v.remove(0)
            t=s._replace(word=w[:-1],q=q^(1<<head),photons=tuple(v))
            yield target(s,t),F(-occ[0],4),'assisted_pop'
        parent=w[-2];np=bit(p,parent)+bit(m,parent)
        if np:
            t=s._replace(word=w[:-1])
            yield target(s,t),F(-np,8),'bare_pop'

seed=State((0,1),3,0,0,3,(0,),0)
alpha=State((0,1),3,3,3,0,(0,),0)
beta=State((0,1,2),7,3,3,0,(0,0),0)
targets={'identity_alpha':alpha,'identity_beta':beta,
         'exchange_alpha':alpha._replace(parity=1),
         'exchange_beta':beta._replace(parity=1)}
check('seed_and_all_target_states_physical',valid(seed) and all(valid(x) for x in targets.values()))
dist={seed:0};level={seed:F(1)};coeff={};layers=[];edge_counts=Counter()
for depth in range(9):
    layers.append({'depth':depth,'states':len(level)})
    for name,t in targets.items():
        if t in level:coeff[name]=level[t]
    if depth==8:break
    following={}
    for s,amplitude in level.items():
        assert valid(s)
        for t,weight,kind in neighbors(s):
            assert weight<0 and valid(t)
            edge_counts[kind]+=1
            if t not in dist:dist[t]=depth+1
            if dist[t]==depth+1:
                following[t]=following.get(t,F(0))+amplitude*weight
    level=following

expected_dist={'identity_alpha':2,'identity_beta':3,'exchange_alpha':8,'exchange_beta':7}
expected_coeff={'identity_alpha':2*J**2,'identity_beta':-F(3,2)*J**2,
                'exchange_beta':-F(3,4)*J**2*JP**4,
                'exchange_alpha':F(13,16)*J**2*JP**4}
for name,t in targets.items():
    check('minimal_distance_'+name,dist.get(t)==expected_dist[name])
    check('coherent_shortest_coefficient_'+name,coeff.get(name)==expected_coeff[name])

# The seven-step witness, then exact reverse links, includes all physical
# modes in the generator. The graph enumeration above includes other words.
steps=[('assisted_append',lambda s:s.word==(0,1,2)),
       ('conversion',lambda s:s.plus==s.minus==1 and s.fuel==2),
       ('plus',lambda s:s.plus==2),('plus',lambda s:s.plus==4),
       ('conversion',lambda s:s.plus==6 and s.minus==3 and s.fuel==0),
       ('plus',lambda s:s.plus==5),('plus',lambda s:s.plus==3)]
s=seed;path=[];path_product=F(1)
for kind,predicate in steps:
    matches=[(t,v,k) for t,v,k in neighbors(s) if k==kind and predicate(t)]
    assert len(matches)==1
    t,v,k=matches[0]
    reverse=sum(value for back,value,_ in neighbors(t) if back==s)
    check('weighted_Hermiticity_witness_'+str(len(path)+1),v*metric(t)==reverse*metric(s))
    path.append({'kind':k,'word':list(t.word),'q':t.q,'plus':t.plus,
                 'minus':t.minus,'fuel':t.fuel,'photons':list(t.photons),
                 'parity':t.parity,'coefficient':ratio(v)})
    path_product*=v;s=t
check('odd_witness_beta_target',s==targets['exchange_beta'])
check('odd_witness_product',path_product==-F(1,4)*J**2*JP**4)
pop=[(t,v) for t,v,k in neighbors(s) if k=='assisted_pop' and t==targets['exchange_alpha']]
check('assisted_pop_witness_alpha',len(pop)==1 and pop[0][1]==-F(1,2))
check('source_diagonal_retained',diagonal(seed)==2*OMEGA)
check('normalized_Fock_metric',metric(alpha)==1 and metric(beta)==2)
norm_bound=48+192+6+3+4+4+1+4*A+2+F(5,48)+2*J
check('full_Hamiltonian_norm_bound',norm_bound<512)

# Exact Choi-state controls for complete finite-dimensional CP operations.
# Gaussian rationals; these illustrative channels are not B2 observations.
def z(a=0,b=0):return F(a),F(b)
def add(a,b):return a[0]+b[0],a[1]+b[1]
def mul(a,b):return a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]
def conj(a):return a[0],-a[1]
def scale(a,s):return a[0]*s,a[1]*s
def abs2(a):return a[0]**2+a[1]**2
def total(xs):
    out=z()
    for x in xs:out=add(out,x)
    return out
def sqrt_rat(x):
    a=isqrt(x.numerator);b=isqrt(x.denominator)
    assert a*a==x.numerator and b*b==x.denominator
    return F(a,b)
examples={
 'identity':[(z(F(1,2)),z(F(1,2)))],
 'unitary_SWAP':[(z(F(1,2)),z(F(-1,2)))],
 'parity_dephasing':[(z(F(1,2)),z()),(z(),z(F(1,2)))],
 'coherent_phase':[(z(F(1,2)),z(F(3,10),F(2,5)))],
 'parity_filter':[(z(F(1,2)),z(F(1,4)))],
 'rare_identity':[(z(F(1,10)),z(F(1,10)))]}
channel_results={}
for name,kraus in examples.items():
    u=sum(abs2(a) for a,b in kraus);v=sum(abs2(b) for a,b in kraus)
    c=total(mul(a,conj(b)) for a,b in kraus);cabs=sqrt_rat(abs2(c))
    check('trace_nonincrease_'+name,0<=u<=1 and 0<=v<=1 and abs2(c)<=u*v)
    phase=scale(c,1/cabs) if cabs else z(1)
    corrected=[(a,mul(phase,b)) for a,b in kraus]
    # Full16x16 output+reference Choi matrix, retaining off-diagonal phases.
    choi=[[z() for _ in range(16)] for _ in range(16)]
    for a,b in corrected:
        d=[a,a,a,b]
        for i in range(4):
            for j in range(4):
                choi[5*i][5*j]=add(choi[5*i][5*j],scale(mul(d[i],conj(d[j])),F(1,4)))
    q=sum(choi[i][i][0] for i in range(16))
    fidelity=scale(total(choi[5*i][5*j] for i in range(4) for j in range(4)),F(1,4)/q)
    formula=(9*u+v+6*cabs)/(4*(3*u+v))
    check('reference_Choi_optimum_'+name,fidelity==(formula,F(0)) and q==(3*u+v)/4)
    # Rank3 entangled input: all9 nonzero Choi entries equal u/3.
    qutrit=[[total(scale(mul(a,conj(a)),F(1,3)) for a,b in kraus)
             for j in range(3)] for i in range(3)]
    check('qutrit_reference_identity_'+name,all(x==(u/3,F(0)) for row in qutrit for x in row))
    check('qutrit_success_independent_'+name,sum(qutrit[i][i][0] for i in range(3))==u)
    perfect=(u==v==cabs and u>0)
    check('perfect_recovery_iff_'+name,(formula==1)==perfect)
    channel_results[name]={'u':ratio(u),'v':ratio(v),
                          'c_real':ratio(c[0]),'c_imag':ratio(c[1]),
                          'Fmax':ratio(formula),'q_Choi':ratio(q),
                          'illustrative_not_B2_data':True}
check('dephasing_Choi_fidelity_ceiling',channel_results['parity_dephasing']['Fmax']==ratio(F(5,8)))
check('filter_Choi_fidelity_ceiling',channel_results['parity_filter']['Fmax']==ratio(F(49,52)))
check('rare_perfect_not_operating',F(1,100)<F(1,10))
check('two_necessary_1percent_conditions_balanced',F(49,100)*2==F(98,100)
      and (F(36,75)+F(37,75))<F(98,100))

for source in json.loads((HERE/'SOURCE_PINS.json').read_text())['files']:
    check('pin_'+source['path'].split('/')[-2]+'_'+source['path'].split('/')[-1],
          hashlib.sha256((ROOT/source['path']).read_bytes()).hexdigest()==source['sha256'])

c0=J**2*abs(coeff['exchange_beta'])/factorial(7)
t0=min(F(1,1024),c0*factorial(10)/(12*1024**10))
check('strict_Gram_determinant_interval',c0>0 and t0>0 and t0<=F(1,1024)
      and 6*1024**10*t0/F(factorial(10))<=c0/2)
result={'status':'pass' if all(x['pass'] for x in checks) else 'frozen_check_failure',
        'checks':checks,'diagnostic_N':N,'physical_receiver_modes':MODES,
        'graph_radius':8,'states_in_ball':len(dist),'layers':layers,
        'inspected_directed_edges':dict(edge_counts),
        'targets':{name:{'distance':dist.get(t),'H_power_coefficient':ratio(coeff[name]),
                         'Fock_norm_squared':metric(t)} for name,t in targets.items()},
        'witness':path,'C0':ratio(c0),'strict_noise_interval_upper_t0':ratio(t0),
        'Gram_determinant_lower_coefficient':ratio(c0**2),
        'Gram_determinant_lower_power':18,'full_H_norm_upper':512,
        'illustrative_channels':channel_results,'physical_time_trajectories':0,
        'actual_S9_receiver_used':False,'original_operating_gate':'unresolved',
        'all_size_formation_or_delivery_bound':False,
        'arithmetic':'exact Fraction; no floating-point decision',
        'external_validation':'none'}
print(json.dumps(result,indent=2,sort_keys=True))
raise SystemExit(0 if result['status']=='pass' else 1)
