"""Exact primitive-neighborhood checks, not a physical evolution solver."""
from collections import defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from typing import NamedTuple
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
checks=[]


def check(name,ok):
    checks.append({'name':name,'pass':bool(ok)})


def squarefree(n):
    a=1;d=2
    while d*d<=n:
        while n%(d*d)==0:
            a*=d;n//=d*d
        d+=1
    return a,n


class Radical:
    """Finite exact sums of rational multiples of square roots of integers."""
    def __init__(self,x=0,terms=None):
        self.t={d:F(v) for d,v in (terms or {1:F(x)}).items() if v}
    def __add__(self,other):
        o=other if isinstance(other,Radical) else Radical(other)
        t=dict(self.t)
        for d,v in o.t.items(): t[d]=t.get(d,F(0))+v
        return Radical(terms=t)
    __radd__=__add__
    def __neg__(self): return Radical(terms={d:-v for d,v in self.t.items()})
    def __sub__(self,o): return self+-asrad(o)
    def __mul__(self,other):
        o=asrad(other);t=defaultdict(F)
        for d,a in self.t.items():
            for e,b in o.t.items():
                factor,key=squarefree(d*e);t[key]+=a*b*factor
        return Radical(terms=t)
    __rmul__=__mul__
    def __eq__(self,o): return self.t==asrad(o).t
    def __bool__(self): return bool(self.t)
    def rational(self):
        assert all(d==1 for d in self.t)
        return self.t.get(1,F(0))


def asrad(x): return x if isinstance(x,Radical) else Radical(x)
def sqrtint(n):
    if not n: return Radical()
    a,d=squarefree(n)
    return Radical(terms={d:F(a)})


class State(NamedTuple):
    word:tuple
    q:int
    plus:tuple
    minus:tuple
    fuels:tuple
    contact:int


jp,jm,J,Ad,eps=F(1,64),F(1,96),F(1,3840),F(13,46080),F(1,4)
g=F(1);N=4


def fields(s):
    charge=0;out=[]
    for x in s.word:
        charge+=int(x in s.plus)-int(x in s.minus)
        out.append(charge)
    assert charge==0
    return tuple(out[:-1])


def valid(s):
    if len(set(s.word))!=len(s.word) or s.word[0]!=0: return False
    if not s.q&1 or s.q>>N: return False
    if not (0<=s.contact<=s.q.bit_count()-1): return False
    if any(len(set(t))!=len(t) or not set(t)<=set(s.word)
           for t in (s.plus,s.minus,s.fuels)): return False
    if len(s.plus)!=len(s.minus) or len(s.plus)+len(s.fuels)!=2: return False
    return max(map(abs,fields(s)),default=0)<=1


def primitives(s,sites,geometry=True,matter=True):
    """Actual off-diagonal primitive maps; untouched receiver modes retained.

    Receiver occupations are contact=c and spectator=(Q-1-c), all remaining
    physical modes vacuum. No receiver hopping or diagonal is truncated
    into a dynamics model: these maps check boundary action/minimal words.
    Those length-preserving terms cannot enter a minimum-length target.
    """
    out={}
    def put(t,value):
        assert valid(t),(s,t)
        out[t]=out.get(t,Radical())+value
    w,q,p,m,f,c=s;h=w[-1]
    nh=int(h in p)+int(h in m)
    if geometry:
        for v in range(N):
            if v in w: continue
            # The physical support constraints make this free child empty.
            if not q>>v&1:
                put(State(w+(v,),q|(1<<v),p,m,f,c+1),-g/F(N)*sqrtint(c+1))
            if nh:
                put(State(w+(v,),q,p,m,f,c),-eps*F(nh,N)*sqrtint(N))
        if len(w)>1 and h not in p+m+f and fields(s)[-1]==0:
            if q>>h&1 and c:
                put(State(w[:-1],q^(1<<h),p,m,f,c-1),-g/F(N)*sqrtint(c))
            parent=w[-2];np=int(parent in p)+int(parent in m)
            if np:
                put(State(w[:-1],q,p,m,f,c),-eps*F(np,N)*sqrtint(N))
    if matter:
        for slots,coupling,species in ((p,jp,0),(m,jm,1)):
            for a,b in zip(w,w[1:]):
                for old,new in ((a,b),(b,a)):
                    if old not in slots or new in slots: continue
                    moved=tuple(sorted(new if x==old else x for x in slots))
                    pp,mm=(moved,m) if species==0 else (p,moved)
                    t=State(w,q,pp,mm,f,c)
                    if max(map(abs,fields(t)),default=0)<=1: put(t,-coupling)
        for x in sites:
            if x not in w: continue
            if x in f:
                if x in p or x in m: continue
                t=State(w,q,tuple(sorted(p+(x,))),tuple(sorted(m+(x,))),
                        tuple(y for y in f if y!=x),c)
            else:
                if x not in p or x not in m: continue
                t=State(w,q,tuple(y for y in p if y!=x),tuple(y for y in m if y!=x),
                        tuple(sorted(f+(x,))),c)
            put(t,-J)
    return {t:a for t,a in out.items() if a}


def fixed_basis(word,sites):
    out=[];q=sum(1<<x for x in word)
    for fn in range(3):
        for f in combinations(sites,fn):
            k=2-fn
            for p in combinations(word,k):
                for m in combinations(word,k):
                    for c in range(len(word)):
                        s=State(word,q,p,m,f,c)
                        if valid(s): out.append(s)
    return out


def form_expected(s):
    w,q,p,m,f,c=s;h=w[-1];parent=w[-2];r=N-len(w)
    nh=int(h in p)+int(h in m);np=int(parent in p)+int(parent in m)
    gh=int(h not in p+m+f and fields(s)[-1]==0)
    return g*g/F(N*N)*(r*(c+1)+c*gh)+eps*eps/F(N)*(r*nh*nh+gh*np*np)


word=(0,1,2);sites=(0,2)
basis=fixed_basis(word,sites)
images={s:primitives(s,sites,matter=False) for s in basis}
inverse=defaultdict(list)
for s,col in images.items():
    for t,a in col.items(): inverse[t].append((s,a))
gram=defaultdict(Radical)
for entries in inverse.values():
    for s,a in entries:
        for t,b in entries: gram[s,t]=gram[s,t]+a*b
check('boundary_quadratic_form_all_entries',all(
      gram[s,t]==(form_expected(s) if s==t else 0) for s in basis for t in basis))
check('boundary_primitives_exact_adjoints',all(
      primitives(t,sites,matter=False).get(s,Radical())==a
      for s,col in images.items() for t,a in col.items()))
check('boundary_images_are_outside_fixed_word_sector',all(t.word!=word
      for col in images.values() for t in col))
check('strict_boundary_injectivity_small_control',all(
      form_expected(s)>=g*g*F(N-len(word),N*N)>0 for s in basis))
seed=State(word,7,(),(),sites,0)
seed_edges=primitives(seed,sites,matter=False)
check('two_fuel_seed_only_assisted_growth',len(seed_edges)==1 and all(
      t.word==(0,1,2,3) and t.q==15 and t.contact==1 for t in seed_edges))
seed_rates=[]
for c in range(3):
    s=seed._replace(contact=c)
    rate=sum((a*a).rational() for a in primitives(s,sites,matter=False).values())
    seed_rates.append(str(rate))
    check('initial_contact_rate_'+str(c),rate==g*g*F(N-3,N*N)*(c+1))


def power_from(source,sites,k):
    v={source:Radical(1)};counts=[1]
    for _ in range(k):
        out=defaultdict(Radical)
        for s,x in v.items():
            for t,a in primitives(s,sites).items(): out[t]=out[t]+a*x
        v={s:x for s,x in out.items() if x};counts.append(len(v))
    return v,counts


third,counts3=power_from(seed,sites,3)
assist_target=State((0,1,2,3),15,(3,),(2,),(0,),1)
bare_target=State((0,1,2,3),7,(3,),(2,),(0,),0)
check('seed_assisted_growth_conversion_outward_H3',
      third.get(assist_target,Radical())==-2*g*J*jp/F(N))
check('seed_bare_growth_conversion_outward_H3',
      third.get(bare_target,Radical())==-2*eps*J*jp/F(N)*sqrtint(N))
check('outward_targets_obey_Gauss_capacity',fields(assist_target)==(0,0,-1)
      and fields(bare_target)==(0,0,-1))
full_seed=State((0,1,2,3),15,(),(),(0,3),0)
check('fully_assembled_seed_initial_leakage_zero',
      not primitives(full_seed,(0,3),matter=False))
fourth,counts4=power_from(full_seed,(0,3),4)
retreat_target=State((0,1,2),15,(2,),(2,),(0,),0)
check('full_record_minimal_bare_retreat_H4',
      fourth.get(retreat_target,Radical())==J*jp*jm/F(N)*sqrtint(N))
check('full_record_retreat_storage_error_retained',
      valid(retreat_target) and (retreat_target.q>>3&1) and 3 not in retreat_target.word)

def matter_diagonal(s):
    density=sum(int(a in s.plus and b in s.minus)+int(a in s.minus and b in s.plus)
                for a,b in zip(s.word,s.word[1:]))
    degree={x:int(i>0)+int(i+1<len(s.word)) for i,x in enumerate(s.word)}
    counter=sum(2-degree[x] for x in s.plus if x in s.minus)
    return (2*len(s.plus)+(2-2*Ad)*len(s.fuels)+sum(e*e for e in fields(s))
            +F(1,4)*density-Ad*counter)


def storage_diagonal(s):
    return 16*(s.q.bit_count()-len(s.word))+64*sum(
           int(bool(s.q>>i&1)!=(i in s.word)) for i in range(N))


old_co=State(word,7,(2,),(2,),(0,),0)
grown_co=State((0,1,2,3),15,(2,),(2,),(0,),1)
drive=defaultdict(Radical)
drive[grown_co]=Radical(matter_diagonal(grown_co)-matter_diagonal(old_co))
for t,a in primitives(grown_co,sites,geometry=False).items(): drive[t]=drive[t]+a
for t,a in primitives(old_co,sites,geometry=False).items():
    embedded=t._replace(word=(0,1,2,3),q=15,contact=1)
    drive[embedded]=drive[embedded]-a
drive={s:a for s,a in drive.items() if a}
defect_sq=sum((a*a).rational() for a in drive.values())
check('post_append_source_defect_three_orthogonal_outputs',len(drive)==3
      and defect_sq==Ad*Ad+jp*jp+jm*jm and drive[grown_co]==Ad)
bare_co=grown_co._replace(q=7,contact=0)
check('bare_append_full_diagonal_ledger',storage_diagonal(bare_co)
      +matter_diagonal(bare_co)-matter_diagonal(old_co)==48+Ad)
before_retreat=retreat_target._replace(word=(0,1,2,3))
check('bare_pop_full_diagonal_ledger',storage_diagonal(retreat_target)
      +matter_diagonal(retreat_target)-matter_diagonal(before_retreat)==80-Ad)
check('assisted_append_storage_consistent',storage_diagonal(grown_co)==0)

# A permutation-symmetric two-color Gram control uses exact vectors and
# verifies that success probability need not be input independent.
chi_plus=(F(3,5),F(0));chi_minus=(F(0),F(4,5))
pplus=sum(x*x for x in chi_plus);pminus=sum(x*x for x in chi_minus)
coherence=sum(x*y for x,y in zip(chi_plus,chi_minus))
check('Gram_control_distinguishes_probability_and_coherence',
      pplus==F(9,25) and pminus==F(16,25) and coherence==0)
state_ceiling=F(717,202000)
check('information_state_error_ceiling',state_ceiling/(F(717,2000)-state_ceiling)==F(1,100))
check('ceiling_also_preserves_output', (F(717,2000)-state_ceiling)**2>F(1,10))
check('ceiling_also_preserves_retention',(F(659,1000)-state_ceiling)**2>F(1,4))
check('parent_comparison_target_is_stricter',F(1,1000)<state_ceiling)
check('absolute_Duhamel_certificate_obstruction_constant',
      F(1)/(6*J)*F(1,2)*F(9,10)==288)

pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
for item in pins['files']:
    check('pinned_'+item['path'].split('/')[-2]+'_'+item['path'].split('/')[-1],
          hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()==item['sha256'])
result={'status':'all_frozen_controls_pass' if all(x['pass'] for x in checks)
        else 'frozen_control_failure','checks':checks,
        'local_fixed_sector_dimension':len(basis),'primitive_image_states':len(inverse),
        'H3_neighborhood_counts':counts3,'H4_neighborhood_counts':counts4,
        'seed_contact_squared_rates':seed_rates,
        'assisted_H3_outward_plus':str(third[assist_target].rational()),
        'bare_H3_outward_plus':str(third[bare_target].rational()),
        'full_record_bare_H4_retreat':str(fourth[retreat_target].rational()),
        'post_append_defect_squared':str(defect_sq),
        'role':'exact primitive algebra and minimum-word coefficients, not time evolution',
        'dynamics_runs':0,'operating_length_simulations':0,
        'receiver_replaced_in_theorem':False,'late_time_error_bound_achieved':False,
        'verdict':'exact_endpoint_embedding_excluded_actual_operation_unresolved'}
print(json.dumps(result,indent=2,sort_keys=True))
raise SystemExit(0 if result['status']=='all_frozen_controls_pass' else 1)
