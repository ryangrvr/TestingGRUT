"""Exact LOCAL algebra and conditional margin checks, not a dynamics run."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
checks=[]


def check(name,value):
    checks.append({'name':name,'pass':bool(value)})


def rational(x):
    x=F(x)
    return {'numerator':str(x.numerator),'denominator':str(x.denominator),
            'display':float(x)}


jp,jm,J,A,kappa=F(1,64),F(1,96),F(1,3840),F(13,46080),F(1,4)
Omega=2-2*A


def field(P,M,length):
    result=[];charge=0
    for x in range(length):
        charge+=int(x in P)-int(x in M)
        if x<length-1: result.append(charge)
    assert charge==0
    return tuple(result)


def states(length,sites):
    out=[]
    for fnum in range(3):
        for ff in combinations(sites,fnum):
            k=2-fnum
            for p in combinations(range(length),k):
                for m in combinations(range(length),k):
                    if max(map(abs,field(p,m,length)),default=0)<=1:
                        out.append((p,m,ff))
    return out


def diagonal(s,length):
    p,m,f=s
    density=sum(int(x in p and x+1 in m)+int(x in m and x+1 in p)
                for x in range(length-1))
    endpoint=sum(int(x in p and x in m) for x in (0,length-1))
    return (2*len(p)+Omega*len(f)+sum(e*e for e in field(p,m,length))
            +kappa*density-A*endpoint)


def neighbors(state,length,sites):
    p,m,f=state
    out={state:diagonal(state,length)}
    def add(s,a): out[s]=out.get(s,F(0))+a
    for species,slots,c in [(0,p,jp),(1,m,jm)]:
        for old in slots:
            for new in (old-1,old+1):
                if new<0 or new>=length or new in slots: continue
                moved=tuple(sorted(new if x==old else x for x in slots))
                pp,mm=(moved,m) if species==0 else (p,moved)
                if max(map(abs,field(pp,mm,length)),default=0)<=1:
                    add((pp,mm,f),-c)
    for x in sites:
        if x in f:
            if x in p or x in m: continue
            target=(tuple(sorted(p+(x,))),tuple(sorted(m+(x,))),
                    tuple(y for y in f if y!=x))
        else:
            if x not in p or x not in m: continue
            target=(tuple(y for y in p if y!=x),tuple(y for y in m if y!=x),
                    tuple(sorted(f+(x,))))
        add(target,-J)
    return out


# A four-event local host and a two-event retained prefix. This is an
# algebraic block witness; neither size is an operating benchmark.
length,sites=4,(0,1)
basis=states(length,sites)
catalog=set(basis)
matrix={s:neighbors(s,length,sites) for s in basis}
check('exact_local_sector_closed',all(t in catalog for col in matrix.values() for t in col))
check('exact_local_matrix_Hermitian',all(matrix[t].get(s,0)==a
      for s,col in matrix.items() for t,a in col.items()))
check('resources_conserved_all_local_states',all(len(p)==len(m)==2-len(f) for p,m,f in basis))
check('all_local_fields_capacity_one',all(max(map(abs,field(p,m,length)),default=0)<=1
      for p,m,f in basis))
source=((),(),(0,1))
plus_target=((2,),(1,),(0,))
minus_target=((1,),(2,),(0,))


def squared_entry(target):
    return sum(matrix[mid].get(target,F(0))*a for mid,a in matrix[source].items())


amp_plus=squared_entry(plus_target)
amp_minus=squared_entry(minus_target)
check('outward_plus_unique_H_squared_word',amp_plus==J*jp==F(1,245760))
check('outward_minus_unique_H_squared_word',amp_minus==J*jm==F(1,368640))
check('outward_plus_Gauss_field',field(plus_target[0],plus_target[1],4)==(0,-1,0))
check('outward_minus_Gauss_field',field(minus_target[0],minus_target[1],4)==(0,1,0))
check('targets_outside_reference_prefix',any(x>=2 for x in plus_target[0]+plus_target[1])
      and any(x>=2 for x in minus_target[0]+minus_target[1]))
co_right=((1,),(1,),(0,))
degree_delta=diagonal(co_right,4)-diagonal(co_right,2)
check('actual_interior_degree_discrepancy',degree_delta==A)
check('initial_mean_energy',matrix[source][source]==2*Omega)
variance=sum(a*a for t,a in matrix[source].items() if t!=source)
check('initial_energy_variance',variance==2*J*J)

# Forward append reads the new child, not old-head emptiness. An admissible
# receiver occupation is contact=0, another mode=2 on a length-three word;
# every other one of the 64N^2 physical modes has occupation zero.
N=5;word=(0,1,2);fuels=(0,2);child=3
P=M=();q=set(word);contact=0;g=F(1)
forward_allowed=(child not in word and child not in P and child not in M
                 and child not in fuels and child not in q)
pop_old_allowed=(word[-1] not in P and word[-1] not in M and word[-1] not in fuels)
new_word=word+(child,)
check('fueled_head_allows_forward_assisted_append',forward_allowed)
check('fueled_head_blocks_its_pop',not pop_old_allowed)
check('same_fuel_becomes_interior',fuels[-1]==2 and new_word[-1]!=2 and new_word.index(2)==len(new_word)-2)
append_coefficient=-g/F(N)  # sqrt(contact+1)=1 exactly
check('assisted_append_nonzero_at_fueled_head',append_coefficient==F(-1,5))
check('receiver_charge_bookkeeping',2==len(word)-1 and 3==len(new_word)-1)

# Time mismatch and complete conditional error budget; no error estimate
# for moving geometry is claimed by these margin calculations.
Tmin=F(2016)/(6*J)
check('B2_window_not_old_B1_window',Tmin==1290240 and Tmin/F(15360)==84)
check('illustrative_size_scaled_window_ratio',F(4096)/(24*J)==655360)
eps=F(1,1000)
p0amp,q0amp=F(717,2000),F(659,1000)
p_after=(p0amp-eps)**2
q_after=(q0amp-eps)**2
conditional_error=eps/(p0amp-eps)
check('conditional_composition_output_margin',p_after==F(20449,160000)>F(1,10))
check('conditional_composition_retention_margin',q_after==F(108241,250000)>F(1,4))
check('conditional_composition_reference_channel_margin',conditional_error==F(2,715)<F(1,100))
check('fixed_path_original_joint_margin',p0amp*p0amp==F(514089,4000000))
check('fixed_path_original_retention_margin',q0amp*q0amp==F(434281,1000000))

pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
pin_results=[]
for x in pins['files']:
    ok=hashlib.sha256((ROOT/x['path']).read_bytes()).hexdigest()==x['sha256']
    check('pinned_input_'+x['path'].split('/')[-2]+'_'+x['path'].split('/')[-1],ok)
    pin_results.append({'path':x['path'],'matches':ok})

ok=all(c['pass'] for c in checks)
result={'status':'all_frozen_controls_pass' if ok else 'frozen_control_failure',
        'role':'exact local algebra, archive identity and hypothetical error-margin checks',
        'dynamics_runs':0,'operating_length_simulations':0,
        'local_scalar_dimension':len(basis),'checks':checks,
        'bounds':{'H2_outward_plus':rational(amp_plus),'H2_outward_minus':rational(amp_minus),
                  'interior_degree_difference':rational(degree_delta),
                  'assisted_append_example':rational(append_coefficient),
                  'hypothetical_epsilon_geo':rational(eps),
                  'hypothetical_output_lower':rational(p_after),
                  'hypothetical_retention_lower':rational(q_after),
                  'hypothetical_conditional_trace_error':rational(conditional_error)},
        'pinned_inputs':pin_results,
        'fixed_path_status':'archived_author_derived_pass_independent_review_pending',
        'moving_geometry_status':'existing_hypotheses_incompatible_operating_claim_unresolved',
        'epsilon_geo_achieved':False}
print(json.dumps(result,indent=2,sort_keys=True))
raise SystemExit(0 if ok else 1)
