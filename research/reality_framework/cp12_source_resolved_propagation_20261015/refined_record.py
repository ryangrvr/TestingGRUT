"""Complete core-and-host instrument: probabilities, work and record entropy.

Only Q=N+1 is covered. The preparation remains the original coherent
neutral sea; all core-host outcomes and the failure outcome are included.
"""
from math import sqrt,lcm
import argparse,json,time
import numpy as np
from shared import *


def entropy_term_interval(lo,hi,count=1):
    def term(p):return -p*np.log2(p/count) if p>0 else 0.
    values=[term(lo),term(hi)];critical=count/np.e
    if lo<=critical<=hi:values.append(term(critical))
    return [min(term(lo),term(hi)),max(values)]


def calculate(n,q,budget,wait):
    start=time.monotonic();assert q==n+1
    m=quotient_model(n,q,budget);sub=m['sub'];H=m['H']
    z=np.load(cache_root()/f'N{n}_Q{q}_B{budget}_W{wait:g}.npz')
    w=z['w'];mask=z['mask'];cores=z['core'];e=float(z['state_error'])
    classes={r:np.zeros(len(w),bool) for r in range(2,n//2+1)}
    for col in np.flatnonzero(mask):
        s=m['states'][col];extra=s[0]^int(cores[col]);assert extra.bit_count()==1
        i,j=sub.es[extra.bit_length()-1];r=distances(n,int(cores[col]),i)[j]
        classes[r][col]=True
    assert np.array_equal(sum(a.astype(int) for a in classes.values()),mask)
    p=float(np.linalg.norm(w*mask)**2);pI=probability_interval(p,e)
    X=w*mask;binary=-2*float(np.vdot(X,H@((1-mask)*w)).real)
    by={m['states'][j]:int(j) for j in np.flatnonzero(mask)}
    # Raw rewiring changes g and leaves the complete electric pattern (and
    # hence C) fixed. Canonical orbit labels alone cannot detect this part:
    # it can contribute to a quotient diagonal after relabelling.
    extra_work=0j;count_transitions=0;sc=lcm(8,10*sub.M)
    for col in np.flatnonzero(mask):
        s=m['states'][col]
        for t,c in transitions(n,s,sc):
            if t[0]==s[0]:continue
            assert t[1]==s[1] and t[2]==s[2]
            key,sign,size=sub.canonical(t)
            if not sign:continue
            row=by[key]
            amplitude=c/sc*sign*sqrt(m['sizes'][col]/size)
            extra_work-=np.conj(w[row])*amplitude*w[col];count_transitions+=1
    assert abs(extra_work.imag)<1e-11
    extra_work=float(extra_work.real);total=binary+extra_work
    knorm=sqrt(2*q*min(n,q));rnorm=q*(sub.M-q)/(10*sub.M)
    shape=2*e+e*e;binary_error=2*knorm*shape+ALLOWANCE
    extra_error=rnorm*shape+ALLOWANCE
    records=[];entropy=entropy_term_interval(1-p,1-p)[0]
    failI=entropy_term_interval(1-pI[1],1-pI[0]);entropyI=failI[:]
    cycle_count=factorial(n-1)//2
    for r,selected in classes.items():
        pr=float(np.linalg.norm(w*selected)**2);prI=probability_interval(pr,e)
        count=cycle_count*(n//2 if 2*r==n else n)
        term=entropy_term_interval(prI[0],prI[1],count)
        entropy+=entropy_term_interval(pr,pr,count)[0]
        entropyI=[entropyI[0]+term[0],entropyI[1]+term[1]]
        records.append(dict(chord_core_distance=r,class_probability=pr,class_probability_interval=prI,
            equivalent_records=count,probability_one_labelled_record=pr/count,
            probability_one_labelled_record_interval=[a/count for a in prI],
            fixed_source_is_chord_endpoint_probability=2*pr/n,
            fixed_source_is_chord_endpoint_probability_interval=[2*a/n for a in prI]))
    result=dict(status='PASS',N=n,Q=q,budget=budget,wait=wait,state_error=e,
        all_cycle_probability=p,all_cycle_probability_interval=pI,classes=records,
        total_successful_records=sum(r['equivalent_records'] for r in records),
        failure_probability=1-p,complete_outcomes=1+sum(r['equivalent_records'] for r in records),
        core_only_record_work=binary,additional_host_record_work=extra_work,
        additional_host_record_work_interval=[extra_work-extra_error,extra_work+extra_error],
        complete_record_work=total,complete_record_work_interval=[total-binary_error-extra_error,total+binary_error+extra_error],
        complete_record_entropy_bits=entropy,complete_record_entropy_interval_bits=entropyI,
        raw_rewiring_transition_count=count_transitions,
        instrument='All projectors P_(C,g) for exact spanning-cycle cores and their complete initial active hosts, plus one aggregate failure projector. No matter or electric-state outcome is selected.',
        work_derivation='Core-only work plus minus the original expectation of rewiring between different hosts within successful cores. K preserves g, R preserves C, and distinct exact spanning cycle masks have no one-step coupling. Raw transitions retain rewiring terms even when a permutation quotient maps them to the same orbit.',
        entropy_derivation='-(1-p) log2(1-p) - sum_r p_r log2(p_r/M_r); M_r counts all equivalent core-host records. This is record entropy, not a dissipated-heat claim.',
        error_scope='Analytic full-state norm bound, hopping norm and complete rewiring degree bound; ordinary floating point with5e-10 allowance for each work contribution.',
        seconds=time.monotonic()-start)
    (P/f'REFINED_RECORD_N{n}_Q{q}_B{budget}_W{wait:g}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,required=True)
    ap.add_argument('--q',type=int,required=True);ap.add_argument('--budget',type=int,required=True)
    ap.add_argument('--wait',type=float,default=5)
    a=ap.parse_args();calculate(a.n,a.q,a.budget,a.wait)
