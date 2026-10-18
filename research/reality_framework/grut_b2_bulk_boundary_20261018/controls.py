#!/usr/bin/env python3
"""Exact actual virtual-path and bulk-boundary coefficient checks; no dynamics."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
JP,JM=F(1,64),F(1,96)
A,J=F(13,46080),F(1,3840)
OMEGA=2-2*A
S,P=JP*JP+JM*JM,JP*JM


def virtual_words(n,steps,E,quadrant):
    result={}

    def walk(i,j,left,weight):
        if left==0:
            if i==j:
                result[i]=result.get(i,F(0))+weight
            return
        for di,dj,c in [(1,0,JP),(-1,0,JP),(0,1,JM),(0,-1,JM)]:
            a,b=i+di,j+dj
            if quadrant and min(a,b)<0:
                continue
            w=weight*c
            if left>1:
                if a==b:
                    continue
                q=2+abs(a-b)+(F(1,4) if abs(a-b)==1 else 0)
                w/=q-E
            elif a!=b:
                continue
            walk(a,b,left-1,w)

    walk(n,n,steps,F(1))
    return result


def main():
    checks=[]

    def check(name,condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    c_lower=F(13,4)-F(5,96)
    gap=c_lower-F(5,2)
    b=F(5,96)
    check('closed_operator_lower_bound',c_lower==F(307,96))
    check('closed_gap',gap==F(67,96))
    check('Sigma_norm',b*b/gap==F(25,6432))
    check('Sigma_derivative_norm',b*b/gap**2==F(25,4489))
    check('new_cut_distance',F(5,2)-OMEGA==F(11533,23040))
    new_tail=J*J/(F(5,2)-OMEGA)**2
    old_tail=J*J/(c_lower-OMEGA)**2
    check('new_tail_exact',new_tail==F(36,133010089))
    check('old_tail_exact',old_tail==F(36,27613**2))
    check('retain_stronger_old_tail',old_tail<new_tail)

    tables=[]
    energies=[('reference_two',F(2)),('interval_left',2-5*J)]
    for name,E in energies:
        a1,a2=F(13,4)-E,4-E
        for steps in (2,4):
            bulk={n:virtual_words(n,steps,E,False) for n in range(3)}
            half={n:virtual_words(n,steps,E,True) for n in range(3)}
            for n in range(3):
                expected={0:2*S/a1,1:2*P/a1,-1:2*P/a1} if steps==2 else {
                    0:2*(S*S+2*P*P)/(a1*a1*a2),
                    1:4*P*S/(a1*a1*a2),-1:4*P*S/(a1*a1*a2),
                    2:2*P*P/(a1*a1*a2),-2:2*P*P/(a1*a1*a2)}
                actual={m-n:v for m,v in bulk[n].items() if v}
                check(f'bulk_CF_coefficients_{name}_{steps}_n{n}',actual==expected)
                for m in range(3):
                    difference=bulk[n].get(m,F(0))-half[n].get(m,F(0))
                    if steps==2:
                        predicted=S/a1 if n==m==0 else F(0)
                    else:
                        numerator={(0,0):S*S+6*P*P,(1,1):S*S-2*P*P,
                                   (0,1):2*P*S,(1,0):2*P*S}.get((n,m),F(0))
                        predicted=numerator/(a1*a1*a2)
                    check(f'boundary_word_{name}_{steps}_{n}_{m}',difference==predicted)
            tables.append({'energy_name':name,'E':str(E),'primitive_steps':steps,
                           'bulk_positive_words':{str(n):{str(m):str(v) for m,v in sorted(vals.items())} for n,vals in bulk.items()},
                           'quadrant_positive_words':{str(n):{str(m):str(v) for m,v in sorted(vals.items())} for n,vals in half.items()}})

    check('W2_PSD_diagonal',S*S-2*P*P>0)
    determinant=(S*S+6*P*P)*(S*S-2*P*P)-(2*P*S)**2
    check('W2_PSD_determinant',determinant==S**4-12*P**4 and determinant>0)
    q=F(5,72);C=F(75,20636)
    check('generic_even_tail_prefactor',b*b/F(3,4)/(1-q*q)==C)
    qb=F(1,24);t=qb*qb;Cb=F(1,460)
    check('band_even_tail_prefactor',b*b/F(5,4)/(1-t)==Cb)
    S1=9*t*t+t**3*(F(7)/(1-t)+2*t/(1-t)**2)
    beta_upper=4*A*J
    trace_W2_upper=(2*S*S+4*P*P)/(F(5,4)**2*2)
    trace_remainder_upper=Cb*S1
    total=beta_upper+trace_W2_upper+trace_remainder_upper
    check('frozen_trace_norm_target',total<J/500)
    check('reference_root_cancellation',S/F(5,4)==A)
    check('reference_pair_hopping',2*P/F(5,4)==J)
    check('band_lower_bracket',OMEGA-F(25,12)*J==2-F(17,4)*J)
    check('band_upper_bracket',OMEGA+F(25,12)*J==2-J/12)
    check('band_bracket_inside_declared_interval',2-5*J<OMEGA-F(25,12)*J and OMEGA+F(25,12)*J<2)

    result={'status':'all_frozen_controls_pass','checks':checks,
            'arithmetic':'exact Fraction; no time evolution, spectral quadrature or operating probabilities',
            'word_tables':tables,
            'bounds':{'Sigma_norm':str(F(25,6432)),'Sigma_derivative':str(F(25,4489)),
                      'new_tail':str(new_tail),'retained_stronger_tail':str(old_tail),
                      'root_beta_upper':str(beta_upper),'W2_trace_upper':str(trace_W2_upper),
                      'higher_boundary_trace_upper':str(trace_remainder_upper),
                      'total_boundary_trace_upper':str(total),'total_boundary_trace_upper_over_J':str(total/J),
                      'target':'J/500','target_pass':True},
            'operating_gate':'unresolved','actual_source_violation':None}
    out=ROOT/'CONTROLS.json'
    out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'checks':len(checks),
                      'word_tables':len(tables),'output_sha256':hashlib.sha256(out.read_bytes()).hexdigest()}))


if __name__=='__main__':
    main()
