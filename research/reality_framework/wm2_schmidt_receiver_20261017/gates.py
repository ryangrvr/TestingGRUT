"""Complete charge-block gates and finite operator-error audits."""
from itertools import combinations_with_replacement
from math import ceil,cos,sin,sqrt,factorial
import numpy as np
import bounds as E
from bounds import S5
import rotation as S7ROT

A=1.3512071919596578
B=1-2*A


def words_and_charges(counts):
    words=S5.type_words(counts);M=sum(counts)
    return words,np.array([M-len(u) for u in words],int)


def boundary_matrix(counts,values,g,local_charge):
    words,charge=words_and_charges(counts);N=sum(counts)+1
    pairs=[(i,local_charge-int(q)) for i,q in enumerate(charge) if q<=local_charge]
    lookup={words[i]:j for j,(i,n) in enumerate(pairs)}
    out=np.zeros((len(pairs),len(pairs)))
    for j,(i,n) in enumerate(pairs):
        u=words[i];out[j,j]=sum(values[a] for a in u)
        for a,capacity in enumerate(counts):
            left=capacity-u.count(a)
            if left and n<sum(counts):
                v=u+(a,)
                if v in lookup:
                    k=lookup[v];coefficient=-g/N*sqrt(left)*sqrt(n+1)
                    out[k,j]=out[j,k]=coefficient
    return pairs,out


def taylor_gate(h,time):
    n=len(h);row=float(np.max(np.sum(abs(h),axis=1))) if n else 0
    term=np.eye(n,dtype=complex);result=term.copy()
    termnorm=1.;terrors=0.;error=0.;total=1.
    for k in range(1,33):
        term=(-1j*time/k)*(h@term);result+=term
        termnorm=E.up(termnorm*abs(time)*row/k)
        terrors=E.up(terrors*abs(time)*row/k+E.gamma(512*n+128)*termnorm)
        total=E.up(total+termnorm);error=E.up(error+terrors+E.gamma(16)*total)
    tail=E.up((row*abs(time))**33/factorial(33))
    primitive=E.gamma(1024)*(1+row)
    delta=E.up(sqrt(max(1,n))*error+tail+abs(time)*primitive)
    return np.ascontiguousarray(result),delta


def boundary_gates(counts,values,g,dt):
    result={};maximum=0.
    for q in range(sum(counts)+1):
        pairs,h=boundary_matrix(counts,values,g,q)
        matrix,error=taylor_gate(h,dt/2)
        result[q]=(pairs,matrix);maximum=max(maximum,error)
    return result,maximum


def field_gates(M,time):
    t=np.array([[cos(time),1j*sin(time)],[1j*sin(time),cos(time)]])
    lifts,delta,_,_=S7ROT.fast_matrices(t,M)
    co=E.OLD.trig_interval(E.OLD.point(time),True);si=E.OLD.trig_interval(E.OLD.point(time),False)
    phase=E.up(max(abs(t[0,0].real-co[0]),abs(co[1]-t[0,0].real))+
        max(abs(t[0,1].imag-si[0]),abs(si[1]-t[0,1].imag)))
    delta=E.up(delta+M*phase*(1+phase)**max(0,M-1))
    result={0:([(0,0)],np.ones((1,1),complex))}
    for q in range(1,M+1):
        pairs=[(p,q-p) for p in range(q+1)]
        result[q]=(pairs,lifts[q])
    return result,delta,t,phase


def local_operators(counts,values,g,spectator):
    words,_=words_and_charges(counts);M=sum(counts);N=M+1
    occ=[list(combinations_with_replacement(range(3),k)) for k in range(M+1)]
    states=[(u,o) for u in words if len(u)>=spectator for o in occ[len(u)-spectator]]
    index={v:i for i,v in enumerate(states)};size=len(states)
    f=np.zeros((size,size));gmat=np.zeros((size,size))
    for j,(u,o) in enumerate(states):
        gmat[j,j]=sum(values[a] for a in u)
        for r in set(o):
            for s in [r-1,r+1]:
                if 0<=s<3:
                    oo=list(o);oo.remove(r);oo.append(s);oo=tuple(sorted(oo))
                    f[index[u,oo],j]-=sqrt(o.count(r)*(o.count(s)+1))
        for a,capacity in enumerate(counts):
            left=capacity-u.count(a)
            if left:
                i=index[u+(a,),(0,)+o]
                coefficient=-g/N*sqrt(left)*sqrt(o.count(0)+1)
                gmat[i,j]=gmat[j,i]=coefficient
    return f,gmat,states


def commutator_certificate(counts,values,g):
    rows=[];gg=ff=0.
    for spectator in range(sum(counts)+1):
        f,h,states=local_operators(counts,values,g,spectator)
        c=h@f-f@h
        first=h@c-c@h;second=-f@c+c@f
        n=len(states);scale=1+np.max(np.sum(abs(f),axis=1))+np.max(np.sum(abs(h),axis=1))
        construction=E.up(E.gamma(4096*n+4096)*scale**3)
        first_norm,first_row=E.hermitian_norm(first,construction)
        second_norm,second_row=E.hermitian_norm(second,construction)
        gg=max(gg,first_norm);ff=max(ff,second_norm)
        rows.append(dict(spectator_photons=spectator,dimension=n,GGF=first_row,FFG=second_row))
    constant=E.up(E.up(gg/24)+E.up(ff/12))
    return dict(GGF_upper=gg,FFG_upper=ff,Strang_constant_upper=constant,blocks=rows)


def one_particle_audit(R,dt):
    """Enclose the entire finite one-particle factorization, not a vector."""
    identity=np.eye(R,dtype=complex);term=identity.copy();exact=identity.copy()
    termnorm=1.;terr=error=0.;total=1.
    for k in range(1,33):
        moved=np.zeros_like(term);moved[1:]-=term[:-1];moved[:-1]-=term[1:]
        term=(-1j*dt/k)*moved;exact+=term
        termnorm=E.up(termnorm*2*abs(dt)/k)
        terr=E.up(terr*2*abs(dt)/k+E.gamma(256)*termnorm)
        total=E.up(total+termnorm);error=E.up(error+terr+E.gamma(16)*total)
    reference_error=E.up(sqrt(R)*error+(2*abs(dt))**33/factorial(33))
    product=identity.copy();absolute=np.ones(R);rounding=np.zeros(R)
    phase_total=0.
    for coefficient in [A,B,A]:
        tau=coefficient*dt/2
        t=np.array([[cos(tau),1j*sin(tau)],[1j*sin(tau),cos(tau)]])
        _,_,_,phase=field_gates(1,tau)
        phase_total=E.up(phase_total+2*(R-1)*phase)
        ta=abs(t.real)+abs(t.imag)
        for sequence in [range(R-1),range(R-2,-1,-1)]:
            for i in sequence:
                pair=[i,i+1];old=product[pair].copy();product[pair]=t@old
                a=ta@absolute[pair];r=ta@rounding[pair]
                absolute[pair]=np.nextafter(a*(1+E.gamma(16)),np.inf)
                rounding[pair]=np.nextafter(r*(1+E.gamma(16))+E.gamma(128)*absolute[pair],np.inf)
    difference=exact-product
    row=float(np.max(np.sum(abs(difference.real)+abs(difference.imag),axis=1)))
    col=float(np.max(np.sum(abs(difference.real)+abs(difference.imag),axis=0)))
    delta=E.up(sqrt(E.up(row*col))*(1+E.gamma(16*R+16)))
    subtraction=E.gamma(16)*sqrt(R)*(total+float(np.max(absolute)))
    product_round=E.up(sqrt(R)*float(np.max(rounding)))
    unitary_correction=E.up(phase_total/(1-phase_total))
    full=E.up(delta+reference_error+product_round+subtraction+unitary_correction)
    return dict(one_particle_error_upper=full,computed_difference_bound=delta,
        reference_arithmetic_upper=reference_error,product_arithmetic_upper=product_round,
        subtraction_arithmetic_upper=subtraction,ideal_unitary_phase_upper=unitary_correction,R=R,dt=dt)


def prepare(counts,values,g,physical_modes,T,strang_budget=.20,cut=True):
    M=sum(counts)
    R=min(physical_modes,ceil(12*T/5)+32) if cut else physical_modes
    certificate=commutator_certificate(counts,values,g)
    constant=certificate['Strang_constant_upper']
    n=max(1,ceil(sqrt(constant*T**3/strang_budget)))
    refinements=0
    while True:
        dt=T/n;audit=one_particle_audit(R,dt);delta=audit['one_particle_error_upper']
        field=E.up(1.01*n*M*delta*(1+delta)**max(0,M-1))
        if field<=.02:break
        n=ceil(1.1*n);refinements+=1
        if refinements>60:raise RuntimeError('field preprocessing budget cannot be certified')
    boundary,boundary_error=boundary_gates(counts,values,g,dt)
    fields=[]
    for a in [A,B,A]:
        gates,error,t,phase=field_gates(M,a*dt/2);fields.append((gates,error))
    return dict(R=R,physical_modes=physical_modes,steps=n,dt=dt,
        receiver_cut_error=(E.up(sqrt(M)/2**32) if R<physical_modes else 0),
        commutator_certificate=certificate,field_audit=audit,
        field_step_error=E.up(M*delta*(1+delta)**max(0,M-1)),
        field_total_allocated_upper=field,preprocessing_refinements=refinements,
        boundary=boundary,boundary_error=boundary_error,fields=fields)
