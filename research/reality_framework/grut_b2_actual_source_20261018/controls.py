#!/usr/bin/env python3
"""B2-AS exact finite actual-operator resolvent verification, no dynamics."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
JP, JM, A, J = F(1,64), F(1,96), F(13,46080), F(1,3840)
OMEGA = 2-2*A


class Q:
    """Gaussian rationals: every resolvent equality is checked exactly."""
    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        if isinstance(re, Q):
            self.re, self.im = re.re, re.im
        else:
            self.re, self.im = F(re), F(im)

    @staticmethod
    def of(x):
        return x if isinstance(x, Q) else Q(x)

    def __add__(self, other):
        y = Q.of(other)
        return Q(self.re+y.re, self.im+y.im)

    __radd__ = __add__

    def __neg__(self):
        return Q(-self.re, -self.im)

    def __sub__(self, other):
        return self + (-Q.of(other))

    def __rsub__(self, other):
        return Q.of(other) + (-self)

    def __mul__(self, other):
        y = Q.of(other)
        return Q(self.re*y.re-self.im*y.im, self.re*y.im+self.im*y.re)

    __rmul__ = __mul__

    def __truediv__(self, other):
        y = Q.of(other)
        d = y.re*y.re+y.im*y.im
        return Q((self.re*y.re+self.im*y.im)/d,
                 (self.im*y.re-self.re*y.im)/d)

    def __rtruediv__(self, other):
        return Q.of(other)/self

    def __bool__(self):
        return bool(self.re or self.im)

    def __eq__(self, other):
        y = Q.of(other)
        return self.re == y.re and self.im == y.im

    def conjugate(self):
        return Q(self.re, -self.im)

    def json(self):
        return {"re": str(self.re), "im": str(self.im)}


def zeros(n, m=None):
    return [[Q() for _ in range(n if m is None else m)] for _ in range(n)]


def dagger(a):
    return [[a[j][i].conjugate() for j in range(len(a))] for i in range(len(a[0]))]


def mul(a,b):
    c=zeros(len(a),len(b[0]))
    for i in range(len(a)):
        for k in range(len(b)):
            if a[i][k]:
                for j in range(len(b[0])):
                    if b[k][j]:
                        c[i][j] = c[i][j] + a[i][k]*b[k][j]
    return c


def sub(a,b):
    return [[a[i][j]-b[i][j] for j in range(len(a[0]))] for i in range(len(a))]


def shifted(z,h):
    return [[(z if i==j else Q())-h[i][j] for j in range(len(h))] for i in range(len(h))]


def inverse(a):
    n=len(a)
    v=[list(row)+[Q(int(i==j)) for j in range(n)] for i,row in enumerate(a)]
    for k in range(n):
        p=next(i for i in range(k,n) if v[i][k])
        v[k],v[p]=v[p],v[k]
        pivot=v[k][k]
        for j in range(k,2*n):
            v[k][j]=v[k][j]/pivot
        for i in range(n):
            if i!=k and v[i][k]:
                coefficient=v[i][k]
                v[i][k]=Q()
                for j in range(k+1,2*n):
                    v[i][j]=v[i][j]-coefficient*v[k][j]
    return [row[n:] for row in v]


def coord(n,r):
    return (n+r,n) if r>=0 else (n,n-r)


def coordinate_h(N,R):
    states=[coord(n,r) for n in range(N+1) for r in range(-R,R+1)]
    index={ij:k+1 for k,ij in enumerate(states)}
    h=zeros(len(states)+1)
    h[0][0]=Q(OMEGA)
    for (i,j),k in index.items():
        h[k][k]=Q(2+abs(i-j)+(F(1,4) if abs(i-j)==1 else 0)-(A if i==j==0 else 0))
        for target,coupling in [((i+1,j),JP),((i,j+1),JM)]:
            if target in index:
                q=index[target]
                h[k][q]=h[q][k]=Q(-coupling)
    b=index[(0,0)]
    h[0][b]=h[b][0]=Q(-J)
    return h


def blocks(R):
    rs=list(range(-R,R+1)); index={r:k for k,r in enumerate(rs)}
    d=zeros(len(rs));t=zeros(len(rs))
    for r,k in index.items():
        d[k][k]=Q(2+abs(r)+(F(1,4) if abs(r)==1 else 0))
        if r>=0 and r+1 in index:
            a=index[r+1]
            d[a][k]=d[k][a]=Q(-JP)
            t[a][k]=Q(-JM)
        if r<=0 and r-1 in index:
            a=index[r-1]
            d[a][k]=d[k][a]=Q(-JM)
            t[a][k]=Q(-JP)
    return d,t


def block_h(N,R,d,t):
    w=len(d);h=zeros(1+(N+1)*w)
    h[0][0]=Q(OMEGA)
    for n in range(N+1):
        for a in range(w):
            for b in range(w):
                h[1+n*w+a][1+n*w+b]=d[a][b]
                if n<N:
                    h[1+n*w+a][1+(n+1)*w+b]=t[a][b]
                    h[1+(n+1)*w+b][1+n*w+a]=t[a][b].conjugate()
    h[1+R][1+R]=h[1+R][1+R]-A
    h[0][1+R]=h[1+R][0]=Q(-J)
    return h


def main():
    checks=[];cases=[]

    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    R=2
    d,t=blocks(R);ts=dagger(t);w=len(d)
    for N in (0,1,2):
        h=coordinate_h(N,R)
        check(f"coordinate_equals_block_N{N}", h==block_h(N,R,d,t))
        check(f"selfadjoint_N{N}", h==dagger(h))
        bare=[row[1:] for row in h[1:]]
        bare[R][R]=bare[R][R]+A
        for label,z in [("zero",Q(0)),("i",Q(0,1))]:
            full_inv=inverse(shifted(z,h))
            bare_inv=inverse(shifted(z,bare))
            ms=[None]*(N+1)
            ms[N]=inverse(shifted(z,d))
            for n in range(N-1,-1,-1):
                ms[n]=inverse(sub(shifted(z,d),mul(mul(t,ms[n+1]),ts)))
            direct_block=[[bare_inv[a][b] for b in range(w)] for a in range(w)]
            check(f"weyl_root_N{N}_{label}", direct_block==ms[0])
            m0=ms[0][R][R]
            denominator=(z-OMEGA)*(1+A*m0)-J*J*m0
            fuel=(1+A*m0)/denominator
            check(f"fuel_green_N{N}_{label}", fuel==full_inv[0][0])
            y=[[ms[0][a][R]] for a in range(w)]
            stationary_y=[row[:] for row in y]
            wrong_stationary_detected=False
            for n in range(N+1):
                if n:
                    y=mul(mul(ms[n],ts),y)
                    stationary_y=mul(mul(ms[0],ts),stationary_y)
                    wrong_stationary_detected |= y!=stationary_y
                predicted=-J*y[R][0]/denominator
                check(f"site_green_N{N}_{label}_n{n}", predicted==full_inv[1+n*w+R][0])
            if N:
                check(f"finite_tail_not_stationary_N{N}_{label}", wrong_stationary_detected)
            cases.append({"N":N,"R":R,"dimension":len(h),"probe":label,
                          "fuel_resolvent_exact":fuel.json(),"all_exact_defects":"0",
                          "finite_tail_stationarity_used":False})

    # Final probability composition and all-L comparison envelope.
    # These are error-margin checks, NOT source probability evaluations.
    def lower_exp(x):
        total=F(1);term=F(1)
        for k in range(1,65):
            term=term*x/k;total+=term
        return total
    Lmin=2016
    delta_upper=F(12,27613)+F(1500)/lower_exp(F(Lmin,105))
    tail_upper=F(101,100)/lower_exp(F(Lmin,105)-F(3,10))
    sigma_upper=F(6,27613)
    check("comparison_delta_all_L", delta_upper<F(1,1000))
    check("tail_squared_all_L", tail_upper**2<F(1,10**6))
    check("high_weight", sigma_upper**2<F(1,10**6))
    check("high_source_amplitude", 2*sigma_upper<F(1,2000))
    joint=(F(17,50)-F(1,2000)-F(1,1000))**2
    retained=(1-F(2,10**6)-(F(2,3)+F(1,10**6))**2-F(1,1000))**2
    check("correct_squared_joint", joint==F(458329,4000000) and joint>F(1,10))
    check("retention_margin", retained==F(24909781540810162540000081,81000000000000000000000000) and retained>F(1,4))
    check("separated_fraction", F(1,19)**2/(1+F(1,19)**2)==F(1,362))
    check("remaining_capture_budget", F(1,25)-F(1,362)==F(337,9050))
    norm_ratio_squared=(JP*JP+JM*JM)/(J*J)
    useful_eta_ratio=F(6,1000*Lmin)
    check("interblock_norm_ratio", norm_ratio_squared==5200 and norm_ratio_squared>72**2)
    check("late_resolution_ratio", useful_eta_ratio==F(1,336000))
    check("generic_contraction_bound_uninformative_at_resolution", norm_ratio_squared/useful_eta_ratio**2>1)

    result={"status":"all_frozen_controls_pass","checks":checks,"finite_compression_cases":cases,
            "arithmetic":"exact Fraction and Gaussian rational; no floating point probability",
            "scope":"verify actual-source resolvent formula and probability composition only; no source trajectory or spectral classification",
            "margin":{"joint_if_source_hypotheses":str(joint),"retention_if_source_hypotheses":str(retained),
                      "remaining_central_capture_budget":str(F(337,9050)),
                      "delta_upper_exact":str(delta_upper),"tail_upper_exact":str(tail_upper)},
            "method_limit":{"norm_T_over_J_squared":str(norm_ratio_squared),
                            "sufficient_eta_over_J_for_001_generic_smoothing_error_at_L2016":str(useful_eta_ratio),
                            "is_impossibility_for_all_methods":False},
            "actual_operating_gate":"unresolved","actual_violating_point":None}
    output=ROOT/'CONTROLS.json'
    output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({"status":result['status'],"checks":len(checks),"cases":len(cases),
                      "sha256":hashlib.sha256(output.read_bytes()).hexdigest()}))


if __name__=='__main__':
    main()
