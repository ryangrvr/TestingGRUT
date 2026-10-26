"""Gauge-invariant scalar census and explicitly scoped one-loop checks.

Uses exact rational polynomial differentiation for the kappa-only Hessian.
No two-loop matching, pole mass or new physical exclusion is calculated.
"""
from fractions import Fraction as F
from itertools import product
import math
import numpy as np
from common import write,read_previous


def add(a,b):
    c=dict(a)
    for p,v in b.items():c[p]=c.get(p,F(0))+v
    return {p:v for p,v in c.items() if v}


def mul(a,b):
    c={}
    for p,v in a.items():
        for q,w in b.items():
            power=tuple(x+y for x,y in zip(p,q));c[power]=c.get(power,F(0))+v*w
    return {p:v for p,v in c.items() if v}


def diff(a,i):
    out={}
    for p,v in a.items():
        if p[i]:
            q=list(p);q[i]-=1;out[tuple(q)]=v*p[i]
    return out


def evaluate(a,x):return sum(float(v)*math.prod(xx**n for xx,n in zip(x,p)) for p,v in a.items())


def norm_poly(i,j):
    p=[0]*4;q=[0]*4;p[i]=2;q[j]=2
    return {tuple(p):F(1),tuple(q):F(1)}


def census():
    # Exponents (X,X*,Y,Y*,Sigma,Sigma*,phi), fields all dimension one.
    charges=[1,-1,3,-3,6,-6,0]
    actual={n for n in product(range(5),repeat=7) if sum(n)<=4 and sum(q*p for q,p in zip(charges,n))==0}
    expected=set()
    def include(n):expected.add(tuple(n))
    for j in range(5):include([0,0,0,0,0,0,j])
    norms=[]
    for i in range(3):
        p=[0]*7;p[2*i]=p[2*i+1]=1;norms.append(p)
        for f in range(3):q=p.copy();q[-1]=f;include(q)
    for i in range(3):
        for j in range(i,3):include([x+y for x,y in zip(norms[i],norms[j])])
    cubic=[0,0,2,0,0,1,0]
    quartic=[3,0,0,1,0,0,0]
    for p in (cubic,quartic):
        conjugate=[p[1],p[0],p[3],p[2],p[5],p[4],p[6]]
        for base in (p,conjugate):
            for f in range(5-sum(base)):
                q=base.copy();q[-1]=f;include(q)
    even={p for p in actual if p[-1]%2==0}
    return actual,expected,even


def main():
    actual,expected,even=census()
    # V_kappa/kappa = 1/2 [u*r^3-3u*r*i^2+3v*r^2*i-v*i^3],
    # X=(r+i*i)/sqrt2, Y=(u+i*v)/sqrt2.
    potential={(3,0,1,0):F(1,2),(1,2,1,0):F(-3,2),(2,1,0,1):F(3,2),(0,3,0,1):F(-1,2)}
    hess=[[diff(diff(potential,i),j) for j in range(4)] for i in range(4)]
    trace={}
    for i in range(4):
        for j in range(4):trace=add(trace,mul(hess[i][j],hess[j][i]))
    nX=norm_poly(0,1);nY=norm_poly(2,3)
    reference={p:9*v for p,v in mul(nX,nX).items()}
    reference=add(reference,{p:18*v for p,v in mul(nX,nY).items()})
    assert trace==reference
    rng=np.random.default_rng(261026)
    trace_err=[];hess_err=[]
    for _ in range(40):
        q=rng.normal(size=4)
        H=np.array([[evaluate(hess[i][j],q) for j in range(4)] for i in range(4)])
        exact=evaluate(trace,q)
        trace_err.append(abs(np.trace(H@H)-exact)/max(1.,abs(exact)))
        # Independent central difference of the gradient, not the Hessian.
        grad=[diff(potential,i) for i in range(4)];fd=np.zeros((4,4));step=2e-5
        for j in range(4):
            delta=np.eye(4)[j]*step
            fd[:,j]=[(evaluate(g,q+delta)-evaluate(g,q-delta))/(2*step) for g in grad]
        hess_err.append(np.linalg.norm(fd-H)/max(1.,np.linalg.norm(H)))
    previous=read_previous();p=previous['normalized_parameters'];mx=previous['light_sector']['mX_GeV']
    kap=p['kappa'];h=p['h_GeV']
    # The isolated kappa^2 contribution inferred from the MS effective
    # potential's renormalization-scale derivative, in this normalization.
    beta=9*kap*kap/(4*math.pi**2)
    gauge_rows=[];fd_errors=[]
    # A dimensionless differentiation avoids subtraction of ~m_A^4 values.
    # V=3 A^2/(64 pi^2)*(1+y)^2[ln(1+y)-5/6],
    # y=2e^2 |X|^2/A, A=m_A^2, at Q_R=m_A.
    f=lambda y:(1+y)**2*(math.log1p(y)-5/6)
    analytic=-2/3
    for step in (1e-3,1e-4,1e-5):
        fd=(f(step)-f(-step))/(2*step)
        fd_errors.append(abs(fd-analytic))
    for e in (1e-7,1e-6,3e-6,1e-5):
        mass=6*e*h
        correction=-e*e*mass*mass/(16*math.pi**2)
        gauge_rows.append(dict(e=e,mA_GeV=mass,delta_mX_squared_GeV2=correction,
                               magnitude_over_target_mX_squared=abs(correction)/(mx*mx)))
    # A specific covariant representative of the supplied -cH*H^2 mass:
    # V_curv=xi_X R |X|^2. In ideal de Sitter R=12H^2, xi=-cH/12.
    # This defines an added matching choice, not a consequence of closure.
    xi=-p['cH']/12
    curv=[];curv_err=[]
    for label,w in (('de_Sitter',-1.),('matter',0.),('radiation',1/3)):
        Hdot_over_H2=-1.5*(1+w)
        r1=6*(2+Hdot_over_H2);r2=3*(1-3*w)
        curv_err.append(abs(r1-r2))
        curv.append(dict(background=label,w=w,R_over_H2=r1,mass_curvature_squared_over_H2=xi*r1))
    xinit=abs(complex(*previous['prepared_state']['values_H_units'][0]))*p['H_GeV']
    planck=2.435e18
    planck_shift=-2*xi*xinit*xinit/(planck*planck)
    groups={
        'complete_scalar_operator_census':dict(passed=actual==expected and len(actual)==26 and len(even)==19,
            gauge_neutral_monomials_degree_le4=len(actual),phi_even_subset=len(even),unexplained=[] if actual==expected else [list(p) for p in actual-expected]),
        'rational_scalar_Hessian_and_counterterm':dict(passed=trace==reference and max(trace_err)<1e-12 and max(hess_err)<1e-8,
            exact_polynomial_identity=True,samples=len(trace_err),max_trace_error=max(trace_err),max_Hessian_FD_error=max(hess_err)),
        'gauge_potential_derivative':dict(passed=fd_errors[-1]<1e-8 and abs(gauge_rows[1]['magnitude_over_target_mX_squared']/gauge_rows[0]['magnitude_over_target_mX_squared']-1e4)<1e-9,
            normalized_derivative_errors=fd_errors,derivative_exact=-2/3,
            e_power_at_fixed_h=4),
        'covariant_curvature_mass_representative':dict(passed=max(curv_err)<1e-14 and abs(curv[0]['mass_curvature_squared_over_H2']+p['cH'])<1e-14 and curv[2]['mass_curvature_squared_over_H2']==0 and 1+planck_shift>0,
            xi_X=xi,FLRW_trace_identity_error=max(curv_err),backgrounds=curv,
            relative_effective_Planck_mass_squared_shift_at_prepared_X=planck_shift,
            scope='One added covariant representative of cH; ideal backgrounds, spectator-size check only. No inflaton, phase holder, reheating trajectory or Einstein-frame perturbation spectrum solved.')
    }
    result=dict(groups=groups,scalar_monomials=[dict(exponents=list(n),degree=sum(n),phi_even=n[-1]%2==0) for n in sorted(actual)],
                scalar_trace_kappa_squared='Tr[(M_scalar^2)^2]|kappa^2 = 36 kappa^2 |X|^4 + 72 kappa^2 |X|^2 |Y|^2',
                isolated_beta_lambda_XY=beta,beta_convention='mu_R*d lambda_XY/d mu_R; V includes lambda_XY |X|^2|Y|^2; isolated kappa^2 term only.',
                gauge_potential='3 m_A^4(X)/(64 pi^2) [ln(m_A^2(X)/mu_R^2)-5/6], Landau gauge MS-bar, Sigma held fixed.',
                gauge_rows=gauge_rows,scope='Operator/matching audit. Gauge derivative is one scheme/gauge-dependent potential contribution, not a pole mass, full threshold matching, tuning criterion or exclusion. Zero portals are matching choices, not identities preserved by running.')
    write('OPERATOR_RESULTS.json',result)
    print('OPERATORS',sum(bool(g['passed']) for g in groups.values()),'/',len(groups),'monomials',len(actual),'betaXY',beta)
    print('GAUGE_MASS_RATIOS',[r['magnitude_over_target_mX_squared'] for r in gauge_rows])
    if not all(g['passed'] for g in groups.values()):raise SystemExit(1)


if __name__=='__main__':main()
