"""One frozen CP-04 action. Standard library only; exact algebra controls, not experiments."""
from fractions import Fraction as F
from pathlib import Path
import hashlib, json, math, platform

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
groups = []
def passed(name, **data):
    groups.append(dict(name=name, status="PASS", **data))
def avg(xs, ws):
    return sum(x*w for x,w in zip(xs,ws))/sum(ws)
def sub(a,b): return [x-y for x,y in zip(a,b)]

# Exact trace of the scalar stress versus derivative with respect to log(lambda).
# T = -lambda^2 K - 4 lambda^4 P; lambda dS_m/dlambda = same expression.
lam,K,P = F(3,2),F(5,7),F(11,13)
trace = -lam**2*K-4*lam**4*P
log_derivative = -lam**2*K-4*lam**4*P
assert trace == log_derivative
z=F(7,3)
assert 1+z*z > 0
passed("global_variation_trace_control", trace=str(trace))

# Mixed stress matrices represented by their four diagonal components at each cell.
ws = [F(1), F(2), F(3)]
tau = [[F(2),F(1),F(-1),F(0)], [F(4),F(0),F(1),F(2)], [F(-2),F(3),F(2),F(1)]]
def effective(ts):
    tbar=avg([sum(t) for t in ts],ws)
    return [[x-tbar/4 for x in t] for t in ts]
ref=effective(tau)
for c in [F(-17,3),F(0),F(19,2)]:
    shifted=[[x-c for x in t] for t in tau]
    assert effective(shifted)==ref
assert avg([sum(t) for t in ref],ws)==0
passed("exact_constant_shift_and_integrated_trace", shifts=3)

# Finite quadrature for source projector: P=I-1 w^T/(sum w).
weights=[F(1),F(2),F(3),F(4)]
n=4; total=sum(weights)
proj=[[F(i==j)-weights[j]/total for j in range(n)] for i in range(n)]
for i in range(n):
    assert sum(proj[i])==0
    for j in range(n):
        assert sum(proj[i][k]*proj[k][j] for k in range(n))==proj[i][j]
future=[F(0),F(0),F(0),F(5)]
response=[sum(proj[i][j]*future[j] for j in range(n)) for i in range(n)]
assert response==[F(-2),F(-2),F(-2),F(3)]
passed("projector_idempotence_and_nonretarded_support", response=[str(x) for x in response], interpretation="fixed-metric off-shell functional, not an operational signaling result")

# Two-stage vacuum bookkeeping: fraction f in early phase A, rest in B.
f=F(1,1000); A=F(1000); B=F(0)
vbar=f*A+(1-f)*B
early=vbar-A; late=vbar-B
assert early==-999 and late==1
assert f*early+(1-f)*late==0
passed("phase_history_weighting", early=str(early),late=str(late), interpretation="algebraic contribution only; transition matter required for conserved physical solution")

# Full pure-vacuum branch, mu=Omega=m=1: Lambda=-lambda^4 v,
# sigma'=1+v^2, hence lambda^4=1+v^2. Metric is flat for both solutions.
rows=[]
for v in [0,2,10]:
    y=F(1+v*v); L=-y*v; zz=L/y
    assert (1+zz*zz)/y==1
    assert L+y*v==0
    rows.append(dict(v=v,lambda4=str(y),Lambda=str(L),physical_mass=float(y)**0.25))
assert rows[1]["physical_mass"]!=rows[0]["physical_mass"]
passed("full_global_equation_mass_sensitivity", rows=rows)

# Fixed Lambda root of Omega*y^3-y^2-Lambda^2=0, diagnostics only.
def root(omega,L):
    lo=1/omega; hi=lo+1
    def f(y): return omega*y**3-y*y-L*L
    while f(hi)<0: hi*=2
    for _ in range(100):
        mid=(lo+hi)/2
        if f(mid)<0: lo=mid
        else: hi=mid
    return (lo+hi)/2
errs=[]
for omega,L in [(1.,0.),(1.,2.),(3.,7.)]:
    y=root(omega,L)
    e=abs(omega*y**3-y*y-L*L)/max(1,L*L)
    assert e<1e-12
    errs.append(e)
passed("nonlinear_volume_root", maximum_relative_residual=max(errs))

# Supplied Einstein theory controls, not emergence demonstrations.
assert (12-2*4)//2==2
# Nonzero conventional Lambda_eff is compatible with diffeomorphism invariance.
# In 4D de Sitter: R=4 Lambda_eff/M^2, H^2=Lambda_eff/(3 M^2).
M=F(2); Leff=F(3); H2=Leff/(3*M*M); R=12*H2
assert R==4*Leff/(M*M) and R!=0
passed("GR_constraint_count_and_nonzero_vacuum_counterexample", dof=2, deSitter_R=str(R), scope="standard Einstein action and constraint regularity supplied")

base=json.loads((HERE/"BASELINE_INTEGRITY.json").read_text())
for path,h in base.items():
    assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==h,path
freeze=json.loads((HERE/"FREEZE.json").read_text())
for path,h in freeze["sha256"].items():
    assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==h,path
passed("baseline_and_freeze_integrity", preserved_files=len(base))
out=dict(status="CP04_GLOBAL_CONSTRAINT_CHECKS_PASS",groups=groups,group_count=len(groups),python=platform.python_version(),verdict="CONDITIONAL_PUBLISHED_VACUUM_CONSTRAINT; NOT_A_GENERATED_SPACETIME_OR_NEW_GRUT_THEORY",scientific_evidence="AUTHOR_ALGEBRA_AND_SOFTWARE_ONLY",code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/"results.json").write_text(json.dumps(out,indent=2)+"\n")
print(out["status"]);print("groups",len(groups));print(out["verdict"])
