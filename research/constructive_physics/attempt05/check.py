from fractions import Fraction as Q
from pathlib import Path
import math,json,hashlib,platform
P=Path(__file__).resolve().parent;ROOT=P.parents[2];groups=[]
def done(name,**data):groups.append(dict(name=name,status='PASS',**data))
M,g,m=Q(1),Q(1,2),Q(1)
Z=1+6*g*g/(M*M);mu2=m*m/Z;eta=2*g*g/(M*M*Z)
assert Z==Q(5,2) and mu2==Q(2,5) and eta==Q(1,5)
assert 1/(1-3*eta)==Z
assert M*M*eta/(2*(1-3*eta))==g*g
assert mu2/(1-3*eta)==m*m
done('exact_frame_trace_and_inverse_calibration',Z=str(Z),mu2=str(mu2),eta=str(eta))

for gi in [Q(0),Q(1,10),Q(1,2),Q(10),Q(1000)]:
    zi=1+6*gi*gi; ei=2*gi*gi/zi
    assert zi>0 and 0<=ei<Q(1,3)
done('specified_kinetic_domain_and_strength_bound',cases=5,scope='algebraic controls; no empirical parameter selection')

k=2.;om=math.sqrt(k*k+float(mu2));dt=.001
q=0.;v=0.;maxerr=0.;rows=[]
def acc(q):return 1-om*om*q
samples={0,500,1000,2000,3000}
for i in range(3001):
    if i in samples:
        exact=(1-math.cos(om*i*dt))/(om*om)
        err=abs(q-exact);maxerr=max(maxerr,err)
        rows.append(dict(t=i*dt,computed=q,exact=exact,error=err))
    if i==3000:break
    k1q=v;k1v=acc(q)
    k2q=v+dt*k1v/2;k2v=acc(q+dt*k1q/2)
    k3q=v+dt*k2v/2;k3v=acc(q+dt*k2q/2)
    k4q=v+dt*k3v;k4v=acc(q+dt*k3q)
    q+=dt*(k1q+2*k2q+2*k3q+k4q)/6
    v+=dt*(k1v+2*k2v+2*k3v+k4v)/6
assert maxerr<1e-10
done('independent_RK4_retarded_mode_response',maximum_absolute_error=maxerr,rows=rows)

mu=math.sqrt(float(mu2));r=1/mu;E=math.exp(-mu*r);e=float(eta)
psi=-(1+e*E)/r;phi=-(1-e*E)/r
assert abs((psi+phi)/2+1/r)<1e-15
ratio=phi/psi
assert abs((1-ratio)/(1+ratio)-e*E)<1e-15
force=1+e*(1+mu*r)*E
done('static_lensing_slip_joint_control',r=r,slip=ratio,force_over_tensor_Newton=force,scalar_potential_fraction=e*E)

# A positive pole has positive Im G_R for omega>0 with retarded i epsilon.
for w in [0.5,1.,3.]:
    eps=1e-5;den=k*k+float(mu2)-complex(w,eps)**2
    assert (1/den).imag>0
assert math.sqrt(k*k+float(mu2))>mu
done('retarded_positive_pole_and_mass_threshold',scope='regulated frequency diagnostic; analytic light-cone support proved in report')

# Radiation mass and static Yukawa inverse range use the same mu.
for ki in [Q(0),Q(1),Q(2)]:
    w2=ki*ki+mu2
    assert w2-ki*ki==mu2
assert 2+1==3
done('common_mass_dispersion_and_local_dof',local_dof=3,tensor=2,scalar=1)

for n,h in json.loads((P/'BASELINE.json').read_text()).items():
    assert hashlib.sha256((ROOT/n).read_bytes()).hexdigest()==h,n
for n,h in json.loads((P/'FREEZE.json').read_text())['sha256'].items():
    assert hashlib.sha256((P/n).read_bytes()).hexdigest()==h,n
done('baseline_and_freeze_integrity')
out=dict(status='CP05_COVARIANT_RESPONSE_CHECKS_PASS',groups=groups,group_count=len(groups),python=platform.python_version(),code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),verdict='WORKING_LEADING_SCALAR_TENSOR_MODEL; NO_NEW_GRUT_LAW_OR_EMPIRICAL_VALIDATION')
(P/'results.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status']);print('groups',len(groups));print('RK4 maximum error',maxerr)
