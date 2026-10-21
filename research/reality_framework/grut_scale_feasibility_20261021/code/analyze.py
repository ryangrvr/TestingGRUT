"""Retrospective feasibility diagnostics. No fit, population simulation or likelihood.

Run from any directory: python3 path/to/code/analyze.py
Only Python's standard library is required. All products stay in this packet.
"""
from pathlib import Path
from fractions import Fraction
import csv
import hashlib
import io
import json
import math
import platform
import sys
import tarfile

PACKET = Path(__file__).resolve().parents[1]
FRAMEWORK = PACKET.parent
G = 0.1
LAMBDA_PHI = 0.01
LAMBDA_X = 0.001
MPL = 2.435323e18  # reduced Planck mass, GeV; external normalization
HBARC = 1.973269804e-14  # GeV cm, conventional conversion
GEV_G = 1.782661922e-24  # grams / GeV, conventional conversion
C_KMS = 299792.458
checks = []


def check(name, ok, evidence):
    checks.append(dict(name=name, passed=bool(ok), evidence=evidence))
    if not ok:
        raise AssertionError((name, evidence))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


profile_path = FRAMEWORK / 'grut_candidate_1_20261018/code/results.json'
cal_path = FRAMEWORK / 'grut_3d_empirical_20261020/calibration_data.json'
archive_path = FRAMEWORK / 'grut_3d_empirical_20261020/raw_diagnostics.tar.gz'
exit_archive_path = FRAMEWORK / 'grut_expanding_charge_20261021/raw_diagnostics.tar.gz'
profile = json.loads(profile_path.read_text())['convergence'][-1]
cal = json.loads(cal_path.read_text())
q = profile['Qbar'] / G**2
e = profile['Ebar'] / profile['Qbar']
rbar = profile['R90_bar']
omega_c = 0.120  # Planck 2018 base-Lambda-CDM retrospective calibration
omega_c_sigma = 0.001
rho_charge_product = omega_c / cal['omega_b'] * cal['mean_mass_per_baryon_GeV_convention']
stat_fraction = math.hypot(omega_c_sigma / omega_c, cal['omega_b_sigma'] / cal['omega_b'])


def props(v, charge=q, energy_ratio=e, radius_bar=rbar):
    mass = energy_ratio * G * v * charge
    radius = radius_bar / (G * v)
    geom = math.pi * (2 * radius * HBARC)**2 / (mass * GEV_G)
    return dict(v_GeV=v, Q=charge, e=energy_ratio, M_GeV=mass,
                R90_cm=radius * HBARC, contact_proxy_cm2_g=geom)


def crossing(v, energy_ratio=e, gst=106.75, gs=106.75):
    """Illustrative constant-g* crossing; NOT a precision SM thermal history."""
    hidden_yield_ratio = rho_charge_product / (energy_ratio * G * v)
    yb = cal['omega_b'] * 1.05375e-5 / (
        cal['mean_mass_per_baryon_GeV_convention'] * 2891.2)
    nc = math.sqrt(LAMBDA_X) * (LAMBDA_PHI / G**2)**1.5 * v**3
    temperature = (nc / (hidden_yield_ratio * yb * (2 * math.pi**2 / 45) * gs))**(1 / 3)
    h_over_m = math.sqrt(math.pi**2 * gst / 90) * temperature**2 / (MPL * G * v)
    return dict(v_GeV=v, e_eff=energy_ratio, r_all_DM=hidden_yield_ratio,
                gstar=gst, gstar_entropy=gs, T_cross_GeV=temperature,
                H_over_mX=h_over_m)


scale_rows = [dict(**props(v), crossing_constant_g=crossing(v))
              for v in [1000.0, 100.0, 1.0, 0.35, 0.1, 0.07]]
base, low = props(1000.0), props(0.1)
check('fixed-shape mass scales as v', abs(low['M_GeV'] / base['M_GeV'] - 1e-4) < 1e-18, low)
check('fixed-shape contact proxy scales as v^-3',
      abs(low['contact_proxy_cm2_g'] / base['contact_proxy_cm2_g'] / 1e12 - 1) < 1e-14, low)
cb, cl = crossing(1000.0), crossing(0.1)
check('all-DM constant-g crossing scales as v^(4/3)',
      abs(cl['T_cross_GeV'] / cb['T_cross_GeV'] / (1e-4)**(4/3) - 1) < 1e-14, cl)
check('all-DM constant-g H/m scales as v^(5/3)',
      abs(cl['H_over_mX'] / cb['H_over_mX'] / (1e-4)**(5/3) - 1) < 1e-14, cl)

# These are contact-PROXY contour labels, not measured SIDM likelihood thresholds.
contours = [dict(contact_proxy_cm2_g=s,
                 v_GeV=1000 * (base['contact_proxy_cm2_g'] / s)**(1/3))
            for s in [0.01, 0.1, 1.0]]
for row in contours:
    got = props(row['v_GeV'])['contact_proxy_cm2_g']
    check('contact-proxy contour inversion %g' % row['contact_proxy_cm2_g'],
          abs(got / row['contact_proxy_cm2_g'] - 1) < 1e-14, row)

# Exact rank/null-direction calculation in the optimistic thin-wall limit.
# Coordinates are (log v, log Q, log r), with e constant.
jacobian = [[Fraction(1), Fraction(0), Fraction(1)],
            [Fraction(-3), Fraction(-1, 3), Fraction(0)]]
null = [Fraction(1), Fraction(-9), Fraction(-1)]
products = [sum(a*b for a,b in zip(row, null)) for row in jacobian]
check('exact calibration null direction', products == [0, 0], [str(x) for x in products])
minor = jacobian[0][0]*jacobian[1][1]-jacobian[0][1]*jacobian[1][0]
check('calibration matrix rank is two', minor != 0, str(minor))
z = 1.2
before = dict(v=0.1, Q=1e6, r=100.0)
after = dict(v=before['v']*z, Q=before['Q']/z**9, r=before['r']/z)
for fn, label in [(lambda x: x['v']*x['r'], 'density'),
                  (lambda x: x['v']**-3*x['Q']**(-1/3), 'contact')]:
    check('finite null curve preserves '+label,
          abs(fn(after)/fn(before)-1) < 1e-14, dict(before=fn(before), after=fn(after)))
degeneracy = dict(approximation='thin-wall, fixed e and uniform density, monodisperse',
                 log_parameter_order=['v','Q','r'], jacobian=[[str(x) for x in row] for row in jacobian],
                 null=[str(x) for x in null], rank=2, before=before, after=after,
                 mass_ratio=z**-8, radius_ratio=z**-4)


def population_moments(charges, probabilities=None):
    """Dimensionless collision-count diagnostics under R=C Q^(1/3), M=A Q.

    Omits common pi*C^2/A. These are NOT drag/heat transport coefficients.
    """
    p = probabilities or [1/len(charges)]*len(charges)
    def moment(power):
        return sum(w*x**power for w,x in zip(p,charges))
    number = 2*(moment(2/3)+moment(1/3)**2)/moment(1)
    mass = (moment(5/3)+2*moment(4/3)*moment(1/3)+moment(1)*moment(2/3))/moment(1)**2
    cw = moment(2)/moment(1)
    mono_cw = 4*cw**(-1/3)
    pair_n = sum(pi*pj*(a**(1/3)+b**(1/3))**2
                 for pi,a in zip(p,charges) for pj,b in zip(p,charges))/moment(1)
    pair_m = sum(pi*pj*a*(a**(1/3)+b**(1/3))**2
                 for pi,a in zip(p,charges) for pj,b in zip(p,charges))/moment(1)**2
    check('number moment vs explicit pair sum', abs(number/pair_n-1) < 3e-13, len(charges))
    check('mass moment vs explicit pair sum', abs(mass/pair_m-1) < 3e-13, len(charges))
    return dict(n=len(charges), Q_number_mean=moment(1), Q_charge_weighted_mean=cw,
                number_collision_coefficient=number, mass_collision_coefficient=mass,
                number_ratio_to_mono_at_charge_weighted_Q=number/mono_cw,
                mass_ratio_to_mono_at_charge_weighted_Q=mass/mono_cw)


target, a, b = 1e6, 2e4, 1e7
p_high = a*(target-a)/((b-a)*(b+a-target))
same_mean = dict(mono=population_moments([target]),
                 mixture=population_moments([a,b],[1-p_high,p_high]),
                 mixture_parameters=dict(Q_low=a,Q_high=b,p_high=p_high))
check('equal charge-weighted mean in different populations',
      abs(same_mean['mixture']['Q_charge_weighted_mean']/target-1) < 1e-14, same_mean)
check('equal charge mean does not fix collision rate',
      same_mean['mixture']['number_ratio_to_mono_at_charge_weighted_Q'] > 2, same_mean)

raw_rows = []
qbind = 17660.0  # quoted a929f94 verifier: classical E=m_X Q crossing, one spherical branch
with tarfile.open(archive_path) as archive:
    for grid in ['coarse', 'medium', 'fine']:
        member = 'raw/production_%s.csv.components.csv' % grid
        rows = list(csv.DictReader(io.StringIO(archive.extractfile(member).read().decode())))
        for threshold in ['0.02','0.05','0.1']:
            allq = [100*float(r['Q']) for r in rows if r['threshold']==threshold]
            selected = [x for x in allq if x>=qbind]
            raw_rows.append(dict(grid=grid,threshold=float(threshold),all_components=len(allq),
                                 above_quoted_binding=len(selected),
                                 selected_charge_fraction_of_identified=sum(selected)/sum(allq),
                                 diagnostic=population_moments(selected)))
check('archived fine 10pct has 39 components and 15 above quoted threshold',
      raw_rows[-1]['all_components']==39 and raw_rows[-1]['above_quoted_binding']==15, raw_rows[-1])

# COM semiclassical parameter for equal masses, k=Mu/2; independent of physical v.
velocities = []
for u_kms in [10,30,200,1000,3000]:
    u = u_kms / C_KMS
    kr = e*q*rbar*u/2
    direct = []
    for v in [1000.0,0.1]:
        m=e*G*v*q
        r=rbar/(G*v)
        direct.append(m*u*r/2)
    check('kR scale cancellation at %g km/s' % u_kms,
          max(abs(x/kr-1) for x in direct)<1e-14, direct)
    velocities.append(dict(relative_speed_kms=u_kms,k_R90=kr,
                           contact_time_R_over_u_in_mX_units=rbar/u))

# Nonzero-gradient pressure control for the dated erratum, plus homogeneous kick.
pressure_controls = []
for K,GG,V in [(0.4,0.7,0.2),(1.0,0.0,0.0)]:
    rho=K+GG+V; p=K-GG/3-V
    corrected=(2*K-2*GG-4*V)/3
    old=(2*K-4*V)/3
    check('pressure identity with G=%g' % GG, abs(p-rho/3-corrected)<1e-15,
          dict(K=K,G=GG,V=V,direct=p-rho/3,corrected=corrected,old=old))
    pressure_controls.append(dict(K=K,G=GG,V=V,direct=p-rho/3,corrected=corrected,old=old))

exit_diagnostics=[]
with tarfile.open(exit_archive_path) as archive:
    for name in sorted(n for n in archive.getnames() if n.endswith('.csv') and
                       not n.startswith(('raw/v_','raw/g_v_'))):
        rows=list(csv.DictReader(io.StringIO(archive.extractfile(name).read().decode())))
        if not rows or 'H_conf' not in rows[0]:
            continue
        h0=float(rows[0]['H_conf'])
        residual=[]; pressure=[]
        for row in rows:
            H=float(row['H_conf']); W=float(row['work_conf'])
            err=abs(H-h0-W)
            residual.append((err/max(abs(h0),abs(H),abs(W)),err/abs(h0)))
            rho=float(row['rho']); p=float(row['p'])
            K=float(row['K_phys']); GG=float(row['G_phys']); U=float(row['U_phys'])
            pressure.append(abs(p-rho/3-(2*K-2*GG-4*U)/3)/max(rho,1e-300))
        peak=max(rows,key=lambda r:min(float(r['genuine_plus']),float(r['genuine_minus'])))
        exit_diagnostics.append(dict(name=name,samples=len(rows),
            max_work_residual_declared=max(x[0] for x in residual),
            max_work_residual_relative_initial=max(x[1] for x in residual),
            max_pressure_identity_relative=max(pressure),
            peak_minority=dict(a=float(peak['a']),plus=float(peak['genuine_plus']),
                               minus=float(peak['genuine_minus']))))
check('raw 3D records satisfy corrected pressure identity',
      len(exit_diagnostics)==12 and max(r['max_pressure_identity_relative'] for r in exit_diagnostics)<1e-14,
      dict(runs=len(exit_diagnostics),max_error=max(r['max_pressure_identity_relative'] for r in exit_diagnostics)))
main101=next(r for r in exit_diagnostics if r['name']=='raw/main_s101.csv')
check('s101 contains a transient minority below formation criterion',
      .07<min(main101['peak_minority']['plus'],main101['peak_minority']['minus'])<.1,main101)

# A fixed-background quantum energy-origin check, not a quantum-gravity theorem.
H=[[Fraction(1),Fraction(2)],[Fraction(2),Fraction(3)]]
rho=[[Fraction(2,3),Fraction(1,6)],[Fraction(1,6),Fraction(1,3)]]
c=Fraction(17,3)
shifted=[[H[i][j]+(c if i==j else 0) for j in range(2)] for i in range(2)]
def commutator(a,b):
    return [[sum(a[i][k]*b[k][j]-b[i][k]*a[k][j] for k in range(2))
             for j in range(2)] for i in range(2)]
check('constant Hamiltonian shift leaves density generator identical',
      commutator(H,rho)==commutator(shifted,rho),str(commutator(H,rho)))
gap2=lambda a:(a[0][0]-a[1][1])**2+4*a[0][1]*a[1][0]
check('constant Hamiltonian shift leaves energy gaps identical',gap2(H)==gap2(shifted),str(gap2(H)))
# For a covariant potential shift C, delta T_mu_nu=-C g_mu_nu, signature -+++.
C=Fraction(5,7)
metric=[-1,1,1,1]
stress=[-C*x for x in metric]
check('vacuum stress shift has positive density and negative pressure',
      stress[0]==C and sum(stress[1:])/3==-C,[str(x) for x in stress])
matter_vac,gravity_counterterm=Fraction(11,5),Fraction(-3,2)
check('compensated renormalization split leaves total vacuum energy fixed',
      matter_vac+gravity_counterterm==(matter_vac+C)+(gravity_counterterm-C),
      str(matter_vac+gravity_counterterm))
vacuum_shift=dict(scope='fixed-background matter dynamics; covariant stress shift separately',
                  energy_shift=str(c),gap_squared=str(gap2(H)),
                  vacuum_density_shift=str(C),stress_diagonal=[str(x) for x in stress],
                  conclusion='closure constraints invariant under H+cI cannot alone fix absolute vacuum energy')

inputs = {str(p.relative_to(FRAMEWORK)):dict(sha256=digest(p),bytes=p.stat().st_size)
          for p in [profile_path,cal_path,archive_path,exit_archive_path]}
result=dict(status='conditional feasibility only; no empirical prediction frozen',
            action_family='fixed dimensionless couplings, m_D0=g_X v; fixed physical M_Pl',
            input_provenance=inputs,
            calibration=dict(omega_c_h2=omega_c,omega_c_sigma=omega_c_sigma,
                             K_GeV=rho_charge_product,calibration_stat_fraction=stat_fraction,
                             baryon_mass_convention_fraction=cal['mean_mass_systematic_fraction'],
                             role='retrospective all-hidden-DM calibration, not abundance prediction'),
            prototype=dict(Q=q,e=e,R90_bar=rbar),scale_rows=scale_rows,
            constant_g_sensitivity_at_v0p1=[crossing(.1,gst=x,gs=x) for x in [10.75,61.75,106.75]],
            contact_proxy_contours=contours,calibration_rank=degeneracy,
            same_mean_counterexample=same_mean,archived_component_diagnostics=raw_rows,
            relative_speed_controls=velocities,pressure_controls=pressure_controls,
            exit_archive_recheck=exit_diagnostics,
            vacuum_shift_controls=vacuum_shift,
            internal_checks=dict(passed=len(checks),failed=0,records=checks))
if '--check' in sys.argv:
    expected=json.loads((PACKET/'RESULTS.json').read_text())
    def compare(a,b,path='root'):
        if isinstance(a,dict):
            if not isinstance(b,dict) or a.keys()!=b.keys():
                raise AssertionError(('keys',path))
            for k in a:
                compare(a[k],b[k],path+'.'+k)
        elif isinstance(a,list):
            if not isinstance(b,list) or len(a)!=len(b):
                raise AssertionError(('list',path))
            for i,(aa,bb) in enumerate(zip(a,b)):
                compare(aa,bb,path+'.'+str(i))
        elif isinstance(a,float):
            if not math.isfinite(a) or not math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-25):
                raise AssertionError(('number',path,a,b))
        elif a!=b:
            raise AssertionError(('value',path,a,b))
    compare(result,expected)
else:
    (PACKET/'RESULTS.json').write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')
    (PACKET/'ENVIRONMENT.json').write_text(json.dumps(dict(python=sys.version,platform=platform.platform(),
        dependencies='standard library only',command='python3 code/analyze.py'),indent=2)+'\n')
print(json.dumps(dict(checks=len(checks),failures=0,K_GeV=rho_charge_product,
                     contours=contours,low_scale=scale_rows[-2],
                     same_mean=same_mean),indent=2))
