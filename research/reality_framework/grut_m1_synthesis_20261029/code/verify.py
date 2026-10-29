"""Same-author consistency and scientific controls for the synthesis packet.

These checks do not establish the cosmological formation history, state
selection, UV completion, isocurvature spectrum, or matter-power prediction.
"""
from pathlib import Path
import ast, json, math
import numpy as np
import phase_holder as ph
ROOT=Path(__file__).resolve().parents[1]
read=lambda name:json.loads((ROOT/name).read_text())
s=read('scale_checks.json');p=read('phase_holder.json');c=read('phase_cubic.json');q=read('holder_evolution.json');h=read('holder_long.json')['cases'][0];r=read('reference_long.json')['cases'][0];sens=read('holder_sensitivity.json');fl=read('phase_fluctuations.json');fs=read('free_stream.json');cr=read('collision_rate.json');cl=read('class_results.json');pa=read('psd_audit.json')
checks=[]
def check(label,condition,detail=''):
 checks.append((label,bool(condition),detail))
def close(a,b,atol=0.,rtol=1e-8):return abs(a-b)<=atol+rtol*abs(b)

check('radiation H at 20 MeV',close(s['H_at_T_GeV'],math.sqrt(math.pi**2*10.75/90)*.02**2/2.435e18))
check('pulse expansion time',close(s['delta_ln_a_constant_H'],s['H_over_mX_at_20MeV']*20000))
check('nonnegative scalar polynomial',1.3**2-4*.75*.95<0 and .75>0)
for K,G,V in [(0.,1.,0.),(3.,5.,2.),(0.,0.,4.),(9.,0.,0.)]:
 check('pressure bounded by energy '+str((K,G,V)), abs(K-G/3-V)<=K+G+V)
check('absolute energy bound',close(s['canonical_energy_gronwall_abs_work_upper_mX'],s['initial_total_field_energy_mX']*math.expm1(3*s['delta_ln_a_constant_H'])))
n=s['density']['ncore_over_mX3'];check('density/spacing/volume',close(s['density']['d_mean_separation_mX_inverse']**3,1/n) and close(s['density']['volume_per_core_mX_inverse_cubed'],1/n))
check('ballistic cylinder',close(s['density']['v_0p7_swept_probability'],math.pi*60**2*.7*20000*n))
check('causal-volume event distinct from target hit',s['density']['poisson_probability_one_center_within_causal_distance']>.04 and s['density']['v_0p7_swept_probability']<2e-7)
a=s['review_arithmetic_correction'];check('heat-retained minus cooled equals checkpoint heat',close(a['hypothetical_if_current_gas_heat_is_retained_percent']-a['hypothetical_if_current_gas_heat_also_fully_redshifts_percent'],a['current_gas_kinetic_percent'],rtol=1e-12))
check('fixed-history calibration',close(s['calibration']['scenarios']['hypothetical_kept_gas_heat']['H_ratio_if_H_three_halves_yield_scaling'],(1+a['hypothetical_if_current_gas_heat_is_retained_percent']/100)**(-2/3)))
check('operator canonical mass dimensions',(2+1+2)-1==4 and (2+1+6)-5==4)
check('operator U(1) charges',-6+3*2==0 and -6+6==0)
st=next(v for v in p['states'] if v['xi']==4e8);z=np.array([complex(*v) for v in st['fields_over_H']]);xi=4e8
V,F,_,_=ph.pot_force(z,xi);rot=np.exp(1j*np.array([1,3,6])*.371);Vr,Fr,_,_=ph.pot_force(z*rot,xi)
check('holder potential gauge invariant',close(V,Vr,rtol=1e-12));check('holder force gauge covariant',np.linalg.norm(Fr-F*rot)<1e-8)
# Independent finite-difference verification of the complex force, at a point
# displaced from the stationary solution so the gradient is nonzero.
zt=z+np.array([.22+.19j,-.13+.25j,.1-.07j]);_,Ft,_,_=ph.pot_force(zt,xi)
err=[]
for j in range(3):
 for component in (1.,1j):
  d=np.zeros(3,dtype=complex);d[j]=component*1e-5
  derivative=(ph.pot_force(zt+d,xi)[0]-ph.pot_force(zt-d,xi)[0])/(2e-5)
  expected=2*(Ft[j].real if component==1 else Ft[j].imag);err.append(abs(derivative-expected)/max(1,abs(expected)))
check('holder force finite-difference gradient',max(err)<2e-6,f'max relative-scaled error {max(err):.3g}')
eig=np.linalg.eigvalsh(ph.hessian6(z,xi));check('local minimum and single gauge direction',abs(eig[0])<1e-4 and eig[1]>9 and st['force_norm_H3']<1e-6)
check('phase shifted by other torques',abs(st['delta_sigma_minus_6x']+math.pi/2)>.5)
check('cubic local root and physical positive modes',all(row['force_norm_over_H3']<1e-5 and row['smallest_eigenvalues_m2_over_H2'][1]>0 for row in c['rows']))
check('mediator hierarchy checked',p['initial_X_over_mediator']>2.5)
old=q['cases'][0];check('parent prepared charge at t=20',abs(old['sampled_NX_at_20']-(-8117.88157476))<2e-5)
check('zero-charge aligned no-holder control',q['cases'][1]['late_mean_NX']==0)
check('long charge ratio',close(h['late_mean_NX']/r['late_mean_NX'],1.0992736203428084,rtol=2e-4))
check('long gauge charge conserved',max(h['max_abs_gauss_charge'],r['max_abs_gauss_charge'])<1e-5)
check('long energy ledger',max(h['energy_ledger_rel_drift'],r['energy_ledger_rel_drift'])<1e-9)
check('Ricci changes sign in switch',h['negative_ricci_sample_fraction']>.1)
check('late charge window control',max(h['relative_late_window_change'],r['relative_late_window_change'])<1e-3)
check('phase-response half-step control',abs(sens['dln_abs_NX_dtheta_sigma_step_0p001']-sens['dln_abs_NX_dtheta_sigma_step_0p0005'])<1e-5)
check('de Sitter massless solution',fl['massless_analytic_relative_error']<1e-6)
check('relative phase removes gauge orbit',fl['goldstone_projection_squared']<1e-12)
check('de Sitter vacuum depth control',abs(fl['rows'][3]['start_depth_1000_vs_2000_fraction'])<1e-6)
x=fs['populations']['96'];y=fs['populations']['64'];check('gas number/charge accounting',close(x['number']/200000,x['charge_fraction']))
check('gas charge and mass fractions distinct',x['mass_fraction_of_cold_plus_held']>x['charge_fraction'])
check('quadrature refinement',abs(fs['checks']['mean_path_4500_vs_9000_fraction'])<1e-6)
check('saved source grid sensitivity < 4%',abs(x['mean_free_streaming_Mpc']/y['mean_free_streaming_Mpc']-1)<.04)
check('finite energy grid floor recorded',.02<x['fraction_floor_bin']<.04 and .08<y['fraction_floor_bin']<.09)
check('equality ballistic scale',3e-5<x['mean_free_streaming_Mpc']<4e-5)
check('power toy bound obeyed',all(abs(v['delta_power_fraction'])<=v['bound_on_abs_fraction_small_perturbation']*1.0001 for v in fs['toy_ballistic_at_equality'].values()))
check('XX collision rate matches parent operator',abs(cr['reference_relative_error'])<1e-12)
rows=cr['rows'];check('collision angular convergence',abs(rows[2]['total_collisions_per_particle_per_Hubble']/rows[3]['total_collisions_per_particle_per_Hubble']-1)<1e-3)
check('collisions are not negligible in inherited model',rows[3]['viscosity_weighted_collisions_per_particle_per_Hubble']>10)
check('control PSD q2 matched',close(fs['psd']['maxwell96']['mean_q2_from_psd'],fs['psd']['maxwell96']['target_node_mean_q2'],rtol=1e-9))
vv=[cl['runs'][k]['percent_difference_from_baseline']['1.0'] for k in ['warm96rk','warm96rkq60','warm96rkq120']]
check('narrow PSD CLASS nonconvergence explicitly detected',max(vv)-min(vv)>.5)
check('PSD derivative quadrature diagnostic detects problem',max(abs(v['adiabatic_derivative_identity']-1) for v in pa['discrete_smoothed_96'][:3])>.05)
check('smooth PSD derivative control',abs(pa['moment_matched_maxwell'][0]['adiabatic_derivative_identity']-1)<.001)
check('CLASS saved spectra numerical and positive',all(np.all(np.isfinite(v)) and np.all(v[:,1]>0) for v in np.load(ROOT/'class_spectra.npz').values()))
for file in ROOT.glob('*.json'):
 def finite(v):
  if isinstance(v,float):return math.isfinite(v)
  if isinstance(v,dict):return all(finite(x) for x in v.values())
  if isinstance(v,list):return all(finite(x) for x in v)
  return True
 check('finite JSON '+file.name,finite(read(file.name)))
for file in ROOT.glob('*.npz'):
 with np.load(file) as ar:check('finite NPZ '+file.name,all(np.all(np.isfinite(ar[k])) for k in ar.files))
for file in (ROOT/'code').glob('*.py'):
 try:ast.parse(file.read_text());ok=True
 except SyntaxError:ok=False
 check('Python syntax '+file.name,ok)
failed=[(a,d) for a,b,d in checks if not b]
for a,b,d in checks:print(('PASS' if b else 'FAIL')+' '+a+(' — '+d if d else ''))
print(f'{len(checks)-len(failed)}/{len(checks)} checks passed. Same-author consistency controls; no independent validation.')
if failed:raise SystemExit(1)
