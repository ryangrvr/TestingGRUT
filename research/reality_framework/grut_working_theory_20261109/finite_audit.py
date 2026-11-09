"""Post-run W1 finite-model arithmetic/integrity audit; no density evolution."""
import argparse
import copy
import hashlib
import itertools
import json
import math
from pathlib import Path

FREEZE='6461c320b87bea70c97916b625e2b14511a3a1a5'
AREA=Path('research/reality_framework/grut_working_theory_20261109')
GRID={'theta1':[-.71,0.,.43,math.pi/3,2.27],
      'theta2':[-1.13,0.,.82,math.pi/3,2.91],
      'beta':[-2.3,0.,.37,2*math.pi/3,math.pi],'eta':[0.,.17,.63,1.]}
DURATIONS={'prepare_hadamards':20.,'prepare_bell_cnot':100.,'writer':100.,
    'hold_phase':20.,'first_pulse':40.,'dephasing_slot':60.,'parity_phase':80.,
    'second_pulse':40.,'actuator_first':100.,'actuator_second':100.}


def p(t1,t2,b,e):
    population=(1+math.cos(t1)*math.cos(t2))/2
    interference=-e*math.sin(t1)*math.sin(t2)*math.cos(b)/2
    return population+interference


def tree_compare(left,right,path=''):
    issues=[];maximum=0.
    if type(left) is not type(right):return [path+'/type'],maximum
    if isinstance(right,dict):
        if left.keys()!=right.keys():issues.append(path+'/keys')
        for key in left.keys()&right.keys():
            bad,delta=tree_compare(left[key],right[key],path+'/'+key);issues+=bad;maximum=max(maximum,delta)
    elif isinstance(right,list):
        if len(left)!=len(right):issues.append(path+'/length')
        for index,(a,b) in enumerate(zip(left,right)):
            bad,delta=tree_compare(a,b,path+'/'+str(index));issues+=bad;maximum=max(maximum,delta)
    elif isinstance(right,(int,float)) and not isinstance(right,bool):
        if not math.isfinite(left) or not math.isfinite(right):issues.append(path+'/nonfinite')
        delta=abs(left-right);maximum=max(maximum,delta)
        if not math.isclose(left,right,rel_tol=1e-12,abs_tol=1e-12):issues.append(path)
    elif left!=right:issues.append(path)
    return issues,maximum


def maximum_trig(center,radius,cosine=False):
    f=math.cos if cosine else math.sin
    origin=0. if cosine else math.pi/2
    lo,hi=center-radius,center+radius
    if math.ceil((lo-origin)/math.pi)<=math.floor((hi-origin)/math.pi):return 1.
    return max(abs(f(lo)),abs(f(hi)))


def scalar_bound(theta1,theta2,beta,eta,radii,remainder):
    e1,e2,eb,ee=radii;sines=abs(math.sin(theta1)*math.sin(theta2))
    return min(1.,(e1+e2+min(1.,eta+ee)*sines*eb+sines*abs(math.cos(beta))*ee)/2+remainder)


def inspect(data):
    issues=[];checks={};arithmetic={}
    def check(name,condition):
        checks[name]=bool(condition)
        if not condition:issues.append(name)
    def close(name,actual,expected,tolerance=1e-12):
        error=abs(actual-expected);arithmetic[name]=error;check(name,error<=tolerance)
    check('theoretical_status',data['status']=='theoretical_internal_checks_passed' and data['empirical_status']=='not_run')
    check('fixed_grid',data['grid']==GRID)
    expected=set(itertools.product(*(GRID[key] for key in ('theta1','theta2','beta','eta'))))
    actual=[tuple(row[key] for key in ('theta1','theta2','beta','eta')) for row in data['grid_records']]
    check('all_500_unique_cases',len(actual)==500 and len(set(actual))==500 and set(actual)==expected
        and data['density_matrix_grid_points']==500)
    density_error=0.;analytic_error=0.
    for row in data['grid_records']:
        target=p(*(row[key] for key in ('theta1','theta2','beta','eta')))
        density_error=max(density_error,abs(row['density_agreement']-target))
        analytic_error=max(analytic_error,abs(row['analytic_agreement']-target))
    check('grid_formula',max(density_error,analytic_error)<1e-12)
    close('reported_grid_error',data['maximum_identity_errors']['agreement_formula'],density_error)
    check('all_reported_identity_checks',all(math.isfinite(v) and 0<=v<1e-12 for v in data['maximum_identity_errors'].values()))
    check('finite_durations',data['finite_segment_durations_ns']==DURATIONS and all(t>0 for t in DURATIONS.values()))
    close('total_finite_duration',data['total_duration_ns'],sum(DURATIONS.values()))
    for index,row in enumerate(data['calibration_indistinguishable_countermodels']):
        for program,value in row['program_agreement'].items():
            close(f'countermodel_{index}_{program}',value,p(math.pi/3*int(program[0]),math.pi/3*int(program[1]),row['beta'],row['eta']))
        for program,target in [('00',1.),('10',.75),('01',.75)]:close(f'same_calibration_{index}_{program}',row['program_agreement'][program],target)
    check('countermodels_have_different_holdout',len({round(row['program_agreement']['11'],10) for row in data['calibration_indistinguishable_countermodels']})>=5)
    for index,row in enumerate(data['local_M2_channel_countermodels']):
        close(f'local_M2_agreement_{index}',row['local_M2_agreement'],row['parity_agreement'])
        close(f'parity_Bell_XX_{index}',row['parity_Bell_XX'],1.)
        close(f'local_M2_Bell_XX_{index}',row['local_M2_Bell_XX'],row['eta'])
    close('wrong_local_M1_agreement',data['negative_control_local_M1_dephaser_agreement'],13/16)
    close('wrong_local_M1_error',data['negative_control_error'],3/16)
    row=data['uncertainty_corner_test'];nominal=row['nominal'];radii=row['radii']
    prep,readout=row['preparation_mixture'],row['readout_flip']
    bound=scalar_bound(*nominal,radii,prep+readout)
    close('uncertainty_allowance',row['derived_upper_bound'],bound)
    corner_error=0.
    for signs in itertools.product((-1,1),repeat=4):
        values=[value+sign*radius for value,sign,radius in zip(nominal,signs,radii)]
        # Every declared channel is unital: the mixed preparation contributes1/2.
        mixed=(1-prep)*p(*values)+prep/2
        reported=readout+(1-2*readout)*mixed
        corner_error=max(corner_error,abs(reported-p(*nominal)))
    close('uncertainty_actual_corner_maximum',row['maximum_actual_error'],corner_error)
    check('uncertainty_corners_within_bound',row['corners']==16 and corner_error<=bound)
    detuning=data['finite_pulse_detuning_test']
    detuning_bound=max((40*abs(a)+40*abs(b))/2 for a,b in ((.003,-.002),(-.003,.002),(.001,.004)))
    close('detuning_allowance',detuning['maximum_residual_channel_bound'],detuning_bound)
    check('reported_detuning_within_allowance',0<=detuning['maximum_actual_error']<=detuning_bound)
    cal=data['illustrative_calibration_design_only']
    radius=math.sqrt(2*math.log(2*cal['number_of_expectations']/cal['family_failure_probability'])/cal['shots_per_Pauli_mean'])
    angle=2*math.asin(radius/math.sqrt(2))
    close('calibration_Pauli_radius',cal['Pauli_mean_radius'],radius)
    close('calibration_angular_radius',cal['ideal_centered_angle_radius_radians'],angle)
    close('calibration_attempt_count',cal['structured_parameter_calibration_attempts'],
        cal['shots_per_Pauli_mean']*cal['number_of_expectations'])
    th,beta=math.pi/3,2*math.pi/3
    residual=sum(cal['assumed_separately_certified_prep_readout_channel_budgets'])
    close('calibration_probability_allowance',cal['probability_uncertainty_bound'],
        scalar_bound(th,th,beta,1.,(angle,angle,angle,radius),residual))
    s1,s2,sb=[maximum_trig(center,angle) for center in (th,th,beta)]
    c1,c2,cb=[maximum_trig(center,angle,True) for center in (th,th,beta)]
    contrast_bound=(c1*s2*cb*angle+s1*c2*cb*angle+s1*s2*sb*angle+s1*s2*cb*(2*radius))/2+2*residual
    close('shared_contrast_allowance',cal['shared_parameter_contrast_uncertainty_bound'],contrast_bound)
    contrast_error=0.
    for signs in itertools.product((-1,1),repeat=5):
        actual_angles=[center+sign*angle for center,sign in zip((th,th,beta),signs[:3])]
        ec=min(1.,max(0.,1+signs[3]*radius));ep=min(1.,max(0.,signs[4]*radius))
        contrast_error=max(contrast_error,abs(p(*actual_angles,ec)-p(*actual_angles,ep)-3/16))
    close('shared_contrast_corner_maximum',cal['contrast_parameter_corner_maximum_error'],contrast_error)
    close('contrast_lower',cal['contrast_lower_before_validation_sampling'],3/16-contrast_bound)
    check('shared_contrast_corners_within_bound',contrast_error<=contrast_bound-2*residual)
    check('calibration_is_not_apparatus_certificate','Not measured calibration' in cal['warning'] and 'EX1' in cal['warning'])
    return {'issues':issues,'checks':checks,'arithmetic_maximum_difference':max(arithmetic.values()),
        'maximum_recomputed_grid_error':density_error,'uncertainty_corner_maximum':corner_error,
        'uncertainty_bound':bound,'calibration_probability_bound':cal['probability_uncertainty_bound'],
        'shared_contrast_bound':contrast_bound,'contrast_lower_before_sampling':3/16-contrast_bound}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--root',type=Path,required=True)
    parser.add_argument('--production',type=Path,required=True);parser.add_argument('--known-regression',type=Path,required=True)
    args=parser.parse_args();root=args.root.resolve();out=args.production.resolve()
    data=json.loads((out/'finite_model.json').read_text());known=json.loads(args.known_regression.read_text())
    provenance=json.loads((out/'provenance.json').read_text());cfg=json.loads((root/AREA/'protocol.json').read_text())
    assessment=json.loads((out/'assessment.json').read_text());execution=json.loads((out/'execution.json').read_text())
    check=inspect(data);issues=list(check['issues'])
    regression_issues,regression_delta=tree_compare(data,known)
    issues+=['known_regression'+value for value in regression_issues]
    source_checks={}
    for name,expected in provenance['source_sha256'].items():
        digest=hashlib.sha256((root/name).read_bytes()).hexdigest();source_checks[name]=digest==expected
        if digest!=expected:issues.append('source_hash:'+name)
    if provenance['freeze_commit']!=FREEZE:issues.append('wrong_freeze_commit')
    if provenance['actual_experimental_data'] or provenance['new_large_lattice_evolution']:issues.append('provenance_scope')
    resource=json.loads((out/'finite_model_resources.json').read_text())
    if (resource['returncode']!=0 or resource['error'] or resource['peak_private_bytes']>cfg['resources']['process_memory_bytes']
            or resource['wall_seconds']>cfg['resources']['point_wall_seconds']):issues.append('finite_resource_contract')
    if not execution['execution_complete'] or not assessment['finite_parametric_model']:issues.append('assessment_status')
    forged=copy.deepcopy(data);forged['grid_records'][0]['density_agreement']+=.01
    omitted=copy.deepcopy(data);omitted['grid_records'].pop()
    fake_empirical=copy.deepcopy(data);fake_empirical['empirical_status']='passed'
    negatives={'tampered_probability_rejected':bool(inspect(forged)['issues']),
        'omitted_grid_point_rejected':bool(inspect(omitted)['issues']),
        'false_empirical_claim_rejected':bool(inspect(fake_empirical)['issues'])}
    result={'passed':not issues and all(negatives.values()),'issues':issues,'freeze_commit':FREEZE,
        'source_hashes_checked':len(source_checks),'source_hashes_all_match':all(source_checks.values()),
        'known_regression_maximum_difference':regression_delta,'grid_cases_checked':len(data['grid_records']),
        'arithmetic_audit':check,'negative_controls':negatives,
        'finite_runtime_seconds':resource['wall_seconds'],'finite_peak_private_bytes':resource['peak_private_bytes'],
        'suite_runtime_seconds':execution['seconds'],
        'audited_finite_json_sha256':hashlib.sha256((out/'finite_model.json').read_bytes()).hexdigest(),
        'scope':'Post-run arithmetic, metadata, known-regression and source-integrity audit only. No density evolution rerun. Finite-pulse outputs are theoretical, not measurements; the illustrative calibration is not an EX1 apparatus or budget certificate.'}
    with (out/'finite_audit.json').open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({key:result[key] for key in ('passed','issues','source_hashes_checked','known_regression_maximum_difference','finite_runtime_seconds')}))


if __name__=='__main__':main()
