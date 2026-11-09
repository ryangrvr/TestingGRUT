"""W1 theoretical finite-pulse Construction 1 model, never hardware data.

Full 16x16 density matrices are evolved by matrix exponentials of declared
finite-duration Hamiltonians. The analytic prediction is used only afterward.
No shot generator, fitting, hardware API, or empirical pass is included.
"""
import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import expm

I2 = np.eye(2,dtype=complex)
X = np.array([[0,1],[1,0]],complex)
Y = np.array([[0,-1j],[1j,0]],complex)
Z = np.diag([1.,-1.]).astype(complex)
I = np.eye(16,dtype=complex)
S,M1,M2,A = range(4)


def local(operator,site):
    factors = [operator if q==site else I2 for q in range(4)]
    return np.kron(np.kron(factors[0],factors[1]),np.kron(factors[2],factors[3]))


XS,ZS = local(X,S),local(Z,S)
X1,Z1 = local(X,M1),local(Z,M1)
X2,Y2,Z2 = local(X,M2),local(Y,M2),local(Z,M2)
XA,ZA = local(X,A),local(Z,A)
P = Z1@Z2
AGREEMENT = (I+ZS@ZA)/2

DURATIONS_NS = {'prepare_hadamards':20.,'prepare_bell_cnot':100.,'writer':100.,
    'hold_phase':20.,'first_pulse':40.,'dephasing_slot':60.,'parity_phase':80.,
    'second_pulse':40.,'actuator_first':100.,'actuator_second':100.}


def apply_hamiltonian(rho,hamiltonian,duration):
    if duration<=0:raise ValueError('Every declared pulse/slot must have positive duration')
    u = expm(-1j*duration*hamiltonian)
    return u@rho@u.conj().T


def cnot_hamiltonian(control,target,duration):
    # exp[-i pi |1><1|_c (I-X_t)/2] = CNOT exactly.
    return math.pi/(2*duration)*((I-local(Z,control))/2)@(I-local(X,target))


def density_model(theta1,theta2,beta,eta,hold_phase=math.pi/2,
                  dephaser='parity',detunings=(0.,0.),prep_mixture=0.,readout_flip=0.):
    """Exact finite-pulse model in ns, angular-frequency Hamiltonians in rad/ns.

    eta is coherence survival, distinct from Construction 1's old hold phase.
    Random finite-duration I/P pulses implement the dephaser; no instantaneous
    operation or infinite Markov decay rate is required, including eta=0.
    Additional detunings are explicit countermodels, not absorbed into eta.
    """
    if not 0<=eta<=1:raise ValueError('eta must lie in [0,1]')
    if not 0<=prep_mixture<=1 or not 0<=readout_flip<=1:raise ValueError('Invalid noise probability')
    rho = np.zeros((16,16),complex);rho[0,0] = 1.
    rho = (1-prep_mixture)*rho+prep_mixture*I/16
    t = DURATIONS_NS
    h = math.pi/(2*t['prepare_hadamards'])*((XS+ZS)+(X1+Z1))/math.sqrt(2)
    rho = apply_hamiltonian(rho,h,t['prepare_hadamards'])
    for name,c,target in (('prepare_bell_cnot',M1,M2),('writer',S,M2)):
        rho = apply_hamiltonian(rho,cnot_hamiltonian(c,target,t[name]),t[name])
    rho = apply_hamiltonian(rho,hold_phase*Z1/(2*t['hold_phase']),t['hold_phase'])
    rho = apply_hamiltonian(rho,theta1*X2/(2*t['first_pulse'])+detunings[0]*Z2/2,t['first_pulse'])
    before_slot = rho.copy()
    operator = {'parity':P,'local_M1':Z1,'local_M2':Z2}[dephaser]
    pulse_branch = apply_hamiltonian(rho,math.pi*operator/(2*t['dephasing_slot']),t['dephasing_slot'])
    idle_branch = apply_hamiltonian(rho,np.zeros((16,16),complex),t['dephasing_slot'])
    rho = (1+eta)*idle_branch/2+(1-eta)*pulse_branch/2
    immediate = rho.copy()
    rho = apply_hamiltonian(rho,beta*P/(2*t['parity_phase']),t['parity_phase'])
    rho = apply_hamiltonian(rho,theta2*X2/(2*t['second_pulse'])+detunings[1]*Z2/2,t['second_pulse'])
    for name,c in (('actuator_first',M1),('actuator_second',M2)):
        rho = apply_hamiltonian(rho,cnot_hamiltonian(c,A,t[name]),t[name])
    probability = float(np.trace(AGREEMENT@rho).real)
    reported = readout_flip+(1-2*readout_flip)*probability
    return {'rho':rho,'before_slot':before_slot,'immediate':immediate,
        'agreement':probability,'reported_agreement':reported}


def prediction(theta1,theta2,beta,eta):
    return (1+math.cos(theta1)*math.cos(theta2)
            -eta*math.sin(theta1)*math.sin(theta2)*math.cos(beta))/2


def independent_flip_prediction(theta1,theta2):
    q1,q2 = math.sin(theta1/2)**2,math.sin(theta2/2)**2
    return 1-q1-q2+2*q1*q2


def expected_final_populations(probability):
    probabilities = np.zeros(16)
    for index in range(16):
        signal,m1,m2,action = [(index>>(3-q))&1 for q in range(4)]
        parity = m1^m2
        if action==parity:
            probabilities[index] = (probability if signal==parity else 1-probability)/4
    return probabilities


def parameter_probability_bound(theta1,theta2,beta,eta,errors,
                                prep_trace_distance=0.,readout_effect_error=0.,channel_half_diamond_sum=0.):
    """Uniform telescoping bound; parameter and remainder errors are disjoint.

    errors are absolute radii (theta1,theta2,beta,eta). All angles are radians.
    Channel budget is ONLY the residual beyond this parameterized model.
    """
    e1,e2,eb,ee = errors
    if min(errors)<0:raise ValueError('Negative uncertainty radius')
    eta_max = min(1.,eta+ee)
    product = abs(math.sin(theta1)*math.sin(theta2))
    structured = .5*(e1+e2+eta_max*product*eb+product*abs(math.cos(beta))*ee)
    remainder = prep_trace_distance+readout_effect_error+channel_half_diamond_sum
    return min(1.,structured+remainder)


def pulse_area_error(omega,duration,omega_error,duration_error):
    return abs(duration)*omega_error+abs(omega)*duration_error+omega_error*duration_error


def maximum_absolute_trig(center,radius,cosine=False):
    lo,hi=center-radius,center+radius
    f=math.cos if cosine else math.sin
    values=[abs(f(lo)),abs(f(hi))]
    shift=0. if cosine else math.pi/2
    first=math.ceil((lo-shift)/math.pi)
    last=math.floor((hi-shift)/math.pi)
    return 1. if first<=last else max(values)


def shared_parameter_contrast_bound(theta1,theta2,beta,eta_c,eta_p,errors,residual_arm_sum=0.):
    """Contrast interval bound when the two arms share calibrated pulse angles.

    errors=(theta1_radius,theta2_radius,beta_radius,eta_c_radius,eta_p_radius).
    Any differential pulse drift belongs in residual_arm_sum, not silently in
    the shared-angle assumption. No independence between calibration errors is
    needed once a joint confidence box has been established.
    """
    e1,e2,eb,ec,ep=errors
    s1,s2,sb=(maximum_absolute_trig(x,e) for x,e in ((theta1,e1),(theta2,e2),(beta,eb)))
    c1,c2,cb=(maximum_absolute_trig(x,e,True) for x,e in ((theta1,e1),(theta2,e2),(beta,eb)))
    difference_max=min(1.,abs(eta_c-eta_p)+ec+ep)
    return min(2.,.5*(difference_max*(c1*s2*cb*e1+s1*c2*cb*e2+s1*s2*sb*eb)
        +s1*s2*cb*(ec+ep))+residual_arm_sum)


def pauli_mean_radius(shots,number_of_expectations=12,family_alpha=.005):
    return math.sqrt(2*math.log(2*number_of_expectations/family_alpha)/shots)


def angle_disk_radius(x_estimate,y_estimate,x_radius,y_radius):
    """Certified angular radius for a unit quadrature vector in an enclosing disk.

    Returns no estimate when the measured uncertainty disk misses the unit
    circle. Branch choice must be fixed by calibration/pulse-area prior.
    """
    magnitude = math.hypot(x_estimate,y_estimate)
    delta = math.hypot(x_radius,y_radius)
    if abs(magnitude-1)>delta+1e-14:raise ValueError('No unit-vector parameter matches calibration confidence region')
    if magnitude==0:return math.pi
    cosine = (1+magnitude*magnitude-delta*delta)/(2*magnitude)
    return math.acos(min(1.,max(-1.,cosine)))


def bell_memory_density():
    v = np.zeros(16,complex);v[0]=v[6]=1/math.sqrt(2)
    return np.outer(v,v.conj())


def theoretical_checks():
    theta1_values = [-.71,0.,.43,math.pi/3,2.27]
    theta2_values = [-1.13,0.,.82,math.pi/3,2.91]
    beta_values = [-2.3,0.,.37,2*math.pi/3,math.pi]
    eta_values = [0.,.17,.63,1.]
    maxima = {'agreement_formula':0.,'final_population_formula':0.,'immediate_population_invariance':0.,
        'trace_error':0.,'hermiticity_error':0.,'negative_eigenvalue':0.,'flip_comparator_at_eta_zero':0.}
    records = []
    for th1,th2,beta,eta in itertools.product(theta1_values,theta2_values,beta_values,eta_values):
        model = density_model(th1,th2,beta,eta)
        probability = prediction(th1,th2,beta,eta)
        maxima['agreement_formula'] = max(maxima['agreement_formula'],abs(model['agreement']-probability))
        maxima['final_population_formula'] = max(maxima['final_population_formula'],float(np.max(abs(np.diag(model['rho']).real-expected_final_populations(probability)))))
        maxima['immediate_population_invariance'] = max(maxima['immediate_population_invariance'],float(np.max(abs(np.diag(model['before_slot'])-np.diag(model['immediate'])))))
        maxima['trace_error'] = max(maxima['trace_error'],float(abs(np.trace(model['rho'])-1)))
        maxima['hermiticity_error'] = max(maxima['hermiticity_error'],float(np.max(abs(model['rho']-model['rho'].conj().T))))
        maxima['negative_eigenvalue'] = max(maxima['negative_eigenvalue'],float(max(0.,-np.linalg.eigvalsh(model['rho'])[0])))
        if eta==0:
            maxima['flip_comparator_at_eta_zero'] = max(maxima['flip_comparator_at_eta_zero'],abs(model['agreement']-independent_flip_prediction(th1,th2)))
        records.append({'theta1':th1,'theta2':th2,'beta':beta,'eta':eta,
            'density_agreement':model['agreement'],'analytic_agreement':probability})
    # Independent phase-sensitive channel calibration identities, using simple states.
    basis0 = np.array([1,0],complex);plus = np.array([1,1],complex)/math.sqrt(2)
    phase_errors = []
    for beta in beta_values:
        for m1 in (0,1):
            v = np.kron(np.kron(basis0,np.eye(2)[:,m1]),np.kron(plus,basis0))
            rho = apply_hamiltonian(np.outer(v,v.conj()),beta*P/(2*DURATIONS_NS['parity_phase']),DURATIONS_NS['parity_phase'])
            phase_errors.extend([abs(np.trace(X2@rho).real-math.cos(beta)),
                abs(np.trace(Y2@rho).real-(-1)**m1*math.sin(beta))])
    maxima['conditional_phase_calibration'] = float(max(phase_errors))
    hold_errors = [abs(density_model(.61,1.27,-.82,.43,hold_phase=phase)['agreement']-prediction(.61,1.27,-.82,.43))
        for phase in (-2.,0.,.48,math.pi/2,2.8)]
    maxima['hold_phase_does_not_change_registered_agreement'] = max(hold_errors)
    th = math.pi/3
    counterexamples = []
    for beta,eta in ((0.,1.),(math.pi/2,1.),(2*math.pi/3,1.),(math.pi,1.),(2*math.pi/3,0.),(2*math.pi/3,.4)):
        values = {}
        for code in ('00','10','01','11'):
            values[code] = density_model(th*int(code[0]),th*int(code[1]),beta,eta)['agreement']
        counterexamples.append({'beta':beta,'eta':eta,'program_agreement':values})
    # Correct parity noise and local M2 noise share this scalar response but differ on Bell coherence.
    bell = bell_memory_density()
    local_counterexamples = []
    for eta in (0.,.37,1.):
        parity_state=(1+eta)*bell/2+(1-eta)*P@bell@P/2
        local_state=(1+eta)*bell/2+(1-eta)*Z2@bell@Z2/2
        correct=density_model(th,th,2*math.pi/3,eta)['agreement']
        local=density_model(th,th,2*math.pi/3,eta,dephaser='local_M2')['agreement']
        local_counterexamples.append({'eta':eta,'parity_agreement':correct,'local_M2_agreement':local,
            'parity_Bell_XX':float(np.trace(X1@X2@parity_state).real),
            'local_M2_Bell_XX':float(np.trace(X1@X2@local_state).real)})
        maxima['agreement_formula']=max(maxima['agreement_formula'],abs(correct-local))
    wrong=density_model(th,th,2*math.pi/3,0.,dephaser='local_M1')['agreement']
    wrong_channel_error=abs(wrong-prediction(th,th,2*math.pi/3,0.))
    # Fixed deterministic perturbation corners; no fabricated shot counts.
    nominal=(.77,1.31,2.02,.68);radii=(.03,.02,.04,.05)
    epsilon_prep,epsilon_readout=.007,.009
    bound=parameter_probability_bound(*nominal,radii,epsilon_prep,epsilon_readout,0.)
    perturbation_error=0.
    for signs in itertools.product((-1,1),repeat=4):
        actual=tuple(value+sign*radius for value,sign,radius in zip(nominal,signs,radii))
        reported=density_model(*actual,prep_mixture=epsilon_prep,readout_flip=epsilon_readout)['reported_agreement']
        perturbation_error=max(perturbation_error,abs(reported-prediction(*nominal)))
    detuning_error=0.;detuning_bound=0.
    for detunings in ((.003,-.002),(-.003,.002),(.001,.004)):
        actual=density_model(*nominal,detunings=detunings)['agreement']
        error=abs(actual-prediction(*nominal))
        bound_one=(DURATIONS_NS['first_pulse']*abs(detunings[0])+DURATIONS_NS['second_pulse']*abs(detunings[1]))/2
        if error>bound_one+1e-12:raise AssertionError('Finite-pulse residual bound violated')
        detuning_error=max(detuning_error,error);detuning_bound=max(detuning_bound,bound_one)
    r=pauli_mean_radius(80000)
    angle_radius=angle_disk_radius(math.cos(th),math.sin(th),r,r)
    example_bound=parameter_probability_bound(th,th,2*math.pi/3,1.,(angle_radius,angle_radius,angle_radius,r),
        prep_trace_distance=.002,readout_effect_error=.002,channel_half_diamond_sum=.004)
    contrast_bound=shared_parameter_contrast_bound(th,th,2*math.pi/3,1.,0.,
        (angle_radius,angle_radius,angle_radius,r,r),residual_arm_sum=.016)
    contrast_maximum_error=0.
    for signs in itertools.product((-1,1),repeat=5):
        actual_angles=[value+sign*angle_radius for value,sign in zip((th,th,2*math.pi/3),signs[:3])]
        actual_eta_c=min(1.,max(0.,1+signs[3]*r))
        actual_eta_p=min(1.,max(0.,signs[4]*r))
        actual_contrast=prediction(*actual_angles,actual_eta_c)-prediction(*actual_angles,actual_eta_p)
        contrast_maximum_error=max(contrast_maximum_error,abs(actual_contrast-3/16))
    if contrast_maximum_error>contrast_bound-.016+1e-12:raise AssertionError('Contrast interval bound violated')
    tests_pass=(max(maxima.values())<1e-12 and wrong_channel_error>.1 and perturbation_error<=bound+1e-12)
    return {'status':'theoretical_internal_checks_passed' if tests_pass else 'theoretical_internal_checks_failed',
        'empirical_status':'not_run','scope':'Prespecified mathematical models and deterministic theoretical checks only. No real or simulated shots; no fit; no hardware realization claimed.',
        'four_work_qubits':['S','M1','M2','A'],'units':'time ns; Hamiltonian angular frequency rad/ns; hbar=1',
        'finite_segment_durations_ns':DURATIONS_NS,'total_duration_ns':sum(DURATIONS_NS.values()),
        'grid':{'theta1':theta1_values,'theta2':theta2_values,'beta':beta_values,'eta':eta_values},
        'density_matrix_grid_points':len(records),'maximum_identity_errors':maxima,'grid_records':records,
        'calibration_indistinguishable_countermodels':counterexamples,
        'local_M2_channel_countermodels':local_counterexamples,
        'negative_control_local_M1_dephaser_agreement':wrong,'negative_control_error':wrong_channel_error,
        'uncertainty_corner_test':{'nominal':nominal,'radii':radii,'corners':16,'preparation_mixture':epsilon_prep,
            'readout_flip':epsilon_readout,'maximum_actual_error':perturbation_error,'derived_upper_bound':bound},
        'finite_pulse_detuning_test':{'maximum_actual_error':detuning_error,'maximum_residual_channel_bound':detuning_bound},
        'illustrative_calibration_design_only':{'shots_per_Pauli_mean':80000,'number_of_expectations':12,
            'structured_parameter_calibration_attempts':960000,
            'family_failure_probability':.005,'Pauli_mean_radius':r,'ideal_centered_angle_radius_radians':angle_radius,
            'assumed_separately_certified_prep_readout_channel_budgets':[.002,.002,.004],
            'probability_uncertainty_bound':example_bound,
            'shared_parameter_contrast_uncertainty_bound':contrast_bound,
            'contrast_parameter_corner_maximum_error':contrast_maximum_error,
            'contrast_lower_before_validation_sampling':3/16-contrast_bound,
            'warning':'Not measured calibration. Residual channel certificates are assumptions, not consequences of these few Pauli probes; this count does not certify the full EX1 apparatus or fit within a separately justified qualification budget.'}}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
    result=theoretical_checks()
    with args.output.open('x',encoding='utf-8') as stream:json.dump(result,stream,indent=2,allow_nan=False);stream.write('\n')
    print(json.dumps({key:result[key] for key in ('status','empirical_status','density_matrix_grid_points','maximum_identity_errors')}))
    raise SystemExit(0 if result['status']=='theoretical_internal_checks_passed' else 1)
