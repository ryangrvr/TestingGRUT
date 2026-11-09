"""Independent finite-operator checks; synthetic N6 data, no histories/outcomes.

The only imported research module is the archived Schwinger Hamiltonian. The
comparison spin operators, hopping matrices and source families are assembled
independently below. This script does not import the D1 implementation.
"""
import argparse,hashlib,importlib.util,json,platform
from pathlib import Path
import numpy as np
import scipy

def maxnorm(matrix):
    return float(np.max(np.abs(matrix)))

def check(source):
    spec=importlib.util.spec_from_file_location('archived_schwinger_for_d1_check',source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    n=6;a=.35;g=1.3;mass_R=.27;mass_lat=mass_R-g*g*a/8
    x=(g*a)**-2
    archived,electric,states=module.hamiltonian(n,x,mass_lat/g)
    physical=archived.toarray()*g/(2*np.sqrt(x));dimension=len(states);eye=np.eye(dimension)
    z=np.array([[1. if int(state)&(1<<j) else -1. for j in range(n)] for state in states])
    epsilon=(-1.)**np.arange(n);links=epsilon[:-1]
    L=np.cumsum((z+epsilon)/2,axis=1)[:,:-1]
    y=.25+np.cumsum(z/2,axis=1)[:,:-1]
    lookup={int(state):i for i,state in enumerate(states)};hops=[];raising=[]
    for j in range(n-1):
        T=np.zeros((dimension,dimension))
        for col,state in enumerate(states):
            state=int(state)
            if not (state&(1<<j)) and (state&(1<<(j+1))):
                row=lookup[state^((1<<j)|(1<<(j+1)))];T[row,col]=1.
        raising.append(T);hops.append((T+T.T)/(2*a))
    hopping=sum(hops)
    shifted_mass=mass_lat/2*np.sum(1+epsilon*z,axis=1)
    renormalized_mass=mass_R/2*np.sum(1+epsilon*z,axis=1)
    built=hopping+np.diag(g*g*a/2*np.sum(L*L,axis=1)+shifted_mass)
    scalar=-g*g*a*(n-1)/32
    symmetric=hopping+np.diag(g*g*a/2*np.sum(y*y,axis=1)+renormalized_mass)+scalar*eye
    results={
        'neutral_sector':maxnorm(z.sum(axis=1)),
        'independent_gauss_fields':maxnorm(electric-L),
        'staggered_flux_identity':maxnorm(L-links/4-y),
        'independent_full_hamiltonian':maxnorm(physical-built),
        'full_zero_source_operator_identity':maxnorm(physical-symmetric)}
    q=np.array([.03,.01,.02,.04,-.02])
    def calibrated_source(q):
        return physical+np.diag(g*g*a*(L@q)+g*g*a/2*np.dot(q,q))-g*g*a/4*np.dot(links,q)*eye
    source_direct=hopping+np.diag(g*g*a/2*np.sum((y+q)**2,axis=1)+renormalized_mass)+scalar*eye
    results['full_finite_source_operator_identity']=maxnorm(calibrated_source(q)-source_direct)
    step=1/64;derivative_error=0.
    for j in range(n-1):
        direction=np.zeros(n-1);direction[j]=step
        derivative=(calibrated_source(direction)-calibrated_source(-direction))/(2*step*g*g*a)
        derivative_error=max(derivative_error,maxnorm(derivative-np.diag(y[:,j])))
    results['central_source_derivative_vs_y']=derivative_error

    # Fixed synthetic coefficients include a negative entry and a nonzero offset.
    # They are unrelated to B2 geometry weights or production results.
    weights=np.array([1.,2.,-1.,3.,2.])/16;other=np.array([2.,-1.,1.,1.,3.])/16
    c_B=3/32;F=np.diag(y@weights+c_B);G=np.diag(y@other-.07)
    control=np.diag([0.,1.]);eye2=np.eye(2);joint_eye=np.eye(2*dimension)
    original_electric=np.diag(g*g*a/2*np.sum(L*L,axis=1))
    remainder=physical-original_electric
    controlled_error=0.
    for u in (3/8,-5/8):
        sourced=np.kron(eye2,remainder)
        for j in range(n-1):
            link=np.kron(eye2,np.diag(L[:,j]))-u*weights[j]/(g*g*a)*np.kron(control,eye)
            sourced+=g*g*a/2*(link@link)
        apparatus=(-u*c_B-u*u*np.dot(weights,weights)/(2*g*g*a)+u*np.dot(links,weights)/4)
        sourced+=apparatus*np.kron(control,eye)
        target=np.kron(eye2,physical)-u*np.kron(control,F)
        controlled_error=max(controlled_error,maxnorm(sourced-target))
    results['controlled_qubit_calibration_identity']=controlled_error
    commutator_error=0.;weyl_error=0.;kappa=.7;weyl=np.diag(np.exp(1j*kappa*np.diag(F)))
    for j,T in enumerate(raising):
        commutator_error=max(commutator_error,maxnorm(F@T-T@F-weights[j]*T))
        weyl_error=max(weyl_error,maxnorm(weyl.conj().T@T@weyl-np.exp(-1j*kappa*weights[j])*T))
    results['hopping_commutator']=commutator_error
    results['hopping_weyl_conjugation']=weyl_error
    inner=physical@G-G@physical
    expected=-sum(weights[j]*other[j]*hops[j] for j in range(n-1))
    results['mixed_double_commutator']=maxnorm(F@inner-inner@F-expected)
    # A missing staggered calibration must produce a detectable wrong derivative.
    uncalibrated=(physical+np.diag(g*g*a*(L@q)+g*g*a/2*np.dot(q,q)))
    missing_calibration_error=maxnorm(uncalibrated-source_direct)
    wrong_source_rejected=missing_calibration_error>1e-4
    return {'scope':'Independent N6 full-operator algebra and synthetic probe checks; no production probabilities or detector histories.',
        'parameters':{'N':n,'a':a,'g':g,'mass_R':mass_R,'mass_lat':mass_lat,'central_step':step,
          'synthetic_weights':weights.tolist(),'synthetic_offset':c_B,'finite_source':q.tolist()},
        'checks':results,'threshold':1e-11,'negative_control_missing_source_calibration_rejected':wrong_source_rejected,
        'negative_control_operator_error':missing_calibration_error,
        'passed':bool(max(results.values())<1e-11 and wrong_source_rejected),
        'archived_source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
        'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__}

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--schwinger',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();result=check(args.schwinger.resolve())
    with args.output.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(result,indent=2))
    raise SystemExit(0 if result['passed'] else 1)
