"""Finite gauge-projected Gaussian boundary kernels and an internal clock.

This is a charged-mode state-kernel test, not a quantization of the full
gravitational field or the physical cosmological condensate.
"""
from itertools import product
from math import factorial,sqrt,log
import numpy as np
from common import write


def entropy(prob):return float(-sum(p*np.log(p) for p in prob if p>1e-15))


def model(cuts):
    basis=[n for n in product(*(range(c+1) for c in cuts)) if n[0]+3*n[1]-6*n[2]==0]
    lookup={n:i for i,n in enumerate(basis)};d=len(basis)
    A=np.zeros((d,d));B=np.zeros((d,d))
    for j,(nx,ny,ns) in enumerate(basis):
        target=(nx-3,ny+1,ns)
        if nx>=3 and target in lookup:A[lookup[target],j]=sqrt(nx*(nx-1)*(nx-2)*(ny+1))
        target=(nx,ny-2,ns-1)
        if ny>=2 and ns>=1 and target in lookup:B[lookup[target],j]=sqrt(ns*ny*(ny-1))
    numbers=np.array(basis);light=numbers[:,0]+3*numbers[:,1]
    H=np.diag(numbers@np.array([1.,1.7,2.3]))+.025*(A+A.T)+.04*(B+B.T)
    # Direct check of every retained transition, also before Q=0 restriction.
    charge_errors=[]
    for n in product(*(range(c+1) for c in cuts)):
        for delta in ((-3,1,0),(0,-2,-1)):
            target=tuple(a+b for a,b in zip(n,delta))
            if all(0<=v<=c for v,c in zip(target,cuts)):
                charge_errors.append(abs(np.dot([1,3,-6],np.array(target)-n)))
    return basis,numbers,light,H,A,B,max(charge_errors)


def coherent(basis,cuts,delta):
    means=np.array([6.,2.,2.]);alpha=np.sqrt(means)*np.exp(1j*np.array([0.,0.,delta]))
    vec=np.array([np.prod([a**n/sqrt(factorial(n)) for a,n in zip(alpha,ns)]) for ns in basis])
    norm=np.vdot(vec,vec).real
    full_norm=np.prod([sum(mean**n/factorial(n) for n in range(cut+1)) for mean,cut in zip(means,cuts)])
    return vec/sqrt(norm),norm/full_norm


def run(cuts,clock=True):
    basis,numbers,light,H,A,B,gauge_error=model(cuts);d=len(basis)
    E,U=np.linalg.eigh(H)
    step=(U*np.exp(-.1j*E))@U.conj().T
    N=5;Ws=[np.eye(d,dtype=complex)]
    for t in range(1,N):Ws.append(step@Ws[-1])
    V=np.concatenate(Ws,axis=0)/sqrt(N);P=V@V.conj().T
    C=np.zeros((N*d,N*d),complex)
    for t in range(N-1):
        sl=slice(t*d,(t+1)*d);sn=slice((t+1)*d,(t+2)*d)
        C[sl,sl]+=.5*np.eye(d);C[sn,sn]+=.5*np.eye(d)
        C[sn,sl]-=.5*step;C[sl,sn]-=.5*step.conj().T
    rows=[]
    for delta in (0.,np.pi/2,-np.pi/2,np.pi):
        vec,prob=coherent(basis,cuts,delta);weights=abs(vec)**2
        mean=weights@numbers;energy=float(np.real(np.vdot(vec,H@vec)))
        torque=float(np.real(1j*np.vdot(vec,(H@np.diag(light)-np.diag(light)@H)@vec)))
        mode_ent=[]
        for i,cut in enumerate(cuts):
            # In the exact neutral sector, different occupation numbers of
            # one charged mode require distinct complement charge sectors.
            mode_ent.append(entropy(np.bincount(numbers[:,i],weights=weights,minlength=cut+1)))
        history=V@vec;rho=np.outer(history,history.conj())
        boundary_hidden=np.outer(vec,vec.conj());SM=np.diag([.2,.8])
        joint=np.kron(boundary_hidden,SM)
        hidden_SM_MI=entropy(np.linalg.eigvalsh(boundary_hidden))+entropy([.2,.8])-entropy(np.linalg.eigvalsh(joint))
        effects=[]
        # Physical event: X occupation >=6 at internal clock reading2.
        t=2;sl=slice(t*d,(t+1)*d)
        event=np.diag((numbers[:,0]>=6).astype(float));Aclock=np.zeros_like(C);Aclock[sl,sl]=event
        effect=N*P@Aclock@P
        born=float(np.real(np.vdot(history,effect@history)))
        direct=float(np.real(np.vdot(Ws[t]@vec,event@(Ws[t]@vec))))
        rows.append(dict(delta=delta,neutral_projection_probability=prob,mean_occupations=mean.tolist(),
                         mean_gauge_charge=float(weights@(numbers@np.array([1,3,-6]))),
                         energy=energy,light_charge=float(weights@light),light_charge_torque=torque,
                         mode_entropies_nats=mode_ent,constraint_residual=float(np.linalg.norm(C@history)),
                         initial_hidden_SM_mutual_information=hidden_SM_MI,
                         state_norm=float(np.vdot(history,history).real),state_purity=float(np.real(np.trace(rho@rho))),
                         physical_effect_commutator=float(np.linalg.norm(C@effect-effect@C)),
                         clock_event_probability=born,direct_clock_event_probability=direct))
    plus,minus=rows[1],rows[2]
    # A product with its charge-conjugate copy has zero initial signed light
    # charge; opposite initial phases yield a nonzero signed torque. No full
    # continuum Fock-mode correspondence is assumed in this elementary test.
    doubled=dict(initial_signed_light_charge=plus['light_charge']-minus['light_charge'],
                 signed_torque=plus['light_charge_torque']-minus['light_charge_torque'],
                 reversed_phase_signed_torque=minus['light_charge_torque']-plus['light_charge_torque'],
                 energy=plus['energy']+minus['energy'])
    return dict(cutoffs=cuts,neutral_dimension=d,clock_readings=N,rows=rows,doubled_conjugate_example=doubled,
                checks=dict(gauge_transition_error=gauge_error,isometry=float(np.linalg.norm(V.conj().T@V-np.eye(d))),
                            constraint_isometry=float(np.linalg.norm(C@V)),
                            CP_energy_difference=abs(plus['energy']-minus['energy']),
                            CP_entropy_difference=max(abs(a-b) for a,b in zip(plus['mode_entropies_nats'],minus['mode_entropies_nats'])),
                            CP_occupation_difference=max(abs(a-b) for a,b in zip(plus['mean_occupations'],minus['mean_occupations'])),
                            CP_torque_sum=abs(plus['light_charge_torque']+minus['light_charge_torque'])))


def main():
    rows=[run((18,8,6)),run((24,12,9))]
    changes={name:abs(rows[0]['rows'][1][name]/rows[1]['rows'][1][name]-1) for name in ('energy','light_charge','light_charge_torque')}
    # A real, diagonal Lorentzian boundary phase changes configuration
    # phases but preserves the initial probabilities in that basis.
    psi=np.sqrt([.2,.8]);after=np.exp(1j*np.array([0.,10.]))*psi
    phase_error=float(np.max(abs(abs(after)**2-abs(psi)**2)))
    groups={
        'gauge_projected_boundary_state':dict(passed=all(r['checks']['gauge_transition_error']==0 and max(abs(s['mean_gauge_charge']) for s in r['rows'])<1e-12 for r in rows)),
        'clock_constraint_and_physical_probability':dict(passed=all(r['checks']['constraint_isometry']<1e-11 and r['checks']['isometry']<1e-11 and max(abs(s['clock_event_probability']-s['direct_clock_event_probability']) for s in r['rows'])<1e-11 for r in rows)),
        'phase_not_selected_by_energy_entropy_or_charge':dict(passed=all(r['checks']['CP_energy_difference']<1e-11 and r['checks']['CP_entropy_difference']<1e-11 and r['checks']['CP_torque_sum']<1e-11 and max(abs(s['initial_hidden_SM_mutual_information']) for s in r['rows'])<1e-11 and abs(r['doubled_conjugate_example']['signed_torque'])>.01 for r in rows)),
        'cutoff_and_real_boundary_phase_controls':dict(passed=max(changes.values())<.02 and phase_error<1e-14,cutoff_relative_changes=changes,real_boundary_phase_probability_error=phase_error)
    }
    write('BOUNDARY_RESULTS.json',dict(models=rows,groups=groups,scope=__doc__,actual_cosmological_state_derived=False))
    print('BOUNDARY',[(r['neutral_dimension'],r['doubled_conjugate_example']) for r in rows],changes,flush=True)
    if not all(g['passed'] for g in groups.values()):raise SystemExit(1)


if __name__=='__main__':main()
