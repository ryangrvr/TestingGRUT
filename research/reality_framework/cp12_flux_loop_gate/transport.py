"""Density response and loop retention in the SAME diagnostic flux branch.

This is explicitly conditioned on a loop filter, with its energy/entropy
backaction audited elsewhere. It is not an autonomous preparation or cone.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply
from loops import P,UP,model,fixed_triangle
from references import energy_basis,inclusion
from persistence import compressed_filter,positive_interval
from evolution import Evolution


def product_interval(a,b):
    values=[x*y for x in a for y in b];return [min(values),max(values)]


def quotient_interval(a,b):
    assert b[0]>0
    values=[x/y for x in a for y in b];return [min(values),max(values)]


def run(snapshot_dir):
    start=time.monotonic();full,energies,V,weights,beta=energy_basis()
    source=model(4,4,256,source=True);J,checks=inclusion(full,source)
    sub,H,B,states,sizes,psi=source;D,cycle=fixed_triangle(4,states)
    source_number=np.array([s[1]&1 for s in states],float)
    density=np.array([sum(s[1]>>i&1 for i in [1,2,3])/3 for s in states])
    trajectory=json.loads((UP/'N4_Q4_B256_T20.json').read_text())
    v=np.load(Path(snapshot_dir)/'N4_Q4_B256_T10.npz')['wave']
    ew=min(trajectory['rows'],key=lambda r:abs(r['T']-10))['state_error_upper']
    inputs=[('neutral_sea_wait_10',(J@v)[:,None],ew)]
    for label,w in weights.items():inputs.append((label,J@(V*np.sqrt(w)[None,:]),0.))
    rows=[];lag=.5;E=Evolution(H,B,2.,4)
    for label,X,ewait in inputs:
        Y,right,res=compressed_filter(X,D['occupied_flux'])
        p=float(np.linalg.norm(Y)**2);Y/=np.sqrt(p)
        conditional_input_error=2*ewait/np.sqrt(p)
        pre_energy=float(np.vdot(Y,H@Y).real)
        paths=[]
        for sign in [1,-1]:
            kick=np.exp(-1j*sign*np.pi/2*source_number)
            Z=kick[:,None]*Y
            work=float(np.vdot(Z,H@Z).real-pre_energy)
            energy_after_kick=pre_energy+work
            eps=conditional_input_error
            for _ in range(5):Z,_,de=E(Z,.1);eps+=de
            pfinal=np.sum(abs(Z)**2,axis=1)
            val=float(pfinal@density)
            retain=float(pfinal@D['protected']);flux=float(pfinal@D['occupied_flux'])
            circulation=float(pfinal@D['electric_circulation'])
            circulating_density=float(pfinal@(density*D['electric_circulation']))
            paths.append(dict(source_phase_sign=sign,mean_triangle_density=val,
                mean_triangle_density_interval=positive_interval(val,eps),
                same_triangle_protection_probability=retain,
                same_triangle_protection_interval=positive_interval(retain,eps),
                same_triangle_occupied_flux_probability=flux,
                same_triangle_electric_circulation_probability=circulation,
                same_triangle_electric_circulation_interval=positive_interval(circulation,eps),
                circulation_weighted_density=circulating_density,
                circulation_weighted_density_interval=positive_interval(circulating_density,eps),
                full_path_error_upper=eps,source_kick_work=work,
                source_kick_work_error_upper=2*np.sqrt(32.)*(2*conditional_input_error+conditional_input_error**2),
                projected_mean_energy_after_kick=energy_after_kick))
        plus,minus=paths
        chi=(plus['mean_triangle_density']-minus['mean_triangle_density'])/2
        interval=[(plus['mean_triangle_density_interval'][0]-minus['mean_triangle_density_interval'][1])/2,
                  (plus['mean_triangle_density_interval'][1]-minus['mean_triangle_density_interval'][0])/2]
        # Quotient derivative at eta=0 separates within-circulation density
        # response from merely changing the probability of the circulation bin.
        Z=Y.copy();e0=conditional_input_error
        for _ in range(5):Z,_,de=E(Z,.1);e0+=de
        pzero=np.sum(abs(Z)**2,axis=1)
        q0=float(pzero@D['electric_circulation']);f0=float(pzero@(density*D['electric_circulation']))
        qI=positive_interval(q0,e0);fI=positive_interval(f0,e0)
        qprime=(plus['same_triangle_electric_circulation_probability']-minus['same_triangle_electric_circulation_probability'])/2
        fprime=(plus['circulation_weighted_density']-minus['circulation_weighted_density'])/2
        qpI=[(plus['same_triangle_electric_circulation_interval'][0]-minus['same_triangle_electric_circulation_interval'][1])/2,
             (plus['same_triangle_electric_circulation_interval'][1]-minus['same_triangle_electric_circulation_interval'][0])/2]
        fpI=[(plus['circulation_weighted_density_interval'][0]-minus['circulation_weighted_density_interval'][1])/2,
             (plus['circulation_weighted_density_interval'][1]-minus['circulation_weighted_density_interval'][0])/2]
        aI=product_interval(fpI,qI);bI=product_interval(fI,qpI)
        numerator=[aI[0]-bI[1],aI[1]-bI[0]]
        conditional_circulation_response=(fprime*q0-f0*qprime)/(q0*q0)
        conditional_circulation_interval=quotient_interval(numerator,[qI[0]**2,qI[1]**2])
        # Independent projected commutator-vector calculation checks the sign
        # and exact +/-pi/2 identity, without a numerical eta derivative.
        width=Y.shape[1];initial=np.column_stack([Y,source_number[:,None]*Y])
        A=H-2*eye(H.shape[0],format='csr')
        linear=expm_multiply(-1j*lag*A,initial,traceA=-1j*lag*float(A.diagonal().sum()))
        a,b=linear[:,:width],linear[:,width:]
        analytic=float(2*np.vdot(a,density[:,None]*b).imag)
        assert abs(analytic-chi)<1e-10
        # The vacuum-flux diagnostic is exactly inert under a density phase.
        vacuum=D['vacuum_flux'][:,None]*X
        assert np.linalg.norm(source_number[:,None]*vacuum)==0
        rows.append(dict(preparation=label,initial_flux_triangle_probability=p,
            initial_conditional_state_error_upper=conditional_input_error,
            conditional_mean_energy_before_source=pre_energy,
            lag=lag,source_vertex=0,density_detector_vertices=[1,2,3],
            density_susceptibility=chi,density_susceptibility_interval=interval,
            commutator_vector_susceptibility=analytic,
            circulation_population_susceptibility=qprime,
            circulation_weighted_density_susceptibility=fprime,
            zero_source_circulation_probability=q0,
            density_susceptibility_conditioned_on_final_circulation=conditional_circulation_response,
            density_susceptibility_conditioned_on_final_circulation_interval=conditional_circulation_interval,
            paths=paths,vacuum_flux_density_phase_response_exact=0.,
            no_hopping_after_same_preparation_response_exact=0.))
        print(json.dumps(dict(stage='conditional_density',preparation=label,chi=chi,
            interval=interval,retention=[a['same_triangle_protection_probability'] for a in paths])),flush=True)
    out=dict(N=4,Q=4,budget=256,named_triangle=[1,2,3],rows=rows,seconds=time.monotonic()-start,
        response_identity='For a projector source n0, the derivative at eta=0 of <D> under exp(-i eta n0) equals [<D>_(+pi/2)-<D>_(-pi/2)]/2 exactly; no finite-difference bias.',
        conditional_circulation_response_identity='For f=<C D>, g=<C>, the derivative of final-circulation-conditioned density at eta=0 is (fprime*g-f*gprime)/g^2, using exact phase-source derivatives and zero-source expectations.',
        scope='Same normalized flux-triangle diagnostic branch for both the source response and final geometry detector. The density detector does not use graph distance or multiply density by a geometry bin, so a nonzero value is a genuine matter-density response. Filtering and phase kicks are external diagnostics with explicit resource costs, not autonomous selection. Vacuum branches have n0 rho=rho n0=0 and exactly zero phase response. With hopping removed after the same prepared state, vertex density commutes with all remaining terms and the response is exactly zero. No distance hierarchy or causal cone is measured. Ordinary floating-point evaluation of full-law bounds.')
    (P/'CONDITIONAL_DENSITY.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--snapshot-dir',required=True)
    run(ap.parse_args().snapshot_dir)
