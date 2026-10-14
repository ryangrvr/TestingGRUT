"""Same named triangle: quantum correlation, conversion, and true loop loss.

The main state remains unmeasured. Sequential probabilities explicitly refer
to a diagnostic projective filter; its energy backaction is reported. The
separate two-time correlators are not called classical survival probabilities.
"""
from pathlib import Path
import argparse,json,time
import numpy as np
from scipy.linalg import svd
from scipy.sparse import eye
from scipy.sparse.linalg import expm_multiply
from loops import P,UP,model,fixed_triangle
from references import energy_basis,inclusion
from evolution import Evolution


def h2(p):
    return float(-sum(x*np.log2(x) for x in [p,1-p] if 0<x<1))


def positive_interval(p,e):
    return [max(0.,np.sqrt(max(0.,p))-e)**2,min(1.,(np.sqrt(max(0.,p))+e)**2)]


def ratio_interval(numerator,denominator):
    a,b=numerator;c,d=denominator
    return [max(0.,a/d) if d else 0.,min(1.,b/c) if c>0 else 1.]


def compressed_filter(X,mask):
    active=np.flatnonzero(mask);Y=X[active]
    U,s,Vh=svd(Y,full_matrices=False,check_finite=False)
    # Retain all singular values. There is no rank-threshold approximation.
    Z=np.zeros((X.shape[0],len(s)),complex);Z[active]=U*s[None,:]
    error=float(np.linalg.norm(Z@Vh-mask[:,None]*X))
    assert error<1e-12
    return Z,Vh.conj().T,error


def case(label,X,source,detectors,wait_error,trajectory,wait,energies=None,output_steps=None,
         initial_channels=None,additional_detectors=None):
    src,H,B,states,sizes,psi=source
    if output_steps is None:output_steps=[0,1,2,5,10,20]
    maxsteps=max(output_steps);step=.1;A=H-2*eye(H.shape[0],format='csr')
    E=Evolution(H,B,2.,4);hopping_norm=np.sqrt(32.)
    # A stationary retained-H reference has an exact simple full-law bound:
    # ||U_H(t)X-U_P(t)X||_F <= t ||BX||_F, since its columns are eigenvectors.
    ref_leak=float(np.linalg.norm(B@X)) if energies is not None else None
    out=[]
    if initial_channels is None:initial_channels=['hole','occupied_flux','vacuum_flux']
    if additional_detectors is None:additional_detectors=[]
    for initial in initial_channels:
        event=detectors[initial];Y,right,compression_error=compressed_filter(X,event)
        p0=float(np.linalg.norm(Y)**2)
        p0interval=positive_interval(p0,wait_error)
        projected_initial=event[:,None]*X
        other=X-projected_initial
        work=float(-2*np.vdot(projected_initial,H@other).real)
        work_error=2*hopping_norm*(2*wait_error+wait_error**2)
        reference_entropy_gain_upper=h2(p0)
        if energies is not None:
            norm_forward_right=float(np.sqrt(np.sum(np.sum(abs(X)**2,axis=0)*np.sum(abs(right)**2,axis=1))))
        rows=[];aux_error=0.
        for j in range(maxsteps+1):
            lag=j*step
            if j in output_steps:
                if energies is None:
                    future=expm_multiply(-1j*lag*A,X,traceA=-1j*lag*float(A.diagonal().sum())) if lag else X.copy()
                    e_forward=min(trajectory['rows'],key=lambda a:abs(a['T']-(wait+lag)))['state_error_upper']
                    forward_right=future@right;norm_forward_right=float(np.linalg.norm(forward_right))
                else:
                    phases=np.exp(-1j*(energies-2)*lag)
                    future=X*phases[None,:]
                    e_forward=lag*ref_leak
                    forward_right=None
                e_filtered=wait_error+aux_error
                values={};correlations={}
                for final in ['hole','occupied_flux','protected',initial]+additional_detectors:
                    if final in values:continue
                    mask=detectors[final]
                    score=float(np.linalg.norm(mask[:,None]*Y)**2)
                    jointinterval=positive_interval(score,e_filtered)
                    cond=score/p0
                    values[final]=dict(sequential_joint_probability=score,
                        sequential_joint_interval=jointinterval,
                        conditional_probability=cond,
                        conditional_interval=ratio_interval(jointinterval,p0interval))
                    if energies is None:
                        corr=np.vdot(forward_right,mask[:,None]*Y)
                    else:
                        active=np.flatnonzero(mask)
                        cross=X[active].conj().T@Y[active]
                        corr=np.sum(np.exp(1j*(energies-2)*lag)[:,None]*right.conj()*cross)
                    corr_error=e_forward*np.linalg.norm(Y)+e_filtered*norm_forward_right+e_forward*e_filtered
                    pfuture=float(np.linalg.norm(mask[:,None]*future)**2)
                    pfutureinterval=positive_interval(pfuture,e_forward)
                    correlations[final]=dict(real=float(corr.real),imaginary=float(corr.imag),
                        complex_absolute_error_upper=float(corr_error),
                        unmeasured_final_probability=pfuture,
                        unmeasured_final_probability_interval=pfutureinterval)
                auto=correlations[initial]
                change=p0+auto['unmeasured_final_probability']-2*auto['real']
                change_error=max(abs(p0-p0interval[0]),abs(p0-p0interval[1]))+max(abs(auto['unmeasured_final_probability']-x) for x in auto['unmeasured_final_probability_interval'])+2*auto['complex_absolute_error_upper']
                assert -1e-10<=change<=1+1e-10
                assert abs(values['protected']['sequential_joint_probability']-values['hole']['sequential_joint_probability']-values['occupied_flux']['sequential_joint_probability'])<1e-11
                rows.append(dict(lag=lag,filtered_evolution_error_upper=aux_error,
                    full_filtered_path_error_upper=e_filtered,unfiltered_path_error_upper=e_forward,
                    destinations=values,two_time_correlations=correlations,
                    noninterventional_projector_change=float(change),
                    noninterventional_projector_change_interval=[max(0.,change-change_error),min(1.,change+change_error)]))
            if j<maxsteps:
                Y,_,increment=E(Y,step);aux_error+=increment
        out.append(dict(initial_channel=initial,initial_probability=p0,
            initial_probability_interval=p0interval,compressed_columns=Y.shape[1],
            filter_compression_residual=compression_error,
            hypothetical_binary_recording_work=work,
            binary_recording_work_error_upper=work_error,
            binary_recording_entropy_increase_bits=h2(p0) if energies is None else None,
            reference_binary_recording_entropy_increase_upper_bits=reference_entropy_gain_upper,
            rows=rows))
        print(json.dumps(dict(stage='persistence',preparation=label,initial=initial,
            p=p0,rank=Y.shape[1],lag2_same_channel=rows[-1]['destinations'][initial]['conditional_probability'],
            lag2_geometry=rows[-1]['destinations']['protected']['conditional_probability'])),flush=True)
    return dict(preparation=label,wait=wait,initial_wait_error_upper=wait_error,
        stationary_reference_residual_Frobenius=ref_leak,channels=out)


def run(snapshot_dir,quick=False):
    started=time.monotonic();full,energies,V,weights,beta=energy_basis()
    source=model(4,4,256,source=True);J,checks=inclusion(full,source)
    D,cycle=fixed_triangle(4,source[3]);trajectory=json.loads((UP/'N4_Q4_B256_T20.json').read_text())
    cases=[]
    for wait in ([10.] if quick else [5.,10.]):
        v=np.load(Path(snapshot_dir)/f'N4_Q4_B256_T{wait:g}.npz')['wave']
        X=(J@v)[:,None]
        err=min(trajectory['rows'],key=lambda a:abs(a['T']-wait))['state_error_upper']
        cases.append(case(f'neutral_sea_wait_{wait:g}',X,source,D,err,trajectory,wait,
                          output_steps=[0,1,5] if quick else None))
    if not quick:
        for label,w in weights.items():
            X=J@(V*np.sqrt(w)[None,:])
            cases.append(case(label,X,source,D,0.,trajectory,0.,energies=energies))
    out=dict(N=4,Q=4,budget=256,named_triangle=[1,2,3],cycle_mask=cycle,
        embedding_checks=checks,cases=cases,quick=quick,seconds=time.monotonic()-started,
        scope='The physical Hamiltonian is unchanged. Initial projectors refer to one fixed triangle and final detectors to that same triangle. Main evolution is unmeasured. Sequential conversion/survival probabilities belong to a hypothetical projective diagnostic; its binary recording work and entropy effect are explicit. Separately reported unforced two-time correlators and projector change are not classical joint probabilities. Reference mean energies are matched but resources and purity are not. Full-law outgoing residuals control subsequent evolution. Ordinary floating point, not formal roundoff intervals.')
    filename='PERSISTENCE_QUICK.json' if quick else 'PERSISTENCE.json'
    (P/filename).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(saved=filename,seconds=out['seconds'])),flush=True);return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--snapshot-dir',required=True);ap.add_argument('--quick',action='store_true')
    a=ap.parse_args();run(a.snapshot_dir,a.quick)
