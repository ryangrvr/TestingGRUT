"""Independent source-response spectra and finite-source Galerkin diagnostics.

Only three sine probes enter the primary whitening/Galerkin diagnostic.
The asymmetric beta probe is retained as an unwhitened scalar witness.
No physical production observable is evaluated by the synthetic precheck.
"""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm,eigh

PROBE_NAMES=('sine1','sine2','sine3','beta_left')
DEFAULT_TIMES=(0.,.29,.83,1.57)

def encode_matrix(matrix):
    x=np.asarray(matrix,dtype=complex)
    return np.stack((x.real,x.imag),axis=-1).tolist()

def decode_matrix(value):
    x=np.asarray(value)
    if x.ndim==3 and x.shape[-1]==2:return x[...,0].astype(float)+1j*x[...,1].astype(float)
    return np.asarray(value,dtype=complex)

def encode_vector(vector):
    x=np.asarray(vector,dtype=complex)
    return np.stack((x.real,x.imag),axis=-1).tolist()

def probe_weights(length,n):
    """Exact source-cell integrals, with explicit identity-valued edge offsets."""
    a=length/(n+1);edges=np.arange(1,n+1)*a
    weights=[];boundaries=[];integrals=[]
    for mode in (1,2,3):
        k=mode*np.pi/length
        weights.append((np.cos(k*edges[:-1])-np.cos(k*edges[1:]))/(2*mode))
        boundaries.append([(1-np.cos(k*a))/(2*mode),
                           (np.cos(k*(length-a))-(-1)**mode)/(2*mode)])
        integrals.append((1-(-1)**mode)/(2*mode))
    def beta_primitive(z):
        x=np.asarray(z)/length
        return 15*x*x-40*x**3+45*x**4-24*x**5+5*x**6
    weights.append(np.diff(beta_primitive(edges)))
    boundaries.append([beta_primitive(a),1-beta_primitive(length-a)])
    integrals.append(1.)
    W=np.array(weights,dtype=float);boundary=np.array(boundaries,dtype=float)
    return {'weights':W,'boundary_integrals':boundary,'constants':boundary.sum(axis=1)/4,
       'signed_integrals':np.array(integrals),'partition_errors':np.abs(W.sum(axis=1)+boundary.sum(axis=1)-integrals)}

def matrix_diagnostics(covariance,first_moment,responses,minimum_eigenvalue=1e-10,condition_limit=1e6):
    """Canonical positive-Gram whitening and a frozen finite-source predictor.

    responses is a dictionary time-string -> complex matrix or encoded pairs.
    Never clips or drops a Gram eigenvalue. The predictor is the source-span
    Galerkin projection; its error is not assumed to be a discretization error.
    """
    V=decode_matrix(covariance);K=decode_matrix(first_moment)
    gram_hermiticity=float(np.max(abs(V-V.conj().T)))
    moment_hermiticity=float(np.max(abs(K-K.conj().T)))
    if max(gram_hermiticity,moment_hermiticity)>1e-9:raise ValueError('Non-Hermitian source moments')
    V=(V+V.conj().T)/2;K=(K+K.conj().T)/2
    eigenvalues,U=np.linalg.eigh(V)
    if not np.all(np.isfinite(eigenvalues)) or eigenvalues[0]<=minimum_eigenvalue:
        raise ValueError('Primary source Gram is singular or below the frozen floor')
    condition=float(eigenvalues[-1]/eigenvalues[0])
    if condition>=condition_limit:raise ValueError('Primary source Gram exceeds frozen condition limit')
    inverse_sqrt=(U*(1/np.sqrt(eigenvalues)))@U.conj().T
    sqrt=(U*np.sqrt(eigenvalues))@U.conj().T
    generator=inverse_sqrt@K@inverse_sqrt;generator=(generator+generator.conj().T)/2
    centroid_spectrum=np.linalg.eigvalsh(generator)
    white={};cos_eigenvalues={};prediction={};residual={};white_prediction={};contraction_excess=0.
    for key,encoded in responses.items():
        t=float(key);M=decode_matrix(encoded);R=inverse_sqrt@M@inverse_sqrt
        predicted_white=expm(-1j*t*generator);predicted=sqrt@predicted_white@sqrt
        white[key]=encode_matrix(R);cos_eigenvalues[key]=np.linalg.eigvalsh((R+R.conj().T)/2).tolist()
        prediction[key]=encode_matrix(predicted);white_prediction[key]=encode_matrix(predicted_white)
        E=R-predicted_white
        residual[key]={'whitened_operator_norm':float(np.linalg.norm(E,2)),
            'whitened_matrix':encode_matrix(E),
            'normalized_diagonal_absolute':(np.abs(np.diag(M-predicted))/np.real(np.diag(V))).tolist()}
        contraction_excess=max(contraction_excess,float(np.linalg.norm(R,2))-1)
    return {'gram_eigenvalues':eigenvalues.tolist(),'gram_condition':condition,
      'gram_hermiticity_error':gram_hermiticity,'first_moment_hermiticity_error':moment_hermiticity,
      'projected_generator':encode_matrix(generator),'centroid_spectrum':centroid_spectrum.tolist(),
      'whitened_responses':white,'whitened_cos_eigenvalues':cos_eigenvalues,
      'galerkin_responses':prediction,'whitened_galerkin_responses':white_prediction,
      'galerkin_residuals':residual,'maximum_response_contraction_excess':max(0.,contraction_excess),
      'scope':'Finite-source projected generator; invariants remove invertible source-basis coordinate changes, not changes of source subspace.'}

def scalar_diagnostics(covariance,first_moment,responses):
    V=decode_matrix(covariance);K=decode_matrix(first_moment);variance=np.real(np.diag(V))
    if np.any(variance<=0):raise ValueError('Nonpositive auto-probe variance')
    normalized={key:encode_vector(np.diag(decode_matrix(M))/variance) for key,M in responses.items()}
    return {'variances':variance.tolist(),'centroids':(np.real(np.diag(K))/variance).tolist(),
      'normalized_auto_responses':normalized}

def reference(card,n,modes):
    """Massless matched-source means and vacuum response, without fitting."""
    length=float(card['length'])
    if card.get('mass_R',0)!=0:raise ValueError('This reference is massless only')
    times=card.get('response_times',DEFAULT_TIMES);a=length/(n+1);probes=probe_weights(length,n)
    W=probes['weights'];heights=W/a;edges=np.arange(1,n+1)*a
    k=np.arange(1,modes+1)*np.pi/length;omega=np.sqrt(k*k+1/np.pi)
    sine_integrals=np.sqrt(2/length)*(np.cos(edges[:-1,None]*k)-np.cos(edges[1:,None]*k))/k
    coefficients=heights@sine_integrals;factor=1/(2*np.pi*omega)
    covariance=(coefficients*factor)@coefficients.T
    first_moment=(W@W.T)/(2*np.pi*a)
    first_truncated=(coefficients@coefficients.T)/(2*np.pi)
    mu=1/np.sqrt(np.pi)
    mean_integrals=(np.sinh(mu*(edges[1:]-length/2))-np.sinh(mu*(edges[:-1]-length/2)))/(4*mu*np.cosh(mu*length/2))
    means=probes['constants']+heights@mean_integrals
    variation=np.abs(heights[:,0])+np.sum(np.abs(np.diff(heights,axis=1)),axis=1)+np.abs(heights[:,-1])
    tails=variation**2*length**2/(2*np.pi**4*modes**2)
    entry_tail=np.sqrt(tails[:,None]*tails[None,:])
    first_tail=length*np.outer(variation,variation)/(np.pi**3*modes)
    responses={f'{float(t):g}':encode_matrix((coefficients*(factor*np.exp(-1j*omega*t)))@coefficients.T) for t in times}
    primary_responses={key:encode_matrix(decode_matrix(M)[:3,:3]) for key,M in responses.items()}
    diagnostics=matrix_diagnostics(covariance[:3,:3],first_moment[:3,:3],primary_responses,
        minimum_eigenvalue=card['gram_minimum_eigenvalue'],condition_limit=card['gram_condition_limit'])
    E=float(np.linalg.norm(entry_tail[:3,:3],2));cmin=diagnostics['gram_eigenvalues'][0];cmax=diagnostics['gram_eigenvalues'][-1]
    # Canonical square-root whitening may rotate under a noncommuting Gram tail.
    # The condition-dependent estimate covers that coordinate effect as well.
    white_tail=E/cmin+E*cmax/(cmin*cmin)
    variance=np.diag(covariance);scalar_tail=2*tails/variance
    centroid_low=np.diag(first_moment)/(variance+tails);centroid_high=np.diag(first_moment)/variance
    return {'N':n,'modes':modes,'a':a,'probe_names':list(PROBE_NAMES),
      'probe_weights':W.tolist(),'probe_boundary_integrals':probes['boundary_integrals'].tolist(),
      'probe_constants':probes['constants'].tolist(),'probe_signed_integrals':probes['signed_integrals'].tolist(),
      'probe_partition_errors':probes['partition_errors'].tolist(),'means':means.tolist(),
      'covariance':covariance.tolist(),'first_moment':first_moment.tolist(),'fsum':(2*first_moment).tolist(),
      'first_moment_truncated':first_truncated.tolist(),'first_moment_tail_bounds':first_tail.tolist(),
      'responses':responses,'variance_tail_bounds':tails.tolist(),'entry_tail_bounds':entry_tail.tolist(),
      'source_variations':variation.tolist(),'primary_diagnostics':diagnostics,
      'scalar_diagnostics':scalar_diagnostics(covariance,first_moment,responses),
      'primary_whitened_response_tail_bound':white_tail,'primary_covariance_tail_norm':E,
      'normalized_auto_response_tail_bounds':scalar_tail.tolist(),
      'centroid_intervals':np.column_stack((centroid_low,centroid_high)).tolist(),
      'scope':'Exact cell-source f-sum by Parseval; truncated Gaussian response with analytic tails evaluated in floating point. Galerkin residual is a model-space diagnostic, not a proven continuum error.'}

def synthetic_precheck():
    """Independent finite-matrix tests, including a causal-interpretation control."""
    frequencies=np.array([.7,1.3,2.1]);times=np.array(DEFAULT_TIMES)
    source=np.array([[.3,.08,-.03],[.05,.27,.04],[-.02,.06,.22]])
    cutoff=2;dim=cutoff+1;eye=np.eye(dim);number=np.diag(np.arange(dim));ann=np.diag(np.sqrt(np.arange(1,dim)),1);position=ann+ann.T
    def local(op,index):
        return np.kron(np.kron(op if index==0 else eye,op if index==1 else eye),op if index==2 else eye)
    H=sum(w*local(number,j) for j,w in enumerate(frequencies));vac=np.zeros(dim**3);vac[0]=1
    operators=[(.1*(i+1))*np.eye(dim**3)+sum(source[i,j]*local(position,j) for j in range(3)) for i in range(3)]
    vectors=np.column_stack([(G-float(vac@G@vac)*np.eye(dim**3))@vac for G in operators])
    V=vectors.conj().T@vectors;K=vectors.conj().T@H@vectors
    responses={f'{t:g}':encode_matrix(vectors.conj().T@expm(-1j*H*t)@vectors) for t in times}
    expected={f'{t:g}':source@np.diag(np.exp(-1j*frequencies*t))@source.T for t in times}
    checks={'dense_vs_spectral_response':max(float(np.max(abs(decode_matrix(responses[k])-v))) for k,v in expected.items()),
      'static_gram':float(np.max(abs(V-source@source.T))),
      'first_moment':float(np.max(abs(K-source@np.diag(frequencies)@source.T)))}
    diagnostics=matrix_diagnostics(V,K,responses)
    checks['invariant_source_span_galerkin_exact']=max(v['whitened_operator_norm'] for v in diagnostics['galerkin_residuals'].values())
    checks['centroid_eigenvalues']=float(np.max(abs(np.array(diagnostics['centroid_spectrum'])-frequencies)))
    mixing=np.array([[1.3,.2,-.1],[.1,.8,.3],[-.2,.1,1.1]])
    mixed={k:encode_matrix(mixing@decode_matrix(M)@mixing.T) for k,M in responses.items()}
    md=matrix_diagnostics(mixing@V@mixing.T,mixing@K@mixing.T,mixed)
    checks['basis_mix_centroid_invariance']=float(np.max(abs(np.array(md['centroid_spectrum'])-diagnostics['centroid_spectrum'])))
    checks['basis_mix_cosine_spectrum_invariance']=max(float(np.max(abs(np.array(md['whitened_cos_eigenvalues'][k])-diagnostics['whitened_cos_eigenvalues'][k]))) for k in responses)
    checks['independent_generalized_eigensolver']=float(np.max(abs(eigh(K,V,eigvals_only=True)-frequencies)))
    weights=np.array([.04,.02,.01]);delta=.005*np.array([frequencies[1]-frequencies[2],frequencies[2]-frequencies[0],frequencies[0]-frequencies[1]])
    other=weights+delta
    checks['reweighted_source_same_variance']=float(abs(other.sum()-weights.sum()))
    checks['reweighted_source_same_centroid']=float(abs(other@frequencies-weights@frequencies))
    response=np.array([np.sum(weights*np.exp(-1j*frequencies*t))/weights.sum() for t in times])
    changed=np.array([np.sum(other*np.exp(-1j*frequencies*t))/other.sum() for t in times])
    shape_difference=float(np.max(abs(response-changed)))
    bad_phase=float(np.max(abs(response-np.conj(response))))
    # Hidden fourth mode makes the projected three-source subspace noninvariant.
    extra=np.column_stack([source,np.array([.11,.04,-.07])]);freq4=np.r_[frequencies,2.8]
    hiddenV=extra@extra.T;hiddenK=extra@np.diag(freq4)@extra.T
    hiddenM={f'{t:g}':encode_matrix(extra@np.diag(np.exp(-1j*freq4*t))@extra.T) for t in times}
    hidden=matrix_diagnostics(hiddenV,hiddenK,hiddenM)
    hidden_residual=max(v['whitened_operator_norm'] for v in hidden['galerkin_residuals'].values())
    # Explicitly reject a singular Gram; no eigenvalue clipping or silent truncation.
    singular_rejected=False
    try:matrix_diagnostics([[1,1],[1,1]],[[1,1],[1,1]],{'0':[[1,1],[1,1]]})
    except ValueError:singular_rejected=True
    return {'scope':'Synthetic oscillator/matrix data only; no physical lattice or matched continuum observables evaluated.',
      'frequencies':frequencies.tolist(),'source_matrix':source.tolist(),'times':times.tolist(),'cutoff':cutoff,
      'checks':checks,'maximum_identity_error':max(checks.values()),
      'same_H_variance_centroid_source_change_response_difference':shape_difference,
      'wrong_phase_error':bad_phase,'noninvariant_span_galerkin_residual':hidden_residual,
      'singular_gram_rejected':singular_rejected,
      'passed':bool(max(checks.values())<1e-12 and shape_difference>1e-3 and bad_phase>1e-3 and hidden_residual>1e-3 and singular_rejected)}

if __name__=='__main__':
    result=synthetic_precheck();out=Path('work/r1_precheck.json')
    with out.open('x',encoding='utf-8') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))
