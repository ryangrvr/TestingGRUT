"""Independent free-field reference for source-native probes and detector backaction.

No import of B2 continuum/history code and no lattice production evaluation.
The caller supplies the already-frozen detector weights and edge coefficients.
Probe cell weights are source integrals, not quadratic reconstruction weights.
"""
import json,itertools
import numpy as np
from scipy.linalg import expm

PROBE_NAMES=('sine1','sine2','sine3','beta_left')

def pair(z):
    return [float(np.real(z)),float(np.imag(z))]

def complex_matrix(matrix):
    return [[pair(z) for z in row] for row in matrix]

def probe_weights(length,n):
    """Return exact-integral cell coefficients and outer-interval integrals.

    sine_m(z)=pi/(2L)*sin(m*pi*z/L), each with L1 norm one;
    beta_left(z)=30/L*(z/L)*(1-z/L)^4, integral/L1 norm one.
    The sine integrals are 1,0,1/3; they must not all be normalized to one.
    """
    a=length/(n+1);edges=np.arange(1,n+1,dtype=float)*a;out={}
    for m in (1,2,3):
        k=m*np.pi/length
        primitive=lambda z:(1-np.cos(k*np.asarray(z)))/(2*m)
        v=np.diff(primitive(edges))
        boundary=np.array([primitive(a),primitive(length)-primitive(length-a)])
        integral=(1-(-1)**m)/(2*m)
        out[f'sine{m}']={'weights':v.tolist(),'boundary_integrals':boundary.tolist(),
          'constant':float(boundary.sum()/4),'integral':float(integral),
          'partition_error':float(abs(v.sum()+boundary.sum()-integral))}
    # Antiderivative of 30*x*(1-x)^4, x=z/L; exact polynomial coefficients.
    coeff=np.array([0.,0.,15.,-40.,45.,-24.,5.])
    primitive=lambda z:np.polynomial.polynomial.polyval(np.asarray(z)/length,coeff)
    v=np.diff(primitive(edges));boundary=np.array([primitive(a),primitive(length)-primitive(length-a)])
    out['beta_left']={'weights':v.tolist(),'boundary_integrals':boundary.tolist(),
        'constant':float(boundary.sum()/4),'integral':1.,
        'partition_error':float(abs(v.sum()+boundary.sum()-1))}
    return out

def spectrum(length,n,modes,weights,boundary):
    """Step-smearing sine coefficients, exact mean and conservative tail estimate.

    Observable = sum_j weights_j * (1/a integral_Ij Phi) + sum(boundary)/4.
    Boundary values of Phi are fixed to 1/4, and its mean is the cosh profile.
    """
    a=length/(n+1);edges=np.arange(1,n+1,dtype=float)*a
    heights=np.asarray(weights,dtype=float)/a;k=np.arange(1,modes+1,dtype=float)*np.pi/length
    coefficients=np.sqrt(2/length)*np.sum(heights[:,None]*(np.cos(edges[:-1,None]*k)-np.cos(edges[1:,None]*k)),axis=0)/k
    mu=1/np.sqrt(np.pi);omega=np.sqrt(k*k+mu*mu)
    mean=float(np.sum(boundary)/4+np.dot(heights,np.sinh(mu*(edges[1:]-length/2))-np.sinh(mu*(edges[:-1]-length/2)))
               /(4*mu*np.cosh(mu*length/2)))
    variation=float(abs(heights[0])+np.abs(np.diff(heights)).sum()+abs(heights[-1]))
    tail=float(variation**2*length**2/(2*np.pi**4*modes**2))
    return {'omega':omega,'coefficients':coefficients,'mean':mean,
      'variance':float(np.sum(coefficients**2/(2*np.pi*omega))),
      'tail_variance_bound':tail,'smearing_total_variation':variation}

def gaussian_matrices(meanF,varF,meanG,varG,covGF,coupling):
    """Gram, G, G^2 matrices in |Omega>, exp(-i(H-E0)t) W|Omega>.

    C=<delta G(t) delta F(0)>, chi=<exp(i*lambda F)>.
    The kick shifts G's mean by -2*lambda*Im(C), retaining its variance.
    """
    lam=coupling;chi=np.exp(1j*lam*meanF-.5*lam*lam*varF)
    shift=-2*lam*np.imag(covGF);cross_mean=meanG+1j*lam*covGF
    gram=np.array([[1,chi],[np.conj(chi),1]],dtype=complex)
    A=np.array([[meanG,chi*cross_mean],[np.conj(chi*cross_mean),meanG+shift]],dtype=complex)
    cross_second=chi*(cross_mean**2+varG)
    B=np.array([[meanG**2+varG,cross_second],[np.conj(cross_second),(meanG+shift)**2+varG]],dtype=complex)
    return gram,A,B

def moments_from_matrices(gram,A,B,phase=0.):
    """Normalized probe moments in vacuum, Weyl kick and K+/- states."""
    coefficients={'vacuum':np.array([1,0],complex),'kick':np.array([0,1],complex),
      'plus':np.array([1,np.exp(1j*phase)],complex)/2,
      'minus':np.array([1,-np.exp(1j*phase)],complex)/2}
    states={};imaginary=0.
    for name,x in coefficients.items():
        prob=x.conj()@gram@x
        if prob.real<=0:raise ValueError('Nonpositive state normalization')
        mean=(x.conj()@A@x)/prob;second=(x.conj()@B@x)/prob
        imaginary=max(imaginary,abs(prob.imag),abs(mean.imag),abs(second.imag))
        states[name]={'probability':float(prob.real),'mean':float(mean.real),
          'second_moment':float(second.real),'variance':float(second.real-mean.real**2)}
    return states,float(imaginary)

def ordered_backaction_matrices(meanF,meanG,varG,omega,hFF,hGF,times,probe_time,coupling):
    """All U_a^dagger {1,G,G^2} U_b insertions, with U_b ordered in time.

    paths b=(b0,...), U_b=W(t_last)^b_last ... W(t0)^b0.
    hFF/hGF are per-mode covariance contributions (not amplitudes).
    """
    omega=np.asarray(omega);hFF=np.asarray(hFF);hGF=np.asarray(hGF)
    r=len(times);paths=[tuple((index>>j)&1 for j in range(r)) for index in range(2**r)];size=len(paths)
    covF={(i,j):np.sum(hFF*np.exp(-1j*omega*(times[i]-times[j]))) for i in range(r) for j in range(r)}
    covG={i:np.sum(hGF*np.exp(-1j*omega*(probe_time-times[i]))) for i in range(r)}
    gram=np.zeros((size,size),complex);A=gram.copy();B=gram.copy()
    for ia,a in enumerate(paths):
        for ib,b in enumerate(paths):
            operators=[(i,-coupling*a[i]) for i in range(r)]+[(i,coupling*b[i]) for i in reversed(range(r))]
            exponent=1j*meanF*sum(q for _,q in operators)
            for position,(i,q) in enumerate(operators):
                exponent-=.5*q*q*covF[i,i]
                exponent-=sum(q*other*covF[i,j] for j,other in operators[position+1:])
            linear=sum(q*np.conj(covG[i]) for i,q in operators[:r])+sum(q*covG[i] for i,q in operators[r:])
            Z=np.exp(exponent);shifted=meanG+1j*linear
            gram[ia,ib]=Z;A[ia,ib]=Z*shifted;B[ia,ib]=Z*(shifted*shifted+varG)
    return paths,gram,A,B

def conditional_path_moments(paths,gram,A,B,phases):
    """Normalized probe moments for every recorded ± outcome path."""
    count=len(paths[0]);states={};imaginary=0.
    for signs in itertools.product((1,-1),repeat=count):
        x=np.array([np.prod([s**bit for s,bit in zip(signs,b)])*np.exp(1j*np.dot(phases,b))/(2**count) for b in paths])
        prob=x.conj()@gram@x
        if prob.real<=0:raise ValueError('Nonpositive branch probability')
        mean=x.conj()@A@x/prob;second=x.conj()@B@x/prob
        imaginary=max(imaginary,abs(prob.imag),abs(mean.imag),abs(second.imag))
        label=''.join('+' if s==1 else '-' for s in signs)
        states[label]={'probability':float(prob.real),'mean':float(mean.real),'second_moment':float(second.real),
                       'variance':float(second.real-mean.real**2)}
    return states,float(imaginary)

def reference(card,n,modes,weightsF,boundaryF):
    """JSON-serializable matched-continuum predictions; no histories evaluated.

    card: length, coupling, optional phase, and probe_times (default 0,.43,1.31).
    F is the fixed B2 functional supplied by the caller. G probes use native
    cell source weights and fixed edge constants, as returned by probe_weights.
    """
    length=float(card['length']);lam=float(card['coupling']);phase=float(card.get('phase',0))
    times=card.get('probe_times',[0.,.43,1.31]);a=length/(n+1)
    backaction_times=card.get('backaction_times',[0.,.43]);two_probe_time=card.get('two_backaction_probe_time',1.31)
    F=spectrum(length,n,modes,weightsF,boundaryF)
    probes=probe_weights(length,n)
    spectra={name:spectrum(length,n,modes,p['weights'],p['boundary_integrals']) for name,p in probes.items()}
    energy=float(lam*lam/(2*np.pi*a)*np.dot(weightsF,weightsF))
    omega=F['omega'];factor=1/(2*np.pi*omega);result={};maximag=0.
    covariance={name:{other:float(np.sum(spectra[name]['coefficients']*spectra[other]['coefficients']*factor))
      for other in PROBE_NAMES} for name in PROBE_NAMES}
    covariance_tail={name:{other:float(np.sqrt(spectra[name]['tail_variance_bound']*spectra[other]['tail_variance_bound']))
      for other in PROBE_NAMES} for name in PROBE_NAMES}
    for name,G in spectra.items():
        pertime={};cross_coeff=G['coefficients']*F['coefficients']*factor
        for t in times:
            C=np.sum(cross_coeff*np.exp(-1j*omega*t))
            gram,A,B=gaussian_matrices(F['mean'],F['variance'],G['mean'],G['variance'],C,lam)
            states,imag=moments_from_matrices(gram,A,B,phase);maximag=max(maximag,imag)
            for state,m in states.items():
                m['energy_above_vacuum']=0. if state=='vacuum' else (energy if state=='kick' else energy/(4*m['probability']))
            pertime[f'{t:g}']={'covariance_GF':pair(C),'gram':complex_matrix(gram),
              'mean_matrix':complex_matrix(A),'second_matrix':complex_matrix(B),'states':states}
        paths,gram,A,B=ordered_backaction_matrices(F['mean'],G['mean'],G['variance'],omega,
          F['coefficients']**2*factor,cross_coeff,backaction_times,two_probe_time,lam)
        two_states,imag=conditional_path_moments(paths,gram,A,B,[phase]*len(backaction_times));maximag=max(maximag,imag)
        two={'backaction_times':list(backaction_times),'probe_time':two_probe_time,'unitary_paths':[list(p) for p in paths],
          'gram':complex_matrix(gram),'mean_matrix':complex_matrix(A),'second_matrix':complex_matrix(B),'states':two_states}
        result[name]={'weights':probes[name],'mean':G['mean'],'variance':G['variance'],
          'tail_variance_bound':G['tail_variance_bound'],'smearing_total_variation':G['smearing_total_variation'],
          'GF_tail_bound':float(np.sqrt(G['tail_variance_bound']*F['tail_variance_bound'])),
          'times':pertime,'two_backaction':two}
    single={f'{t:g}':{'gram':result[PROBE_NAMES[0]]['times'][f'{t:g}']['gram'],
      'probes':{name:{'first_matrix':result[name]['times'][f'{t:g}']['mean_matrix'],
        'second_matrix':result[name]['times'][f'{t:g}']['second_matrix'],
        'states':result[name]['times'][f'{t:g}']['states'],
        'covariance_GF':result[name]['times'][f'{t:g}']['covariance_GF']} for name in PROBE_NAMES}} for t in times}
    sample=result[PROBE_NAMES[0]]['two_backaction']
    double={k:sample[k] for k in ('gram','unitary_paths','backaction_times','probe_time')}
    double['probes']={name:{'first_matrix':result[name]['two_backaction']['mean_matrix'],
      'second_matrix':result[name]['two_backaction']['second_matrix'],'states':result[name]['two_backaction']['states']}
      for name in PROBE_NAMES}
    metadata={name:{k:result[name][k] for k in ('mean','variance','tail_variance_bound','smearing_total_variation','GF_tail_bound')}
      for name in PROBE_NAMES}
    return {'N':n,'modes':modes,'a':a,'detector':{k:F[k] for k in ('mean','variance','tail_variance_bound','smearing_total_variation')},
      'energy_kick':energy,'probe_weights':probes,'probe_metadata':metadata,'single':single,'double':double,
      'probe_covariance':covariance,'probe_covariance_tail_bounds':covariance_tail,
      'max_imaginary_residual':maximag,
      'scope':'Free matched step-smearing reference; energy from exact Parseval identity, tail estimates evaluated in floating point.'}

def synthetic_oscillator_precheck():
    """Independent Hilbert-space check; no production lattice observable data."""
    cfg={'mF':.17,'mG':-.08,'hF':.21,'hG':.13,'omega':.83,'coupling':2.,
      'times':[0.,.43,1.31],'phases':[0.,.29],'cutoffs':[24,32]}
    records=[];two_records=[];maxerror=0.;wrong_sign_error=0.
    for cutoff in cfg['cutoffs']:
        annihilation=np.diag(np.sqrt(np.arange(1,cutoff+1)),1);q=annihilation+annihilation.T
        F=cfg['mF']*np.eye(cutoff+1)+cfg['hF']*q;G=cfg['mG']*np.eye(cutoff+1)+cfg['hG']*q
        vacuum=np.zeros(cutoff+1,complex);vacuum[0]=1
        W=expm(1j*cfg['coupling']*F)
        for t in cfg['times']:
            basis=np.column_stack([vacuum,np.exp(-1j*cfg['omega']*np.arange(cutoff+1)*t)*(W@vacuum)])
            gram_direct=basis.conj().T@basis;A_direct=basis.conj().T@G@basis;B_direct=basis.conj().T@G@G@basis
            C=cfg['hG']*cfg['hF']*np.exp(-1j*cfg['omega']*t)
            gram,A,B=gaussian_matrices(cfg['mF'],cfg['hF']**2,cfg['mG'],cfg['hG']**2,C,cfg['coupling'])
            matrix_error=max(float(np.max(abs(gram-gram_direct))),float(np.max(abs(A-A_direct))),float(np.max(abs(B-B_direct))))
            moment_error=0.
            for phase in cfg['phases']:
                states,imag=moments_from_matrices(gram,A,B,phase)
                for name,x in {'vacuum':[1,0],'kick':[0,1],'plus':np.array([1,np.exp(1j*phase)])/2,
                               'minus':np.array([1,-np.exp(1j*phase)])/2}.items():
                    v=basis@x;prob=np.vdot(v,v).real
                    direct={'probability':prob,'mean':float(np.vdot(v,G@v).real/prob),'second_moment':float(np.vdot(v,G@G@v).real/prob)}
                    moment_error=max(moment_error,max(abs(direct[k]-states[name][k]) for k in direct))
            maxerror=max(maxerror,matrix_error,moment_error)
            # Wrong sign for the coherent shift must fail at a nonzero time.
            wrong_sign_error=max(wrong_sign_error,float(abs(cfg['mG']+2*cfg['coupling']*C.imag-A_direct[1,1])))
            records.append({'cutoff':cutoff,'time':t,'matrix_error':matrix_error,'state_moment_error':moment_error})
        backaction_times=[0.,.43];probe_time=1.31;energies=cfg['omega']*np.arange(cutoff+1)
        paths,gram,A,B=ordered_backaction_matrices(cfg['mF'],cfg['mG'],cfg['hG']**2,
          [cfg['omega']],[cfg['hF']**2],[cfg['hG']*cfg['hF']],backaction_times,probe_time,cfg['coupling'])
        vectors=[]
        for b in paths:
            v=vacuum.copy()
            if b[0]:v=W@v
            v=np.exp(-1j*energies*.43)*v
            if b[1]:v=W@v
            vectors.append(np.exp(-1j*energies*(probe_time-.43))*v)
        basis=np.column_stack(vectors);gram_direct=basis.conj().T@basis;A_direct=basis.conj().T@G@basis;B_direct=basis.conj().T@G@G@basis
        matrix_error=max(float(np.max(abs(gram-gram_direct))),float(np.max(abs(A-A_direct))),float(np.max(abs(B-B_direct))))
        moment_error=0.
        for phase in cfg['phases']:
            states,imag=conditional_path_moments(paths,gram,A,B,[phase,phase])
            for signs in itertools.product((1,-1),repeat=2):
                x=np.array([np.prod([s**bit for s,bit in zip(signs,b)])*np.exp(1j*phase*sum(b))/4 for b in paths])
                v=basis@x;prob=np.vdot(v,v).real;label=''.join('+' if s==1 else '-' for s in signs)
                direct={'probability':prob,'mean':float(np.vdot(v,G@v).real/prob),'second_moment':float(np.vdot(v,G@G@v).real/prob)}
                moment_error=max(moment_error,max(abs(direct[k]-states[label][k]) for k in direct))
        maxerror=max(maxerror,matrix_error,moment_error)
        two_records.append({'cutoff':cutoff,'matrix_error':matrix_error,'conditional_moment_error':moment_error})
    return {'configuration':cfg,'checks':records,'two_backaction_checks':two_records,'max_error':maxerror,
      'wrong_shift_sign_negative_control_error':wrong_sign_error,'wrong_shift_sign_rejected':wrong_sign_error>1e-3,
      'passed':bool(maxerror<1e-12 and wrong_sign_error>1e-3),
      'scope':'Synthetic single oscillator only; production source observables have not been evaluated.'}

if __name__=='__main__':
    print(json.dumps(synthetic_oscillator_precheck(),indent=2))
