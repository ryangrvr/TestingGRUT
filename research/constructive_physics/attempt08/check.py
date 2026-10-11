"""CP-08 exact finite causal-order dynamics. NumPy/SciPy, no parameter scan."""
import itertools, json, math
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix, csr_matrix, kron, eye, diags
from scipy.sparse.linalg import expm_multiply, eigsh
J=1.;Omega=.25;gamma=math.pi/4

def posets(n):
    pairs=[(i,j) for i in range(n) for j in range(n) if i!=j]
    states=[]
    for mask in range(1<<len(pairs)):
        C=np.zeros((n,n),bool)
        for a,(i,j) in enumerate(pairs):C[i,j]=bool(mask>>a&1)
        if np.any(C&C.T):continue
        if np.all(((C.astype(int)@C.astype(int))>0)<=C):states.append(mask)
    return states,pairs

def pauli(n,i,kind):
    basis=np.arange(1<<n); bits=(basis>>(n-1-i))&1
    if kind=='Z':return diags(1-2*bits,dtype=complex,format='csr')
    return coo_matrix((np.ones(1<<n),(basis^(1<<(n-1-i)),basis)),shape=(1<<n,1<<n)).tocsr()

def build(n,gam=gamma):
    states,pairs=posets(n); index={m:a for a,m in enumerate(states)};M=1<<n
    z=[pauli(n,i,'Z') for i in range(n)];x=[pauli(n,i,'X') for i in range(n)]
    U={(i,j):math.cos(gam)*eye(M,dtype=complex)-1j*math.sin(gam)*(z[i]@x[j]) for i,j in pairs}
    H=kron(eye(len(states)),Omega*sum(x),format='csr')
    transitions=[]
    for a,mask in enumerate(states):
        for b,(i,j) in enumerate(pairs):
            new=mask|(1<<b)
            if new!=mask and new in index:
                dest=index[new];T=coo_matrix(([1.],([dest],[a])),shape=(len(states),len(states)))
                term=-J*kron(T,U[i,j],format='csr')
                H=H+term+term.conj().T;transitions.append((a,dest,i,j))
    return H.tocsr(),states,pairs,index,U,transitions

def probabilities(v,n):return np.sum(abs(v.reshape(-1,1<<n))**2,axis=1)

def geometry(p,states,pairs,n):
    hist={}; rel=0.;height=0.
    for prob,mask in zip(p,states):
        C=np.zeros((n,n),bool)
        for a,(i,j) in enumerate(pairs):C[i,j]=bool(mask>>a&1)
        h=np.ones(n,int)
        for _ in range(n):
            nxt=h.copy()
            for j in range(n):nxt[j]=1+max(h[np.flatnonzero(C[:,j])],default=0)
            if np.array_equal(h,nxt):break
            h=nxt
        ht=int(h.max());hist[ht]=hist.get(ht,0)+float(prob)
        height+=prob*ht;rel+=prob*mask.bit_count()
    return {'mean_height':float(height),'mean_relations':float(rel),'height_distribution':hist}

def permutation_map(states,pairs,n):
    perm=list(range(n));perm[0],perm[1]=perm[1],perm[0]
    pairidx={p:a for a,p in enumerate(pairs)};lookup={m:a for a,m in enumerate(states)}
    mapping=[]
    for mask in states:
        new=0
        for a,(i,j) in enumerate(pairs):
            if mask>>a&1:new|=1<<pairidx[perm[i],perm[j]]
        for matter in range(1<<n):
            v=0
            for i in range(n):v|=((matter>>(n-1-i))&1)<<(n-1-perm[i])
            mapping.append(lookup[new]*(1<<n)+v)
    return mapping

def evolve(H,v,t):return expm_multiply((-1j*t)*H,v)

def test():
    results={'cases':[]};errors={}
    for n,count in [(3,19),(4,219)]:
        H,states,pairs,index,U,trans=build(n);assert len(states)==count
        assert np.max(abs((H-H.conj().T).data),initial=0)<1e-12
        mapping=permutation_map(states,pairs,n)
        delta=H[mapping,:][:,mapping]-H
        covariance=float(np.max(abs(delta.data),initial=0));assert covariance<1e-12
        v=np.zeros(H.shape[0],complex);v[0]=1;vt=evolve(H,v,1.)
        norm=abs(np.vdot(vt,vt)-1); inverse=float(np.linalg.norm(evolve(H,vt,-1.)-v))
        energy=abs(np.vdot(vt,H@vt)-np.vdot(v,H@v));assert max(norm,inverse,energy)<1e-10
        # A scalar link-count phase is a gradient on configuration space.
        D=np.repeat(np.exp(1j*math.pi/3*np.array([m.bit_count() for m in states])),1<<n)
        HD=diags(D)@H@diags(D.conj());vD=evolve(HD,D*v,1.)
        pure_gauge=float(np.linalg.norm(vD-D*vt));assert pure_gauge<1e-10
        # Transport is instead tested through the same-Hamiltonian gamma=0 control.
        H0,*_=build(n,0.)
        v0=evolve(H0,v,1.); tv=float(np.sum(abs(probabilities(vt,n)-probabilities(v0,n)))/2)
        assert tv>1e-3
        groundE,ground=eigsh(H,k=1,which='SA',tol=1e-10,v0=np.ones(H.shape[0]))
        results['cases'].append({'N':n,'posets':len(states),'Hilbert_dimension':H.shape[0],
          'reversible_configuration_edges':len(trans),'permutation_error':covariance,
          'norm_error':float(norm),'inverse_error':inverse,'energy_error':float(energy),
          'pure_gradient_phase_error':pure_gauge,'geometry_TV_vs_gamma0_at_t1':tv,
          't1_geometry':geometry(probabilities(vt,n),states,pairs,n),
          'ground_energy':float(groundE[0]),
          'ground_geometry':geometry(probabilities(ground[:,0],n),states,pairs,n)})
        if n==3:
            # Two different relation-addition paths end at the SAME final poset.
            A=U[0,1].toarray();B=U[1,2].toarray();matter=np.eye(1,8,0).ravel()
            hol=A.conj().T@B.conj().T@A@B
            hol_trace=float(np.trace(hol).real/8);assert abs(hol_trace-.5)<1e-12
            coefficient=float(np.linalg.norm((A@B+B@A)@matter)**2/4)
            assert abs(coefficient-.75)<1e-12
            init=np.zeros(H.shape[0],complex);init[index[1<<pairs.index((0,2))]*8]=1
            target=sum(1<<pairs.index(p) for p in [(0,1),(1,2),(0,2)])
            small=[]
            for t in [.04,.02,.01]:
                prob=probabilities(evolve(H,init,t),n)[index[target]]
                small.append({'t':t,'probability_over_t4':float(prob/t**4)})
            assert abs(small[-1]['probability_over_t4']-.75)<1e-3
            results['same_endpoint_interference']={'normalized_holonomy_trace':hol_trace,
                'exact_short_time_probability_coefficient':coefficient,
                'gamma0_coefficient':1.,'full_Hamiltonian_controls':small}
            # Response and fluctuations of one relation, from the SAME Hamiltonian.
            prepared=evolve(H,v,.3)
            a=np.repeat(np.array([float(bool(m>>(pairs.index((0,1)))&1)) for m in states]),8)
            t=.4;w=evolve(H,prepared,t);z=evolve(H,a*prepared,t)
            corr=np.vdot(w,a*z);mean0=float(np.vdot(prepared,a*prepared).real);meant=float(np.vdot(w,a*w).real)
            retarded=float(-2*corr.imag);noise=float(corr.real-mean0*meant)
            eps=1e-5
            vals=[float(np.vdot(u,a*u).real) for u in [evolve(H,np.exp(1j*eps*a)*prepared,t),evolve(H,np.exp(-1j*eps*a)*prepared,t)]]
            derivative=(vals[0]-vals[1])/(2*eps);assert abs(derivative-retarded)<1e-9
            covariance_matrix=np.array([[mean0*(1-mean0),noise],[noise,meant*(1-meant)]])
            assert np.linalg.eigvalsh(covariance_matrix).min()>-1e-12
            # Recording same relation distinguishes histories; before later evolution
            # it leaves diagonal probabilities unchanged, but afterward matters.
            r0=(1-a)*prepared;r1=a*prepared
            reduced=np.outer(r0,r0.conj())+np.outer(r1,r1.conj())
            assert np.max(abs(np.diag(reduced)-abs(prepared)**2))<1e-12
            recorded_final=abs(evolve(H,r0,t))**2+abs(evolve(H,r1,t))**2
            history_record_effect=float(np.sum(abs(recorded_final-abs(w)**2))/2)
            assert history_record_effect>1e-3
            # Explicit copying and inverse copying are an exact unitary permutation.
            ancillary=np.stack([prepared,np.zeros_like(prepared)],axis=1)
            copied=ancillary.copy();selected=a.astype(bool);copied[selected]=copied[selected,::-1]
            restored=copied.copy();restored[selected]=restored[selected,::-1]
            assert np.linalg.norm(restored-ancillary)==0
            # Affine mixture channel is evaluated by independent pure-state columns.
            s2=v.copy();s2[0]=0;s2[1]=1;out1=evolve(H,v,t);out2=evolve(H,s2,t)
            F=np.stack([out1,out2],axis=1)
            rA=np.array([[.5,.5],[.5,.5]]);rB=np.diag([.2,.8]);mix=.37
            ml=F@(mix*rA+(1-mix)*rB)@F.conj().T
            mr=mix*(F@rA@F.conj().T)+(1-mix)*(F@rB@F.conj().T)
            mixture=float(np.linalg.norm(ml-mr));assert mixture<1e-12
            # A separate stationary ground-state probe avoids the zero response of
            # this particular evolving initial state. The Hamiltonian is unchanged.
            gs=ground[:,0];gst=evolve(H,gs,t);source=evolve(H,a*gs,t)
            cgs=np.vdot(gst,a*source);rg=float(-2*cgs.imag)
            gp=evolve(H,np.exp(1j*eps*a)*gs,t);gm=evolve(H,np.exp(-1j*eps*a)*gs,t)
            rnum=float((np.vdot(gp,a*gp).real-np.vdot(gm,a*gm).real)/(2*eps))
            assert abs(rg-rnum)<1e-9 and abs(rg)>.001
            # Independent spectral reconstruction from the full finite Hamiltonian.
            energies,modes=np.linalg.eigh(H.toarray())
            gs_exact=modes[:,0]; gaps=energies-energies[0]
            spectral_weights=abs(modes.conj().T@(a*gs_exact))**2
            excited=gaps>1e-8
            spectral_response=float(2*np.sum(spectral_weights[excited]*np.sin(gaps[excited]*t)))
            spectral_noise=float(np.sum(spectral_weights[excited]*np.cos(gaps[excited]*t)))
            susceptibility=float(2*np.sum(spectral_weights[excited]/gaps[excited]))
            shifts=[]
            for sign in [1,-1]:
                _,mshift=np.linalg.eigh(H.toarray()-sign*eps*np.diag(a))
                shifts.append(float(np.vdot(mshift[:,0],a*mshift[:,0]).real))
            static_derivative=(shifts[0]-shifts[1])/(2*eps)
            assert abs(spectral_response-rg)<1e-10
            assert abs(static_derivative-susceptibility)<1e-8
            results['spectral_response']={'retarded_from_spectrum':spectral_response,
                'connected_noise_from_spectrum':spectral_noise,'static_susceptibility':susceptibility,
                'independent_ground_state_source_derivative':static_derivative,
                'positive_spectral_weight_total':float(spectral_weights[excited].sum()),
                'gap':float(gaps[excited].min())}
            results['stationary_ground_probe']={'retarded_response':rg,'independent_kick_derivative':rnum,
                'eigenvector_residual':float(np.linalg.norm(H@gs-groundE[0]*gs))}
            # Exact two-layer hypercube variational state: all allowed transports
            # commute on a fixed bipartition, so its hopping energy is -J*a*b.
            partition_edges=[(0,1),(0,2)]
            trial=np.zeros_like(gs)
            initial_matter=np.eye(1,8,0).ravel()
            for bits in itertools.product([0,1],repeat=2):
                mask=0;matter=initial_matter.copy()
                for bit,eij in zip(bits,partition_edges):
                    if bit:mask |= 1<<pairs.index(eij);matter=U[eij]@matter
                trial[index[mask]*8:(index[mask]+1)*8]=matter/2
            HT=H-kron(eye(len(states)),Omega*sum(pauli(n,i,'X') for i in range(n)),format='csr')
            variational=float(np.vdot(trial,HT@trial).real)
            assert abs(variational+2)<1e-12
            results['layered_variational_control']={'N':3,'partition_sizes':[1,2],
                'exact_hopping_expectation':variational,
                'general_total_energy_upper_bound':'-J*floor(N^2/4)+N*abs(Omega)'}
            results['response_records']={'relation_mean_initial':mean0,'relation_mean_final':meant,
               'retarded_impulse_response':retarded,'independent_kick_derivative':derivative,
               'symmetric_connected_noise':noise,'noise_matrix_min_eigenvalue':float(np.linalg.eigvalsh(covariance_matrix).min()),
               'record_effect_on_later_joint_probabilities_TV':history_record_effect,'record_inverse_error':0.,'mixture_error':mixture}
    results['verdict']='FINITE_RECOMBINING_MATTER_GEOMETRY_DYNAMICS_WORKS; CONTINUUM_GRAVITY_NOT_ESTABLISHED'
    Path(__file__).with_name('results.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps(results,indent=2));print('CP08_REVERSIBLE_CAUSAL_DYNAMICS_CHECKS_PASS')
if __name__=='__main__':test()
