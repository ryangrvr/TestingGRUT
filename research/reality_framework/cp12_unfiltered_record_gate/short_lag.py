"""Independent operator test: the first remote density response is coherence.

This computes exact matrix identities away from regulator boundaries and
tests them with an independently built physical Gauss basis. It does not
replace finite-time evolution by a Taylor extrapolation.
"""
from itertools import combinations
from math import sqrt
import argparse,json,time
import numpy as np
from scipy.sparse import coo_matrix,diags,eye
from scipy.sparse.linalg import expm_multiply
from audit import P,REL,model,inclusion,detectors,energy_basis,divide,positive_interval
from frozen_model import build,hop,gauss


def comm_diagonal(values,A):
    C=A.tocoo()
    return coo_matrix((C.data*(values[C.row]-values[C.col]),(C.row,C.col)),shape=A.shape).tocsr()


def double_response(H,S,D):
    first=-comm_diagonal(D,H)
    second=H@first-first@H
    return -.5j*comm_diagonal(S,second)


def dressed_transfer(n,states,sub=None,sizes=None):
    """sum_j a0^dag U_(0e) U_(je)^dag a_j, including n_graph(e)."""
    es=list(combinations(range(n),2));index={s:i for i,s in enumerate(states)}
    rr=[];cc=[];vv=[];outside=0
    for col,(g,b,E) in enumerate(states):
        for j in range(1,n):
            e=es.index((0,j))
            if not(g>>e&1):continue
            a=hop(b,0,j)
            if a is None:continue
            bits,fermion_sign=a;F=list(E);F[2*e]+=1;F[2*e+1]-=1
            target=(g,bits,tuple(F));sign=fermion_sign;factor=1.
            if sub is not None:
                target,symmetry_sign,size=sub.canonical(target)
                sign*=symmetry_sign
                if not sign:continue
                factor=sqrt(sizes[col]/size)
            if target not in index:outside+=1;continue
            rr.append(index[target]);cc.append(col);vv.append(sign*factor)
    W=coo_matrix((vv,(rr,cc)),shape=(len(states),len(states))).tocsr()
    K=.5j/(n-1)*(W-W.T)
    return W,K,outside


def independent_operator_control():
    start=time.monotonic();n=q=4;L=3
    H,states,lookup,degrees=build(n,q,L)
    O,events=detectors(n,states);S=O['source_density'];D=O['density'];C=O['circulation']
    W,K,outside=dressed_transfer(n,states)
    safe=np.array([max(map(abs,s[2]))<=L-2 for s in states])
    K2=double_response(H,S,D)
    delta=(K2-K)@diags(safe.astype(float))
    error=float(np.linalg.norm(delta.data));assert error<1e-12
    first=comm_diagonal(S,-comm_diagonal(D,H))
    first_error=float(np.linalg.norm(first.data));assert first_error<1e-12
    K2C=double_response(H,S,C)
    c_error=float(np.linalg.norm(K2C.data));assert c_error<1e-12
    K2CD=double_response(H,S,C*D)
    cd_delta=(K2CD-diags(C)@K)@diags(safe.astype(float))
    cd_error=float(np.linalg.norm(cd_delta.data));assert cd_error<1e-12
    pinch_errors={label:float(np.linalg.norm(comm_diagonal(mask,K).data)) for label,mask in events.items()}
    assert max(pinch_errors.values())<1e-12
    # Independently inspect the compressed physical matrix spectrum. The
    # unregulated bound follows from a two-mode fermionic star in angle space.
    from scipy.sparse.linalg import eigsh
    actual_norm=float(abs(eigsh(K,k=1,which='LM',return_eigenvectors=False)[0]))
    norm_bound=1/(2*sqrt(n-1))
    assert actual_norm<=norm_bound+1e-10
    es=list(combinations(range(n),2))
    degree=np.array([sum(bool(s[0]>>e&1) for e,(i,j) in enumerate(es) if i==0) for s in states])
    capacity_excess=K@K-diags(degree/(4*(n-1)**2))
    capacity_violation=float(eigsh(capacity_excess,k=1,which='LA',return_eigenvectors=False)[0])
    assert capacity_violation<1e-10
    edge=H.tocoo();unchanged=np.array([states[i][1:]==states[j][1:] for i,j in zip(edge.row,edge.col)])
    nonhop=coo_matrix((edge.data[unchanged],(edge.row[unchanged],edge.col[unchanged])),shape=H.shape).tocsr()
    inert=float(np.linalg.norm(comm_diagonal(S,nonhop).data));assert inert==0
    # A separate physical-configuration unfolding checks the signed quotient
    # and all new readouts on a genuinely different representation.
    sub,HP,BP,ss,sz,psi=model(4,4,8,source=True)
    index={s:i for i,s in enumerate(ss)};rr=[];cc=[];vv=[]
    inside=np.array([sum(x*x for x in s[2])<=8 for s in states])
    for i,s in enumerate(states):
        if not inside[i]:continue
        key,sign,size=sub.canonical(s)
        if sign and key in index:rr.append(i);cc.append(index[key]);vv.append(sign/sqrt(size))
    J=coo_matrix((vv,(rr,cc)),shape=(len(states),len(ss))).tocsr()
    physical_projection=diags(inside.astype(float));HR=physical_projection@H@physical_projection
    iso=float(np.linalg.norm((J.T@J-eye(len(ss))).data))
    intertwine=float(np.linalg.norm((HR@J-J@HP).data))
    outsideH=diags((~inside).astype(float))@H@J
    gram=float(np.linalg.norm((outsideH.T@outsideH-BP.T@BP).data))
    OP,EP=detectors(4,ss)
    observer=max(float(np.linalg.norm((diags(O[k])@J-J@diags(OP[k])).data)) for k in O)
    assert max(iso,intertwine,gram,observer)<1e-10
    wait=.7;lag=.4
    xp=expm_multiply(-1j*wait*HP,psi);xr=expm_multiply(-1j*wait*HR,J@psi)
    maxdiff=0.;cases=0
    for channel,mask in [('unfiltered',np.ones(len(ss)))]+list(EP.items()):
        physical_mask=np.ones(len(states)) if channel=='unfiltered' else events[channel]
        xp0=mask*xp;xr0=physical_mask*xr
        zp=expm_multiply(-1j*lag*HP,np.column_stack([xp0,OP['source_density']*xp0]))
        zr=expm_multiply(-1j*lag*HR,np.column_stack([xr0,O['source_density']*xr0]))
        for name in O:
            a=float(2*np.vdot(zp[:,0],OP[name]*zp[:,1]).imag)
            b=float(2*np.vdot(zr[:,0],O[name]*zr[:,1]).imag)
            maxdiff=max(maxdiff,abs(a-b));cases+=1
    assert maxdiff<1e-10
    return dict(physical_dimension=len(states),safe_interior_columns=int(safe.sum()),
        first_order_remote_response_operator_norm=first_error,
        dressed_transfer_second_order_residual=error,
        circulation_second_order_response_operator_norm=c_error,
        circulation_weighted_second_order_residual=cd_error,
        leading_response_record_commutator_norms=pinch_errors,
        dressed_transfer_operator_norm=actual_norm,
        unregulated_global_coefficient_norm_bound=norm_bound,
        operator_square_degree_bound_maximum_violation=capacity_violation,
        no_hopping_source_commutator_norm=inert,
        quotient_dimension=len(ss),isometry_error=iso,H_intertwiner_error=intertwine,
        outgoing_Gram_error=gram,readout_intertwiner_error=observer,
        independent_response_comparisons=cases,maximum_response_difference=maxdiff,
        discarded_transfer_targets_at_regulator_boundary=outside,
        seconds=time.monotonic()-start)


def main_coefficients(snapshot_dir):
    rows=[]
    for n,budget in [(4,256),(5,128)]:
        full=model(n,n,budget);source=model(n,n,budget,source=True)
        J,checks=inclusion(full,source);sub,H,B,states,sizes,psi=source
        O,events=detectors(n,states);W,K,outside=dressed_transfer(n,states,sub,sizes)
        x=J@np.load(str(snapshot_dir)+f'/N{n}_Q{n}_B{budget}_T10.npz')['wave']
        plain=float(np.vdot(x,K@x).real)
        assert abs(plain)<1e-12
        weighted=float(np.vdot(x,O['circulation']*(K@x)).real)
        conjugate=float(np.vdot(x.conj(),O['circulation']*(K@x.conj())).real)
        assert abs(weighted+conjugate)<1e-12
        p=float(np.linalg.norm(O['circulation']*x)**2)
        old=json.loads((REL/f'RELATIONAL_N{n}_Q{n}_B{budget}.json').read_text())
        ew=min(old['rows'],key=lambda r:abs(r['T']-10))['state_error_upper']
        knorm=1/(2*sqrt(n-1));allowance=knorm*(2*ew+ew**2)
        weighted_interval=[weighted-allowance,weighted+allowance]
        p_interval=positive_interval(p,ew)
        conditional_interval=divide(weighted_interval,p_interval)
        record_differences={label:float(np.vdot(x,K@x).real-
            np.vdot(mask*x,K@(mask*x)).real-
            np.vdot((1-mask)*x,K@((1-mask)*x)).real) for label,mask in events.items()}
        assert max(map(abs,record_differences.values()))<1e-12
        samples=[]
        if n==4:
            for lag in [.005,.01,.02,.05,.1,.2,.5]:
                z=expm_multiply(-1j*lag*H,np.column_stack([x,O['source_density']*x]))
                chi={name:float(2*np.vdot(z[:,0],d*z[:,1]).imag) for name,d in O.items()}
                samples.append(dict(lag=lag,response=chi,density_over_lag_squared=chi['density']/lag**2,
                    circulation_density_over_lag_squared=chi['circulation_density']/lag**2))
        reference_checks=None
        if n==4:
            _,energies,V,weights,beta=energy_basis()
            reference_checks={}
            for label,w in weights.items():
                X=J@(V*np.sqrt(w)[None,:])
                value=float(np.vdot(X,O['circulation'][:,None]*(K@X)).real)
                assert abs(value)<1e-12
                reference_checks[label]=value
        rows.append(dict(N=n,Q=n,budget=budget,wait=10,projected_density_T2_coefficient=plain,
            projected_circulation_weighted_density_T2_coefficient=weighted,
            full_law_circulation_weighted_density_T2_coefficient_interval=weighted_interval,
            coefficient_global_norm_bound=knorm,input_state_error_upper=ew,
            exchangeable_fixed_Q_weighted_coefficient_bound=float(sqrt(2*n/n)/(2*(n-1))),
            conjugated_state_weighted_coefficient=conjugate,
            coefficient_absolute_error_upper=allowance,
            projected_initial_circulation_probability=p,
            projected_conditional_density_T2_coefficient=weighted/p,
            full_law_conditional_density_T2_coefficient_interval=conditional_interval,
            real_energy_reference_weighted_T2_coefficients=reference_checks,
            leading_record_dephasing_differences=record_differences,finite_lag_samples=samples))
    return rows


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--snapshot-dir',required=True)
    a=ap.parse_args();result=dict(status='CP12_SHORT_LAG_CONTROLS_PASS',
        independent_physical_control=independent_operator_control(),
        main_state_coefficients=main_coefficients(a.snapshot_dir),
        scope='Short-lag coefficients are separate from the finite-lag calculations. Exact operator identities are analytically unregulated; independent numerical matrix checks restrict input columns away from the separate electric regulator boundary. Weighted and conditional leading coefficients carry explicit full-law input-error intervals from a global bounded operator norm. The additional finite-lag samples are projected diagnostics for local Taylor behavior, not uniform-in-size or long-time predictions. Full permutation symmetry makes the unweighted quadratic coefficient exactly zero.')
    (P/'SHORT_LAG_CONTROLS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
