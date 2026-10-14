"""Independent directed-graph, physical-basis and dual-bound controls."""
from itertools import combinations,product
from math import sqrt
import json,time
import numpy as np
from scipy.linalg import expm
from scipy.sparse import coo_matrix,diags,eye
from scipy.sparse.linalg import expm_multiply
from loops import P,model,masks,incidence_masks,directed_support,fixed_triangle,circuits
from persistence import compressed_filter
from evolution import Evolution
from frozen_model import build,gauss,wilson_vector


def oriented_graphs():
    n=4;es=list(combinations(range(n),2));checked=0
    for codes in product([0,1,-1],repeat=len(es)):
        A=np.zeros((n,n),int);flux=positive=0
        for e,((a,b),s) in enumerate(zip(es,codes)):
            if not s:continue
            flux|=1<<e
            if s>0:A[a,b]=1;positive|=1<<e
            else:A[b,a]=1
        reach=np.eye(n,dtype=int);power=np.eye(n,dtype=int)
        for _ in range(1,n):power=power@A;reach+=power
        expected=sum(1<<e for e,((a,b),s) in enumerate(zip(es,codes))
                     if (s>0 and reach[b,a]) or (s<0 and reach[a,b]))
        assert directed_support(n,flux,positive)==expected
        checked+=1
    return checked


def physical():
    HF,physical_states,_,_=build(4,4,3)
    sub,H,B,states,sizes,psi=model(4,4,8,source=True)
    index={s:j for j,s in enumerate(states)};rr=[];cc=[];vv=[]
    e2=np.array([sum(x*x for x in s[2]) for s in physical_states]);inside=e2<=8
    for i,s in enumerate(physical_states):
        if not inside[i]:continue
        key,sign,size=sub.canonical(s)
        if sign and key in index:rr.append(i);cc.append(index[key]);vv.append(sign/sqrt(size))
    J=coo_matrix((vv,(rr,cc)),shape=(len(physical_states),len(states))).tocsr()
    iso=float(np.linalg.norm((J.T@J-eye(len(states))).data))
    HP=diags(inside.astype(float))@HF@diags(inside.astype(float))
    hi=float(np.linalg.norm((HP@J-J@H).data))
    outside=diags((~inside).astype(float))@HF@J
    gram=float(np.linalg.norm((outside.T@outside-B.T@B).data))
    assert max(iso,hi,gram)<1e-10
    D,_=fixed_triangle(4,states);DF,_=fixed_triangle(4,physical_states)
    edge=HF.tocoo()
    transfer=(DF['protected'][edge.row]==0)&(DF['occupied_flux'][edge.col]>0)
    transfer_residual=float(np.linalg.norm(edge.data[transfer]))
    vacuum_two=diags(1-DF['protected'])@HF@HF@diags(DF['vacuum_flux'])
    vacuum_two_residual=float(np.linalg.norm(vacuum_two.data))
    assert transfer_residual<1e-12 and vacuum_two_residual<1e-12
    es=list(combinations(range(4),2));ring_rows=0;circulations=0;charged_rings=0
    for s in physical_states:
        pp,ff,hh,sign=masks(4,s);direct=directed_support(4,ff,sign)
        deg=[sum(bool(direct>>e&1) for e,(a,b) in enumerate(es) if v in (a,b)) for v in range(4)]
        if all(deg):
            # At Q=N, a spanning cyclic support must be the full simple ring.
            assert s[1]&15==0 and all(d==2 for d in deg) and direct==s[0]
            ring_rows+=1
        # Independent full bipartite incidence graph, including charged edges.
        A=np.zeros((10,10),bool)
        for e,(a,b) in enumerate(es):
            for leg,v in enumerate([a,b]):
                value=s[2][2*e+leg]
                if value>0:A[v,4+e]=True
                if value<0:A[4+e,v]=True
        reach=A|np.eye(10,dtype=bool)
        for v in range(10):reach|=reach[:,v,None]&reach[v,None,:]
        independent=sum(1<<e for e,(a,b) in enumerate(es) if reach[a,4+e] and reach[4+e,a])
        through,positive=incidence_masks(4,s);actual=directed_support(4,through,positive)
        assert independent==actual
        if actual:
            circulations+=1
            for length,mask,pattern,order in circuits(4):
                if through&mask==mask and positive&mask in (pattern,mask^pattern):
                    flow=order if positive&mask==pattern else tuple(reversed(order))
                    extracted=wilson_vector(4,tuple(reversed(flow)),s)
                    assert not np.any(gauss(4,extracted))
                    assert sum(map(abs,s[2]))-sum(map(abs,extracted[2]))==2*length
                    break
        deg=[sum(bool(actual>>e&1) for e,(a,b) in enumerate(es) if v in (a,b)) for v in range(4)]
        if all(deg):
            nv=(s[1]&15).bit_count()
            assert actual==s[0] and all(d==2 for d in deg)
            assert sum(x*x for x in s[2])>=8+3*nv
            charged_rings+=int(nv>0)
    wait=.7;lag=.4
    v=expm_multiply(-1j*wait*H,psi);vf=expm_multiply(-1j*wait*HP,J@psi)
    discrepancy=float(np.linalg.norm(J@v-vf))
    assert discrepancy<1e-10
    endpoint=expm_multiply(-1j*lag*H,v)
    full_endpoint=expm_multiply(-1j*lag*HP,vf)
    tests=0;maxdiff=0.
    for initial in ['hole','occupied_flux','vacuum_flux','electric_circulation']:
        z=expm_multiply(-1j*lag*H,D[initial]*v)
        zf=expm_multiply(-1j*lag*HP,DF[initial]*vf)
        for final in ['hole','occupied_flux','protected','electric_circulation',initial]:
            a=float(np.linalg.norm(D[final]*z)**2)
            b=float(np.linalg.norm(DF[final]*zf)**2)
            ca=np.vdot(endpoint,D[final]*z)
            cb=np.vdot(full_endpoint,DF[final]*zf)
            maxdiff=max(maxdiff,abs(a-b),abs(ca-cb));tests+=1
        # A full physical finite unitary dilation tests the actual error bounds.
        E=Evolution(H,B,2.,16)
        z0=D[initial]*v
        zp,_,ez=E(z0,lag)
        true=expm_multiply(-1j*lag*(HF-2*eye(HF.shape[0],format='csr')),J@z0)
        assert np.linalg.norm(true-J@zp)<=ez+1e-12
        for final in ['hole','occupied_flux','protected','electric_circulation']:
            a=float(np.linalg.norm(D[final]*zp)**2)
            b=float(np.linalg.norm(DF[final]*true)**2)
            assert max(0.,sqrt(a)-ez)**2-1e-12<=b<=min(1.,(sqrt(a)+ez)**2)+1e-12
    # Independently check compression and the stationary-reference correlator.
    eigen,V=np.linalg.eigh(H.toarray());weights=np.exp(-.7*(eigen-eigen[0]));weights/=weights.sum()
    X=V*np.sqrt(weights)[None,:];maxcorr=0.
    for initial in ['hole','occupied_flux','vacuum_flux','electric_circulation']:
        Y,right,res=compressed_filter(X,D[initial]);Yt=expm_multiply(-1j*lag*H,Y)
        Xt=X*np.exp(-1j*eigen*lag)[None,:]
        direct=expm_multiply(-1j*lag*H,D[initial][:,None]*X)
        for final in ['hole','occupied_flux','protected','electric_circulation',initial]:
            a=np.vdot(Xt,D[final][:,None]*direct)
            cross=X.conj().T@(D[final][:,None]*Yt)
            b=np.sum(np.exp(1j*eigen*lag)[:,None]*right.conj()*cross)
            maxcorr=max(maxcorr,float(abs(a-b)))
            assert abs(np.linalg.norm(D[final][:,None]*direct)**2-np.linalg.norm(D[final][:,None]*Yt)**2)<1e-11
    assert maxdiff<1e-10 and maxcorr<1e-10
    return dict(raw_physical_dimension=len(physical_states),source_dimension=len(states),
        isometry_error=iso,H_intertwiner_error=hi,outgoing_Gram_error=gram,
        independent_same_triangle_checks=tests,maximum_probability_or_correlator_difference=maxdiff,
        compression_correlation_error=maxcorr,Q_equals_N_spanning_occupied_ring_configurations_checked=ring_rows,
        occupied_flux_one_hop_geometry_loss_residual=transfer_residual,
        vacuum_flux_two_step_geometry_loss_residual=vacuum_two_residual,
        independent_incidence_circulation_configurations_checked=circulations,
        charged_spanning_circulation_configurations_checked=charged_rings)


def positive_duals():
    rng=np.random.default_rng(72417);checks=0;max_violation=0.
    for j in range(24):
        full=8;kept=4
        A=rng.normal(size=(full,full));HF=(A+A.T)/4+np.diag(np.arange(full)/3)
        H=coo_matrix(HF[:kept,:kept]).tocsr();B=coo_matrix(HF[kept:,:kept]).tocsr()
        x=rng.normal(size=kept)+1j*rng.normal(size=kept);x/=np.linalg.norm(x)
        t=.1+.025*j;E=Evolution(H,B,0.,8)
        v,_,e=E(x,t)
        O=np.diag(rng.integers(0,2,size=full).astype(float));O[0,0]=1
        dvec=np.diag(O)[:kept];_,_,d=E(dvec*v,t,back=True)
        embed=np.r_[v,np.zeros(full-kept)];true=expm(-1j*t*HF)@np.r_[x,np.zeros(full-kept)]
        p=float(np.vdot(v,dvec*v).real);actual=float(np.vdot(true,O@true).real)
        cross=abs(np.vdot(O@embed,true-embed))
        assert cross<=e*d+1e-11
        lower=max(0.,p-e*d)**2/p if p else 0.
        upper=min(1.,p+2*e*d+e*e)
        max_violation=max(max_violation,lower-actual,actual-upper)
        assert lower-1e-11<=actual<=upper+1e-11
        checks+=1
    return dict(dense_unitary_dilation_positive_dual_cases=checks,maximum_bound_violation=max_violation)


def run():
    start=time.monotonic();out=dict(status='CP12_FLUX_LOOP_CONTROLS_PASS',
        oriented_graphs_checked=oriented_graphs(),physical=physical(),positive_duals=positive_duals(),
        seconds=time.monotonic()-start)
    (P/'CONTROLS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True);return out


if __name__=='__main__':run()
