"""Exact coherent action covariance and Krylov-support controls, no trajectories."""
from math import sqrt
from pathlib import Path
import json
import numpy as np
import fock as F
import conditioned as C
import controls


def annihilation_columns(frame,x,k):
    f=F.space(frame.shape[1],k)
    out=np.zeros((F.space(frame.shape[1],k-1).size,frame.shape[1]),complex)
    for mode,(src,dst,factor) in enumerate(f.annihilation):out[dst,mode]=factor*x[src]
    return out


def covariance_factor(components,k):
    """Return Z,M with gamma=Z M Z^dagger, including all coherent phases."""
    frames=[w for w,x in components]
    columns=[annihilation_columns(w,x,k) for w,x in components]
    offsets=np.r_[0,np.cumsum([w.shape[1] for w in frames])]
    z=np.column_stack(frames)
    matrix=np.zeros((int(offsets[-1]),int(offsets[-1])),complex)
    for a,wa in enumerate(frames):
        for b,wb in enumerate(frames):
            contraction=F.Map(wb.conj().T@wa,k-1)
            transformed=np.column_stack([contraction(columns[a][:,j]) for j in range(wa.shape[1])])
            block=(columns[b].conj().T@transformed).T
            matrix[offsets[a]:offsets[a+1],offsets[b]:offsets[b+1]]=block
    return z,matrix


def leading_missing_action(components,k,current):
    z,matrix=covariance_factor(components,k)
    z=z-current@(current.conj().T@z)
    q,r=np.linalg.qr(z,mode='reduced')
    small=r@matrix@r.conj().T
    small=(small+small.conj().T)/2
    values,vectors=np.linalg.eigh(small)
    return q@vectors[:,-1],float(values[-1])


def action_components(e,u):
    k=len(u); output=[]
    for p,v,a in e.edges:
        if v==u:
            w=e.frames[p]
            q=e.contact-w@(w.conj().T@e.contact)
            norm=float(np.linalg.norm(q))
            extended=np.column_stack((w,q/norm)) if norm>1e-13 else w
            vector=F.Map(extended.conj().T@w,k-1)(e.x[p])
            vector=F.space(extended.shape[1],k-1).create(extended.conj().T@e.contact,vector)
            output.append((extended,-a*vector))
        if p==u:
            w=e.frames[v]
            vector=F.space(w.shape[1],k+1).ann(w.conj().T@e.contact,e.x[v])
            output.append((w,-a*vector))
    return output


def project_complete(e,model,y):
    out=np.zeros_like(y)
    if not hasattr(e,'_control_lifts'):
        e._control_lifts={}
        for u in e.words:
            w=controls.dst(e.frames[u],type=1,norm='ortho',axis=0)
            # This contraction has separately been checked against the
            # independent permanent lift. The physical B matrix is S4's.
            e._control_lifts[u]=F.symmetric_power(w,len(u))
    for u in e.words:
        iso=e._control_lifts[u]
        indices=[model['index'][u,o] for o in F.space(e.modes,len(u)).occ]
        out[indices]=iso@(iso.conj().T@y[indices])
    return out


def run():
    rng=np.random.default_rng(623020)
    covariance=[]
    for k in (1,2,3,4):
        components=[]
        for rank in (2,3,2):
            w=controls.frame(rng,6,rank)
            x=rng.normal(size=F.space(rank,k).size)+1j*rng.normal(size=F.space(rank,k).size)
            components.append((w,x/5))
        z,m=covariance_factor(components,k)
        physical=sum(controls.S5.field_isometry(w,k)@x for w,x in components)
        exact=F.space(6,k).density(physical)
        error=float(np.linalg.norm(exact-z@m@z.conj().T))
        assert error<5e-11
        covariance.append(dict(photons=k,coherent_components=3,physical_density_error=error,
                               trace_identity_error=abs(float(np.trace(exact).real)-k*float(np.vdot(physical,physical).real))))
    w=controls.frame(rng,6,3)
    x=rng.normal(size=F.space(3,3).size)+1j*rng.normal(size=F.space(3,3).size)
    z,m=covariance_factor([(w,x),(w,-x)],3)
    cancellation=float(np.linalg.norm(z@m@z.conj().T))
    assert cancellation<5e-11
    # Explicit physical B action in every one-photon target.
    e=C.Evolution([1,1],[-.25,.25],6)
    for u in e.words:
        e.frames[u]=controls.frame(rng,6,1)
        e.x[u]=rng.normal(size=1)+1j*rng.normal(size=1)
    e.normalize()
    model=controls.M.vacuum(3,6,[0,-.25,.25],.9,classes=True)
    physical_action=model['H_emission']@controls.lift(e,model)
    action=[]
    for u in e.words:
        if len(u)!=1:continue
        components=action_components(e,u)
        z,m=covariance_factor(components,1)
        indices=[model['index'][u,o] for o in F.space(6,1).occ]
        field=controls.dst(physical_action[indices],type=1,norm='ortho')
        covariance_error=float(np.linalg.norm(z@m@z.conj().T-np.outer(field,field.conj())))
        direction,value=leading_missing_action(components,1,e.frames[u])
        old_missing=float(np.linalg.norm(field-e.frames[u]@(e.frames[u].conj().T@field)))
        enlarged=np.column_stack((e.frames[u],direction))
        after=float(np.linalg.norm(field-enlarged@(enlarged.conj().T@field)))
        assert covariance_error<5e-11 and after<5e-11
        action.append(dict(word=list(u),action_density_error=covariance_error,
                           missing_before=old_missing,missing_after_one_direction=after,largest_outside_eigenvalue=value))
    # Complete B^l Phi support, all append/pop paths, arbitrary complex packets.
    e=C.Evolution([2,1],[-.25,.25],16)
    for u in e.words:
        e.frames[u]=controls.frame(rng,16,1)
        e.x[u]=rng.normal(size=1)+1j*rng.normal(size=1)
    e.normalize()
    model=controls.M.vacuum(4,16,[0,-.25,-.25,.25],.9,classes=True)
    phi=controls.lift(e,model)
    adjacency={u:set() for u in e.words}
    for p,v,a in e.edges:adjacency[p].add(v);adjacency[v].add(p)
    krylov=[]
    for q in (1,2,3):
        enriched=C.Evolution([2,1],[-.25,.25],16,rank_cap=16)
        for u in e.words:
            reached={u}; frontier={u}
            for _ in range(q):
                frontier=set().union(*(adjacency[v] for v in frontier))-reached
                reached.update(frontier)
            basis=np.column_stack([e.contact[:,None]]+[e.frames[v] for v in e.words if v in reached])
            w=np.linalg.qr(basis,mode='reduced')[0]
            enriched.frames[u]=w
            enriched.x[u]=F.Map(w.conj().T@e.frames[u],len(u))(e.x[u])
        errors=[]; y=phi.copy()
        for order in range(q+1):
            if order:y=model['H_emission']@y
            errors.append(float(np.linalg.norm(y-project_complete(enriched,model,y))))
        assert max(errors)<5e-11
        krylov.append(dict(order=q,all_path_support_errors=errors,
                           maximum_frame_rank=max(w.shape[1] for w in enriched.frames.values()),
                           physical_modes=16,complete_dimension=len(phi)))
    out=dict(coherent_density_controls=covariance,exact_cancellation_density_norm=cancellation,
             one_photon_action_capture=action,exact_Krylov_enclosure=krylov,
             used_in_primary_trajectories=False,new_formation_result=False,
             external_validation=False,roundoff_certified=False)
    Path(__file__).with_name('ACTION_CONTROLS.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


if __name__=='__main__':print(json.dumps(run(),indent=2))
