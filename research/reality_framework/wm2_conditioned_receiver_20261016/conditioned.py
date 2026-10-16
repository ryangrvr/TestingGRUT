"""History-conditioned receiver frames and full-action residual accounting."""
from math import ceil, exp, factorial, sqrt
from pathlib import Path
import importlib.util
import time
import numpy as np
from fock import space, Map

_path=Path(__file__).resolve().parent.parent/'wm2_history_formation_20261016/active.py'
_spec=importlib.util.spec_from_file_location('frozen_s5_active',_path)
S5=importlib.util.module_from_spec(_spec); _spec.loader.exec_module(S5)


class Limit(RuntimeError):
    pass


def norm2(x):
    return sum(float(np.vdot(v,v).real) for v in x.values())


class Operator:
    def __init__(self,words,frames,edges,contact,g,cache_cap=8_000_000,dense_limit=65536):
        self.words,self.frames,self.edges,self.contact,self.g=words,frames,edges,contact,g
        self.c={u:w.conj().T@contact for u,w in frames.items()}
        self.maps={}; self.dense_entries=0; self.cache_cap=cache_cap
        self.dense_limit=dense_limit; self.maximum_negative_residual_squared=0.0
        self.children={u:[] for u in words}; self.parent={}
        for p,v,a in edges:
            self.children[p].append((v,a)); self.parent[v]=(p,a)

    def transfer(self,dest,source,k,x):
        if k==0:
            return x.copy()
        if dest==source:
            return x.copy()
        key=(dest,source,k)
        if key not in self.maps:
            g=self.frames[dest].conj().T@self.frames[source]
            m=Map(g,k,self.dense_limit)
            if self.dense_entries+m.dense_entries>self.cache_cap:
                raise Limit('cached overlap-entry ceiling')
            self.maps[key]=m; self.dense_entries+=m.dense_entries
        return self.maps[key](x)

    def hint(self,x):
        y={u:np.zeros_like(v) for u,v in x.items()}
        for p,v,a in self.edges:
            k=len(p); rp=self.frames[p].shape[1]; rv=self.frames[v].shape[1]
            z=self.transfer(v,p,k,x[p])
            y[v] -= a*space(rv,k).create(self.c[v],z)
            z=space(rv,k+1).ann(self.c[v],x[v])
            y[p] -= a*self.transfer(p,v,k,z)
        return y

    def full_action_norm2(self,x,return_blocks=False):
        value=0.0; blocks={}
        for u in self.words:
            k=len(u); total=0.0
            if u in self.parent:
                p,a=self.parent[u]; rp=self.frames[p].shape[1]
                ap=(space(rp,k-1).ann(self.c[p],x[p]) if k>1 else np.zeros(0,complex))
                total += a*a*(float(np.vdot(x[p],x[p]).real)+float(np.vdot(ap,ap).real))
            z=[]
            for v,a in self.children[u]:
                rv=self.frames[v].shape[1]
                reduced=space(rv,k+1).ann(self.c[v],x[v])
                total += a*a*float(np.vdot(reduced,reduced).real)
                for w,b,other in z:
                    total += 2*a*b*float(np.vdot(other,self.transfer(w,v,k,reduced)).real)
                z.append((v,a,reduced))
                if k:
                    p,b=self.parent[u]
                    twice=space(rv,k).ann(self.c[v],reduced)
                    total += 2*a*b*float(np.vdot(x[p],self.transfer(p,v,k-1,twice)).real)
            blocks[u]=total; value+=total
        return (value,blocks) if return_blocks else value

    def residual(self,x,cx=None):
        if cx is None:
            cx=self.hint(x)
        squared=self.full_action_norm2(x)-norm2(cx)
        self.maximum_negative_residual_squared=max(self.maximum_negative_residual_squared,-squared)
        if squared < -2e-10*max(1,norm2(x)):
            raise RuntimeError('cross-history residual lost positivity beyond floating-point tolerance')
        return sqrt(max(0.0,squared))

    def polynomial(self,x,h,degree=8,residual_terms=3):
        term={u:v.copy() for u,v in x.items()}
        out={u:v.copy() for u,v in x.items()}
        leakage=0.0
        for j in range(degree):
            image=self.hint(term)
            if j<residual_terms:
                leakage += abs(h)/(j+1)*self.residual(term,image)
            term={u:(-1j*h/(j+1))*v for u,v in image.items()}
            for u in out:
                out[u] += term[u]
        z=self.g*abs(h)
        leakage += exp(z)*z**(residual_terms+1)/factorial(residual_terms+1)
        tail=exp(z)*z**(degree+1)/factorial(degree+1)
        return out,leakage,tail


class Evolution:
    def __init__(self,counts,values,modes,g=.9,rank_cap=10,state_cap=2_000_000,
                 cache_cap=8_000_000,dense_limit=65536,borrow=True):
        self.counts,self.values=tuple(counts),tuple(values)
        self.N=sum(counts)+1; self.modes=modes; self.g=g
        self.rank_cap,self.state_cap,self.cache_cap=rank_cap,state_cap,cache_cap
        self.dense_limit,self.borrow= dense_limit,borrow
        self.words=S5.type_words(counts)
        j=np.arange(1,modes+1)
        self.frequency=-2*np.cos(np.pi*j/(modes+1))
        self.contact=sqrt(2/(modes+1))*np.sin(np.pi*j/(modes+1))
        self.frames={u:self.contact[:,None].astype(complex) for u in self.words}
        self.x={u:np.zeros(1,complex) for u in self.words}; self.x[()][0]=1
        self.edges=[]; self.neighbors={u:[] for u in self.words}
        for u in self.words:
            for a,n in enumerate(counts):
                if u.count(a)<n:
                    v=u+(a,); coefficient=g/self.N*sqrt(n-u.count(a))
                    self.edges.append((u,v,coefficient))
                    self.neighbors[u].append((v,coefficient)); self.neighbors[v].append((u,coefficient))
        self.error=dict(orbital_projection=0.0,projected_emission=0.0,exponential_tail=0.0,
                        normalization=0.0,frame_change=0.0)
        self.projections=0; self.maximum_rank=1; self.maximum_size=len(self.words)
        self.maximum_frame_entries=modes*len(self.words)
        self.maximum_overlap_entries=0; self.maximum_gram_error=0.0
        self.maximum_negative_residual_squared=0.0

    def normalize(self):
        z=sqrt(norm2(self.x))
        self.error['normalization']+=abs(z-1)
        for u in self.words:
            self.x[u]/=z

    def free(self,t):
        phase=np.exp(-1j*self.frequency*t)
        for u in self.words:
            self.frames[u]*=phase[:,None]
            self.x[u]*=np.exp(-1j*t*sum(self.values[a] for a in u))

    def add(self,u,q):
        if not u:
            return False
        w=self.frames[u]; rank=w.shape[1]
        q=q-w@(w.conj().T@q); q=q-w@(w.conj().T@q)
        z=float(np.linalg.norm(q))
        if z<1e-13 or rank==self.modes:
            return False  # The complete physical residual still pays this exclusion.
        if rank>=self.rank_cap:
            raise Limit('word-conditioned rank ceiling')
        old,new=space(rank,len(u)),space(rank+1,len(u))
        size=sum(len(v) for v in self.x.values())-len(self.x[u])+new.size
        if size>self.state_cap:
            raise Limit('total coefficient ceiling')
        out=np.zeros(new.size,complex)
        out[[new.index[o] for o in old.occ]]=self.x[u]
        self.x[u]=out; self.frames[u]=np.column_stack((w,q/z))
        self.maximum_rank=max(self.maximum_rank,rank+1)
        self.maximum_size=max(self.maximum_size,size)
        self.maximum_frame_entries=max(self.maximum_frame_entries,sum(w.size for w in self.frames.values()))
        return True

    def compress(self,u,allowed):
        w=self.frames[u]; rank=w.shape[1]; k=len(u)
        if not k or rank<=1 or allowed<=0:
            return False
        f=space(rank,k)
        values,vectors=np.linalg.eigh(f.density(self.x[u]))
        if max(0.0,float(values[0]))>allowed*allowed/2:
            return False
        direction=vectors[:,0].copy(); last=rank-1; c=self.x[u]
        for i in range(last):
            a,b=direction[i],direction[last]; z=sqrt(abs(a)**2+abs(b)**2)
            if z==0:
                continue
            t=np.array([[b,-a],[a.conjugate(),b.conjugate()]])/z
            c=f.beam(c,i,last,t)
            w[:,[i,last]]=w[:,[i,last]]@t.conj().T
            direction[i],direction[last]=0,z
        self.x[u]=c; self.frames[u]=w
        self.normalize()
        discard=[i for i,o in enumerate(f.occ) if last in o]
        q=float(np.sum(abs(self.x[u][discard])**2))
        change=sqrt(2*q/(1+sqrt(max(0.0,1-q))))
        if change>allowed+1e-9:
            raise RuntimeError('history-conditioned projection exceeded its allowance')
        smaller=space(rank-1,k)
        self.x[u]=self.x[u][[f.index[o] for o in smaller.occ]]
        self.frames[u]=w[:,:-1]
        for v in self.words:
            self.x[v]/=sqrt(max(1e-300,1-q))
        self.error['orbital_projection']+=change; self.projections+=1
        return True

    def enrich(self,projection_budget,remaining):
        for u in self.words:
            if not u:
                continue
            unspent=max(0.0,projection_budget-self.error['orbital_projection'])
            if self.frames[u].shape[1]>=self.rank_cap:
                self.compress(u,unspent)
            self.add(u,self.contact)
        if not self.borrow:
            return
        densities={u:space(w.shape[1],len(u)).density(self.x[u])
                   for u,w in self.frames.items() if u}
        source_frames={u:w.copy() for u,w in self.frames.items()}
        for u in self.words:
            if not u or self.frames[u].shape[1]>=self.rank_cap:
                continue
            candidates=[]
            for v,a in self.neighbors[u]:
                if not v:
                    continue
                values,vectors=np.linalg.eigh(densities[v])
                candidates.append(a*source_frames[v]@(vectors*np.sqrt(np.maximum(values,0))))
            if not candidates:
                continue
            z=np.column_stack(candidates); w=self.frames[u]
            z-=w@(w.conj().T@z); z-=w@(w.conj().T@z)
            gram=z.conj().T@z
            values,vectors=np.linalg.eigh(gram)
            if values[-1]>1e-26:
                self.add(u,z@vectors[:,-1]/sqrt(values[-1]))

    def orthogonalize(self):
        for u,w in self.frames.items():
            if not u:
                continue
            gram=w.conj().T@w
            eta=float(np.linalg.norm(gram-np.eye(w.shape[1]),2))
            self.maximum_gram_error=max(self.maximum_gram_error,eta)
            if eta>1e-14:
                values,vectors=np.linalg.eigh(gram)
                if values[0]<=0:
                    raise Limit('frame lost linear independence')
                self.frames[u]=w@((vectors/np.sqrt(values))@vectors.conj().T)
                k=len(u)
                self.error['frame_change']+=k*eta*exp(k*eta)*sqrt(float(np.vdot(self.x[u],self.x[u]).real))
        self.normalize()

    def populations(self):
        return [sum(float(np.vdot(self.x[u],self.x[u]).real) for u in self.words if len(u)==k)
                for k in range(self.N)]

    def run(self,duration,strang_budget=.15,projection_budget=.20,residual_budget=.20,
            max_seconds=180,progress=None,steps_override=None,strict_step_residual=True):
        constant=S5.commutator_bounds(self.g)['strang_constant']
        steps=steps_override or max(1,ceil(sqrt(constant*duration**3/strang_budget)))
        dt=duration/steps; completed=0; start=time.monotonic()
        status,reason='COMPLETE',None
        rejected_residual=rejected_allowance=None
        self.free(dt/2)
        for step in range(steps):
            if time.monotonic()-start>max_seconds:
                status,reason='UNRESOLVED','per-instance wall-time budget'; break
            try:
                self.enrich(projection_budget,steps-step)
                op=Operator(self.words,self.frames,self.edges,self.contact,self.g,
                            self.cache_cap,self.dense_limit)
                out,residual,tail=op.polynomial(self.x,dt)
                self.maximum_overlap_entries=max(self.maximum_overlap_entries,op.dense_entries)
                self.maximum_negative_residual_squared=max(self.maximum_negative_residual_squared,
                                                          op.maximum_negative_residual_squared)
                available=max(0.0,residual_budget-self.error['projected_emission'])
                allowance=available/(steps-step) if strict_step_residual else available
                if residual>allowance+1e-12:
                    rejected_residual,rejected_allowance=residual,allowance
                    raise Limit('projected-emission residual step allowance' if strict_step_residual
                                else 'cumulative projected-emission residual budget')
                self.x=out
                self.error['projected_emission']+=residual; self.error['exponential_tail']+=tail
                self.normalize()
                for u in self.words:
                    unspent=max(0.0,projection_budget-self.error['orbital_projection'])
                    ceiling=self.frames[u].shape[1]>=self.rank_cap
                    allowed=unspent if ceiling else unspent/((steps-step)*max(1,len(self.words)-1))
                    self.compress(u,allowed)
                self.free(dt/2 if step==steps-1 else dt)
                if step%50==0 or step==steps-1:
                    self.orthogonalize()
                completed=step+1
                if progress and (completed%100==0 or completed==steps):
                    progress(dict(step=completed,steps=steps,coefficients=sum(len(v) for v in self.x.values()),
                                  max_rank=max(w.shape[1] for w in self.frames.values()),
                                  residual_error=self.error['projected_emission'],
                                  projection_error=self.error['orbital_projection'],seconds=time.monotonic()-start))
            except Limit as error:
                status,reason='UNRESOLVED',str(error); break
        if completed<steps:
            # The merged free steps otherwise leave a stopped state at the
            # next Strang midpoint. Return it to the last completed endpoint.
            self.free(-dt/2)
        error=dict(strang=constant*duration**3/steps**2,**self.error)
        e=sum(error.values()); population=self.populations()
        p=sum(population[max(2,ceil(self.N/2))-1:]); full=population[-1]
        def lower(value):
            return max(0.0,max(0.0,sqrt(min(1,value))-e)**2,
                       1-(sqrt(max(0.0,1-value))+e)**2) if status=='COMPLETE' else None
        low,full_low=lower(p),lower(full)
        gates=dict(strang=error['strang']<=strang_budget+1e-10,
                   projection=error['orbital_projection']<=projection_budget+1e-9,
                   residual=error['projected_emission']<=residual_budget+1e-10,
                   other=sum(v for k,v in error.items() if k not in
                             ('strang','orbital_projection','projected_emission'))<=.05,
                   total=e<=.60)
        return dict(status=status,stop_reason=reason,N=self.N,receiver_modes=self.modes,
             nonroot_type_counts=self.counts,detuning_values=self.values,g=self.g,
             observation_time=duration,steps=steps,completed_steps=completed,dt=dt,
             completed_steps_time=completed*dt,elapsed_seconds=time.monotonic()-start,
             populations=population,macro_probability_estimate=p,full_assembly_probability_estimate=full,
             error_components=error,analytic_state_error_account=e,formation_lower=low,
             full_assembly_lower=full_low,gate_pass=bool(status=='COMPLETE' and all(gates.values()) and low>=.5),
             numerical_gates=gates,residual_policy='local quota' if strict_step_residual else 'cumulative budget',
             rejected_candidate_residual=rejected_residual,rejected_candidate_allowance=rejected_allowance,
             maximum_word_rank=self.maximum_rank,maximum_coefficient_dimension=self.maximum_size,
             maximum_frame_entries=self.maximum_frame_entries,maximum_cached_overlap_entries=self.maximum_overlap_entries,
             maximum_gram_error=self.maximum_gram_error,
             maximum_negative_residual_squared=self.maximum_negative_residual_squared,
             history_count=len(self.words),final_ranks=[self.frames[u].shape[1] for u in self.words],
             projections=self.projections,exact_ordered_histories=True,all_reverse_transitions=True,
             conditional_frames=True,extra_photon_cutoff=False,residual_terms=3,polynomial_degree=8,
             physical_contact_exclusions_paid_in_full_residual=True,roundoff_certified=False,
             incomplete_population_notice=None if status=='COMPLETE' else 'Intermediate only; target-time probability unknown')
