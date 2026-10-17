"""S7 shared-frame action capture with a positive, roundoff-aware error ledger.

Arithmetic bounds use the explicitly stated IEEE reduction model, not a
formally verified BLAS/compiler. Eigenvectors propose a basis; actual positive
vacuum probabilities decide whether its projection is allowed.
"""
from collections import OrderedDict
from fractions import Fraction as Q
from math import ceil,sqrt
import resource
import time
import numpy as np
import enclosures as E
from layout import Layout,Limit,embed,restrict,S5
import native
import rotation


class Evolution:
    def __init__(self,counts,values,modes,g=.9,rank_cap=28,state_cap=3_000_000):
        self.counts,self.values=tuple(counts),tuple(values)
        self.N=sum(counts)+1;self.M=self.N-1
        self.modes,self.g,self.rank_cap,self.state_cap=modes,g,min(modes,rank_cap),state_cap
        self.contact,self.contact_error,self.frequency=E.receiver(modes)
        self.cache=OrderedDict();self.layout=self.get_layout(0)
        self.c=np.ones(1,complex);self.w=np.zeros((modes,0),complex)
        self.error={k:0.0 for k in ('strang','projection','emission_residual','polynomial',
            'arithmetic','basis','normalization','phase','time_grid','final_normalization')}
        self.maximum_rank=self.maximum_size=self.projections=self.enrichments=self.rotations=0
        self.maximum_gram=0.0;self.maximum_contact_defect=0.0;self.paid_skip_steps=0
        self.phase_cache={};self.metrics={'contact_checks':0,'natural_density_calls':0,'direct_drops':0}

    def get_layout(self,rank):
        if rank not in self.cache:
            x=Layout(self.counts,self.values,rank,self.g,self.state_cap)
            # Each edge is a rational g/N times two square roots. Bound the
            # relative discrepancy of every precomputed positive factor.
            relative=0.0
            for k in range(1,self.M+1):
                a,b=E.sqrt_bounds(float(k))
                v=sqrt(k);relative=max(relative,E.upper(max(v-a,b-v)/a))
            for capacity in self.counts:
                for n in range(1,capacity+1):
                    a,b=E.sqrt_bounds(float(n));lo,hi=E.rational_bounds(Q(str(self.g))/self.N)
                    lo=E.lower(lo*a);hi=E.upper(hi*b)
                    v=self.g/self.N*sqrt(n)
                    relative=max(relative,E.upper(max(v-lo,hi-v)/lo))
            x.relative_factor_error=E.add_up(2*relative,E.gamma(8))
            self.cache[rank]=x
            while len(self.cache)>3:self.cache.popitem(last=False)
        self.cache.move_to_end(rank)
        return self.cache[rank]

    def charge(self,key,value):
        assert value>=0 and np.isfinite(value),(key,value)
        self.error[key]=E.add_up(self.error[key],value)

    def other_error(self):return sum(v for k,v in self.error.items() if k not in ('strang','projection','emission_residual'))

    def gates(self):
        if self.error['strang']>.200000000001:raise Limit('Strang error gate')
        if self.error['projection']>.35:raise Limit('projection error gate')
        if self.error['emission_residual']>.02:raise Limit('emission residual gate')
        if self.other_error()>.03:raise Limit('remaining error gate')

    def gram(self):
        eta=E.gram_upper(self.w);self.maximum_gram=max(self.maximum_gram,eta)
        if eta>E.ETA:raise Limit('receiver Gram enclosure exceeds 1e-9')
        return eta

    def normalize(self):
        lo,hi=E.norm_interval(self.c)
        if not lo>0:raise Limit('zero or unenclosed coefficient norm')
        n=sqrt(native.raw_norm(self.c));scale=1/n
        # This pays the actual transformation of the current numerical vector,
        # without assuming its norm was exactly one on entry.
        change=E.mul_up(E.upper(abs(scale-1)),hi)
        self.c*=scale
        self.charge('normalization',E.add_up(change,E.mul_up(E.gamma(8),E.mul_up(abs(scale),hi))))

    def free(self,dt):
        rank=self.layout.rank
        key=float(dt)
        if key not in self.phase_cache:
            field,fe=E.phases(self.frequency,dt)
            # All declared detunings are exact binary fractions. The finite
            # words here have at most five entries, so their sums are exact.
            word,we=E.phases(E.point(self.layout.word_energy),dt)
            self.phase_cache[key]=(field,fe,word,we)
        field,fe,word,we=self.phase_cache[key]
        eta=self.gram();n=E.norm_up(self.c)
        for wi,k in enumerate(self.layout.lengths):
            start=self.layout.offset[wi];count=self.layout.nocc[k]
            self.c[start:start+count]*=word[wi]
        self.w*=field[:,None]
        # Relative to polar-isometric interpretations before/after the update.
        neweta=self.gram()
        frame_round=E.mul_up(E.gamma(16),1+eta)
        change=E.mul_up(2*self.M,E.add_up(E.add_up(eta,neweta),E.add_up(fe,frame_round)))
        self.charge('phase',E.mul_up(n,E.add_up(we,E.add_up(E.gamma(16),change))))

    def contact_data(self):
        eta=self.gram();z=self.w.conj().T@self.contact
        sigma,coordinate=E.contact_bounds(self.w,z,self.contact,self.contact_error,eta)
        self.metrics['contact_checks']+=1
        return z,sigma,coordinate,eta

    def include_contact(self,dt,remaining):
        z,sigma,coordinate,eta=self.contact_data()
        pref=E.upper(self.g*sqrt(self.M)/self.N)
        cost=E.mul_up(abs(dt)*pref,E.mul_up(sigma,E.norm_up(self.c)))
        unspent=max(0.0,.02-self.error['emission_residual'])
        if cost>unspent/remaining and self.layout.rank<self.modes:
            if self.layout.rank==self.rank_cap:
                self.compress(max(0.0,.35-self.error['projection']))
                z,sigma,coordinate,eta=self.contact_data()
            if self.layout.rank<self.rank_cap:
                q=self.contact-self.w@z
                q-=self.w@(self.w.conj().T@q)
                norm=sqrt(native.raw_norm(q))
                if norm>1e-12:
                    new=self.get_layout(self.layout.rank+1)
                    self.c=embed(self.c,self.layout,new);self.layout=new
                    self.w=np.column_stack((self.w,q/norm))
                    neweta=self.gram()
                    self.charge('basis',E.mul_up(2*self.M*E.norm_up(self.c),E.add_up(eta,neweta)))
                    self.maximum_rank=max(self.maximum_rank,new.rank)
                    self.maximum_size=max(self.maximum_size,new.size)
                    self.enrichments+=1
                    z,sigma,coordinate,eta=self.contact_data()
        self.maximum_contact_defect=max(self.maximum_contact_defect,sigma)
        cost=E.mul_up(abs(dt)*pref,E.mul_up(sigma,E.norm_up(self.c)))
        if cost>unspent:raise Limit('contact action cannot fit remaining residual allowance')
        if sigma>1e-10:self.paid_skip_steps+=1
        self.charge('emission_residual',cost)
        return z,coordinate

    def interact(self,dt,z,coordinate):
        norm=E.norm_up(self.c);znorm=E.norm_up(z)
        relative=self.layout.relative_factor_error
        gnum=E.upper(self.g*znorm*(1+relative))
        perturb=E.add_up(self.g*coordinate,self.g*znorm*relative)
        self.charge('arithmetic',E.mul_up(abs(dt)*norm,perturb))
        term=self.c.copy();out=term.copy()
        for j in range(1,4):
            term=(-1j*dt/j)*self.layout.hint(term,z)
            out+=term
        self.c=out
        self.charge('polynomial',E.upper(norm*(gnum*abs(dt))**4/24))
        # Explicit scalar complex products and bounded sparse row sums.
        ops=4096*(self.layout.rank*(len(self.counts)+1)+self.M+10)
        self.charge('arithmetic',E.upper(E.gamma(ops)*norm*(1+gnum*abs(dt))**4))
        self.normalize()

    def drop_costs(self):
        raw=self.layout.vacancy(self.c);nlo,nhi=E.norm2_bounds(self.c)
        result=[]
        for v in raw:
            _,b=E.positive_raw_bounds(float(v),self.c.size,native.LONG_DIGITS)
            q=min(1.0,E.upper(b/nlo))
            value=E.sqrt_up(E.upper(2*q/E.lower(1+E.sqrt_down(max(0.0,E.lower(1-q))))))
            result.append(value)
        return result

    def drop(self,index,cost):
        old=self.layout;new=self.get_layout(old.rank-1)
        eta=self.gram();lo,hi=E.norm_interval(self.c)
        self.c=restrict(self.c,old,new,index)
        self.w=np.delete(self.w,index,axis=1);self.layout=new
        # Projection formula is for a normalized vector. Norm deviation and
        # division roundoff are explicitly paid here, not hidden in the formula.
        self.charge('projection',cost)
        self.charge('normalization',max(abs(lo-1),abs(hi-1)))
        n=sqrt(native.raw_norm(self.c));self.c/=n
        self.charge('normalization',E.mul_up(E.gamma(12),max(1,hi/n)))
        neweta=self.gram()
        self.charge('basis',E.mul_up(2*self.M,E.add_up(eta,neweta)))
        self.projections+=1

    def compress(self,allowed):
        rank=self.layout.rank
        if rank<=1 or allowed<=0:return False
        costs=self.drop_costs();j=int(np.argmin(costs))
        if costs[j]<=allowed:
            self.drop(j,costs[j]);self.metrics['direct_drops']+=1;return True
        rho=self.layout.density(self.c);self.metrics['natural_density_calls']+=1
        values,vectors=np.linalg.eigh(rho)
        # A failed proposal is safe; its eigenvalue is never an error bound.
        if float(values[0])>allowed*allowed:return False
        u=vectors[:,0].copy();last=rank-1
        for i in range(last):
            a,b=u[i],u[last];n=sqrt(abs(a)**2+abs(b)**2)
            if n<1e-15:continue
            t=np.array([[b,-a],[a.conjugate(),b.conjugate()]])/n
            matrices,delta,rows,te=rotation.fast_matrices(t,self.M)
            norm=E.norm_up(self.c);eta=self.gram()
            self.layout.beam(self.c,i,last,matrices)
            self.w[:,[i,last]]=self.w[:,[i,last]]@t.conj().T
            u[i],u[last]=0,n
            neweta=self.gram()
            # Gamma(T) vs exact unitary polar(T), native matrix products, and
            # physical frame polar interpretations. No cancellation subtraction.
            fee=E.add_up(delta,E.mul_up(E.gamma(512),rows))
            fee=E.add_up(fee,2*self.M*(te+eta+neweta+E.gamma(64)))
            self.charge('basis',E.mul_up(norm,fee));self.rotations+=1
        self.normalize()
        cost=self.drop_costs()[last]
        if cost<=allowed:self.drop(last,cost);return True
        return False

    def populations(self):return [native.raw_norm(self.layout.block(self.c,k)) for k in range(self.N)]

    def events(self):
        k=max(2,ceil(self.N/2))-1
        total=native.raw_norm(self.c)
        macro=native.raw_norm(self.c[self.layout.koff[k]:])
        full=native.raw_norm(self.layout.block(self.c,self.M))
        return {name:E.probability_interval(raw,total,self.c.size) for name,raw in [('macro',macro),('full',full)]}

    def run(self,duration,strang_budget=.20,projection_budget=.35,max_seconds=600,
            started=None,progress=None,steps_override=None,rss_limit=2048):
        started=time.monotonic() if started is None else started
        # Bound declared g and sqrt(2) outward, and round the scalar expression
        # conservatively. Scalar construction is negligible next to .20.
        constant=E.upper(S5.commutator_bounds(self.g)['strang_constant']*(1+E.gamma(128)))
        steps=steps_override or max(1,ceil(sqrt(constant*duration**3/strang_budget)))
        while not steps_override and E.upper(constant*duration**3/steps**2)>strang_budget:steps+=1
        dt=duration/steps
        mismatch=abs(Q.from_float(dt)*steps-Q(str(duration)))
        self.charge('time_grid',E.rational_bounds(mismatch*Q(str(self.M*(max(map(abs,self.values))+2)+self.g)))[1])
        status,reason,done='COMPLETE',None,0
        self.free(dt/2)
        try:
            for step in range(steps):
                if time.monotonic()-started>max_seconds:raise Limit('whole-worker wall-time ceiling')
                if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024>rss_limit:raise Limit('soft peak-RSS ceiling')
                z,coordinate=self.include_contact(dt,steps-step)
                self.interact(dt,z,coordinate)
                unspent=max(0.0,projection_budget-self.error['projection'])
                allowed=unspent if self.layout.rank==self.rank_cap else unspent/(steps-step)
                self.compress(allowed)
                self.free(dt/2 if step==steps-1 else dt)
                done=step+1;self.gates()
                if progress and (done%200==0 or done==steps):
                    progress(dict(step=done,steps=steps,rank=self.layout.rank,coefficients=self.layout.size,
                        error=self.error.copy(),elapsed_seconds=time.monotonic()-started))
        except Limit as error:
            status,reason='UNRESOLVED',str(error)
        if status=='UNRESOLVED':
            # Each completed nonfinal step has advanced A by an extra half step.
            # Undo it for the reported intermediate state, paying the operation.
            self.free(-dt/2)
        self.charge('strang',E.upper(constant*done*dt**3))
        lo,hi=E.norm_interval(self.c)
        self.charge('final_normalization',max(abs(lo-1),abs(hi-1)))
        error=E.upper(sum(self.error.values())*(1+E.gamma(len(self.error)+4)))
        event=self.events();paid={key:E.physical_probability(value,error) for key,value in event.items()}
        if status=='COMPLETE':
            try:self.gates()
            except Limit as failure:status,reason='UNRESOLVED',str(failure)
        return dict(status=status,stop_reason=reason,N=self.N,receiver_modes=self.modes,
            nonroot_type_counts=self.counts,detuning_values=self.values,g=self.g,
            observation_time=duration,completed_steps=done,steps=steps,dt=dt,
            completed_physical_steps_time=done*dt,exponential_degree=3,
            maximum_active_rank=self.maximum_rank,final_active_rank=self.layout.rank,
            maximum_coefficient_dimension=self.maximum_size,projections=self.projections,
            enrichments=self.enrichments,rotations=self.rotations,metrics=self.metrics,
            positive_contact_residual_paid_steps=self.paid_skip_steps,
            maximum_contact_defect_upper=self.maximum_contact_defect,
            maximum_gram_enclosure=self.maximum_gram,populations=self.populations(),
            computed_normalized_event_intervals=event,
            analytic_state_error_account=error,error_components=self.error,
            physical_probability_intervals=(paid if status=='COMPLETE' else None),
            finite_instance_gate=bool(status=='COMPLETE' and error<=.60 and paid['macro'][0]>=.5),
            elapsed_seconds=time.monotonic()-started,
            arithmetic_status='roundoff-aware under declared IEEE/reduction model; not formal software certification',
            exact_history_retained=True,extra_photon_cutoff=False,
            incomplete_population_notice=(None if status=='COMPLETE' else 'Intermediate state only; no target-time formation bound'))
