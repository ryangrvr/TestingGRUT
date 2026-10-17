"""S9: complete coherent MPS with monotone orthogonal projection batches.

The history register is one physical site. Its words are coherent labels,
not individually normalized trajectories. Every receiver occupation allowed
by the exact total charge M is available at every cut.
"""
from pathlib import Path
import sys
sys.path.append(str(Path(__file__).resolve().parent.parent/'wm2_schmidt_receiver_20261017'))
from math import ceil,sqrt
from fractions import Fraction as Q
from functools import lru_cache
import resource
import time
import numpy as np
from scipy.linalg.lapack import zgesdd
import bounds as E
import gates as G
from allocation import global_keep
from certificate import residual_enclosure


class Limit(RuntimeError):pass


@lru_cache(maxsize=2048)
def partitions(left,right,c1,c2,M):
    result=[]
    for middle in range(M+1):
        rows=[];cols=[];nr=nc=0
        for q,dl in left:
            for p,charge in enumerate(c1):
                if q+charge==middle:rows.append((q,p,nr,dl));nr+=dl
        for q,dr in right:
            for n,charge in enumerate(c2):
                if q-charge==middle:cols.append((n,nc,dr));nc+=dr
        if nr and nc:result.append((middle,tuple(rows),tuple(cols),nr,nc))
    return tuple(result)


class Receiver:
    def __init__(self,counts,values,modes,g=.9,rank_cap=96,total_rank_cap=96):
        self.counts,self.values=tuple(counts),tuple(values)
        self.M=sum(counts);self.N=self.M+1;self.physical_modes=modes;self.g=g
        self.rank_cap,self.total_rank_cap=rank_cap,total_rank_cap
        self.words,self.word_charge=G.words_and_charges(counts)
        self.word_index={u:i for i,u in enumerate(self.words)}
        self.error={k:0. for k in ['strang','schmidt_truncation','field_factorization',
            'gate_matrix','arithmetic','svd_reconstruction','normalization',
            'receiver_cut','time_grid','final_normalization']}
        self.gate_count=self.svd_count=self.truncations=0
        self.maximum_bond_rank=self.maximum_charge_rank=1;self.maximum_tensor_entries=1
        self.maximum_gauge_sum=self.maximum_svd_residual=0.
        self.maximum_active_receiver=0;self.stop_detail=None
        self.error['in_range_defect']=0.
        self.batch=None;self.completed_batches=0;self.sweep_checks=0
        self.maximum_in_range_defect=0.;self.maximum_projection_loss=0.
        self.maximum_residual_rounding=0.;self.maximum_batch_gates=0

    def initialize(self,R):
        self.R=R
        self.charges=[self.word_charge]+[np.arange(self.M+1)]*R
        self.bonds=[{0:1}]+[{self.M:1} for _ in range(R+1)]
        self.cores=[{(0,self.word_index[()]):np.ones((1,1),complex)}]
        self.cores += [{(self.M,0):np.ones((1,1),complex)} for _ in range(R)]
        self.eta=[0.]*(R+1);self.eta_sum=0.;self.center=0
        self.entries_current=R+1;self.maximum_tensor_entries=R+1

    def charge(self,name,value):
        if not np.isfinite(value) or value<0:raise Limit('nonfinite error account')
        self.error[name]=E.up(self.error[name]+value)

    def set_eta(self,i,value):
        self.eta_sum=max(0.,E.up(E.up(self.eta_sum-self.eta[i])+value))
        self.eta[i]=value;self.maximum_gauge_sum=max(self.maximum_gauge_sum,self.eta_sum)
        if self.eta_sum>=1e-6:raise Limit('summed canonical-isometry enclosure')

    def environment(self):return E.up(1/(1-self.eta_sum))

    def center_norm(self):
        lo=hi=0.
        for a in self.cores[self.center].values():
            x,y=E.norm2(a);lo=max(0.,E.down(lo+x));hi=E.up(hi+y)
        return lo,hi

    def normalize(self):
        lo,hi=self.center_norm();center=sqrt((lo+hi)/2)
        if not center>0:raise Limit('zero center norm')
        norm=E.up(sqrt(hi)*self.environment());scale=1/center
        change=E.up(abs(scale-1)*norm+E.gamma(32)*max(1,norm*abs(scale)))
        for key in self.cores[self.center]:self.cores[self.center][key]*=scale
        self.charge('normalization',change)

    def begin_batch(self,direction,maximum):
        if self.batch is not None:raise RuntimeError('unclosed projection batch')
        self.batch=dict(direction=direction,maximum=max(1,maximum),last=None,
            squares=0.,norm=0.,gates=0)

    def end_batch(self):
        if self.batch is None:raise RuntimeError('missing projection batch')
        self.maximum_batch_gates=max(self.maximum_batch_gates,self.batch['gates'])
        self.completed_batches+=1;self.batch=None

    def apply(self,i,gate,gate_error,direction,remaining_steps,max_gates):
        singleton=self.batch is None
        if singleton:self.begin_batch(direction,1)
        batch=self.batch
        if batch['direction']!=direction:raise RuntimeError('direction reversal inside batch')
        if batch['last'] is not None:
            wanted=batch['last']+(1 if direction=='right' else -1)
            if i!=wanted:raise RuntimeError('nonmonotone eligible sweep')
        batch['last']=i;batch['gates']+=1
        if batch['gates']>batch['maximum']:raise RuntimeError('batch gate ceiling')
        self.sweep_checks+=1
        if self.center not in [i,i+1]:raise RuntimeError('two-site center mismatch')
        left,right=self.cores[i],self.cores[i+1]
        c1,c2=self.charges[i],self.charges[i+1]
        theta={}
        right_index={}
        for (q,n),b in right.items():right_index.setdefault(q,[]).append((n,b))
        for (q,p),a in left.items():
            middle=q+int(c1[p])
            for n,b in right_index.get(middle,()):theta[q,p,n]=a@b
        # Apply each complete local physical charge block coherently. Outer
        # charge labels are environment labels, not selected histories.
        out={}
        for q,dl in self.bonds[i].items():
            for total,(pairs,u) in gate.items():
                dr=self.bonds[i+2].get(q+total)
                if dr is None:continue
                x=np.zeros((len(pairs),dl*dr),complex);found=False
                for row,(p,n) in enumerate(pairs):
                    a=theta.get((q,p,n))
                    if a is not None:x[row]=a.reshape(-1);found=True
                if not found:continue
                y=u@x
                for row,(p,n) in enumerate(pairs):
                    a=y[row].reshape(dl,dr)
                    out[q,p,n]=a
        decomposed={};rmax=max(self.bonds[i].values()) if self.bonds[i] else 1
        recon=[];full_gram=0.;center_squared=0.
        structure=partitions(tuple(self.bonds[i].items()),tuple(self.bonds[i+2].items()),
            tuple(map(int,c1)),tuple(map(int,c2)),self.M)
        for middle,rows,cols,nr,nc in structure:
            x=np.zeros((nr,nc),complex)
            for q,p,r,dl in rows:
                for n,c,dr in cols:
                    value=out.get((q,p,n))
                    if value is not None:x[r:r+dl,c:c+dr]=value
            if not x.any():continue
            u,s,v,info=zgesdd(x,compute_uv=1,full_matrices=0,overwrite_a=0)
            if info:raise Limit('LAPACK SVD did not converge')
            self.svd_count+=1
            if not np.isfinite(s).all():raise Limit('nonfinite SVD')
            eu,ev=E.gram(u),E.gram(v.conj().T)
            if max(eu,ev)>1e-8:raise Limit('SVD isometry enclosure')
            full_gram=max(full_gram,eu+ev)
            residual=x-(u*s)@v
            reconstruction=E.norm(residual)
            recon.append(reconstruction);self.maximum_svd_residual=max(self.maximum_svd_residual,reconstruction)
            center_squared+=float((s*s).sum())
            decomposed[middle]=dict(x=x,u=u,s=s,v=v,rows=rows,cols=cols,eu=eu,ev=ev)
        if not decomposed:raise Limit('zero two-site state')
        env=self.environment();unspent=max(0.,.35-self.error['schmidt_truncation'])
        if env*sqrt(center_squared)*(1+full_gram+E.gamma(8192))>1.01:
            raise Limit('intermediate physical-norm guard')
        allowance=unspent/max(1,8*remaining_steps)/sqrt(batch['maximum'])
        guard=E.up(env*(1+full_gram))
        keep=global_keep({q:d['s'] for q,d in decomposed.items()},self.total_rank_cap)
        # Rank choices are numerical proposals; the accepted tail is positively
        # measured below. First respect all hard caps, then discard only tiny
        # additional entries with a factor-two safety margin to the allocation.
        tail_est=sum(float((d['s'][keep[q]:]**2).sum()) for q,d in decomposed.items())
        candidates=sorted((float(d['s'][j]),q,j) for q,d in decomposed.items() for j in range(keep[q]))
        target=(allowance/max(1,guard)/2)**2
        for value,q,j in candidates:
            if j!=keep[q]-1:continue
            if tail_est+value*value>target:break
            tail_est+=value*value;keep[q]-=1
        # Stored factors define exact nested spans. Residuals are formed
        # directly, with no difference of squared norms.
        anew={};bnew={};bond={};eta=0.;losses=[];defects=[]
        env=self.environment();h=E.up(self.eta_sum/(1-self.eta_sum))
        for q,d in decomposed.items():
            r=keep[q];u=d['u'][:,:r];v=d['v'][:r];sv=d['s'][:r]
            if direction=='right':
                ua=u;vb=sv[:,None]*v;etaq=d['eu'];x=d['x'];K=ua;Z=vb
            else:
                ua=u*sv;vb=v;etaq=d['ev'];x=d['x'].conj().T;K=vb.conj().T;Z=ua.conj().T
            proof=residual_enclosure(x,K,Z,etaq,self.eta_sum,env)
            losses.append(proof['orthogonal_loss_upper'])
            self.maximum_residual_rounding=max(self.maximum_residual_rounding,proof['residual_rounding_upper'])
            if r:
                defects.append(proof['in_range_defect_upper'])
                bond[q]=r;eta=max(eta,etaq)
                for ql,p,start,width in d['rows']:anew[ql,p]=np.ascontiguousarray(ua[start:start+width])
                for n,start,width in d['cols']:bnew[q,n]=np.ascontiguousarray(vb[:,start:start+width])
        loss=E.scalar_norm(losses);defect=E.scalar_norm(defects)
        squares=E.up(batch['squares']+E.up(loss*loss));newnorm=E.up(sqrt(squares))
        increment=E.up(max(0.,newnorm-batch['norm']))
        if increment>unspent:
            self.stop_detail=dict(required_increment=increment,remaining=unspent,
                batch_norm=batch['norm'],candidate_batch_norm=newnorm,
                gate=self.gate_count,center=self.center)
            raise Limit('sweep projection cannot fit remaining truncation allowance')
        dimension=max(max(d['u'].shape+d['v'].shape) for d in decomposed.values())
        arithmetic=E.up(E.gamma(8192*(dimension+len(c1)*len(c2)+self.total_rank_cap+1))
                        *env*(1+sqrt(center_squared)))
        self.charge('arithmetic',arithmetic)
        self.charge('svd_reconstruction',E.up(env*E.scalar_norm(recon)))
        self.charge('gate_matrix',E.up(gate_error*env*(1+sqrt(center_squared))))
        self.charge('in_range_defect',defect)
        self.charge('schmidt_truncation',increment)
        self.maximum_in_range_defect=max(self.maximum_in_range_defect,defect)
        self.maximum_projection_loss=max(self.maximum_projection_loss,loss)
        batch['squares']=squares;batch['norm']=newnorm
        truncation=loss
        old_entries=sum(a.size for a in left.values())+sum(a.size for a in right.values())
        new_entries=sum(a.size for a in anew.values())+sum(a.size for a in bnew.values())
        self.entries_current+=new_entries-old_entries
        self.maximum_tensor_entries=max(self.maximum_tensor_entries,self.entries_current)
        self.cores[i],self.cores[i+1]=anew,bnew;self.bonds[i+1]=bond
        if direction=='right':self.set_eta(i,eta);self.set_eta(i+1,0.);self.center=i+1
        else:self.set_eta(i,0.);self.set_eta(i+1,eta);self.center=i
        self.gate_count+=1;self.truncations+=int(truncation>1e-14)
        self.maximum_charge_rank=max(self.maximum_charge_rank,max(bond.values()))
        self.maximum_bond_rank=max(self.maximum_bond_rank,sum(bond.values()))
        if singleton:self.end_batch()

    def vacuum_tail(self):
        # q=M at the right of this center forbids photons in the suffix.
        return set(self.bonds[self.center])=={self.M}

    def field_sweep(self,gate,error,remaining,max_gates):
        if self.center!=1:raise RuntimeError('field starts at first receiver')
        last=0;self.begin_batch('right',self.R-1)
        for i in range(1,self.R):
            if self.vacuum_tail():break
            self.apply(i,gate,error,'right',remaining,max_gates);last=i
            self.maximum_active_receiver=max(self.maximum_active_receiver,i+1)
        self.end_batch();self.begin_batch('left',max(1,last))
        for i in range(last,0,-1):self.apply(i,gate,error,'left',remaining,max_gates)
        self.end_batch()
        if self.center!=1:raise RuntimeError('field did not return to boundary')

    def state_norm_interval(self):
        lo,hi=self.center_norm();s=self.eta_sum
        return max(0.,E.down(sqrt(max(0,lo*(1-s))))),E.up(sqrt(hi/(1-s)))

    def event_interval(self,minimum):
        if self.center!=0:raise RuntimeError('event needs boundary center')
        eventlo=eventhi=totallo=totalhi=0.
        for (q,p),a in self.cores[0].items():
            lo,hi=E.norm2(a);totallo=max(0.,E.down(totallo+lo));totalhi=E.up(totalhi+hi)
            if len(self.words[p])+1>=minimum:
                eventlo=max(0.,E.down(eventlo+lo));eventhi=E.up(eventhi+hi)
        lower=max(0.,E.down(eventlo*(1-self.eta_sum)**2/totalhi))
        upper=min(1.,E.up(eventhi/(totallo*(1-self.eta_sum)**2)))
        return lower,upper

    def amplitude(self,word,occupation):
        if any(r>=self.R for r in occupation):return 0j
        q=0;p=self.word_index[tuple(word)];a=self.cores[0].get((0,p))
        if a is None:return 0j
        v=a[0];q=int(self.word_charge[p]);numbers={r:occupation.count(r) for r in set(occupation)}
        for r in range(self.R):
            n=numbers.get(r,0);b=self.cores[r+1].get((q,n))
            if b is None:return 0j
            v=v@b;q+=n
        return complex(v[0])

    def other_error(self):
        return sum(v for k,v in self.error.items() if k not in ['strang','schmidt_truncation','field_factorization'])

    def gates(self):
        for name,limit in [('strang',.20),('schmidt_truncation',.35),('field_factorization',.02)]:
            if self.error[name]>limit+1e-12:raise Limit(name+' error gate')
        if self.other_error()>.03:raise Limit('remaining arithmetic/error gate')

    def tensor_entries(self):return self.entries_current

    def run(self,T,max_seconds=900,started=None,progress=None,strang_budget=.20,cut=True):
        started=time.monotonic() if started is None else started
        setup=time.monotonic();configuration=G.prepare(self.counts,self.values,self.g,self.physical_modes,T,strang_budget,cut)
        setup_seconds=time.monotonic()-setup
        self.initialize(configuration['R']);n=configuration['steps'];dt=configuration['dt']
        max_gates=6*(self.R-1)+2;constant=configuration['commutator_certificate']['Strang_constant_upper']
        self.charge('receiver_cut',configuration['receiver_cut_error'])
        mismatch=abs(Q.from_float(dt)*n-Q(str(T)))
        hnorm=self.M*(2+max(map(abs,self.values)))+self.g
        self.charge('time_grid',float(mismatch)*hnorm*(1+E.gamma(64)))
        status,reason,done='COMPLETE',None,0
        try:
            for step in range(n):
                if time.monotonic()-started>max_seconds:raise Limit('whole-worker wall-time ceiling')
                if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024>2048:raise Limit('soft peak-RSS ceiling')
                before_norm=self.state_norm_interval()[1]
                if before_norm>1.01:raise Limit('field reference input-norm guard')
                self.apply(0,configuration['boundary'],configuration['boundary_error'],'right',n-step,max_gates)
                for field,error in configuration['fields']:self.field_sweep(field,error,n-step,max_gates)
                self.apply(0,configuration['boundary'],configuration['boundary_error'],'left',n-step,max_gates)
                self.normalize();done=step+1
                self.charge('strang',constant*abs(dt)**3*before_norm)
                self.charge('field_factorization',configuration['field_step_error']*1.01)
                self.gates();entries=self.tensor_entries();self.maximum_tensor_entries=max(self.maximum_tensor_entries,entries)
                if progress and (done%25==0 or done==n):progress(dict(step=done,steps=n,
                    time=done*dt,tensors=entries,maximum_bond=self.maximum_bond_rank,
                    maximum_charge_bond=self.maximum_charge_rank,gates=self.gate_count,
                    error=self.error.copy(),elapsed_seconds=time.monotonic()-started))
        except Limit as failure:status,reason='UNRESOLVED',str(failure)
        if self.center==0:
            lo,hi=self.state_norm_interval();self.charge('final_normalization',max(abs(lo-1),abs(hi-1)))
            events=dict(macro=self.event_interval(max(2,ceil(self.N/2))),full=self.event_interval(self.N))
        else:events=None
        error=E.up(sum(self.error.values())*(1+E.gamma(64)))
        if status=='COMPLETE':
            try:
                self.gates()
                if error>.60:raise Limit('total state-error envelope')
            except Limit as failure:status,reason='UNRESOLVED',str(failure)
        paid=({name:E.physical_probability(interval,error) for name,interval in events.items()}
            if status=='COMPLETE' and events else None)
        audit={k:v for k,v in configuration.items() if k not in ['boundary','fields']}
        return dict(status=status,stop_reason=reason,stop_detail=self.stop_detail,N=self.N,
            physical_receiver_modes=self.physical_modes,computational_receiver_sites=self.R,
            nonroot_type_counts=self.counts,detuning_values=self.values,g=self.g,
            observation_time=T,steps=n,completed_steps=done,dt=dt,
            completed_physical_steps_time=done*dt,whole_step_state=(self.center==0),
            maximum_bond_rank=self.maximum_bond_rank,maximum_charge_rank=self.maximum_charge_rank,
            maximum_tensor_entries=self.maximum_tensor_entries,maximum_active_receiver_site=self.maximum_active_receiver,
            two_site_gates=self.gate_count,SVDs=self.svd_count,truncations=self.truncations,
            maximum_isometry_sum_upper=self.maximum_gauge_sum,
            maximum_observed_SVD_residual=self.maximum_svd_residual,
            error_components=self.error,analytic_state_error_account=error,
            computed_normalized_event_intervals=events,physical_probability_intervals=paid,
            finite_instance_gate=bool(status=='COMPLETE' and error<=.60 and paid['macro'][0]>=.5),
            operator_audits=audit,operator_preprocessing_seconds=setup_seconds,
            elapsed_seconds=time.monotonic()-started,exact_history_retained=True,
            extra_photon_cutoff=False,emission_action_projection=False,
            arithmetic_status='roundoff-aware under declared IEEE/reduction model; not formal software certification',
            monotone_projection_certificate=True,completed_projection_batches=self.completed_batches,
            eligible_gate_checks=self.sweep_checks,maximum_batch_gates=self.maximum_batch_gates,
            maximum_in_range_defect=self.maximum_in_range_defect,
            maximum_projection_loss=self.maximum_projection_loss,
            maximum_residual_rounding=self.maximum_residual_rounding,
            pending_batch=(None if self.batch is None else self.batch),
            incomplete_population_notice=(None if status=='COMPLETE' else 'No frozen target-time formation probability'))
