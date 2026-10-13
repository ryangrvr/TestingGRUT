"""Exact signed S_N-invariant full-law evolution from unordered pregeometry.

All fixed-Q graphs and physical matter/rotor states at the declared cutoff.
No spatial pool or graph classifier appears in H or the initial amplitudes.
The initial state is the neutral equal superposition of all labelled graphs.
"""
from array import array
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import itertools
import json
import math
import sys
import time

import numpy as np
from scipy.sparse import coo_matrix,eye
from scipy.sparse.linalg import expm_multiply

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P))
assert hashlib.sha256((P/'frozen_model.py').read_bytes()).hexdigest()=='e896fce8a3e61948ff41891576b4c3ac028ad26785cafabb30bddef93c5adb9a'
from frozen_model import graph_basis,gauss,build
from exact_destination import transitions


class Invariant:
    def __init__(self,n,Q,L):
        self.n=n;self.Q=Q;self.L=L;self.es=list(itertools.combinations(range(n),2));self.M=len(self.es)
        self.group=[self.permutation(p) for p in itertools.permutations(range(n))]
        gid={tuple(p['vertex']):i for i,p in enumerate(self.group)}
        self.inverse=[]
        for p in self.group:
            inv=[0]*n
            for a,b in enumerate(p['vertex']):inv[b]=a
            self.inverse.append(gid[tuple(inv)])
        self.gcache={};self.cache={};self.maps={};self.width=n+self.M
        self.base=2*L+1;self.powers=[self.base**i for i in range(2*self.M)]
        self.baseline=L*sum(self.powers)
        masks=[sum(1 << e for e in es) for es in itertools.combinations(range(self.M),Q)]
        reps=sorted({self.graph_info(g)[0] for g in masks})
        orbits={};counts={}
        for g in reps:
            local={};count=0
            for s in graph_basis(n,g,L):
                count+=1;key,sign,size=self.canonical(s)
                if sign:orbits[key]=size
                local[self.pack(s)]=(self.pack(key),sign,size)
            self.maps[g]=local;counts[g]=count
        self.states=sorted(orbits)
        self.sizes=np.array([orbits[s] for s in self.states],np.int64)
        self.index={s:i for i,s in enumerate(self.states)}
        self.packed={}
        for j,s in enumerate(self.states):self.packed.setdefault(s[0],{})[self.pack(s)]=j
        self.full_dimension=sum(count*len(self.group)//len(self.graph_info(g)[1]) for g,count in counts.items())
        self.graph_representatives=len(reps)
        self.cache.clear()

    def pack(self,s):
        _,bits,E=s
        field=self.baseline+sum(x*self.powers[j] for j,x in enumerate(E) if x)
        return field << self.width | bits

    def permutation(self,p):
        p=list(p)
        emap=[self.es.index(tuple(sorted((p[a],p[b])))) for a,b in self.es]
        flip=[p[a]>p[b] for a,b in self.es]
        parity=sum(p[a]>p[b] for a in range(self.n) for b in range(a+1,self.n))%2
        return dict(vertex=p,edge=emap,mode=p+[self.n+e for e in emap],
                    leg=[2*emap[e]+(side ^ int(flip[e])) for e in range(self.M) for side in [0,1]],
                    flipmask=sum(1 << (self.n+e) for e,f in enumerate(flip) if f),
                    background=(-1)**(parity*(self.n-1)))

    def graph_info(self,g):
        if g in self.gcache:return self.gcache[g]
        occupied=[e for e in range(self.M) if g >> e & 1]
        minimum=min(sum(1 << p['edge'][e] for e in occupied) for p in self.group)
        occ=[e for e in range(self.M) if minimum >> e & 1]
        transport={}
        for i,p in enumerate(self.group):
            target=sum(1 << p['edge'][e] for e in occ)
            transport.setdefault(target,[]).append(self.inverse[i])
        for target,ids in transport.items():self.gcache[target]=(minimum,ids)
        return self.gcache[g]

    def transform(self,s,p,known=None):
        g,bits,E=s;newbits=0;parity=(bits & p['flipmask']).bit_count();remaining=bits
        while remaining:
            bit=remaining & -remaining;j=bit.bit_length()-1;target=p['mode'][j]
            parity+=(newbits >> target).bit_count();newbits|=1 << target;remaining^=bit
        if known is None:known=sum(1 << p['edge'][e] for e in range(self.M) if g >> e & 1)
        field=[0]*(2*self.M)
        for j,x in enumerate(E):
            if x:field[p['leg'][j]]=x
        return (known,newbits,tuple(field)),p['background']*(-1 if parity%2 else 1)

    def canonical(self,s):
        if hasattr(self,'packed'):
            g,ids=self.graph_info(s[0]);t=s;phase=1
            if s[0]!=g:t,phase=self.transform(s,self.group[ids[0]],g)
            packed,sign,size=self.maps[g][self.pack(t)]
            if not sign:return t,0,size
            return self.states[self.packed[g][packed]],phase*sign,size
        if s in self.cache:return self.cache[s]
        g,ids=self.graph_info(s[0]);best=None;sign=0;count=0;conflict=False
        for i in ids:
            t,phase=self.transform(s,self.group[i],g)
            if best is None or t<best:best=t;sign=phase;count=1;conflict=False
            elif t==best:count+=1;conflict=conflict or sign!=phase
        result=best,(0 if conflict else sign),len(self.group)//count
        if len(self.cache)>150000:self.cache.clear()
        self.cache[s]=result
        return result

    def matrix(self,raw=False):
        sc=math.lcm(8,10*self.M);rr=array('i');cc=array('i');vv=array('d')
        if raw:self.raw_columns=[[] for _ in self.states];self.raw_scale=sc
        for col,s in enumerate(self.states):
            assert not np.any(gauss(self.n,s))
            for target,c in transitions(self.n,s,sc):
                if any(abs(x)>self.L for x in target[2]):continue
                key,sign,size=self.canonical(target)
                if sign:
                    rr.append(self.index[key]);cc.append(col)
                    vv.append(c/sc*sign*math.sqrt(self.sizes[col]/size))
                    if raw:self.raw_columns[col].append((self.index[key],c*sign))
        H=coo_matrix((np.frombuffer(vv,dtype=np.float64),(np.frombuffer(rr,dtype=np.int32),np.frombuffer(cc,dtype=np.int32))),
                     shape=(len(self.states),)*2).tocsr()
        assert np.linalg.norm((H-H.T).data)<1e-10
        return H

    def initial(self):
        sea=sum(1 << (self.n+e) for e in range(self.M))
        total=math.comb(self.M,self.Q)
        psi=np.zeros(len(self.states),complex)
        count=0
        for i,(g,bits,E) in enumerate(self.states):
            if bits==sea and not any(E):
                psi[i]=math.sqrt(self.sizes[i]/total);count+=int(self.sizes[i])
        assert count==total and abs(np.linalg.norm(psi)-1)<1e-12
        return psi


def shape(n,g,es):
    adj=[set() for _ in range(n)]
    for e,(a,b) in enumerate(es):
        if g >> e & 1:adj[a].add(b);adj[b].add(a)
    reached={0};todo=[0]
    while todo:
        a=todo.pop()
        for b in adj[a]-reached:reached.add(b);todo.append(b)
    ds=[len(x) for x in adj];connected=len(reached)==n
    return [float(connected),float(connected and max(ds)<=2),sum(d*(d-1)/2 for d in ds),
            max(ds),sum(d==0 for d in ds)]


def intertwiner(n,Q,L):
    sub=Invariant(n,Q,L);HQ=sub.matrix();full,states,_,_=build(n,Q,L)
    rr=[];cc=[];vv=[]
    for i,s in enumerate(states):
        key,sign,size=sub.canonical(s)
        if sign:rr.append(i);cc.append(sub.index[key]);vv.append(sign/math.sqrt(size))
    T=coo_matrix((vv,(rr,cc)),shape=(len(states),len(sub.states))).tocsr()
    a=T.T@T-eye(T.shape[1],format='csr');b=full@T-T@HQ
    assert np.linalg.norm(a.data)<1e-10 and np.linalg.norm(b.data)<1e-9
    psi=sub.initial();physical=T@psi
    sea=sum(1 << (n+e) for e in range(sub.M));expected=np.array([1/math.sqrt(math.comb(sub.M,Q))
             if bits==sea and not any(E) else 0 for g,bits,E in states])
    assert np.linalg.norm(physical-expected)<1e-12
    rows=[]
    O=np.array([shape(n,s[0],sub.es) for s in sub.states],float).T
    Of=np.array([shape(n,s[0],sub.es) for s in states],float).T
    for t in [.4,2.,8.]:
        v=expm_multiply(-1j*t*HQ,psi);w=expm_multiply(-1j*t*full,physical)
        diff=float(np.linalg.norm(T@v-w))
        odiff=float(np.max(abs(O@abs(v)**2-Of@abs(w)**2)))
        assert diff<1e-10 and odiff<1e-10
        rows.append(dict(T=t,state_error=diff,shape_error=odiff))
    return dict(N=n,Q=Q,L=L,full_dimension=len(states),invariant_dimension=len(sub.states),
                isometry_residual=float(np.linalg.norm(a.data)),Hamiltonian_intertwiner_residual=float(np.linalg.norm(b.data)),
                dynamics_controls=rows)


def run(n,L,end):
    start=time.monotonic();Q=n-1
    sub=Invariant(n,Q,L)
    print(json.dumps(dict(stage='basis',N=n,L=L,full_dimension=sub.full_dimension,
                         invariant_dimension=len(sub.states),graph_representatives=sub.graph_representatives)),flush=True)
    H=sub.matrix();psi=sub.initial()
    energy=float(np.vdot(psi,H@psi).real)
    variance=float(np.linalg.norm(H@psi)**2-energy**2)
    assert abs(energy-Q/2)<1e-10 and abs(variance-2*Q)<1e-9
    flags=np.array([shape(n,s[0],sub.es) for s in sub.states],float).T
    labels=['connected','spanning_path','degree_wedges','maximum_degree','isolated_vertices']
    boundary=np.array([max(map(abs,s[2]))==L for s in sub.states],float)
    nv=np.array([(s[1]&((1 << n)-1)).bit_count() for s in sub.states],float)
    electric=np.array([sum(x*x for x in s[2]) for s in sub.states],float)
    initial=dict(zip(labels,(flags@abs(psi)**2).tolist()))
    assert abs(initial['connected']-n**(n-2)/math.comb(sub.M,Q))<1e-12
    assert abs(initial['spanning_path']-math.factorial(n)/2/math.comb(sub.M,Q))<1e-12
    rows=[];last=0.
    for T in [.1,.5,1.,2.,4.,8.,16.]:
        if T>end+1e-10:continue
        steps=max(1,int(np.ceil(T-last)));dt=(T-last)/steps
        for _ in range(steps):psi=expm_multiply(-1j*dt*H,psi)
        last=T;prob=abs(psi)**2
        norm_error=float(abs(prob.sum()-1));energy_error=float(abs(np.vdot(psi,H@psi).real-energy))
        assert norm_error<2e-9 and energy_error<2e-8
        row=dict(T=T,**dict(zip(labels,(flags@prob).tolist())),cutoff_boundary=float(boundary@prob),
                 charge_density=float(nv@prob/n),electric_square_per_vertex=float(electric@prob/n),
                 norm_error=norm_error,energy_error=energy_error)
        rows.append(row);print(json.dumps(dict(N=n,L=L,**row)),flush=True)
    result=dict(status='CP12_UNORDERED_INVARIANT_EVOLUTION_PASS',N=n,Q=Q,L=L,
                full_physical_dimension=sub.full_dimension,invariant_dimension=len(sub.states),
                graph_representatives=sub.graph_representatives,full_group_size=len(sub.group),
                initial_energy=energy,initial_variance=variance,initial_shapes=initial,rows=rows,
                seconds=time.monotonic()-start,
                scope='Complete S_N-invariant dynamics at finite electric cutoff, not the full spectrum or a continuum theorem.')
    (P/f'invariant_N{n}_L{L}.json').write_text(json.dumps(result,indent=2)+'\n')


def exact_shapes(n):
    Q=n-1;sub=Invariant(n,Q,3);sub.matrix(raw=True)
    sea=sum(1 << (n+e) for e in range(sub.M));total=math.comb(sub.M,Q)
    y0=[int(sub.sizes[i]) if bits==sea and not any(E) else 0 for i,(g,bits,E) in enumerate(sub.states)]
    ys=[y0]
    for k in range(6):
        y=[0]*len(sub.states)
        for col,a in enumerate(ys[-1]):
            if not a:continue
            # Through this order the cutoff cannot change a matrix action.
            assert max(map(abs,sub.states[col][2]))<3
            for row,c in sub.raw_columns[col]:y[row]+=c*a
        ys.append(y)
    coefficients={}
    graph_counts={}
    for i,(g,bits,E) in enumerate(sub.states):
        if y0[i]:graph_counts[g]=int(sub.sizes[i])
    for order in [2,4,6]:
        values={g:F(0) for g in graph_counts}
        for a in range(order+1):
            phase=(2*a-order)%4;sign=1 if phase==0 else (-1 if phase==2 else 0)
            if not sign:continue
            denom=total*sub.raw_scale**order*math.factorial(a)*math.factorial(order-a)
            for i,(g,_,_) in enumerate(sub.states):
                if ys[a][i] and ys[order-a][i]:
                    values[g]+=F(sign*ys[a][i]*ys[order-a][i],denom*int(sub.sizes[i]))
        assert sum(values.values())==0
        coefficients[order]=values
    assert all(v==0 for k in [2,4] for v in coefficients[k].values())
    mean=F(n*math.comb(n-1,2)*Q*(Q-1),sub.M*(sub.M-1))
    factor=None;rows=[]
    for g,count in sorted(graph_counts.items()):
        W=int(shape(n,g,sub.es)[2]);per_graph=coefficients[6][g]*total/count
        if W!=mean:
            slope=per_graph/(mean-W)
            if factor is None:factor=slope
            assert slope==factor
        else:assert per_graph==0
        rows.append(dict(graph_mask=g,labelled_graphs=count,degree_wedges=W,
                         probability_T6_coefficient=str(coefficients[6][g]),
                         relative_per_graph_T6_coefficient=str(per_graph)))
    assert factor>0
    return dict(N=n,Q=Q,L=3,invariant_dimension=len(sub.states),
                lower_probability_coefficients_vanish_through_order=4,
                uniform_mean_degree_wedges=str(mean),
                degree_bias_factor=str(factor),graph_class_coefficients=rows,
                scope='Exact integer/rational sixth-order coefficient for the unordered neutral preparation; no weak-hopping approximation.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=6);ap.add_argument('--cutoff',type=int,default=2)
    ap.add_argument('--end',type=float,default=16.);ap.add_argument('--check',action='store_true')
    ap.add_argument('--exact',action='store_true')
    ap.add_argument('--exact-n',type=int)
    a=ap.parse_args()
    if a.exact_n is not None:
        out=dict(status='CP12_UNORDERED_EXACT_DEGREE_BIAS_PASS',controls=[exact_shapes(a.exact_n)])
        (P/f'unordered_exact_N{a.exact_n}.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(dict(status=out['status'],N=a.exact_n,degree_bias_factor=out['controls'][0]['degree_bias_factor'])))
    elif a.exact:
        out=dict(status='CP12_UNORDERED_EXACT_DEGREE_BIAS_PASS',controls=[exact_shapes(n) for n in [4,5,6]])
        (P/'unordered_exact_results.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(out,indent=2))
    elif a.check:
        results=[intertwiner(4,3,3),intertwiner(5,4,2)]
        out=dict(status='CP12_FULL_PERMUTATION_INTERTWINER_PASS',controls=results)
        (P/'invariant_intertwiner_results.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps(out,indent=2))
    else:run(a.n,a.cutoff,a.end)
