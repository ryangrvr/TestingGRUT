"""All S_N representations of frozen CP12, without a singlet-ground assumption.

Young's orthogonal matrices implement ordinary permutation representation
theory. Stabilizer projectors select the multiplicity spaces. All their null
directions are removed before solving for the lowest energy.
Small full matrices independently validate the reduction, including every
eigenvalue with its representation multiplicity.
"""
from pathlib import Path
from array import array
import argparse
import itertools
import json
import math
import os
import pickle
import tempfile
import time
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import LinearOperator, eigsh
from invariant_evolution import Invariant
from exact_destination import transitions
from frozen_model import graph_basis, gauss, build
from ground import topology

P=Path(__file__).resolve().parent
CACHE=Path(os.environ.get('CP12_CACHE_DIR',str(Path(tempfile.gettempdir())/'cp12-representation-cache')))
ALGORITHM='stabilizer-multiplicity-v2'

def partitions(n,maximum=None):
    if not n:yield ();return
    for first in range(min(n,maximum or n),1-1,-1):
        for tail in partitions(n-first,first):yield (first,)+tail

def tableaux(shape):
    n=sum(shape);pos=[None]*n;used=[0]*len(shape);out=[]
    def recurse(k):
        if k==n:out.append(tuple(pos));return
        for row in range(len(shape)):
            col=used[row]
            if col==shape[row] or (row and col>=used[row-1]):continue
            pos[k]=(row,col);used[row]+=1;recurse(k+1);used[row]-=1
    recurse(0);return out

def representation(n,shape,group):
    ts=tableaux(shape);d=len(ts);ix={t:i for i,t in enumerate(ts)};gens=[]
    for j in range(n-1):
        D=np.zeros((d,d))
        for a,t in enumerate(ts):
            r=(t[j+1][1]-t[j+1][0])-(t[j][1]-t[j][0]);D[a,a]=1/r
            swap=list(t);swap[j],swap[j+1]=swap[j+1],swap[j];b=ix.get(tuple(swap))
            if b is not None:D[b,a]=math.sqrt(1-1/(r*r))
        assert np.linalg.norm(D@D-np.eye(d))<1e-12;gens.append(D)
    perms=[tuple(x['vertex']) for x in group];ids={p:i for i,p in enumerate(perms)}
    R=np.empty((len(perms),d,d));identity=ids[tuple(range(n))]
    R[identity]=np.eye(d);known={identity};todo=[identity]
    while todo:
        i=todo.pop();p=perms[i]
        for j,D in enumerate(gens):
            new=tuple(j+1 if v==j else j if v==j+1 else v for v in p);k=ids[new]
            if k in known:continue
            R[k]=D@R[i];known.add(k);todo.append(k)
    assert len(known)==math.factorial(n)
    # Independently check group products, beyond the generation tree.
    rng=np.random.default_rng(919)
    for _ in range(100):
        i,j=rng.integers(0,len(perms),2);k=ids[tuple(perms[i][perms[j][a]] for a in range(n))]
        assert np.linalg.norm(R[k]-R[i]@R[j])<2e-12
    return R

class AllOrbits(Invariant):
    def __init__(self,n,Q,L):
        self.n=n;self.Q=Q;self.L=L;self.es=list(itertools.combinations(range(n),2));self.M=len(self.es)
        self.group=[self.permutation(p) for p in itertools.permutations(range(n))]
        gid={tuple(p['vertex']):i for i,p in enumerate(self.group)};self.inverse=[]
        for p in self.group:
            inv=[0]*n
            for a,b in enumerate(p['vertex']):inv[b]=a
            self.inverse.append(gid[tuple(inv)])
        self.gcache={};self.cache={};self.maps={};self.width=n+self.M
        self.base=2*L+1;self.powers=[self.base**i for i in range(2*self.M)]
        self.baseline=L*sum(self.powers)
        masks=[sum(1<<e for e in es) for es in itertools.combinations(range(self.M),Q)]
        reps=sorted({self.graph_info(g)[0] for g in masks})
        orbits={}
        for g in reps:
            local={}
            for s in graph_basis(n,g,L):
                key,sign,size,pg=self.canonical(s);orbits[key]=size
                local[self.pack(s)]=(self.pack(key),sign,size,pg)
            self.maps[g]=local
        self.states=sorted(orbits);self.sizes=np.array([orbits[s] for s in self.states],np.int64)
        self.index={s:i for i,s in enumerate(self.states)};self.packed={}
        for j,s in enumerate(self.states):self.packed.setdefault(s[0],{})[self.pack(s)]=j
        self.full_dimension=int(self.sizes.sum());self.graph_representatives=len(reps);self.cache.clear()
        self.stabilizers=[];self.patterns={}
        for i,s in enumerate(self.states):
            stable=[]
            for j in self.graph_info(s[0])[1]:
                t,phase=self.transform(s,self.group[j],s[0])
                if t==s:stable.append((j,phase))
            assert len(stable)*self.sizes[i]==len(self.group)
            key=tuple(stable)
            self.stabilizers.append(key);self.patterns.setdefault(key,[]).append(i)
        self.patterns={key:np.array(ids,int) for key,ids in self.patterns.items()}

    def canonical(self,s):
        if hasattr(self,'packed'):
            g,ids=self.graph_info(s[0]);t=s;phase=1;transport=tuple(range(self.n))
            if s[0]!=g:
                p=self.group[ids[0]];t,phase=self.transform(s,p,g);transport=tuple(p['vertex'])
            packed,sign,size,pg=self.maps[g][self.pack(t)]
            if s[0]!=g:
                # Both maps are physical-to-canonical; compose them in that order.
                p=self.group[pg]['vertex'];target=tuple(p[transport[a]] for a in range(self.n))
                pg=self.group_id[target]
            return self.states[self.packed[g][packed]],phase*sign,size,pg
        if s in self.cache:return self.cache[s]
        g,ids=self.graph_info(s[0]);best=None;sign=0;count=0;pg=0
        for i in ids:
            t,phase=self.transform(s,self.group[i],g)
            if best is None or t<best:best=t;sign=phase;count=1;pg=i
            elif t==best:count+=1
        result=best,sign,len(self.group)//count,pg
        if len(self.cache)>150000:self.cache.clear()
        self.cache[s]=result;return result

    def transition_data(self):
        self.group_id={tuple(p['vertex']):i for i,p in enumerate(self.group)}
        sc=math.lcm(8,10*self.M);data={}
        for col,s in enumerate(self.states):
            assert not np.any(gauss(self.n,s))
            for t,c in transitions(self.n,s,sc):
                if any(abs(x)>self.L for x in t[2]):continue
                key,sign,size,pg=self.canonical(t)
                if pg not in data:data[pg]=(array('i'),array('i'),array('d'))
                rr,cc,vv=data[pg];rr.append(self.index[key]);cc.append(col)
                vv.append(c/sc*sign*math.sqrt(self.sizes[col]/size))
        return {pg:(np.array(rr,np.int32),np.array(cc,np.int32),np.array(vv,float))
                for pg,(rr,cc,vv) in data.items()}

def get_orbits(n,L):
    CACHE.mkdir(parents=True,exist_ok=True);path=CACHE/f'all_N{n}_Q{n}_L{L}.pickle'
    if path.exists():
        with path.open('rb') as f:
            class CacheReader(pickle.Unpickler):
                def find_class(self,module,name):
                    if name=='AllOrbits' and module=='__main__':return AllOrbits
                    return super().find_class(module,name)
            return CacheReader(f).load()
    sub=AllOrbits(n,n,L);sub.data=sub.transition_data()
    with path.open('wb') as f:pickle.dump(sub,f,protocol=5)
    print(json.dumps(dict(stage='all_orbits',N=n,L=L,orbits=len(sub.states),
        full_dimension=sub.full_dimension,transition_permutations=len(sub.data))),flush=True)
    return sub

class Reduced:
    def __init__(self,sub,shape):
        self.sub=sub;self.shape=shape;self.R=representation(sub.n,shape,sub.group)
        self.d=self.R.shape[1];self.patterns=[];self.multiplicity=0
        projections=[];ranks=np.zeros(len(sub.states),int)
        for key,ids in sub.patterns.items():
            proj=sum(sign*self.R[pg] for pg,sign in key)/len(key)
            assert np.linalg.norm(proj@proj-proj)<3e-12
            ev,vec=np.linalg.eigh(proj);V=vec[:,ev>.5];rank=V.shape[1]
            assert rank==int(round(np.trace(proj)));ranks[ids]=rank
            projections.append((ids,V));self.multiplicity+=len(ids)*rank
        offsets=np.r_[0,np.cumsum(ranks)]
        for ids,V in projections:
            indices=offsets[ids,None]+np.arange(V.shape[1])[None,:]
            self.patterns.append((ids,V,indices))
        self.shape_operator=(self.multiplicity,)*2
        self.H=LinearOperator(self.shape_operator,matvec=self.action,dtype=float)
        rng=np.random.default_rng(631)
        a=rng.standard_normal(self.multiplicity);b=rng.standard_normal(self.multiplicity)
        ha=self.action(a);hb=self.action(b)
        self.hermiticity_error=float(abs(a@hb-ha@b)/
                                     max(1,np.linalg.norm(a)*np.linalg.norm(b)))
        assert self.hermiticity_error<1e-11,self.hermiticity_error

    def embed(self,v):
        result=np.zeros((len(self.sub.states),self.d))
        for ids,V,indices in self.patterns:result[ids]=v[indices]@V.T
        return result

    def reduce(self,v):
        result=np.empty(self.multiplicity)
        for ids,V,indices in self.patterns:result[indices]=v[ids]@V
        return result

    def action(self,v):
        z=self.embed(v);out=np.zeros_like(z)
        for pg,(rr,cc,c) in self.sub.data.items():
            contribution=(z[cc]@self.R[pg].T)*c[:,None]
            np.add.at(out,rr,contribution)
        return self.reduce(out)

    def matrix(self):
        # Explicit projected matrix is useful only for independent small checks.
        n=self.H.shape[0];eye=np.eye(n)
        return np.column_stack([self.action(eye[:,i]) for i in range(n)])

def full_check(n,L):
    sub=get_orbits(n,L);HF,states,_,_=build(n,n,L);vals=np.linalg.eigvalsh(HF.toarray())
    collected=[];results=[]
    for shape in partitions(n):
        red=Reduced(sub,shape)
        if not red.multiplicity:continue
        H=red.matrix();eig=np.linalg.eigvalsh(H)
        assert len(eig)==red.multiplicity
        collected.extend(np.repeat(eig,red.d))
        results.append(dict(shape=shape,dimension=red.d,multiplicity=red.multiplicity,
                            ground=float(eig.min()),hermiticity_error=red.hermiticity_error))
    collected=np.sort(collected)
    assert len(collected)==len(vals)
    error=float(np.max(abs(collected-vals)));assert error<2e-9,error
    result=dict(status='CP12_ALL_REPRESENTATIONS_FULL_SPECTRUM_PASS',N=n,L=L,
                full_dimension=len(vals),max_full_spectrum_difference=error,irreps=results)
    (P/f'representation_check_N{n}_L{L}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)

def run(n,L,shapes=None,resume=True):
    start=time.monotonic();sub=get_orbits(n,L);flags=np.array([topology(n,s[0],sub.es) for s in sub.states],float).T
    names=['connected','components','largest_component_fraction','isolates','degree_defect',
           'degree_two','bridges_on_connected_graph','maximum_degree','wedges']
    boundary=np.array([max(map(abs,s[2]))==L for s in sub.states])
    output=P/f'representations_N{n}_L{L}.json';rows=[]
    if resume and output.exists():
        previous_result=json.loads(output.read_text())
        if previous_result.get('algorithm')==ALGORITHM:rows=previous_result['irreps']
    previous={tuple(x['shape']) for x in rows}
    for shape in (shapes or list(partitions(n))):
        if shape in previous:continue
        t0=time.monotonic();red=Reduced(sub,shape)
        if not red.multiplicity:continue
        rng=np.random.default_rng(763);initial=rng.standard_normal(red.multiplicity)
        E,v=eigsh(red.H,k=1,which='SA',tol=3e-10,maxiter=15000,v0=initial,ncv=30)
        ground=v[:,0];res=float(np.linalg.norm(red.H@ground-E[0]*ground));assert res<2e-8,res
        probability=np.sum(red.embed(ground)**2,axis=1)
        row=dict(shape=shape,representation_dimension=red.d,multiplicity_dimension=red.multiplicity,
            augmented_dimension=red.H.shape[0],lowest_energy=float(E[0]),residual=res,
            hermiticity_error=red.hermiticity_error,observables=dict(zip(names,(flags@probability).tolist())),
            boundary_probability=float(boundary@probability),seconds=time.monotonic()-t0)
        rows.append(row)
        result=dict(status='PARTIAL' if len(rows)<len(list(partitions(n))) else 'CP12_ALL_REPRESENTATIONS_GROUND_PASS',
                    algorithm=ALGORITHM,
                    N=n,Q=n,L=L,full_dimension=sub.full_dimension,orbits=len(sub.states),irreps=rows,
                    seconds=time.monotonic()-start)
        output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(row),flush=True)
    if len(rows)==len(list(partitions(n))):
        accounted=sum(r['representation_dimension']*r['multiplicity_dimension'] for r in rows)
        assert accounted==sub.full_dimension,(accounted,sub.full_dimension)
        result=json.loads(output.read_text())
        result['dimension_accounted_for']=accounted
        result['all_irrep_multiplicities_accounted_for']=True
        output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--N',type=int,required=True);a.add_argument('--L',type=int,required=True)
    a.add_argument('--check',action='store_true');a.add_argument('--shapes',default='')
    a.add_argument('--fresh',action='store_true',help='Recompute requested sectors rather than resume saved results.')
    a=a.parse_args()
    if a.check:full_check(a.N,a.L)
    else:run(a.N,a.L,[tuple(map(int,s.split(','))) for s in a.shapes.split(';')] if a.shapes else None,not a.fresh)
