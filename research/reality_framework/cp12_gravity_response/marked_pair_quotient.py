"""Exact signed permutation quotient of the frozen CP-12 Hamiltonian.

Keep states invariant under S_2 x S_(N-2), i.e. two marked vertices as an
unordered pair. This contains complete trivial, standard, and pair-pattern
representations, one subgroup-invariant vector per irrep copy. No fixed
graph or charge-pair approximation is made. Finite electric cutoff only.
"""
from pathlib import Path
from collections import Counter
from array import array
import argparse
import itertools
import json
import math
import sys

import numpy as np
from scipy.sparse import coo_matrix, eye
from scipy.sparse.linalg import eigsh, expm_multiply

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parents[1]/"constructive_physics/attempt12"))
from model import graph_basis, gauss
from exact_short_time import transitions


class Quotient:
    def __init__(self,n,q,L):
        self.n=n;self.q=q;self.L=L;self.edges=list(itertools.combinations(range(n),2));self.M=len(self.edges)
        self.group=[self.perm_data(list(front)+list(back)) for front in itertools.permutations(range(2))
                    for back in itertools.permutations(range(2,n))]
        group_index={tuple(p['vertex']):j for j,p in enumerate(self.group)}
        self.inverse_group=[]
        for p in self.group:
            inverse=[0]*n
            for i,j in enumerate(p['vertex']):inverse[j]=i
            self.inverse_group.append(group_index[tuple(inverse)])
        self.graph_cache={};self.state_cache={}
        self.local_maps={};self.mode_width=n+self.M
        self.field_base=2*L+1;self.field_powers=[self.field_base**j for j in range(2*self.M)]
        self.field_baseline=L*sum(self.field_powers)
        masks=[sum(1 << e for e in es) for es in itertools.combinations(range(self.M),q)]
        representatives=sorted({self.graph_info(mask)[0] for mask in masks})
        orbits={}; graph_dimensions={}
        for mask in representatives:
            count=0
            local={}
            for state in graph_basis(n,mask,L):
                count+=1;key,sign,size=self.canonical(state)
                if sign:orbits[key]=size
                local[self.pack(state)]=(self.pack(key),sign,size)
            graph_dimensions[mask]=count
            self.local_maps[mask]=local
        self.states=sorted(orbits);self.orbit_sizes=np.array([orbits[s] for s in self.states],int)
        self.index={s:i for i,s in enumerate(self.states)}
        self.packed_index={}
        for j,state in enumerate(self.states):self.packed_index.setdefault(state[0],{})[self.pack(state)]=j
        # Exact complete physical dimension from graph-orbit representatives.
        self.full_dimension=sum(count*len(self.group)//len(self.graph_info(mask)[1]) for mask,count in graph_dimensions.items())
        self.graph_representatives=len(representatives)

    def pack(self,state):
        _,bits,E=state
        electric=self.field_baseline+sum(value*self.field_powers[j] for j,value in enumerate(E) if value)
        return (electric << self.mode_width) | bits

    def perm_data(self,perm):
        emap=[self.edges.index(tuple(sorted((perm[i],perm[j])))) for i,j in self.edges]
        flip=[perm[i]>perm[j] for i,j in self.edges]
        parity=sum(perm[i]>perm[j] for i in range(self.n) for j in range(i+1,self.n))%2
        return dict(vertex=perm,edge=emap,flip=flip,mode=perm+[self.n+e for e in emap],
                    leg=[2*emap[e]+(side ^ int(flip[e])) for e in range(self.M) for side in [0,1]],
                    flip_fermion_mask=sum(1 << (self.n+e) for e,f in enumerate(flip) if f),
                    background_phase=(-1)**(parity*(self.n-1)))

    def graph_info(self,mask):
        if mask in self.graph_cache:return self.graph_cache[mask]
        occupied=[e for e in range(self.M) if mask >> e & 1]
        mapped=[sum(1 << perm['edge'][e] for e in occupied) for perm in self.group]
        minimum=min(mapped)
        # Fill the entire graph orbit in one pass. Inverse permutations are
        # the complete transporters from each member back to the minimum.
        # This is an exact cache, not a sampled isomorphism/canonicalization.
        occupied_minimum=[e for e in range(self.M) if minimum >> e & 1]
        transports={}
        for j,perm in enumerate(self.group):
            target=sum(1 << perm['edge'][e] for e in occupied_minimum)
            transports.setdefault(target,[]).append(self.inverse_group[j])
        for target,ids in transports.items():self.graph_cache[target]=(minimum,ids)
        return self.graph_cache[mask]

    def transform(self,state,perm,known_mask=None):
        mask,bits,E=state;mode=perm['mode']
        newbits=0;parity=(bits & perm['flip_fermion_mask']).bit_count();remaining=bits
        while remaining:
            bit=remaining & -remaining;j=bit.bit_length()-1;destination=mode[j]
            parity+=(newbits >> destination).bit_count();newbits |= 1 << destination;remaining ^= bit
        newmask=known_mask
        if newmask is None:
            newmask=0;remaining=mask
            while remaining:
                bit=remaining & -remaining;e=bit.bit_length()-1;newmask |= 1 << perm['edge'][e];remaining ^= bit
        newE=[0]*(2*self.M)
        for leg,value in enumerate(E):
            if value:newE[perm['leg'][leg]]=value
        return (newmask,newbits,tuple(newE)),perm['background_phase']*(-1 if parity%2 else 1)

    def canonical(self,state):
        if hasattr(self,'packed_index'):
            minimum,ids=self.graph_info(state[0]);phase=1;target=state
            if state[0]!=minimum:target,phase=self.transform(state,self.group[ids[0]],minimum)
            packed,sign,size=self.local_maps[minimum][self.pack(target)]
            if not sign:return target,0,size
            key=self.states[self.packed_index[minimum][packed]]
            return key,phase*sign,size
        if state in self.state_cache:return self.state_cache[state]
        minimum,ids=self.graph_info(state[0]);best=None;sign=None;count=0;conflict=False
        for j in ids:
            target,phase=self.transform(state,self.group[j],minimum)
            if best is None or target<best:best=target;sign=phase;count=1;conflict=False
            elif target==best:count+=1;conflict=conflict or phase!=sign
        result=(best,0 if conflict else sign,len(self.group)//count)
        # Reproducible recomputation trades CPU for bounded temporary memory;
        # this does not truncate states, operators, or graph orbits.
        if len(self.state_cache)>=150000:self.state_cache.clear()
        self.state_cache[state]=result
        return result

    def matrices(self):
        scale=math.lcm(8,10*self.M)
        rows=array('i');cols=array('i');data=array('d');integer_data=array('q')
        cr=array('i');cc=array('i');cv=array('d');integer_cv=array('q')
        transpositions=[]
        for i,j in self.edges:
            # Inside the marked-pair stabilizer every transposition is the
            # identity. Only cross-group swaps need an explicit action.
            if not ((i<2) != (j<2)):continue
            p=list(range(self.n));p[i],p[j]=p[j],p[i];transpositions.append(self.perm_data(p))
        for col,state in enumerate(self.states):
            assert not np.any(gauss(self.n,state))
            for target,coefficient in transitions(self.n,state,scale):
                if any(abs(x)>self.L for x in target[2]):continue
                key,sign,size=self.canonical(target)
                if not sign:continue
                row=self.index[key]
                rows.append(row);cols.append(col);data.append(coefficient/scale*sign*math.sqrt(self.orbit_sizes[col]/size))
                integer_data.append(coefficient*sign)
            for perm in transpositions:
                target,phase=self.transform(state,perm);key,sign,size=self.canonical(target)
                if not sign:continue
                cr.append(self.index[key]);cc.append(col);cv.append(phase*sign*math.sqrt(self.orbit_sizes[col]/size))
                integer_cv.append(phase*sign)
            diagonal_class=1+(self.n-2)*(self.n-3)//2
            cr.append(col);cc.append(col);cv.append(diagonal_class);integer_cv.append(diagonal_class)
        shape=(len(self.states),)*2
        rows=np.frombuffer(rows,dtype=np.int32);cols=np.frombuffer(cols,dtype=np.int32)
        cr=np.frombuffer(cr,dtype=np.int32);cc=np.frombuffer(cc,dtype=np.int32)
        H=coo_matrix((np.frombuffer(data,dtype=np.float64),(rows,cols)),shape=shape).tocsr()
        C=coo_matrix((np.frombuffer(cv,dtype=np.float64),(cr,cc)),shape=shape).tocsr()
        assert np.linalg.norm((H-H.T).data)<1e-10
        assert np.linalg.norm((C-C.T).data)<1e-10
        complete=shape[0]<=100000
        columns=np.arange(shape[0]) if complete else np.unique(np.linspace(0,shape[0]-1,64,dtype=int))
        comm=H @ C[:,columns]-C @ H[:,columns]
        assert np.linalg.norm(comm.data)<1e-8
        # Both physical Hermitian matrices are a common diagonal similarity
        # transform of exact integer matrices (H also has the declared scale).
        # Thus group closure and H/C commutation can be checked without roundoff.
        HI=coo_matrix((np.frombuffer(integer_data,dtype=np.int64),(rows,cols)),shape=shape).tocsr()
        CI=coo_matrix((np.frombuffer(integer_cv,dtype=np.int64),(cr,cc)),shape=shape).tocsr()
        # Full exact integer products at moderate sizes; deterministic column
        # checks at large sizes avoid forming unnecessarily dense intermediates.
        exact_comm=HI @ CI[:,columns]-CI @ HI[:,columns]
        exact_comm.eliminate_zeros();assert exact_comm.nnz==0
        casimirs=[self.n*(self.n-1)//2,self.n*(self.n-3)//2,(self.n-1)*(self.n-4)//2]
        I=eye(shape[0],format='csr',dtype=np.int64)
        polynomial=(CI-casimirs[0]*I) @ ((CI-casimirs[1]*I) @ ((CI-casimirs[2]*I)[:,columns]))
        polynomial.eliminate_zeros();assert polynomial.nnz==0
        self.integer_controls=dict(H_scale=scale,exact_H_C_commutator_nonzeros=exact_comm.nnz,
                                   exact_C_cubic_polynomial_nonzeros=polynomial.nnz,allowed_C_eigenvalues=casimirs)
        self.integer_controls.update(complete_integer_matrix_check=complete,checked_columns=len(columns))
        return H,C,float(np.linalg.norm(comm.data))


def run(n,q,L,k=32):
    quotient=Quotient(n,q,L)
    print(json.dumps(dict(stage="basis",N=n,Q=q,L=L,full_dimension=quotient.full_dimension,
                         quotient_dimension=len(quotient.states),graph_representatives=quotient.graph_representatives)),flush=True)
    H,C,comm=quotient.matrices()
    print(json.dumps(dict(stage="matrices",N=n,Q=q,L=L,H_nonzeros=H.nnz,
                          exact_integer_checks=quotient.integer_controls)),flush=True)
    count=min(k,len(quotient.states)-2)
    values,vectors=eigsh(H,k=count,which='SA',tol=2e-11,maxiter=15000,
                         v0=np.random.default_rng(412+n+L).normal(size=len(quotient.states)))
    order=np.argsort(values);values=values[order];vectors=vectors[:,order]
    # Resolve any energy degeneracy using the exactly commuting central class sum.
    blocks=[];first=0
    for j in range(1,count+1):
        if j==count or abs(values[j]-values[first])>1e-8:
            frame=vectors[:,first:j];_,rotation=np.linalg.eigh(frame.T @ C @ frame)
            vectors[:,first:j]=frame @ rotation;first=j
    casimir=np.sum(vectors*(C @ vectors),axis=0)
    expected=dict(trivial=n*(n-1)/2,standard=n*(n-3)/2,pair=(n-1)*(n-4)/2)
    rows=[]
    for j,(energy,c) in enumerate(zip(values,casimir)):
        kind=min(expected,key=lambda x:abs(expected[x]-c))
        cres=float(np.linalg.norm(C @ vectors[:,j]-expected[kind]*vectors[:,j]))
        hres=float(np.linalg.norm(H @ vectors[:,j]-energy*vectors[:,j]))
        assert cres<1e-7 and hres<1e-7
        rows.append(dict(energy=float(energy),gap=float(energy-values[0]),representation=kind,
                         casimir=float(c),casimir_residual=cres,eigen_residual=hres))
    assert rows[0]['representation']=='trivial'
    first_pair=next(row for row in rows[1:] if row['representation']=='pair')
    first_standard=next(row for row in rows[1:] if row['representation']=='standard')
    first_invariant=next(row for row in rows[1:] if row['representation']=='trivial')
    pair_dimension=n*(n-3)//2
    ground=vectors[:,0];edge_index=quotient.edges.index((0,1))
    observable=np.array([float(s[0] >> edge_index & 1) for s in quotient.states])
    centered=(observable-float(observable @ ground**2))*ground
    amplitudes=vectors.T @ centered
    for j,row in enumerate(rows):row['marked_bond_transition_weight']=float(amplitudes[j]**2)
    # A permutation-invariant matter/gauge probe has no nontrivial-irrep transitions.
    active=np.array([sum(x!=0 for x in s[2])/(2*q) for s in quotient.states],float)
    invariant_rhs=(active-float(active @ ground**2))*ground
    invariant_coeff=vectors.T @ invariant_rhs
    assert max(abs(invariant_coeff[j]) for j,row in enumerate(rows) if row['representation']!='trivial')<1e-8
    responses=[]
    for T in [.1,.4,1.,2.]:
        evolved=expm_multiply(-1j*T*H,np.column_stack([ground,active*ground]))
        R=-2*float(np.imag(np.vdot(evolved[:,0],active*evolved[:,1])))
        eps=1e-5
        plus=expm_multiply(-1j*T*H,np.exp(1j*eps*active)*ground)
        minus=expm_multiply(-1j*T*H,np.exp(-1j*eps*active)*ground)
        derivative=float((active @ abs(plus)**2-active @ abs(minus)**2)/(2*eps))
        assert abs(derivative-R)<2e-8
        responses.append(dict(T=T,invariant_activity_response=R,impulse_derivative=derivative,
                              impulse_error=abs(derivative-R)))
    out=dict(status="CP12_MARKED_PAIR_QUOTIENT_PASS",N=n,Q=q,L=L,full_physical_dimension=quotient.full_dimension,
         subgroup_invariant_dimension=len(quotient.states),subgroup_size=len(quotient.group),
         graph_representatives=quotient.graph_representatives,H_nonzeros=H.nnz,H_class_sum_commutator_norm=comm,
         exact_integer_symmetry_controls=quotient.integer_controls,
         invariant_ground_energy=float(values[0]),first_pair_shell=first_pair,full_pair_irrep_dimension=pair_dimension,
         first_standard_shell=first_standard,first_invariant_excitation=first_invariant,low_spectrum=rows,
         invariant_activity_controls=responses,
         scope="Exact symmetry-subspace compression with all graphs, fermions and regulated gauge states. N>=6 global ground outside this subspace is not searched. Relabelling-as-gauge remains an additional interpretation.")
    (P/f"marked_pair_N{n}_Q{q}_L{L}.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({k:out[k] for k in ['N','Q','L','full_physical_dimension','subgroup_invariant_dimension',
                        'invariant_ground_energy','first_pair_shell','full_pair_irrep_dimension','first_invariant_excitation']},indent=2),flush=True)
    return out


if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=5);ap.add_argument('--q',type=int)
    ap.add_argument('--cutoff',type=int,default=2);ap.add_argument('--eigenvectors',type=int,default=32)
    args=ap.parse_args();run(args.n,args.q if args.q is not None else args.n-1,args.cutoff,args.eigenvectors)
