"""All ordered histories, compact shared occupation maps, no path removal."""
from pathlib import Path
import importlib.util
from itertools import combinations_with_replacement
from math import comb,sqrt
import numpy as np
import native

_path=Path(__file__).resolve().parent.parent/'wm2_history_formation_20261016/active.py'
_spec=importlib.util.spec_from_file_location('s7_frozen_s5',_path)
S5=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(S5)


class Limit(RuntimeError):pass


class Layout:
    def __init__(self,counts,values,rank,g,state_cap=3_000_000):
        self.counts,self.values,self.rank,self.M=tuple(counts),tuple(values),rank,sum(counts)
        self.words=S5.type_words(counts)
        estimate=sum(comb(rank+len(u)-1,len(u)) if rank else int(not u) for u in self.words)
        if estimate>state_cap:raise Limit('coefficient ceiling')
        self.occ=[list(combinations_with_replacement(range(rank),k)) for k in range(self.M+1)]
        self.index=[{o:i for i,o in enumerate(occ)} for occ in self.occ]
        self.nocc=np.array([len(occ) for occ in self.occ],np.int64)
        self.nwords=[sum(len(u)==k for u in self.words) for k in range(self.M+1)]
        self.koff=np.r_[0,np.cumsum(self.nocc*np.array(self.nwords))].astype(np.int64)
        self.size=int(self.koff[-1])
        lookup={u:i for i,u in enumerate(self.words)}
        self.lengths=np.array([len(u) for u in self.words],np.int32)
        self.offset=np.zeros(len(self.words),np.int64)
        self.auxoff=np.zeros(len(self.words),np.int64)
        self.parent=np.zeros(len(self.words),np.int32)
        self.coupling=np.zeros(len(self.words))
        self.energy=np.empty(self.size)
        self.word_energy=np.array([sum(values[a] for a in u) for u in self.words])
        offset=aux=0
        for i,u in enumerate(self.words):
            k=len(u);self.offset[i]=offset;self.auxoff[i]=aux
            self.energy[offset:offset+self.nocc[k]]=self.word_energy[i]
            offset+=int(self.nocc[k])
            if k:
                p=u[:-1];a=u[-1];self.parent[i]=lookup[p]
                self.coupling[i]=g/(self.M+1)*sqrt(counts[a]-p.count(a))
                aux+=int(self.nocc[k-1])
        self.auxsize=aux
        self.annoff=np.zeros(self.M+1,np.int64)
        modes=[];lower=[];factor=[]
        for k in range(1,self.M+1):
            self.annoff[k]=len(modes)
            for o in self.occ[k]:
                for j,m in enumerate(o):
                    if j and o[j-1]==m:
                        modes.append(-1);lower.append(0);factor.append(0.0)
                    else:
                        modes.append(m);lower.append(self.index[k-1][o[:j]+o[j+1:]])
                        factor.append(sqrt(o.count(m)))
        self.annmode=np.array(modes,np.int32)
        self.annlower=np.array(lower,np.int32)
        self.annfactor=np.array(factor)
        self.rotation_groups={};self.drop_maps={}

    def block(self,c,k):return c[self.koff[k]:self.koff[k+1]].reshape(self.nwords[k],int(self.nocc[k]))
    def hint(self,c,z):return native.action(self,c,np.ascontiguousarray(z))
    def density(self,c):return native.density(self,c)
    def vacancy(self,c):return native.vacancy(self,c)
    def field_groups(self,k,i,j):return S5.Layout.field_groups(self,k,i,j)

    def beam(self,c,i,j,matrices):
        for k in range(1,self.M+1):
            for n,index in self.field_groups(k,i,j).items():
                native.beam(self.block(c,k),self.nwords[k],int(self.nocc[k]),index,matrices[n])

    def kept_map(self,new,removed):
        if removed not in self.drop_maps:
            self.drop_maps[removed]=[
                np.array([self.index[k][tuple(m+int(m>=removed) for m in o)] for o in new.occ[k]],np.int64)
                for k in range(self.M+1)]
        return self.drop_maps[removed]


def embed(c,old,new):
    out=np.zeros(new.size,complex)
    for k in range(old.M+1):
        mapping=np.array([new.index[k][o] for o in old.occ[k]],np.int64)
        new.block(out,k)[:,mapping]=old.block(c,k)
    return out


def restrict(c,old,new,removed):
    out=np.zeros(new.size,complex)
    for k,mapping in enumerate(old.kept_map(new,removed)):
        new.block(out,k)[:]=old.block(c,k)[:,mapping]
    return out
