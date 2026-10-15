"""Exact C++ acceleration of the inherited signed S_N canonicalizer.

The Python canonicalizer is retained as an independent reference. The C++
helper only performs integer permutations and Fock-parity counts.
"""
from __future__ import annotations
import ctypes as ct
from pathlib import Path
import hashlib
import subprocess
import sys
import tempfile
import numpy as np

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'cp12_ground_connectivity'))
from invariant_evolution import Invariant

_LIB=None


def library():
    global _LIB
    if _LIB is None:
        source=P/'symmetry.cpp'
        digest=hashlib.sha256(source.read_bytes()).hexdigest()[:16]
        target=Path(tempfile.gettempdir())/f'grut_signed_symmetry_{digest}.so'
        if not target.exists():
            subprocess.run(['g++','-std=c++17','-O3','-shared','-fPIC',str(source),'-o',str(target)],check=True,capture_output=True)
        lib=ct.CDLL(str(target))
        f=lib.canonical_state
        i32=ct.POINTER(ct.c_int32);u64=ct.POINTER(ct.c_uint64)
        f.argtypes=[ct.c_uint64,i32,ct.c_int,ct.c_int,i32,ct.c_int,
                    i32,i32,u64,i32,ct.c_int,u64,i32,i32,i32]
        f.restype=ct.c_int
        _LIB=lib
    return _LIB


class FastInvariant(Invariant):
    def initialize_engine(self):
        if hasattr(self,'engine_arrays'):return
        self.engine_arrays=(
            np.asarray([p['mode'] for p in self.group],dtype=np.int32),
            np.asarray([p['leg'] for p in self.group],dtype=np.int32),
            np.asarray([p['flipmask'] for p in self.group],dtype=np.uint64),
            np.asarray([p['background'] for p in self.group],dtype=np.int32))
        self.engine_ids={}

    def canonical(self,s):
        if s in self.cache:return self.cache[s]
        self.initialize_engine()
        g,ids=self.graph_info(s[0])
        if s[0] not in self.engine_ids:self.engine_ids[s[0]]=np.asarray(ids,dtype=np.int32)
        choices=self.engine_ids[s[0]]
        modes,legs,flips,backgrounds=self.engine_arrays
        E=(ct.c_int32*(2*self.M))(*s[2]);outE=(ct.c_int32*(2*self.M))()
        outbits=ct.c_uint64();sign=ct.c_int32();size=ct.c_int32()
        p32=lambda x:x.ctypes.data_as(ct.POINTER(ct.c_int32))
        answer=library().canonical_state(s[1],E,self.n+self.M,2*self.M,
            p32(choices),len(choices),p32(modes),p32(legs),
            flips.ctypes.data_as(ct.POINTER(ct.c_uint64)),p32(backgrounds),
            len(self.group),ct.byref(outbits),outE,ct.byref(sign),ct.byref(size))
        if answer:raise RuntimeError(('canonicalization failed',answer))
        result=(g,int(outbits.value),tuple(outE)),int(sign.value),int(size.value)
        if len(self.cache)>300000:self.cache.clear()
        self.cache[s]=result
        return result

    def __getstate__(self):
        return {k:v for k,v in self.__dict__.items() if k not in ['engine_arrays','engine_ids']}
