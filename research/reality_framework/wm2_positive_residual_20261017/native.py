"""Strict C++ kernels; NumPy owns every buffer and all physical data."""
import ctypes as C
from pathlib import Path
import subprocess
import hashlib
import numpy as np

BASE=Path(__file__).resolve().parent
SOURCE=BASE/'native.cpp'
FLAGS=['-O3','-std=c++17','-fPIC','-shared','-fno-fast-math','-ffp-contract=off']
BUILD=BASE/'build'
BUILD.mkdir(exist_ok=True)
LIB=BUILD/('kernels_'+hashlib.sha256(SOURCE.read_bytes()+str(FLAGS).encode()).hexdigest()[:12]+'.so')
if not LIB.exists():subprocess.run(['c++',*FLAGS,str(SOURCE),'-o',str(LIB)],check=True,capture_output=True)
lib=C.CDLL(str(LIB))
P=C.c_void_p;I=C.c_int;L=C.c_int64
lib.round_to_nearest.restype=I
lib.long_double_digits.restype=I
assert lib.round_to_nearest(),'Required IEEE round-to-nearest environment unavailable'
LONG_DIGITS=int(lib.long_double_digits())
lib.positive_norm.argtypes=[L,P];lib.positive_norm.restype=C.c_double
lib.action.argtypes=[L,I]+[P]*12
lib.fill_density.argtypes=[I,I]+[P]*10
lib.vacancy_weights.argtypes=[I,I]+[P]*7
lib.beam.argtypes=[L,L,I,L,P,P,P]


def ptr(a):return C.c_void_p(a.ctypes.data)
def raw_norm(x):
    x=np.ascontiguousarray(x,dtype=np.complex128)
    return float(lib.positive_norm(x.size,ptr(x)))


def action(layout,x,z):
    out=np.zeros_like(x)
    lib.action(layout.size,len(layout.words),*[ptr(a) for a in
        (layout.lengths,layout.offset,layout.parent,layout.coupling,layout.nocc,layout.annoff,
         layout.annmode,layout.annlower,layout.annfactor,z,x,out)])
    return out


def density(layout,x):
    d=np.zeros((layout.auxsize,layout.rank),complex)
    lib.fill_density(len(layout.words),layout.rank,*[ptr(a) for a in
        (layout.lengths,layout.offset,layout.auxoff,layout.nocc,layout.annoff,layout.annmode,
         layout.annlower,layout.annfactor,x,d)])
    return d.T@d.conj()


def vacancy(layout,x):
    out=np.zeros(layout.rank)
    lib.vacancy_weights(len(layout.words),layout.rank,*[ptr(a) for a in
        (layout.lengths,layout.offset,layout.nocc,layout.annoff,layout.annmode,x,out)])
    return out


def beam(x,rows,columns,indices,matrix):
    indices=np.ascontiguousarray(indices,dtype=np.int32)
    matrix=np.ascontiguousarray(matrix,dtype=np.complex128)
    lib.beam(rows,columns,matrix.shape[0],len(indices),ptr(indices),ptr(matrix),ptr(x))
