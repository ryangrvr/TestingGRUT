"""Exact algebra for a NEW neutral-probe transport design, not frozen CP-12.

The local through-field gate preserves geometry formation in probe vacuum,
but an exact charged-boundary control exposes unsuppressed activation.
"""
from collections import defaultdict
from fractions import Fraction as F
from math import lcm
import json
from common import P,transitions,gauss,combinations


def gate(s,e):return bool(s[0]>>e&1 and s[2][2*e]*s[2][2*e+1]<0)


def geom(n,s):
    scale=lcm(8,10*(n*(n-1)//2))
    return [(t,F(c,scale)) for t,c in transitions(n,s,scale)]


def probe(n,s,pos):
    if pos is None:return []
    out=[]
    for e,(i,j) in enumerate(combinations(range(n),2)):
        if not gate(s,e):continue
        if pos==i:out.append(((s,j),F(1)))
        if pos==j:out.append(((s,i),F(1)))
    return out


def apply(n,v,geometry=True,transport=True):
    out=defaultdict(F)
    for (s,pos),amplitude in v.items():
        if geometry:
            for t,c in geom(n,s):gauss(n,t);out[t,pos]+=c*amplitude
        if transport:
            for t,c in probe(n,s,pos):out[t]+=c*amplitude
    return {s:a for s,a in out.items() if a}


def control(n,flux):
    es=list(combinations(range(n),2));M=len(es);j=n//2;g=0;E=[0]*(2*M)
    for u in range(n):
        v=(u+1)%n;e=es.index(tuple(sorted((u,v))));g|=1<<e
        if u<v:E[2*e:2*e+2]=[flux,-flux]
        else:E[2*e:2*e+2]=[-flux,flux]
    f=es.index((0,j));g|=1<<f;E[2*f]=1
    e=es.index((j,j+1));E[2*e]+=1
    bits=sum(1<<(n+k) for k in range(M))^(1<<(n+f))^(1<<(n+e))
    bits|=(1<<0)|(1<<j);s=(g,bits,tuple(E));gauss(n,s);assert not gate(s,f)
    # A single original hop a_j -> b_f creates through-field (1,-1).
    targetE=E.copy();targetE[2*f+1]-=1
    target=(g,bits^(1<<j)^(1<<(n+f)),tuple(targetE));gauss(n,target);assert gate(target,f)
    initial={(s,0):F(1)};first=apply(n,initial);second=apply(n,first)
    amplitude=second.get((target,j),F(0));assert abs(amplitude)==1
    # The target has a changed geometry state, so the ordinary two-edge
    # ring path cannot contribute to this particular amplitude at order 2.
    geometry_then_probe=apply(n,apply(n,initial,transport=False),geometry=False)
    probe_then_geometry=apply(n,apply(n,initial,geometry=False),transport=False)
    comm=geometry_then_probe.get((target,j),0)-probe_then_geometry.get((target,j),0)
    assert abs(comm)==1
    diagonal=lambda ss:next(c for t,c in geom(n,ss) if t==ss)
    detuning=diagonal(target)-diagonal(s);assert detuning==F(-15,8)
    vacuum=apply(n,{(s,None):F(1)})
    original={(t,None):c for t,c in geom(n,s)};assert vacuum==original
    reverse=apply(n,{(target,j):F(1)}).get((s,0),0)
    assert reverse==0  # two-step activation, not a direct first-order edge
    for state,pos in second:gauss(n,state);assert pos is None or 0<=pos<n
    return dict(N=n,ring_flux=flux,source=0,receiver=j,ring_distance=j,
        initial_shortcut_gate=False,one_geometry_hop_activates_gate=True,
        Hnew_squared_target_amplitude=str(amplitude),
        target_transfer_probability_T4_coefficient='J^2/4',
        geometry_probe_commutator_target=str(comm),activation_diagonal_detuning=str(detuning),
        probe_vacuum_preserves_original_geometry_law=True,second_power_states=len(second))


if __name__=='__main__':
    rows=[control(n,l) for n in [4,6,8,12] for l in [1,3,8]]
    data=dict(status='PASS',law_status='NEW PROPOSED LAW, NOT FROZEN CP-12 AND NOT VALIDATED PHYSICS',
        Hamiltonian='Hnew=H_CP12 tensor I + J sum_e g_e 1_(E_ei E_ej<0) (c_i^dagger c_j+c_j^dagger c_i)',
        probe='New electrically neutral fermion species; one-probe sector for these exact controls.',
        inputs=['Hilbert/composition and original gauge structure','new probe species','new through-field coupling rule','free coupling J','supplied control states'],
        rows=rows,conclusion='Hermitian, gauge-preserving gate and probe-vacuum formation identity hold, but a generic allowed one-leg charged shortcut activates in one geometry hop. Its specified transfer probability has leading coefficient J^2/4 at order T^4, independent of ring distance and flux. This candidate alone does not earn locality.',
        limits='Exact supplied-state transfer controls, not original-neutral formation or density-phase experiments. Fixed-N small-time coefficients; no uniform finite-time large-N remainder bound. Original charged a,b excitations remain and universality is unproved.')
    (P/'PATH_B_ALGEBRA.json').write_text(json.dumps(data,indent=2)+'\n');print(json.dumps(data),flush=True)
