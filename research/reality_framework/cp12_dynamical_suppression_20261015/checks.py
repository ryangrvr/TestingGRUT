"""Independent full-basis references, exact finite spectral integrals and exits."""
from pathlib import Path
from math import sqrt,lcm
import json,sys
import numpy as np
from scipy.linalg import eigh
from scipy.sparse import coo_matrix,vstack,eye
from scipy.sparse.linalg import expm_multiply
from boundary import quotient_model,recover_outer,extend,RF
from mechanisms import source_statistics,NAMES,action_operators
from action_integrals import ActionEvolution
from hopping_bound import hopping_bound,endpoint_maximum
from linear_propagator import MatrixFreePropagator
from operators import physical_model,flow,hopping_edge
from integrals import Propagator
sys.path.insert(0,str(RF/'cp12_ground_connectivity'))
from exact_destination import transitions,gauss
from invariant_evolution import Invariant

P=Path(__file__).resolve().parent


def main():
    m=quotient_model(4,5,12);full=physical_model(4,5,12);sub=m['sub']
    old=m['states'];qi={s:i for i,s in enumerate(old)};phys=full['states'];pi={s:i for i,s in enumerate(phys)}
    outer={};rr=[];cc=[];vv=[];br=[];bc=[];bv=[];dref=np.zeros(len(phys));hopcount=0
    filling_counts={}
    for col,s in enumerate(phys):
        filling=(s[1]&15).bit_count();up=down=0
        for t,c in transitions(4,s,120):
            e=hopping_edge(4,s,t)
            if e is None:continue
            if (t[1]&15).bit_count()>filling:up+=1
            else:down+=1
            gauss(4,t);w=1-.5*(bool(flow(4,s)[0]>>e&1)+bool(flow(4,t)[0]>>e&1))
            dref[col]+=(c/120*w)**2;hopcount+=1
            if t in pi:rr.append(pi[t]);cc.append(col);vv.append(c/120*w)
            else:
                row=outer.setdefault(t,len(outer));br.append(row);bc.append(col);bv.append(c/120*w)
        assert up<=endpoint_maximum(4,4-filling,5-filling)
        assert down<=endpoint_maximum(4,filling,filling)
        before=filling_counts.setdefault(filling,[0,0])
        before[0]=max(before[0],up);before[1]=max(before[1],down)
    V=coo_matrix((vv,(rr,cc)),shape=full['H'].shape).tocsr()
    BV=coo_matrix((bv,(br,bc)),shape=(len(outer),len(phys))).tocsr()
    jr=[];jc=[];jv=[];sub.cache.clear()
    for row,s in enumerate(phys):
        key,sign,size=Invariant.canonical(sub,s)
        if sign and key in qi:jr.append(row);jc.append(qi[key]);jv.append(sign/sqrt(size))
    J=coo_matrix((jv,(jr,jc)),shape=(len(phys),len(old))).tocsr()
    stats=np.asarray([source_statistics(4,5,s) for s in old]).T
    reference=np.asarray((J.multiply(J)).T@dref).ravel()
    difference=float(np.max(abs(stats[0]-reference)));assert difference<1e-12
    wrong=np.asarray(m['V'].multiply(m['V']).sum(axis=0)+m['BV'].multiply(m['BV']).sum(axis=0)).ravel()
    distinction=float(np.max(abs(wrong-reference)));assert distinction>.1
    rng=np.random.default_rng(768);v=rng.normal(size=len(old))+1j*rng.normal(size=len(old));v/=np.linalg.norm(v)
    coherent=float(np.linalg.norm(m['V']@v)**2+np.linalg.norm(m['BV']@v)**2)
    coherent_phys=float(np.linalg.norm(V@(J@v))**2+np.linalg.norm(BV@(J@v))**2)
    assert abs(coherent-coherent_phys)<1e-11
    ds=dict(zip(NAMES,stats));bounds=dict(zip(NAMES,[10,10,5,10,10,10,640,1,1,1]))
    actions=action_operators(m,4)
    T=.375
    norm_bound=hopping_bound(4,5)[0]
    action_norms={name:float(np.linalg.eigvalsh((vstack(blocks).T@vstack(blocks)).toarray())[-1]**.5)
        for name,blocks in actions.items()}
    assert max(action_norms.values())<norm_bound
    result=ActionEvolution(m['H'],m['B'],actions,ds,bounds,5,2.5,dt=1/64,
        hopping_norm_bound=norm_bound).run(v,[T])[T]
    energies,Q=eigh(m['H'].toarray());z=Q.T@v
    delta=energies[:,None]-energies[None,:]
    F=T*np.exp(.5j*delta*T)*np.sinc(delta*T/(2*np.pi))
    spectral=[]
    for name in ['outside','inside','dephased_outside_action','neutral_edge_vacancies']:
        if name in actions:
            L=vstack(actions[name],format='csr');G=(Q.T@(L.T@L@Q))
        else:G=Q.T@(ds[name][:,None]*Q)
        exact=float(np.sum(z.conj()[:,None]*G*z[None,:]*F).real)
        I=result['integrals'][name]['retained_quadrature_interval']
        assert I[0]-1e-10<=exact<=I[1]+1e-10
        spectral.append(dict(action=name,exact_spectral_integral=exact,
            computed=result['integrals'][name]['value'],quadrature_interval=I))
    olderr=Propagator(m['H'],{'B':m['B']},2.5,dt=1/64).run(v,[T])[T]
    stream_error=float(np.linalg.norm(olderr[0]-result['state']))
    residual_error=abs(float(olderr[1]['B'][0])-result['lag_state_error'])
    assert stream_error<1e-11 and residual_error<1e-11
    lean=ActionEvolution(m['H'],m['B'],actions,ds,bounds,5,2.5,dt=1/64,
        hopping_norm_bound=norm_bound,matrix_free=True).run(v,[T])[T]
    lean_error=max(float(np.linalg.norm(lean['state']-result['state'])),
        max(abs(lean['integrals'][k]['value']-result['integrals'][k]['value']) for k in result['integrals']),
        abs(lean['lag_state_error']-result['lag_state_error']))
    lean_stream=MatrixFreePropagator(m['H'],{'B':m['B']},2.5,dt=1/64).run(v,[T])[T]
    lean_stream_error=max(float(np.linalg.norm(lean_stream[0]-olderr[0])),
        abs(float(lean_stream[1]['B'][0])-float(olderr[1]['B'][0])))
    assert max(lean_error,lean_stream_error)<1e-11
    boundary=recover_outer(m,4,5,12);m['outer_states']=boundary['states']
    choices=np.sort(rng.choice(len(m['outer_states']),23,replace=False))
    ext=extend(m,4,5,choices,'CHECK_N4_PARTIAL23')
    si={s:i for i,s in enumerate(ext['states'])};oi={s:i for i,s in enumerate(ext['outer_states'])}
    hr=[];hc=[];hv=[];xr=[];xc=[];xv=[];rewiring_exits=0
    for col,s in enumerate(ext['states']):
        for t,c in transitions(4,s,120):
            key,sign,size=Invariant.canonical(sub,t)
            if not sign:continue
            amp=c/120*sign*sqrt(ext['sizes'][col]/size)
            if key in si:hr.append(si[key]);hc.append(col);hv.append(amp)
            else:
                xr.append(oi[key]);xc.append(col);xv.append(amp)
                rewiring_exits+=s[0]!=t[0]
    HH=coo_matrix((hv,(hr,hc)),shape=ext['H'].shape).tocsr()
    BB=coo_matrix((xv,(xr,xc)),shape=ext['B'].shape).tocsr()
    assert np.linalg.norm((HH-ext['H']).data)<1e-10 and np.linalg.norm((BB-ext['B']).data)<1e-10
    assert rewiring_exits>0
    aa=action_operators(ext,4);no=np.array([(s[1]&15).bit_count() for s in ext['outer_states']])
    ni=np.array([(s[1]&15).bit_count() for s in ext['states']]);b=ext['B'].tocoo();mask=no[b.row]==ni[b.col]
    rnorm=float(np.linalg.norm(b.data[mask]));assert rnorm>0
    out=dict(status='PASS',physical_hopping_Gauss_checks=hopcount,
        full_physical_dephased_reference_error=difference,
        orbit_diagonal_reference_would_be_wrong_by=distinction,
        coherent_action_full_physical_difference=abs(coherent-coherent_phys),
        exact_spectral_integral_controls=spectral,
        inherited_propagator_state_difference=stream_error,inherited_residual_bound_difference=residual_error,
        matrix_free_action_and_residual_difference=lean_error,
        matrix_free_inherited_propagator_difference=lean_stream_error,
        adaptive_full_column_checks=len(ext['states']),adaptive_rewiring_exit_transitions=rewiring_exits,
        adaptive_outgoing_nonhopping_norm=rnorm,
        physical_filling_sector_maximum_up_down_counts=filling_counts,
        finite_action_norm_checks=action_norms,
        integer_hopping_bound_controls={str(n):dict(zip(['bound','Jacobi_offdiagonal_squares','positive_leading_minors'],
            hopping_bound(n,q))) for n,q in [(4,5),(5,10),(6,12)]},
        interpretation='Independent labelled physical baseline and signed embedding, finite spectral integration, inherited residual stream, and every unregulated adaptive column. Numerical basis selection does not change the physical law.')
    (P/'CHECKS.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))


if __name__=='__main__':main()
