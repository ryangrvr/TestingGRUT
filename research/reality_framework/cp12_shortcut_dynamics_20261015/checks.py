"""Independent full-basis, explicit-cycle and inherited-propagator controls."""
from pathlib import Path
from math import sqrt
import json,sys
import numpy as np
from scipy.sparse import coo_matrix,eye
from scipy.sparse.linalg import expm_multiply
from operators import quotient_model,physical_model,hopping_edge,outside_weight,flow,RF
from integrals import Propagator
from recorded_transport import ring_for,ring_present,selected,opposite
sys.path.insert(0,str(RF/'charged_density_dynamics_20261015'))
from diagnostics import cycle_edges
sys.path.insert(0,str(RF/'cp12_ground_connectivity'))
from exact_destination import transitions,gauss
from invariant_evolution import Invariant
sys.path.insert(0,str(RF/'cp12_distance_history'))
from evolution import Evolution


def main():
    physical=physical_model(4,5,12);states=physical['states'];index={s:i for i,s in enumerate(states)}
    H=physical['H'];B=physical['B'];rr=[];cc=[];vv=[]
    checked=reverse=unequal=boundary=perm=0
    for col,s in enumerate(states):
        mask,reach=flow(4,s)
        assert mask==sum(1<<e for e,u,v in cycle_edges(4,s))
        for t,c in transitions(4,s,120):
            gauss(4,t);e=hopping_edge(4,s,t)
            if e is None:continue
            weight=1.-.5*(bool(mask>>e&1)+bool(flow(4,t)[0]>>e&1))
            assert outside_weight(4,s,t,e,(mask,reach))==weight
            assert outside_weight(4,t,s,e)==weight
            reverse+=1;checked+=1;boundary+=weight==.5
            unequal+=s[2][2*e]*s[2][2*e+1]<0 and s[2][2*e]!=-s[2][2*e+1]
            if t in index:rr.append(index[t]);cc.append(col);vv.append(c/120*weight)
    V=coo_matrix((vv,(rr,cc)),shape=H.shape).tocsr()
    controls=[]
    for source in [False,True]:
        q=quotient_model(4,5,12,source);sub=q['sub'];qi={s:i for i,s in enumerate(q['states'])}
        jr=[];jc=[];jv=[];sub.cache.clear()
        for row,s in enumerate(states):
            # Independently use the old pure-Python signed canonicalizer.
            key,sign,size=Invariant.canonical(sub,s)
            if sign and key in qi:jr.append(row);jc.append(qi[key]);jv.append(sign/sqrt(size))
        J=coo_matrix((jv,(jr,jc)),shape=(len(states),len(qi))).tocsr()
        iso=float(np.linalg.norm((J.T@J-eye(len(qi))).data))
        hi=float(np.linalg.norm((H@J-J@q['H']).data));vi=float(np.linalg.norm((V@J-J@q['V']).data))
        gram=float(np.linalg.norm((J.T@B.T@B@J-q['B'].T@q['B']).data))
        assert max(iso,hi,vi,gram)<1e-10
        for s in states[::max(1,len(states)//25)]:
            cm=flow(4,s)[0]
            for p in sub.group:
                t,sign=sub.transform(s,p)
                expected=sum(1<<p['edge'][e] for e in range(6) if cm>>e&1)
                assert flow(4,t)[0]==expected;perm+=1
        controls.append(dict(source_stabilizer=source,physical_dimension=len(states),
            quotient_dimension=len(qi),isometry_error=iso,H_intertwiner_error=hi,
            V_intertwiner_error=vi,outgoing_Gram_error=gram))
    # Compare the new independent multiple-integral stream with the inherited
    # full-grid method, and its action integral against dense-time quadrature.
    q=quotient_model(4,4,32);H=q['H'];B=q['B'];V=q['V'];psi=q['psi']
    w,_,e=Evolution(H,B,2.,steps=64)(psi,.5)
    ws,es=Propagator(H,{'electric':B,'Vret':V},2.,dt=1/128).run(psi,[.5])[.5]
    difference=float(np.linalg.norm(w-ws));assert difference<1e-11
    grid=expm_multiply(-1j*(H-2*eye(H.shape[0])),psi,start=0,stop=.5,num=2001)
    numeric=float(np.trapezoid(np.linalg.norm(V@grid.T,axis=0),dx=.5/2000))
    assert numeric<=es['Vret'][0]+1e-10
    # Direct physical graph and record checks, with no drawing-based distance.
    for g in sorted({s[0] for s in states}):
        C,chord=ring_for(g);assert g&C==C;assert C.bit_count()==4
        j=opposite(C);assert not any(0 in pair and j in pair for e,pair in enumerate(__import__('itertools').combinations(range(4),2)) if C>>e&1)
        for s in [s for s in states if s[0]==g]:
            if selected(s):assert flow(4,s)[0]==C and ring_present(s,C)
    result=dict(status='PASS',full_basis_hops_and_Gauss_checks=checked,
        Hermitian_reverse_weight_checks=reverse,half_weight_core_change_hops=boundary,
        unequal_opposite_sign_hole_hops_checked=unequal,
        permutation_core_checks=perm,independent_physical_quotient_controls=controls,
        inherited_full_grid_state_difference=difference,
        action_integral_numerical_quadrature=numeric,action_integral_analytic_upper=float(es['Vret'][0]),
        scope='Independent labelled physical basis, old Python signed orbit embedding, direct source/target cycle indicators, reverse transitions, permutation checks, and inherited exponential evolution. Internal controls, not external review.')
    Path(__file__).with_name('CHECKS.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
