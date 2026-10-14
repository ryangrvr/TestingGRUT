"""Matched-state term controls after the intervention-time record.

Only diagnostic comparator evolutions change terms. The reported main runs
always use the unchanged frozen Hamiltonian.
"""
from pathlib import Path
from math import lcm,sqrt
import json,numpy as np
from scipy.sparse import coo_matrix,diags,csr_matrix
from basis import model,distance,transitions
from evolution import Evolution,dual

P=Path(__file__).resolve().parent

def run():
    n=5;q=4;budget=10;sub,H,B,states,sizes,psi=model(n,q,budget)
    index={s:i for i,s in enumerate(states)};rr=[];cc=[];vv=[]
    scale=lcm(8,10*sub.M)
    for j,s in enumerate(states):
        for t,c in transitions(n,s,scale,hopping=False):
            key,sign,size=sub.canonical(t)
            if not sign:continue
            assert key in index  # rewiring never leaves the budget projection
            rr.append(index[key]);cc.append(j);vv.append(c/scale*sign*sqrt(sizes[j]/size))
    no_hop=coo_matrix((vv,(rr,cc)),shape=H.shape).tocsr()
    diag=np.array([q/2+2*(bits&((1<<n)-1)).bit_count()+sum(x*x for x in E)/8 for g,bits,E in states])
    R=no_hop-diags(diag);no_rewire=H-R
    assert np.linalg.norm(R@psi)<1e-10
    D=distance(n,states,sub.es);F=(D==n).astype(float)
    S=np.array([b&1 for g,b,e in states],float);T=np.array([(b>>1)&1 for g,b,e in states],float)
    source_comm=(no_hop@diags(S)-diags(S)@no_hop).tocsr()
    assert np.linalg.norm(source_comm.data)<1e-12
    distance_comm=(no_rewire@diags(D)-diags(D)@no_rewire).tocsr()
    assert np.linalg.norm(distance_comm.data)<1e-10
    Efull=Evolution(H,B,q/2,16);v,_,ew=Efull(psi,.1);a=F*v
    rows=[]
    for lag in [.1,.3]:
        Enr=Evolution(no_rewire,B,q/2,16)
        Z,ec,_=Enr(np.column_stack([a,S*a]),lag)
        response=dual(Enr,Z[:,0],Z[:,1],T,lambda z:F*z,lambda z:F*S*z,.1,lag,float(ew+ec[0]),float(ew+ec[1]),pre=Efull)
        assert abs(response['value'])<=response['error']+1e-12
        Enh=Evolution(no_hop,csr_matrix(B.shape,dtype=float),q/2,16)
        Z,_,_=Enh(np.column_stack([a,S*a]),lag)
        null=float(2*np.vdot(Z[:,0],T*Z[:,1]).imag)
        assert abs(null)<1e-12
        changed=float(np.linalg.norm((1-F)*Z[:,0])**2)
        main=json.loads((P/f'N5_Q4_B10_W0.1_T{lag:g}.json').read_text())
        full=next(row for row in main['past'] if row['disconnected'])
        rows.append(dict(lag=lag,full_response=full['value'],full_error=full['error'],rewiring_removed_after_same_preparation=response,
                         hopping_removed_response=null,hopping_removed_geometry_change_probability=changed))
    out=dict(status='CP12_DISTANCE_HISTORY_ATTRIBUTION_PASS',N=n,Q=q,budget=budget,
             same_full_law_waiting_state_and_distance_record=True,
             no_hopping_source_commutator_norm=float(np.linalg.norm(source_comm.data)),
             no_rewiring_distance_commutator_norm=float(np.linalg.norm(distance_comm.data)),rows=rows,
             scope='Diagnostic term controls only. Exact full-law no-rewiring disconnected-response null follows from component factorization; numerical global electric budgets can introduce spurious effects that must remain within their bounds.')
    (P/'ATTRIBUTION.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)

if __name__=='__main__':run()
