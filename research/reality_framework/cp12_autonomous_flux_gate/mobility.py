"""Finite, gauge-invariant density probes in the autonomously evolved state.

No initial graph-distance record or geometry postselection is introduced.
The source is an external diagnostic; it does not select a spatial phase.
"""
from pathlib import Path
import json,sys,time
import numpy as np
from scipy.sparse import diags,csr_matrix

P=Path(__file__).resolve().parent
sys.path.insert(0,str(P.parent/'cp12_unmasked_formation'))
sys.path.insert(0,str(P.parent/'cp12_distance_history'))
from run import model
from evolution import Evolution
from dynamics import resources,rewire_matrix


def main():
    started=time.monotonic();n=4;q=4;budget=256
    sub,H,B,states,sizes,psi=model(n,q,budget,source=True)
    r,k,nv,e2=resources(n,q,sub,states)
    R=rewire_matrix(n,q,sub,states,sizes,H=H)
    S=np.array([s[1]&1 for s in states],float)
    O=(nv-S)/(n-1)
    no_hop=R+diags(q/2+2*nv+e2/8)
    comm=no_hop@diags(S)-diags(S)@no_hop
    assert np.linalg.norm(comm.data)<1e-10
    E=Evolution(H,B,q/2,16)
    rows=[];v=psi.copy();wait_error=0.;last=0.
    for wait in [5.,10.]:
        for _ in range(round((wait-last)/.1)):
            v,_,de=E(v,.1);wait_error+=de
        last=wait
        A=np.column_stack([(1+(-1j-1)*S)*v,(1+(1j-1)*S)*v])
        works=[float(np.vdot(A[:,j],H@A[:,j]).real-np.vdot(v,H@v).real) for j in range(2)]
        # Source unitaries commute with R and all diagonal energy terms.
        # Only the bounded hopping difference contributes to source work.
        hop_norm_bound=float(np.sqrt(2*n*q))
        work_error=2*hop_norm_bound*(2*wait_error+wait_error**2)
        for lag in [.1,.5,1.]:
            Z,ec,_=E(A,lag)
            population=np.asarray(O@abs(Z)**2)
            susceptibility=float((population[0]-population[1])/2)
            state_errors=wait_error+ec
            bound=float(sum(2*state_errors+state_errors**2)/2)
            Enull=Evolution(no_hop,csr_matrix(B.shape),q/2,16)
            Zn,_,_=Enull(A,lag)
            null=float((O@abs(Zn[:,0])**2-O@abs(Zn[:,1])**2)/2)
            assert abs(null)<1e-11
            rows.append(dict(wait=wait,lag=lag,
                before_source_removable_fraction=float(abs(v)**2@r/q),
                before_source_flux_protected_fraction=float(abs(v)**2@k/q),
                waiting_state_error_upper=wait_error,
                detector='mean vertex density excluding source, independent of graph-distance bins',
                finite_source_phases=[float(np.pi/2),float(-np.pi/2)],
                detector_populations=population.tolist(),
                susceptibility=susceptibility,
                susceptibility_error_upper=bound,
                susceptibility_interval=[susceptibility-bound,susceptibility+bound],
                pulse_work=works,pulse_work_error_upper=work_error,
                no_hopping_after_same_preparation_response=null,
                after_source_removable_fractions=[float(abs(Z[:,j])**2@r/q) for j in range(2)]))
    out=dict(status='CP12_AUTONOMOUS_FLUX_MOBILITY_PASS',N=n,Q=q,budget=budget,
             dimension=len(states),source_gauge_invariant=True,
             initial_distance_record=False,geometry_postselection=False,
             rows=rows,seconds=time.monotonic()-started,
             scope='Finite-system source diagnostic after unforced partial flux protection. The preparation, source and time remain supplied. No macroscopic speed, intrinsic cone, extended 3D geometry, or autonomous apparatus is derived.')
    (P/'MOBILITY.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)


if __name__=='__main__':main()
