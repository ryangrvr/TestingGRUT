"""Which term changes protected-core topology? Reuse final wave snapshots."""
from pathlib import Path
import argparse,json,time
import numpy as np
from analyze import P,model,resources,setup,rewire_matrix

CASES=[(4,4,256,[1.,3.,5.,10.]),(5,5,128,[1.,3.,5.,10.]),
       (6,6,128,[1.,3.,5.,10.]),(6,12,10,[.4,.6,1.])]


def run(snapshot_dir):
    start=time.monotonic();out=[]
    for n,q,budget,times in CASES:
        sub,H,B,states,sizes,psi=model(n,q,budget)
        r,k,nv,e2=resources(n,q,sub,states);S=setup(n,q,states,sizes,r,nv,e2)
        R=rewire_matrix(n,q,sub,states,sizes,H=H);C=R.tocoo()
        different=S['classes'][C.row]!=S['classes'][C.col]
        residual=float(np.linalg.norm(C.data[different]))
        assert residual<1e-12
        old=json.loads((P/f'RELATIONAL_N{n}_Q{q}_B{budget}.json').read_text())
        hbound=float(np.sqrt(2*n*q))
        rows=[]
        for t in times:
            v=np.load(Path(snapshot_dir)/f'N{n}_Q{q}_B{budget}_T{t:g}.npz')['wave']
            Hv=H@v;Rv=R@v
            rates=2*np.imag(S['actual'].T@(v.conj()*Hv))
            rrates=2*np.imag(S['actual'].T@(v.conj()*Rv))
            ref=min(old['rows'],key=lambda a:abs(a['T']-t))
            assert np.max(abs(rrates))<1e-11
            delta=ref['bounded_operator_error_upper']
            rows.append(dict(T=t,protected_core_geometry_derivatives=dict(zip(S['names'],map(float,rates))),
                largest_rewiring_only_derivative=float(np.max(abs(rrates))),
                generic_core_derivative_error_upper=2*hbound*delta,
                shortcut_rewiring_current_error_upper=.1*(sub.M-q)/sub.M*delta))
        # Every physical diagonal configuration obeys this SHARP gap.
        a=.1*(sub.M-q)/sub.M
        diagonal=q/2+2*nv+e2/8+a*r
        gap_lower=a*q
        assert np.min(diagonal-q/2)>=gap_lower-1e-12
        out.append(dict(N=n,Q=q,budget=budget,rewiring_core_shape_commutator_residual=residual,
            hopping_norm_upper=hbound,configuration_diagonal_mean_energy_lower=q/2+gap_lower,
            neutral_preparation_mean_energy=q/2,configuration_diagonal_energy_gap=gap_lower,
            rows=rows))
        print(json.dumps(dict(stage='mechanism',N=n,Q=q)),flush=True)
    result=dict(models=out,seconds=time.monotonic()-start,
        exact_theorem='Every protected-edge projector p_e and every function of the protected core commute with H_rewire and H_diagonal on the full physical Hilbert space. Only H_hop changes the core or its topology. Neutral rewiring can still affect later core dynamics through the state.',
        energy_theorem='At frozen couplings and 0<Q<M, every full-configuration-diagonal density operator has <H> >= kappa Q + Gamma Q(M-Q)/M > kappa Q. Sharp equality at the neutral diagonal configurations. Thus no fully configuration-incoherent null matches the initial mean energy.',
        scope='Core rates and current errors are ordinary floating-point evaluations of finite-matrix expectations and analytic bounded-operator inequalities. No inference about the lifetime of the same edge or loop is made.')
    (P/'MECHANISM.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(stage='saved',seconds=result['seconds'])),flush=True)
    return result


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--snapshot-dir',required=True)
    run(ap.parse_args().snapshot_dir)
