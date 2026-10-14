"""Exact finite-pulse extraction for a projector source.

For S^2=S, [f(+eta)-f(-eta)]/(2 sin eta) equals the linear response exactly.
This checks the susceptibility without an infinitesimal finite difference.
"""
from pathlib import Path
import json,numpy as np
from scipy.sparse.linalg import expm_multiply
from basis import model,distance

P=Path(__file__).resolve().parent

def main():
    n=5;q=4;budget=10;sub,H,B,states,sizes,psi=model(n,q,budget)
    v=expm_multiply(-.1j*H,psi)
    S=np.array([bits&1 for g,bits,E in states]);T=np.array([(bits>>1)&1 for g,bits,E in states])
    F=(distance(n,states,sub.es)==n);a=F*v;p=float(np.linalg.norm(a)**2)
    plus=np.exp(-.5j*np.pi*S)*a;minus=np.exp(.5j*np.pi*S)*a
    energy=float(np.vdot(a,H@a).real)
    works=[float((np.vdot(z,H@z).real-energy)/p) for z in [plus,minus]]
    rows=[]
    for lag in [.1,.3]:
        Z=expm_multiply(-1j*lag*H,np.column_stack([plus,minus]))
        pop=[float(T@abs(Z[:,j])**2) for j in [0,1]]
        measured=(pop[0]-pop[1])/2
        main=json.loads((P/f'N5_Q4_B10_W0.1_T{lag:g}.json').read_text())
        expected=next(x for x in main['past'] if x['disconnected'])
        diff=abs(measured-expected['value'])
        assert diff<1e-13
        rows.append(dict(lag=lag,positive_pulse_target_population=pop[0],negative_pulse_target_population=pop[1],
                         exact_pulse_extracted_response=measured,history_bilinear_response=expected['value'],absolute_numerical_difference=diff))
    out=dict(status='CP12_FINITE_PROJECTOR_PULSE_PASS',N=n,Q=q,budget=budget,initial_distance='disconnected',class_probability=p,
             pulse_phase_magnitudes='plus/minus pi/2',conditional_energy_changes_for_plus_and_minus=works,rows=rows,
             source_probability_preserved=True,scope='Finite-pulse extraction identity is exact for the projector source. Pulse energies are computed resource diagnostics in the retained space; the controller is external.')
    (P/'PULSE.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out),flush=True)

if __name__=='__main__':main()
