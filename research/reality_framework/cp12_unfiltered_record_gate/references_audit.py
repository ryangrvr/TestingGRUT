"""Apply the no-preselection readout to the existing physical energy references."""
import json,time
import numpy as np
from audit import P,model,energy_basis,inclusion,detectors,evolve_paths


def run():
    start=time.monotonic();full,energies,V,weights,beta=energy_basis()
    source=model(4,4,256,source=True);J,checks=inclusion(full,source)
    sub,H,B,states,sizes,psi=source;O,events=detectors(4,states)
    rows=[]
    for label,w in weights.items():
        X=J@(V*np.sqrt(w)[None,:])
        result,_=evolve_paths(X,0.,O['source_density'],H,B,O,.5,save_waves=False)
        mean=float(np.vdot(X,H@X).real)
        assert abs(mean-2)<1e-11
        second=float(np.linalg.norm(H@X)**2+np.linalg.norm(B@X)**2)
        rows.append(dict(preparation=label,full_mean_energy=mean,full_energy_variance=second-mean**2,
            entropy_bits=float(-sum(x*np.log2(x) for x in w if x>0)),response=result))
        print(json.dumps(dict(stage='energy_reference',preparation=label,
            density=result['susceptibilities']['density'],
            inside_circulation=result['density_response_inside_final_circulation'])),flush=True)
    out=dict(N=4,Q=4,budget=256,lag=.5,inverse_temperature=beta,
        embedding_checks=checks,rows=rows,seconds=time.monotonic()-start,
        scope='The two positive configuration-coherent finite-support physical reference density operators from 9d0f2ef and 7d41496. They match full mean energy, not all resources/purity or full-law equilibrium. No initial loop filter is applied. Their subsequent unbounded-H evolution is controlled using the full outgoing residual. Final circulation conditioning remains readout postselection.')
    (P/'ENERGY_REFERENCE_RESPONSE.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


if __name__=='__main__':run()
