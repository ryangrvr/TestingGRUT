"""Explicit energy backaction and conditional energies of diagnostic filters."""
from pathlib import Path
import argparse,json,time
import numpy as np
from loops import P,UP,model,fixed_triangle
from references import energy_basis,inclusion


def division_interval(numerator,denominator):
    a,b=numerator;c,d=denominator
    assert c>0
    vals=[a/c,a/d,b/c,b/d]
    return [min(vals),max(vals)]


def run(snapshot_dir):
    start=time.monotonic();full,energies,V,weights,beta=energy_basis()
    source=model(4,4,256,source=True);J,checks=inclusion(full,source)
    src,H,B,states,sizes,psi=source;D,_=fixed_triangle(4,states)
    trajectory=json.loads((UP/'N4_Q4_B256_T20.json').read_text())
    cases=[]
    for wait in [5.,10.]:
        v=np.load(Path(snapshot_dir)/f'N4_Q4_B256_T{wait:g}.npz')['wave']
        e=min(trajectory['rows'],key=lambda a:abs(a['T']-wait))['state_error_upper']
        cases.append((f'neutral_sea_wait_{wait:g}',(J@v)[:,None],e,np.sqrt(12.)))
    for label,w in weights.items():
        X=J@(V*np.sqrt(w)[None,:]);cases.append((label,X,0.,None))
    rows=[]
    for label,X,e,exact_Hnorm in cases:
        Hnorm=float(np.sqrt(np.linalg.norm(H@X)**2+np.linalg.norm(B@X)**2))
        entries=[]
        for key in ['hole','occupied_flux','vacuum_flux','electric_circulation']:
            A=D[key][:,None]*X;Z=X-A;p=float(np.linalg.norm(A)**2)
            numerator=float(np.vdot(A,H@A).real)
            pinterval=[max(0.,np.sqrt(p)-e)**2,min(1.,(np.sqrt(p)+e)**2)]
            # Pi commutes with H_D and H_R; [H,Pi]=[Hhop,Pi] is bounded.
            # ||Pi H Pi psi|| <= ||H psi||+2||Hhop||; likewise for X.
            en=e*((exact_Hnorm if exact_Hnorm is not None else Hnorm)+Hnorm+4*np.sqrt(32.))
            cond=division_interval([numerator-en,numerator+en],pinterval)
            work=float(-2*np.vdot(A,H@Z).real)
            workerr=2*np.sqrt(32.)*(2*e+e*e)
            entries.append(dict(initial_channel=key,event_probability=p,event_probability_interval=pinterval,
                conditional_full_mean_energy=numerator/p,conditional_full_mean_energy_interval=cond,
                binary_nonselective_recording_work=work,binary_recording_work_interval=[work-workerr,work+workerr]))
        rows.append(dict(preparation=label,mean_energy_before_filter=2.,channels=entries))
    out=dict(N=4,Q=4,budget=256,rows=rows,seconds=time.monotonic()-start,
        scope='References and main inputs are mean-energy matched BEFORE conditioning. Postselected channel energies are computed and need not match. Main-state unbounded conditional-energy errors use ||H psi(T)||=sqrt(12), conservation of the initial second energy moment, and the bounded [H,Pi]=[Hhop,Pi], rather than a norm-only estimate for an unbounded observable. Counterfactual measurement backaction is explicit; apparatus implementation cost is not derived. Ordinary floating point.')
    (P/'MEASUREMENT_CONTROLS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(stage='measurement_controls',seconds=out['seconds'])),flush=True);return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--snapshot-dir',required=True)
    run(ap.parse_args().snapshot_dir)
