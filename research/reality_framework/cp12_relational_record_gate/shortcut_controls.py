"""Control shortcut statistics for Gauss resources and parent topology."""
from pathlib import Path
import argparse,json,time
import numpy as np
from analyze import P,model,resources,setup
from mechanism import CASES


def run(snapshot_dir):
    start=time.monotonic();models=[]
    for n,q,budget,times in CASES:
        sub,H,B,states,sizes,psi=model(n,q,budget)
        r,k,nv,e2=resources(n,q,sub,states);S=setup(n,q,states,sizes,r,nv,e2)
        references={}
        for name in ['gauss_resources','gauss_parent_resources']:
            inverse=S['controls'][name][1]
            z=np.bincount(inverse,weights=sizes)
            averages=np.column_stack([np.bincount(inverse,weights=sizes*S['shortcut'][:,c])/z for c in range(3)])
            references[name]=averages[inverse]
        source=json.loads((P/f'RELATIONAL_N{n}_Q{q}_B{budget}.json').read_text())
        rows=[]
        for t in times:
            v=np.load(Path(snapshot_dir)/f'N{n}_Q{q}_B{budget}_T{t:g}.npz')['wave'];p=abs(v)**2
            actual=p@S['shortcut']
            ref=min(source['rows'],key=lambda a:abs(a['T']-t))
            rows.append(dict(T=t,actual_neutral_edges_per_Q=actual.tolist(),
                bounded_score_or_contrast_error_upper=ref['bounded_operator_error_upper'],
                controls={name:dict(neutral_edges_per_Q=(p@f).tolist(),
                    contrasts_per_Q=(p@(S['shortcut']-f)).tolist()) for name,f in references.items()}))
        models.append(dict(N=n,Q=q,budget=budget,bins=S['shortcut_bins'],rows=rows))
        print(json.dumps(dict(stage='shortcut_controls',N=n,Q=q)),flush=True)
    out=dict(models=models,seconds=time.monotonic()-start,
        scope='Same restricted physical configuration nulls and identity extension as the core-shape audit, now applied to the joint parent/core neutral-edge statistics. Resources, and for the parent version the parent unlabelled graph, are matched. Full Hamiltonian energy is not matched. Uniform completion is a different, less constrained geometric reference.')
    (P/'SHORTCUT_CONTROLS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(stage='saved',seconds=out['seconds'])),flush=True);return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--snapshot-dir',required=True)
    run(ap.parse_args().snapshot_dir)
