"""Post-run integrity audit; does not change the frozen B2 decision rule."""
import argparse,copy,hashlib,itertools,json,sys
from pathlib import Path
import numpy as np

def transfer_policy(massless_pass,has_result):
    return not has_result or massless_pass

def distribution_errors(points,card):
    errors=[]
    for name,d in points.items():
        # Check both lattice and reference groups separately, including equal names.
        groups={**{'lattice_'+k:v for k,v in d['histories'].items()},
                **{'reference_'+k:v['histories'] for k,v in d['references'].items() if k in ('target','B1','B2')}}
        expected_keys={','.join(f'{t:g}' for t in times) for times in card['protocols']}
        for group,hist in groups.items():
            if set(hist)!=expected_keys:errors.append(name+'/'+group+'/protocols')
            for key,p in hist.items():
                labels={''.join(s) for s in itertools.product('+-',repeat=key.count(',')+1)}
                if set(p)!=labels or not all(np.isfinite(v) and v>=-1e-12 for v in p.values()) or abs(sum(p.values())-1)>1e-9:
                    errors.append(name+'/'+group+'/'+key)
    return errors

def run(root,out):
    here=root/'research/reality_framework/grut_detector_reconstruction_20261108'
    sys.path.insert(0,str(here))
    from model import reconstruction,base_weights
    from lattice import schwinger
    from analyze import analyze,massive
    card=json.loads((here/'protocol_b2_01.json').read_text())
    provenance=json.loads((out/'provenance.json').read_text())
    hashes={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in provenance['source_sha256']}
    source_ok=hashes==provenance['source_sha256']
    points={p.stem:json.loads(p.read_text()) for p in out.glob('N*.json')}
    distribution_failures=distribution_errors(points,card)
    checkpoint_checks={}
    for name,d in points.items():
        n=d['N'];a=d['a'];saved=np.load(out/(name+'.npz'),allow_pickle=False)
        w,electric,states=schwinger.hamiltonian(n,a**-2,-a/8);h=w*a/2
        ground=saved['ground'];_,_,weights,boundary,_=reconstruction(card,n)
        _,_,oldweights,b=base_weights(card,n);stag=(-1.)**np.arange(n-1)
        fields={'B1':electric@oldweights+(b-oldweights@stag)/4,
                'B2':electric@weights+(boundary.sum()-weights@stag)/4}
        first=max(abs(.5*(1+np.vdot(ground,np.exp(1j*card['coupling']*field)*ground).real)-d['histories'][label]['0']['+'])
                  for label,field in fields.items())
        checks={'ground_norm':float(abs(np.vdot(ground,ground).real-1)),
            'ground_residual':float(np.linalg.norm(h@ground-saved['eigenvalues'][0]*ground)),
            'weights':float(np.max(abs(saved['weights']-weights))),
            'boundary':float(np.max(abs(saved['boundary']-boundary))),
            'first_event_independent_formula':float(first),
            'states_identical':bool(np.array_equal(states,saved['states']))}
        checkpoint_checks[name]=checks
    computed=analyze(out);saved_assessment=json.loads((out/'assessment.json').read_text())
    transfer=massive(out);saved_transfer=json.loads((out/'massive_status.json').read_text())
    transfer_allowed=transfer_policy(computed['massless_pass'],(out/'massive_N20.json').exists())
    damaged=copy.deepcopy(points);damaged['N20']['histories']['B2']['0']['+']+=.01
    forged=copy.deepcopy(saved_assessment);forged['massless_pass']=not computed['massless_pass']
    fake_hashes=dict(provenance['source_sha256']);first_key=next(iter(fake_hashes));fake_hashes[first_key]='0'*64
    negatives={'corrupt_probability_rejected':bool(distribution_errors(damaged,card)),
               'forged_verdict_rejected':forged!=computed,'source_hash_change_rejected':fake_hashes!=hashes,
               'bypassed_massive_gate_rejected':not transfer_policy(False,True) and transfer_policy(False,False) and transfer_policy(True,True)}
    checks_ok=all(c['states_identical'] and max(v for k,v in c.items() if k!='states_identical')<1e-9 for c in checkpoint_checks.values())
    result={'source_hashes_match':source_ok,'distribution_failures':distribution_failures,
      'checkpoint_checks':checkpoint_checks,'assessment_recomputed_identically':computed==saved_assessment,
      'massive_status_recomputed_identically':transfer==saved_transfer,'massive_execution_obeys_rule':transfer_allowed,
      'negative_controls':negatives,
      'passed':bool(source_ok and not distribution_failures and checks_ok and computed==saved_assessment and transfer==saved_transfer and transfer_allowed and all(negatives.values())),
      'scope':'Post-run arithmetic, integrity and stored-ground-state audit; not an independent full physical implementation.'}
    with (out/'audit.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();run(args.root.resolve(),args.output.resolve())
