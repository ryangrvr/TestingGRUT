"""Verify B1 archival integrity, complete probabilities and transfer barrier."""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path
import shutil
import tempfile

import numpy as np
from analyze import analyze,analyze_massive

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]


def require(test,message):
    if not test: raise ValueError(message)


def validate(out):
    out=Path(out);provenance=json.loads((out/'provenance.json').read_text())
    require(provenance['source_clean'],'Unclean production source')
    for p,expected in provenance['source_and_input_sha256'].items():
        require(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==expected,'Frozen source/input changed: '+p)
    analysis=json.loads((out/'analysis.json').read_text())
    require(analysis==analyze(out),'Frozen analysis mismatch')
    require(json.loads((out/'massive_stage.json').read_text())==analyze_massive(out),'Massive-stage verdict mismatch')
    if not analysis['massless_pass']:
        require(not (out/'massive_N20.json').exists() and not (out/'massive_N20.npz').exists(),
                'Corrected massive data evaluated after a massless gate miss')
    count=0;maxnorm=0.;minprob=1.;maxclosure=0.
    for point in list(out.glob('N*.json'))+list(out.glob('massive_N20.json')):
        d=json.loads(point.read_text());n=d['N']
        event=json.loads((out/('event_'+point.stem+'.json')).read_text())
        require(event['status']=='complete','Point did not complete')
        checkpoint=np.load(point.with_suffix('.npz'),allow_pickle=False)
        require(d['dimension']==math.comb(n,n//2)==len(checkpoint['ground']), 'Wrong neutral sector dimension')
        require(abs(np.vdot(checkpoint['ground'],checkpoint['ground']).real-1)<1e-10,'Ground not normalized')
        require(np.max(abs(checkpoint['cell_weights']-d['weights']))<1e-14,'Weight checkpoint mismatch')
        require(abs(sum(d['weights'])+d['boundary_weight']-1)<1e-12,'Weights renormalized incorrectly')
        groups=[d['new_histories'],d['raw_replay_histories']]
        if 'reference' in d:
            groups += [d['reference']['target_histories'],d['reference']['cell_histories'],
                       d['gaussian_lattice_summary'],d['gaussian_static_summary']]
        for histories in groups:
            require(set(histories)=={'0','0,1','0,2','0,1,2'},'Missing history')
            for k,p in histories.items():
                expected={''.join(s) for s in itertools.product('+-',repeat=len(k.split(',')))}
                require(set(p)==expected,'Missing or postselected outcome')
                v=np.array(list(p.values()));require(np.isfinite(v).all() and v.min()>=-1e-10,'Invalid probability')
                err=abs(v.sum()-1);require(err<1e-9,'Unnormalized distribution')
                maxnorm=max(maxnorm,float(err));minprob=min(minprob,float(v.min()));count+=1
        if point.stem in analysis.get('diagnostic_decomposition',{}):
            diag=analysis['diagnostic_decomposition'][point.stem]
            for k,p in d['new_histories'].items():
                for s,v in p.items():
                    residual=sum(c[k][s] for c in diag['signed_contrasts'].values())-(v-d['reference']['target_histories'][k][s])
                    maxclosure=max(maxclosure,abs(residual))
    require(maxclosure<1e-12,'Signed decomposition fails closure')
    return {'verification':'preservation_checks_pass','source_and_input_hashes':len(provenance['source_and_input_sha256']),
            'complete_distributions':count,'max_probability_norm_error':maxnorm,'minimum_probability':minprob,
            'max_signed_closure_error':maxclosure,'massless_verdict':analysis['verdict'],
            'massless_gates':analysis.get('gates'),'massive_stage':analyze_massive(out),
            'scope':'Same-code archival checks, not independent physics replication, evaluator certification or interval arithmetic'}


def negatives(out):
    status={}
    with tempfile.TemporaryDirectory(prefix='grut_b1_verify_') as directory:
        trial=Path(directory)/'raw';shutil.copytree(out,trial)
        original=(trial/'analysis.json').read_text();bad=json.loads(original);bad['massless_pass']=not bad['massless_pass']
        (trial/'analysis.json').write_text(json.dumps(bad))
        try: validate(trial);raise AssertionError('False gate accepted')
        except ValueError: status['forged_acceptance_rejected']=True
        (trial/'analysis.json').write_text(original)
        original_stage=(trial/'massive_stage.json').read_text()
        (trial/'massive_stage.json').write_text(json.dumps({'status':'physical_continuum_certified'}))
        try: validate(trial);raise AssertionError('False stage accepted')
        except ValueError: status['forged_massive_stage_rejected']=True
        (trial/'massive_stage.json').write_text(original_stage)
        if not json.loads(original)['massless_pass']:
            (trial/'massive_N20.json').write_text('{}')
            try: validate(trial);raise AssertionError('Bypassed transfer barrier accepted')
            except ValueError: status['bypassed_transfer_barrier_rejected']=True
    return status


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw',type=Path,default=HERE/'raw')
    p.add_argument('--negative-controls',action='store_true');p.add_argument('--output',type=Path);a=p.parse_args()
    result=validate(a.raw)
    if a.negative_controls:result['negative_controls']=negatives(a.raw)
    if a.output:
        with a.output.open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(result,indent=2))
