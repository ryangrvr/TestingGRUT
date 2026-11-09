"""Post-run N22 integrity, arithmetic and checkpoint audit (not production)."""
import argparse,hashlib,itertools,json,math,sys
from pathlib import Path
import numpy as np

def tv(p,q):return .5*sum(abs(p[s]-q[s]) for s in p)

def audit(root,out):
    area=root/'research/reality_framework/grut_n22_refinement_20261109'
    sys.path.insert(0,str(area))
    from shared import B2,card,reconstruction,run_b2_point
    from b3_analyze import analyze
    from lattice import schwinger
    cfg=card();g=cfg['gates'];saved=json.loads((out/'assessment.json').read_text())
    prov=json.loads((out/'provenance.json').read_text());errors=[]
    for name,digest in prov['source_sha256'].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest:errors.append('source_hash:'+name)
    if not saved['complete']:
        result={'passed':False,'scope':'incomplete experiment','errors':errors,'assessment':saved}
    else:
        points={label:json.loads((out/(label+'.json')).read_text()) for label in ('N22','N22_refined')}
        old={f'N{n}':json.loads((B2/'raw'/f'N{n}.json').read_text()) for n in (10,14,18,20)}
        maximum_norm=0.;metric_error=0.;signed_error=0.;checkpoint={};guard={}
        protocol_keys={','.join(f'{t:g}' for t in p) for p in cfg['protocols']}
        for label,d in points.items():
            groups=list(d['histories'].values())+[d['references'][r]['histories'] for r in ('target','B1','B2')]
            for group in groups:
                if set(group)!=protocol_keys:errors.append('protocol_keys:'+label)
                for h,p in group.items():
                    if set(p)!={''.join(s) for s in itertools.product('+-',repeat=h.count(',')+1)}:errors.append('outcomes:'+label+'/'+h)
                    maximum_norm=max(maximum_norm,abs(sum(p.values())-1))
                    if not all(math.isfinite(v) and v>=-g['normalization'] for v in p.values()):errors.append('invalid_probability:'+label+'/'+h)
            comparison=saved['comparisons'][label]
            for h,p in d['histories']['B2'].items():
                matched=d['references']['B2']['histories'][h];target=d['references']['target']['histories'][h]
                independent={'matched_tv':tv(p,matched),'geometry_tv':tv(matched,target),'target_tv':tv(p,target),
                    'triangle_plus_reserve':tv(p,matched)+tv(matched,target)+g['numerical_reserve']}
                for key,value in independent.items():metric_error=max(metric_error,abs(value-comparison[key][h]))
            data=np.load(out/(label+'.npz'),allow_pickle=False);n=d['N'];a=d['a']
            ham,electric,states=schwinger.hamiltonian(n,a**-2,-a/8);hmat=ham*a/2;ground=data['ground']
            _,_,weights,boundary,_=reconstruction(cfg,n);field=electric@weights+(boundary.sum()-weights@((-1.)**np.arange(n-1)))/4
            first=.5*(1+np.vdot(ground,np.exp(1j*cfg['coupling']*field)*ground).real)
            checkpoint[label]={'ground_norm_error':float(abs(np.vdot(ground,ground).real-1)),
              'ground_residual':float(np.linalg.norm(hmat@ground-data['eigenvalues'][0]*ground)),
              'first_event_error':float(abs(first-d['histories']['B2']['0']['+'])),
              'weight_error':float(np.max(abs(data['weights']-weights))),
              'boundary_error':float(np.max(abs(data['boundary']-boundary))),
              'basis_identical':bool(np.array_equal(states,data['states']))}
            if not checkpoint[label]['basis_identical'] or max(v for k,v in checkpoint[label].items() if k!='basis_identical')>1e-9:errors.append('checkpoint:'+label)
            del ham,hmat,electric,states,ground,data
        for label,d in {**old,**points}.items():
            for h,p in d['histories']['B2'].items():
                for s,value in p.items():
                    delta=value-d['references']['B2']['histories'][h][s]
                    signed_error=max(signed_error,abs(delta-saved['signed_vectors'][label]['delta'][h][s]),
                        abs(delta/d['a']**2-saved['signed_vectors'][label]['delta_over_a2'][h][s]))
        for h in protocol_keys:
            guard[h]=g['lattice_to_matched']-max(saved['comparisons'][label]['matched_tv'][h] for label in points)-g['numerical_reserve']
        if guard!=saved['component_margin_after_reserve']:errors.append('guard_arithmetic')
        predicted={};forecast_checks=[]
        for label,d in points.items():
            for h,p in d['histories']['B2'].items():
                source={s:(old['N20']['histories']['B2'][h][s]-old['N20']['references']['B2']['histories'][h][s])/old['N20']['a']**2 for s in p}
                current={s:(p[s]-d['references']['B2']['histories'][h][s])/d['a']**2 for s in p}
                ratio=sum(abs(current[s]-source[s]) for s in p)/sum(abs(v) for v in source.values())
                if abs(ratio-saved['forecast'][label][h]['relative_signed_vector_change'])>1e-13:errors.append('forecast_ratio:'+label+'/'+h)
                forecast_checks.append(ratio<=cfg['forecast']['relative_signed_vector_tolerance'])
        if all(forecast_checks)!=saved['forecast_accepted']:errors.append('forecast_decision')
        recomputed=analyze(out)
        if recomputed!=saved:errors.append('frozen_assessment_replay')
        if maximum_norm>g['normalization'] or metric_error>1e-14 or signed_error>1e-14:errors.append('arithmetic')
        if any('massive' in p.name.lower() for p in out.glob('*.npz')):errors.append('unauthorized_massive_checkpoint')
        for label in points:
            r=json.loads((out/(label+'_resources.json')).read_text())
            if r['returncode']!=0 or r['peak_private_bytes']>cfg['resources']['process_memory_bytes'] or r['wall_seconds']>cfg['resources']['point_wall_seconds']+2:errors.append('resources:'+label)
        damaged=json.loads(json.dumps(saved));damaged['nominal_pass']=not saved['nominal_pass']
        bad_p=dict(points['N22']['histories']['B2']['0']);bad_p['+']+=.01
        negative={'forged_verdict_rejected':damaged!=recomputed,'corrupt_normalization_rejected':abs(sum(bad_p.values())-1)>g['normalization']}
        result={'passed':not errors and all(negative.values()),'errors':errors,'max_normalization_error':maximum_norm,
          'max_independent_metric_error':metric_error,'max_signed_vector_error':signed_error,'source_hash_count':len(prov['source_sha256']),
          'checkpoint_checks':checkpoint,'negative_controls':negative,'assessment_replayed_identically':recomputed==saved,
          'scope':'Independent arithmetic plus stored-state and source audit; Hamiltonian/weights reused. No independent full-size physical solver or rigorous ground-state certificate.'}
    with (out/'audit.json').open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();audit(a.root.resolve(),a.output.resolve())
