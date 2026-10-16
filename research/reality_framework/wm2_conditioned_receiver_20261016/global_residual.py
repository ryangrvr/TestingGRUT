"""S6-G1: preserve original stops, pay the same total residual budget."""
import hashlib
import json
from pathlib import Path
import campaign as P
import conditioned as C

BASE=Path(__file__).resolve().parent


def run():
    original_path=BASE/'RESULTS.json'
    original=json.loads(original_path.read_text())
    digest=hashlib.sha256(original_path.read_bytes()).hexdigest()
    result=dict(charter_commit='23de49378a64c16c405db87c7b1a703952db9813',
                original_S6_results_sha256=digest,original_S6_preserved=True,
                physical_law_changed=False,physical_test_changed=False,global_error_limits_changed=False,
                policy='cumulative residual <=.20; no per-step quota',attempts=[],families=[],
                all_size_formation_proved=False,roundoff_certified=False)
    seen={}
    for family in P.families():
        row=dict(**family,components=[])
        for root in (0,1):
            c=P.component(family,root)
            direct=P.key(c['counts'],c['values'],c['g'])
            mirrored=P.mirror(c['counts'],c['values'],c['g'])
            relation='new G1 trajectory'
            if direct in seen:
                index=seen[direct]; relation='identical component'
            elif mirrored in seen:
                index=seen[mirrored]; relation='S5 proved sign identity'
            else:
                index=len(result['attempts'])
                complete=[r for r in original['attempts'] if r['status']=='COMPLETE' and
                          P.key(r['nonroot_type_counts'],r['detuning_values'],r['g']) in (direct,mirrored)]
                if complete:
                    r=complete[0]; relation='completed original S6 component reused'
                else:
                    print('G1 START',index,family['N'],family['bar_delta'],c['counts'],flush=True)
                    e=C.Evolution(c['counts'],c['values'],64*family['N']**2)
                    r=e.run(4*family['N']**2,strict_step_residual=False,
                        progress=lambda v:print('G1 PROGRESS',index,json.dumps(v),flush=True))
                result['attempts'].append(r); seen[direct]=index
                print('G1 FINISH',index,r['status'],r['stop_reason'],r['formation_lower'],flush=True)
                P.save(BASE/'GLOBAL_RESULTS.json',result)
            row['components'].append(dict(**c,index=index,relation=relation))
        lower=full=missing=0.0
        for c in row['components']:
            r=result['attempts'][c['index']]
            if r['status']=='COMPLETE':
                lower+=c['weight']*r['formation_lower'];full+=c['weight']*r['full_assembly_lower']
            else:missing+=c['weight']
        row.update(status='COMPLETE' if missing==0 else 'UNRESOLVED',formation_lower=lower,
                   full_assembly_lower=full,unresolved_root_weight=missing,
                   gate_pass=bool(missing==0 and lower>=.5 and
                                 all(result['attempts'][c['index']]['gate_pass'] for c in row['components'])))
        result['families'].append(row);P.save(BASE/'GLOBAL_RESULTS.json',result)
        print('G1 FAMILY',family['N'],family['bar_delta'],row['status'],lower,flush=True)
    assert hashlib.sha256(original_path.read_bytes()).hexdigest()==digest
    P.save(BASE/'GLOBAL_RESULTS.json',result)
    return result


if __name__=='__main__':run()
