"""Add the small primordial-vector piece and evaluate the exact central ledger."""
import json
from concurrent.futures import ProcessPoolExecutor
from run import ROOT,write
from relic import run


def job(args):
    index,eps,R,ym=args
    return dict(candidate=index,**run(eps,R,ym,steps=1040))


def main():
    dec=json.loads((ROOT/'DECAY_RESULTS.json').read_text());vec=json.loads((ROOT/'VECTOR_RESULTS.json').read_text())
    thresholds=json.loads((ROOT/'THRESHOLD_RESULTS.json').read_text());results=[];jobs=[]
    for r in dec['rows']:
        v=next(x for x in vec['estimates'] if x['candidate']==r['candidate'])
        rbd=sum(v['corrected_R_range'])/2;pbd=sum(v['corrected_pairs_per_net_range'])/2
        combined={k:value for k,value in r.items() if k!='ensemble_rows'}
        combined.update(R_primordial_vector_range=v['corrected_R_range'],
            R_total_central_including_primordial=r['R_central']+rbd,
            pairs_primordial_vector_range=v['corrected_pairs_per_net_range'],
            pairs_total_central_including_primordial=r['pairs_total_central']+pbd,
            thresholds=[x for x in thresholds if x['candidate']==r['candidate']])
        results.append(combined)
        for ym in [1.,r['ymult_raw']]:jobs.append((r['candidate'],r['epsilon_net'],r['R_central']+rbd,ym))
    with ProcessPoolExecutor(max_workers=4) as pool:relics=list(pool.map(job,jobs))
    write(ROOT/'COMBINED_RESULTS.json',dict(rows=results,full_central_benchmark_relics=relics,
        decision='Both specified low-scale releases fail the supplied case-L benchmark. The conservative Y+Sigma ledger fails even with all heavy-cascade and primordial vectors omitted.',
        scope='Closure of the two published prepared releases under a stated late-number decay ledger and equilibrium lump proxy. Not a theorem for every continuous phase, modified history or full quantum kinetic model.'))
    print('COMBINED_COMPLETE',[(x['candidate'],x['ymult'],x['f_extra']) for x in relics],flush=True)


if __name__=='__main__':main()
