"""Two original-C release controls plus all 96 published endpoint projections."""
import csv,json
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np
from build import build,WORK
from run import ROOT,evolve,write
from modes import project,vacuum


def job(spec):
    label,args=spec;r=evolve(**args)
    write(WORK/'releases'/f'{label}.json',r)
    return label,r


def main():
    build();source=json.loads((ROOT/'inputs/PARENT_CASES.json').read_text())
    rows=[];groups={}
    for r in csv.DictReader((ROOT/'inputs/PARENT_ENDPOINTS.csv').open()):
        index=int(r['ensemble_id'].split('_')[1]);state=[float(r[f'final_state_{i}']) for i in range(12)]
        x=project(float(r['H_GeV']),1000.,state)
        # Window X charge/energy are the published ensemble statistics; mode
        # actions come from each canonical endpoint, not a nominal trajectory.
        x.update(candidate=index,member=int(r['member']),NX_window=float(r['NX']),EX_window=float(r['EX']),eps_window=float(r['epsX']))
        rows.append(x)
    for index in [10,17]:
        a=[r for r in rows if r['candidate']==index]
        keys=['NX_window','EX_window','NYH','NYL','NSIG','JY','ESIG','ESANG']
        g={k:float(np.mean([r[k] for r in a])) for k in keys}
        for k in ['NYH','NYL','NSIG']:g[k+'_SEM']=float(np.std([r[k] for r in a],ddof=1)/np.sqrt(len(a)))
        g.update(source[index.__str__()]);g['member_count']=len(a);groups[str(index)]=g
    write(ROOT/'MODE_RESULTS.json',dict(vacuum=vacuum(),ensemble_means=groups,endpoint_rows=rows,
        scope='Own real-matrix projection of every published canonical endpoint. Parent window X means retained; particle interpretation requires a late quadratic regime.'))
    jobs=[]
    for index in [10,17]:
        s=source[str(index)];base=dict(H=s['H_GeV'],cH=s['cH'],delta=s['delta'])
        jobs.extend([(f'{index}_long',base|dict(end=3000,rtol=2e-12,atol=2e-15,ns=160)),
                     (f'{index}_tight',base|dict(end=1000,rtol=2e-14,atol=2e-16,ns=160))])
    result={}
    with ProcessPoolExecutor(max_workers=4) as pool:
        for f in as_completed([pool.submit(job,j) for j in jobs]):
            label,r=f.result();q=project(r['H_GeV'],r['a_end'],r['final_canonical_state'])
            last=r['history'][-1]
            diffs={k:abs(q[k]-last[k])/max(abs(q[k]),1e-30) for k in ['NYH','NYL','NSIG','ESIG']}
            result[label]=dict(H_GeV=r['H_GeV'],cH=r['cH'],delta=r['delta'],a_end=r['a_end'],rtol=r['rtol'],
                               windows=r['windows'],controls=r['controls'],final_projection=q,
                               C_vs_real_matrix_endpoint_relative_error=diffs)
            write(ROOT/'RELEASE_RESULTS.json',result)
            print('RELEASE',label,'eps',r['windows']['600-1000']['epsX'],'Nsigma',r['windows']['600-1000']['NSIG'],flush=True)
    print('RELEASE_CONTROLS_COMPLETE',flush=True)


if __name__=='__main__':main()
