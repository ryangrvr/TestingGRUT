"""Explicit replay of pinned reviewer code; not an independent relic solver.

argv[1] is the extracted grut_formation_review_20261025 directory at f7d641a.
Never writes inside that archive. Uses the author's computed finite-band R.
"""
from pathlib import Path
import sys,json,time
from multiprocessing import Pool
from background import ROOT,write

sys.dont_write_bytecode=True


def job(args):
    directory,epsilon,injection,kind,steps=args
    sys.path.insert(0,str(Path(directory)/'tracks/ls_threshold/code'))
    import lsrelic as LR
    LR.UC=.8611951308217225
    LR.EC=LR.UC if kind!='original_proxy' else .9233161353726057
    if kind=='original_proxy':
        result=LR.proxy_run(epsilon,injection,1.,N=steps)
    else:
        result=LR.ke_run(epsilon,injection,1.,N=steps,antilumps=True,s_channel=True)
    result['kind']=kind;result['implementation']='Pinned reviewer f7d641a lsrelic.py, replayed unmodified except supplied documented global UC=EC.'
    return result


def main():
    directory=sys.argv[1]
    b=json.loads((ROOT/'BACKGROUND_RESULTS.json').read_text())
    R=json.loads((ROOT/'BUDGET_RESULTS.json').read_text())['rows'][0]['R_vector_band_only']
    eps=abs(b['late']['epsilon'])
    jobs=[(directory,eps,R,'ke_plus_s',600),(directory,eps,R,'ke_plus_s_refined',1200),
          (directory,eps,R,'original_proxy',600)]
    rows=[]
    with Pool(3,maxtasksperchild=1) as pool:
        for result in pool.imap(job,jobs):
            rows.append(result)
            write(ROOT/'RELIC_REPLAY.json',dict(rows=rows,scope=__doc__,review_archive=directory))
            print('RELIC_REPLAY',result['kind'],result['f'],result.get('energy_work_cumulative_residual'),flush=True)


if __name__=='__main__':main()
