"""Serial reproduction: numerical caches are disposable, full-H bounds are not.

Default: verify published outputs. --full: regenerate simulations/controls
and then verify. Large jobs are deliberately serial to fit in ~10 GB RAM.
"""
import argparse,os,subprocess,sys,tempfile
from pathlib import Path
P=Path(__file__).resolve().parent


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true');a=ap.parse_args()
    env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
    env.setdefault('GRUT_SHORTCUT_CACHE',str(Path(tempfile.gettempdir())/'grut_shortcut_cache'))
    env.setdefault('GRUT_SOURCE_CACHE',str(Path(tempfile.gettempdir())/'grut_source_cache'))
    def run(script,*args):
        subprocess.run([sys.executable,str(P/script),*map(str,args)],env=env,check=True)
    if a.full:
        for n,q,b,wait in [(4,5,256,5),(5,6,256,5),(5,6,96,'2,5'),(5,10,64,2)]:
            run('neutral_prepare.py','--n',n,'--q',q,'--budget',b,'--waits',wait,'--dt',1/64)
        run('neutral_krylov.py','--n',6,'--q',7,'--budget',192,'--columns',128,'--waits',5)
        for n,q,b in [(4,5,256),(5,6,256),(6,7,192)]:
            run('geometry_scan.py','--n',n,'--q',q,'--budget',b)
            run('refined_record.py','--n',n,'--q',q,'--budget',b)
        for n,q,b,L in [(4,5,256,12),(5,6,256,10)]:
            run('response_run.py','--n',n,'--q',q,'--budget',b,'--wait',5,'--level',L,'--lags','.0625,.125,.25,.5')
        for n,q,b,L in [(4,5,256,12),(5,6,256,10),(6,7,192,8)]:
            run('response_run.py','--n',n,'--q',q,'--budget',b,'--wait',5,'--level',L,
                '--chord',n//2,'--lags','.0625,.125,.25,.5')
        run('response_run.py','--n',5,'--q',10,'--budget',64,'--wait',2,'--level',6,'--lags','.125,.25')
        run('checks.py');run('refined_checks.py');run('krylov_check.py')
        run('assess.py');run('plot.py')
    if a.full:run('verify.py','--refresh')
    else:run('verify.py')


if __name__=='__main__':main()
