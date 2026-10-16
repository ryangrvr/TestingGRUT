"""Serial full reproduction or saved-evidence verification."""
from pathlib import Path
import argparse,os,subprocess,sys,tempfile
P=Path(__file__).resolve().parent


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--full',action='store_true');a=ap.parse_args()
    env=os.environ.copy();env['OPENBLAS_NUM_THREADS']='1';env['OMP_NUM_THREADS']='1'
    env.setdefault('GRUT_SHORTCUT_CACHE',str(Path(tempfile.gettempdir())/'grut_shortcut_cache'))
    env.setdefault('GRUT_SOURCE_CACHE',str(Path(tempfile.gettempdir())/'grut_host_cache'))
    def run(script,*args):subprocess.run([sys.executable,str(P/script),*map(str,args)],env=env,check=True)
    if a.full:
        for n,q,b in [(4,5,256),(5,6,256),(6,7,192)]:
            run('graphs.py','--n',n,'--q',q)
            run('neutral.py','--n',n,'--q',q,'--budget',b)
            args=['--n',n,'--q',q,'--budget',b,'--level',8 if n==6 else 6]
            run('host_response.py',*args)
            if n==6:run('host_response.py',*args,'--close-initial')
        run('graphs.py','--sample');run('checks.py');run('response_check.py');run('path_b.py');run('report.py');run('plot.py')
        run('verify.py','--refresh')
    else:run('verify.py')
