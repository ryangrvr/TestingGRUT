import sys, json, time
sys.path.insert(0, sys.argv[1])
import ks_common as K
from moment_controls import Control
from concurrent.futures import ProcessPoolExecutor
specs=[('e1e7_band_gas',K.R_VEC['e1e-7_band']+K.R_MED['e1e-7'],'gas',False),
       ('e1e6_band_gas',K.R_VEC['e1e-6_band']+K.R_MED['e1e-6'],'gas',False),
       ('diag_R5_gas',5.0,'gas',False),('diag_R7_gas',7.0,'gas',False)]
def job(s):
    t=time.time(); o=Control(s[1],s[2],recycle=s[3],rtol=2e-9).run(); return s[0],s[1],o.get('final'),o.get('stop'),time.time()-t
out={}
with ProcessPoolExecutor(max_workers=4) as p:
    for r in p.map(job,specs):
        out[r[0]]=dict(R=r[1],final=r[2],stop=r[3],sec=r[4]); print(r[0],r[1],json.dumps(r[2])[:400],r[3],round(r[4]),flush=True)
json.dump(out,open(sys.argv[2],'w'),indent=1,default=str)
