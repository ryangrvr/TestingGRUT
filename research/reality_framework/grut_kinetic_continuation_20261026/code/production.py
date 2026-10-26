"""Predeclared continuation controls; no point or source retuning."""
from concurrent.futures import ProcessPoolExecutor,as_completed
import ks_common as K
from moment_controls import Control,save

def job(spec):
    label,R,formation,recycle=spec
    try:
        out=Control(R,formation,recycle=recycle,rtol=2e-9).run()
        return label,out
    except Exception as e:return label,dict(error=repr(e),complete=False,R=R,formation=formation,recycle=recycle)

def main():
    sources=[('release_only',0.),('mediators_only',K.R_MED['e1e-7']),
        ('e1e7_band',K.R_VEC['e1e-7_band']+K.R_MED['e1e-7']),
        ('e1e7_broad',K.R_VEC['e1e-7_broad']+K.R_MED['e1e-7']),
        ('e1e6_band',K.R_VEC['e1e-6_band']+K.R_MED['e1e-6'])]
    specs=[(name+'_'+form,R,form,form!='gas') for name,R in sources for form in ['gas','instant']]
    # Diagnostic source-strength samples, same formation recipes. R is total
    # extra conformal energy, with no assertion of a cosmological point.
    specs += [('diagnostic_R10_instant',10.,'instant',True),('diagnostic_R20_instant',20.,'instant',True),
              ('diagnostic_R90_instant',90.,'instant',True)]
    out={}
    with ProcessPoolExecutor(max_workers=3) as pool:
        tasks=[pool.submit(job,s) for s in specs]
        for task in as_completed(tasks):
            label,res=task.result();out[label]=res;save('MOMENT_CONTROLS.json',dict(cases=out,full_kinetic_complete=False))
            print('PRODUCTION',label,res.get('final',res.get('error')),res.get('stop'),res.get('seconds'),flush=True)

if __name__=='__main__':main()
