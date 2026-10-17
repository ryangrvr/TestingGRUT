"""One fresh measured worker per frozen primary component or S5 baseline."""
import time
START=time.monotonic();CPU=time.process_time()
import argparse
import json
import resource
from pathlib import Path
import solver
from layout import S5


def main():
    p=argparse.ArgumentParser();p.add_argument('--method',choices=['s7','s5'],required=True)
    p.add_argument('--counts',required=True);p.add_argument('--output',required=True)
    p.add_argument('--seconds',type=float,default=600)
    a=p.parse_args();counts=list(map(int,a.counts.split(',')));n=sum(counts)+1
    def progress(x):print(json.dumps(x),flush=True)
    if a.method=='s7':
        evolution=solver.Evolution(counts,[-.25,.25],64*n*n)
        result=evolution.run(4*n*n,max_seconds=a.seconds,started=START,progress=progress)
    else:
        evolution=S5.Evolution(counts,[-.25,.25],64*n*n,.9,rank_cap=24,state_cap=2_000_000,degree=8)
        remaining=max(.001,a.seconds-(time.monotonic()-START))
        result=evolution.run(4*n*n,strang_budget=.20,compression_budget=.35,max_seconds=remaining,progress=progress)
        error=result['analytic_state_error_account'];estimate=result['macro_probability_estimate']
        from enclosures import physical_probability
        result['physical_probability_intervals']=(
            dict(macro=physical_probability((estimate,estimate),error),
                 full=physical_probability((result['full_assembly_probability_estimate'],)*2,error))
            if result['status']=='COMPLETE' else None)
        result['arithmetic_status']='immutable S5 ordinary double arithmetic, no roundoff enclosure'
    result.update(method=a.method,whole_worker_wall_seconds=time.monotonic()-START,
        whole_worker_cpu_seconds=time.process_time()-CPU,
        peak_resident_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        wall_limit_seconds=a.seconds,wall_overshoot_seconds=max(0,time.monotonic()-START-a.seconds),
        receiver_setup_and_imports_included=True)
    path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)


if __name__=='__main__':main()
