"""One fresh, measured S8 primary component. Setup belongs to its budget."""
import time
START=time.monotonic();CPU=time.process_time()
import argparse
import json
import resource
from pathlib import Path
from mps_global import Receiver


def main():
    p=argparse.ArgumentParser();p.add_argument('--counts',required=True)
    p.add_argument('--output',required=True);p.add_argument('--seconds',type=float,default=900)
    a=p.parse_args();counts=list(map(int,a.counts.split(',')));N=sum(counts)+1
    def progress(x):print(json.dumps(x),flush=True)
    e=Receiver(counts,[-.25,.25],64*N*N)
    r=e.run(4*N*N,max_seconds=a.seconds,started=START,progress=progress)
    r.update(method='S8-G1',whole_worker_wall_seconds=time.monotonic()-START,
        whole_worker_cpu_seconds=time.process_time()-CPU,
        peak_resident_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        wall_limit_seconds=a.seconds,wall_overshoot_seconds=max(0,time.monotonic()-START-a.seconds),
        imports_setup_and_operator_preprocessing_included=True,
        retained_tensor_count_is_not_total_memory=True)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)


if __name__=='__main__':main()
