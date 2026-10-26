"""Optional replay of one pinned reviewer benchmark, explicitly NOT own code.

Pass an archive extraction root of f1aa7df. The reference source is read
only. Main production uses thermo.py/relic.py, not these imported routines.
"""
import argparse,json,sys,time
from pathlib import Path
from run import ROOT,write


def main():
    p=argparse.ArgumentParser();p.add_argument('review_root');a=p.parse_args()
    path=Path(a.review_root)/'research/reality_framework/grut_formation_review_20261025/tracks/sp_yledger/code'
    sys.dont_write_bytecode=True;sys.path.insert(0,str(path))
    import thermal_relic
    data=json.loads((ROOT/'DECAY_RESULTS.json').read_text())['rows'][0]
    t=time.perf_counter()
    result=thermal_relic.run_inj(data['epsilon_net'],'L',data['R_lower'],data['ymult_raw'])
    write(ROOT/'inputs/REVIEW_BENCHMARK_REPLAY.json',dict(
        commit='f1aa7df968cb3dc90ca32f4f32aeb0a819ddfce4',
        scope='Executed reviewer adaptive-integral/260-step benchmark for ONE low-scale lower-ledger case; not an outside replication or a replay of all review tracks.',
        source='tracks/sp_yledger/code/thermal_relic.py, using tracks/gp_leftovers/code/equilibrium.py',
        result=result,runtime_seconds=time.perf_counter()-t))
    print('REVIEW_REFERENCE',result['f_extra'],time.perf_counter()-t,flush=True)


if __name__=='__main__':main()
