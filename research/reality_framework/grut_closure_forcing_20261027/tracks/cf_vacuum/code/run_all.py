"""Run the cf_vacuum track: python3 -B run_all.py  (~8 min; SMDR must be built, see c1_fit_grid.py)."""
import subprocess, sys
for s in ('c1_fit_grid.py', 'c2_sm_relation.py', 'c3_threshold.py', 'c4_hidden_and_costs.py'):
    print('==', s, flush=True)
    subprocess.run([sys.executable, '-B', s], check=True)
