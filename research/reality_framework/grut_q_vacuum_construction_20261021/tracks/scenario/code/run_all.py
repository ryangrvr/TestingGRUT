"""Reproduce every result of the scenario track.

    cd research/reality_framework/grut_q_vacuum_construction_20261021/tracks/scenario
    python3 -I code/run_all.py            # about 6-10 min on 4 cores

Downloads (inputs, never stored in the repository) go to the scratchpad directory DATA (edit if needed):
  Saikawa-Shirai SM EoS tables (arXiv:1803.01038 ancillary file + authors' page), DESI DR2 BAO Gaussian likelihood
  (CobayaSampler/bao_data), Pantheon+ distances and STAT+SYS covariance (PantheonPlusSH0ES/DataRelease).
Outputs: results/symbolic_checks.json, reduction_validation.json, early_relaxation.json, realistic_history.json,
variants_comparison.json.
"""
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, '..', 'results')
DATA = ('/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/qtheory/'
        'scenario_dl')

sys.path.insert(0, HERE)
import sm_eos  # noqa: E402
import s5_variants_data  # noqa: E402

os.makedirs(RES, exist_ok=True)
sm_eos.fetch(os.path.join(DATA, 'ss'))
s5_variants_data.fetch(DATA)
for script, args in (('s1_symbolic_checks.py', [RES]),
                     ('s2_reduction_validation.py', [RES]),
                     ('s3_early_relaxation.py', [RES]),
                     ('s4_realistic_history.py', [RES, os.path.join(DATA, 'ss')]),
                     ('s5_variants_data.py', [RES, DATA, os.path.join(DATA, 'ss')])):
    t0 = time.time()
    print('==>', script, flush=True)
    subprocess.run([sys.executable, '-I', os.path.join(HERE, script)] + args, check=True)
    print('    done in %.0f s' % (time.time() - t0), flush=True)
