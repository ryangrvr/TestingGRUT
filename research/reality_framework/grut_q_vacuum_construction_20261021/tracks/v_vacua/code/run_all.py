"""Reproduce every result of verifier v_vacua.

    cd code && sh v3_build_smdr.sh && python3 run_all.py

v3_build_smdr.sh downloads SMDR v1.3 (Martin-Robertson, GPL) into the scratchpad and compiles smdr_driver.c against
it (about 1 minute). v3a (16 SMDR fits, about 2-3 minutes on 4 cores) is skipped if its cache exists; delete
../results/v3a_smdr_fit_grid.json to redo it. Total about 4 minutes.
"""
import os
import runpy

HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
scripts = ['v1_qsector.py', 'v2_nucleation.py']
if not os.path.exists(os.path.join(HERE, '..', 'results', 'v3a_smdr_fit_grid.json')):
    scripts.append('v3a_smdr_fit_grid.py')
scripts += ['v3b_sm_check.py', 'v3c_g3_qcd.py', 'v4_verdicts.py']
for f in scripts:
    print('=== running', f)
    runpy.run_path(os.path.join(HERE, f), run_name='__main__')
