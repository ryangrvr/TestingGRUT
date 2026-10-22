"""Reproduce every result of track `vacua`:  cd code && python3 run_all.py   (about 1 minute)."""
import runpy, os
HERE = os.path.dirname(os.path.abspath(__file__))
for f in ('vacua_q_sector.py', 'nucleation.py', 'sm_criticality.py'):
    print('=== running', f)
    runpy.run_path(os.path.join(HERE, f), run_name='__main__')
