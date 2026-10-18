#!/usr/bin/env python3
"""Sequential, finite reproduction run; execute in a fresh package copy."""
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
for script in ['qball.py','results_additional_checks.py','reporter.py','verify.py','plot_regime.py']:
    print('Running '+script,flush=True)
    subprocess.run([sys.executable,str(ROOT/script)],check=True,cwd=ROOT.parent)
print('Completed all requested finite calculations.',flush=True)
