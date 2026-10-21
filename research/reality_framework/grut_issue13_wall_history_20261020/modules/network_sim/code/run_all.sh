#!/bin/bash
# One-command reproduction of module network_sim (about 2.5 CPU-hours on 2 single-threaded workers).
# Existing results/runs/*.json are skipped by run_networks.py; delete results/runs to recompute everything.
set -e
cd "$(dirname "$0")"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 validate.py
python3 run_networks.py run            # 128 lattice jobs, 2 worker processes, each job capped at 24 min
python3 bubbles.py                     # isolated-bubble control (single process)
python3 scales.py
python3 analyze.py
