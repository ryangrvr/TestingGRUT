#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=2
export MPLBACKEND=Agg
case "${1:---verify}" in
  --verify)
    python code/check_inputs.py
    python code/verify.py
    python code/manifest.py --check
    ;;
  --all)
    python code/check_inputs.py
    python code/scale_checks.py
    python code/phase_holder.py
    python code/phase_cubic.py
    python code/phase_fluctuations.py
    python code/holder_evolution.py
    python code/holder_evolution.py --end 800 --rtol 1e-8 --case original_prepared_no_holder --out reference_long
    python code/holder_evolution.py --end 800 --rtol 1e-8 --case local_equilibrium_holder_4e8 --out holder_long
    python code/holder_sensitivity.py
    python code/free_stream.py
    python code/collision_rate.py
    python code/psd_audit.py
    python code/figure.py
    python code/verify.py
    ;;
  *) echo 'Usage: ./reproduce.sh [--verify|--all]' >&2; exit 2;;
esac
