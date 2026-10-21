#!/usr/bin/env bash
# Usage: run_case.sh NAME BINARY [key value ...]
# Runs the compiled (unmodified, packet-hash e62701e5...) solver with OMP_NUM_THREADS=2 and --threads 2,
# writing outputs to $SCR/runs/NAME.* ; records wall time and command to $SCR/runs/NAME.cmd.json
set -euo pipefail
SCR=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/rerun
name=$1; shift; bin=$1; shift
mkdir -p "$SCR/runs" "$SCR/states"
start=$(date +%s.%N)
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 "$bin" "$@" --threads 2 --out "$SCR/runs/$name" > "$SCR/runs/$name.stdout" 2>&1
end=$(date +%s.%N)
python3 -I -c "import json,sys; print(json.dumps({'name':sys.argv[1],'command':sys.argv[4:],'wall_seconds':float(sys.argv[3])-float(sys.argv[2])}))" "$name" "$start" "$end" "$bin" "$@" --threads 2 > "$SCR/runs/$name.cmd.json"
cat "$SCR/runs/$name.stdout"
