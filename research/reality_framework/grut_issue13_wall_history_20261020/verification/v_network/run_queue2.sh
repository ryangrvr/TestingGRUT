#!/bin/bash
# Second verifier queue: statistically independent seeds (the verifier's seeds 1,2 turned out to
# reproduce the module's seed-1,2 realizations because both use default_rng(seed).standard_normal).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
R=results/runs
while pgrep -f "run_queue.sh" > /dev/null; do sleep 5; done
( python3 run_network.py --d 2 --N 2048 --dx 0.5 --alpha 2 --beta 0 --tau0 1 --tau_end 256 --s_c 0.75 --std 0.1 --seed 11 --nout 36 --out $R/prs2d_a2_N2048_s11.json > $R/prs2d_a2_N2048_s11.log 2>&1 ) &
( python3 run_network.py --d 2 --N 2048 --dx 0.5 --alpha 2 --beta 0 --tau0 1 --tau_end 256 --s_c 0.75 --std 0.1 --seed 12 --nout 36 --out $R/prs2d_a2_N2048_s12.json > $R/prs2d_a2_N2048_s12.log 2>&1
  python3 run_network.py --d 2 --N 1024 --dx 0.5 --alpha 2 --beta 0 --tau0 1 --tau_end 256 --s_c 0.75 --std 0.1 --seed 11 --p 0.35 --nout 40 --out $R/prs2d_a2_N1024_p035_s11.json > $R/prs2d_a2_N1024_p035_s11.log 2>&1 ) &
wait
echo QUEUE2_DONE
