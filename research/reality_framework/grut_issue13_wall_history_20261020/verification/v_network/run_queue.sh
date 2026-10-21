#!/bin/bash
# Verifier job queue: two single-thread workers (2 CPU threads total).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd "$(dirname "$0")"
R=results/runs
# wait for the 2D 2048^2 runs
while pgrep -f "prs2d_a2_N2048" > /dev/null; do sleep 5; done
# worker A: 3D 160^3 seeds 1,2 (alpha=3, beta=0)
( python3 run_network.py --d 3 --N 160 --dx 0.5 --alpha 3 --beta 0 --tau0 1 --tau_end 40 --s_c 0.75 --std 0.1 --seed 1 --dt_frac 0.25 --nout 32 --out $R/prs3d_a3_N160_s1.json > $R/prs3d_a3_N160_s1.log 2>&1
  python3 run_network.py --d 3 --N 160 --dx 0.5 --alpha 3 --beta 0 --tau0 1 --tau_end 40 --s_c 0.75 --std 0.1 --seed 2 --dt_frac 0.25 --nout 32 --out $R/prs3d_a3_N160_s2.json > $R/prs3d_a3_N160_s2.log 2>&1 ) &
# worker B: biased 2D runs and bubbles
( python3 run_network.py --d 2 --N 1024 --dx 0.5 --alpha 2 --beta 0 --tau0 1 --tau_end 256 --s_c 0.75 --std 0.1 --seed 1 --p 0.35 --nout 40 --out $R/prs2d_a2_N1024_p035_s1.json > $R/prs2d_a2_N1024_p035_s1.log 2>&1
  python3 run_network.py --d 2 --N 1024 --dx 0.5 --alpha 2 --beta 0 --tau0 1 --tau_end 256 --s_c 0.75 --std 0.1 --seed 2 --p 0.35 --nout 40 --out $R/prs2d_a2_N1024_p035_s2.json > $R/prs2d_a2_N1024_p035_s2.log 2>&1
  python3 run_network.py --d 2 --N 1024 --dx 0.5 --alpha 2 --beta 0 --tau0 1 --tau_end 256 --s_c 0.75 --std 0.1 --seed 1 --nout 40 --out $R/prs2d_a2_N1024_p050_s1.json > $R/prs2d_a2_N1024_p050_s1.log 2>&1
  python3 run_network.py --d 2 --N 1024 --dx 0.5 --alpha 2 --beta 0 --tau0 1 --tau_end 256 --s_c 0.75 --std 0.1 --seed 1 --p 0.45 --nout 40 --out $R/prs2d_a2_N1024_p045_s1.json > $R/prs2d_a2_N1024_p045_s1.log 2>&1
  python3 bubble_v.py > results/bubble_v.log 2>&1 ) &
wait
echo QUEUE_DONE
