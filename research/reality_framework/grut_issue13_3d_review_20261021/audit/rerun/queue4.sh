#!/usr/bin/env bash
# Queue 4: large box L=115.2 (N=288, dx=0.4): rerun of archived large_s101 row, then the two rows Codex left
# unexecuted at its cost cap (large_s202, large_s303) -- new evidence, same frozen source and arguments.
set -u
SCR=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/rerun; A=/home/user/TestingGRUT/research/reality_framework/grut_issue13_3d_review_20261021/audit/rerun
E=$SCR/build/expanding
C="--dt 0.06 --T 48 --taui 40.0 --kind charge --noise 1.0 --sample 1.2 --top-sample 4.8 --topology 1 --kappa 2.0"
$A/run_case.sh full_large_s101 $E --N 288 --dx 0.4 $C --modes $SCR/modes/modes_L115.2_s101.txt
$A/run_case.sh new_large_s202 $E --N 288 --dx 0.4 $C --modes $SCR/modes/modes_L115.2_s202.txt
$A/run_case.sh new_large_s303 $E --N 288 --dx 0.4 $C --modes $SCR/modes/modes_L115.2_s303.txt
