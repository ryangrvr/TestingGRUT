#!/usr/bin/env bash
# Queue 3: (i) intermediate-time final states for the independent topology reconstruction (same runs,
# stopped early: bitwise identical evolution up to T), (ii) extra seeds at reduced cost (N=128, dx=0.4, L=51.2).
set -u
SCR=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/rerun; A=/home/user/TestingGRUT/research/reality_framework/grut_issue13_3d_review_20261021/audit/rerun
E=$SCR/build/expanding; G=$SCR/build/gaussian
C="--dt 0.06 --taui 40.0 --kind charge --noise 1.0 --sample 1.2 --top-sample 4.8 --topology 1 --kappa 2.0"
$A/run_case.sh mid_gaussian_s202_T24 $G --N 192 --dx 0.4 --dt 0.06 --T 24 --noise-file $SCR/noise/g192_s202.f64 --mirror 1 --state $SCR/states/mid_gaussian_s202_T24.f64
$A/run_case.sh mid_main_s101_T27.6 $E --N 192 --dx 0.4 --T 27.6 $C --modes $SCR/modes/modes_L76.8_s101.txt --state $SCR/states/mid_main_s101_T27.6.f64
for s in 303 404 505; do
 $A/run_case.sh seed128dx04_gaussian_s$s $G --N 128 --dx 0.4 --dt 0.06 --T 48 --noise-file $SCR/noise/g128_dx04_s$s.f64 --mirror 1
 $A/run_case.sh seed128dx04_main_s$s $E --N 128 --dx 0.4 --T 48 $C --modes $SCR/modes/modes_L51.2_s$s.txt
done
$A/run_case.sh red128dx04_main_s202 $E --N 128 --dx 0.4 --T 48 $C --modes $SCR/modes/modes_L51.2_s202.txt
