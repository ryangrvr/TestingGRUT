#!/usr/bin/env bash
# Queue 1: full-size reruns of two archived rows (sequential, 2 threads each); final states kept in scratch
# for an independent topology reconstruction (state files are not archived in the audit dir: 340 MB each).
set -u
SCR=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/rerun; A=/home/user/TestingGRUT/research/reality_framework/grut_issue13_3d_review_20261021/audit/rerun
$A/run_case.sh full_main_s101 $SCR/build/expanding --N 192 --dx 0.4 --dt 0.06 --T 48 --taui 40.0 --kappa 2.0 --kind charge --noise 1.0 --sample 1.2 --top-sample 4.8 --topology 1 --modes $SCR/modes/modes_L76.8_s101.txt --state $SCR/states/full_main_s101.f64
$A/run_case.sh full_gaussian_s202 $SCR/build/gaussian --N 192 --dx 0.4 --dt 0.06 --T 48 --noise-file $SCR/noise/g192_s202.f64 --mirror 1 --state $SCR/states/full_gaussian_s202.f64
