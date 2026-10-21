#!/usr/bin/env bash
# Queue 2: reduced-cost reruns (what changed is in the run name and stated in NOTES.md), then further
# full-size rows. Sequential, OMP_NUM_THREADS=2 each. Common: dt=0.06, T=48 (a: 1 -> 2.2), rate 1/40, kappa 2.
set -u
SCR=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/rerun; A=/home/user/TestingGRUT/research/reality_framework/grut_issue13_3d_review_20261021/audit/rerun
E=$SCR/build/expanding; G=$SCR/build/gaussian
C="--dt 0.06 --T 48 --taui 40.0 --kind charge --noise 1.0 --sample 1.2 --top-sample 4.8 --topology 1"
# (a) same box L=76.8, coarser grid: N=128 dx=0.6 ; N=96 dx=0.8. Same sparse modes; Gaussian = Fourier restriction of the N=192 draw.
$A/run_case.sh red128dx06_main_s101 $E --N 128 --dx 0.6 $C --kappa 2.0 --modes $SCR/modes/modes_L76.8_s101.txt --state $SCR/states/red128dx06_main_s101.f64
$A/run_case.sh red128dx06_gaussian_s202 $G --N 128 --dx 0.6 --dt 0.06 --T 48 --noise-file $SCR/noise/g128r_s202.f64 --mirror 1 --state $SCR/states/red128dx06_gaussian_s202.f64
# (b) same dx=0.4 (same L/N), smaller box N=128 L=51.2: new sparse mode draw (L-dependent), native N=128 Gaussian draw.
$A/run_case.sh red128dx04_main_s101 $E --N 128 --dx 0.4 $C --kappa 2.0 --modes $SCR/modes/modes_L51.2_s101.txt --state $SCR/states/red128dx04_main_s101.f64
$A/run_case.sh red128dx04_gaussian_s202 $G --N 128 --dx 0.4 --dt 0.06 --T 48 --noise-file $SCR/noise/g128_dx04_s202.f64 --mirror 1 --state $SCR/states/red128dx04_gaussian_s202.f64
$A/run_case.sh red96dx08_main_s101 $E --N 96 --dx 0.8 $C --kappa 2.0 --modes $SCR/modes/modes_L76.8_s101.txt
$A/run_case.sh red96dx08_gaussian_s202 $G --N 96 --dx 0.8 --dt 0.06 --T 48 --noise-file $SCR/noise/g96r_s202.f64 --mirror 1
# (c) further full-size rows of the archived table
$A/run_case.sh full_gaussian_s101 $G --N 192 --dx 0.4 --dt 0.06 --T 48 --noise-file $SCR/noise/g192_s101.f64 --mirror 1 --state $SCR/states/full_gaussian_s101.f64
$A/run_case.sh full_main_s202 $E --N 192 --dx 0.4 $C --kappa 2.0 --modes $SCR/modes/modes_L76.8_s202.txt
$A/run_case.sh full_cold_s101 $E --N 192 --dx 0.4 --dt 0.06 --T 48 --taui 40.0 --kind cold --noise 1.0 --sample 1.2 --top-sample 4.8 --topology 1 --kappa 0.9 --modes $SCR/modes/modes_L76.8_s101.txt
$A/run_case.sh full_kappa4_s101 $E --N 192 --dx 0.4 $C --kappa 4 --modes $SCR/modes/modes_L76.8_s101.txt
$A/run_case.sh full_main_s303 $E --N 192 --dx 0.4 $C --kappa 2.0 --modes $SCR/modes/modes_L76.8_s303.txt
$A/run_case.sh full_gaussian_large_s101 $G --N 288 --dx 0.4 --dt 0.06 --T 48 --noise-file $SCR/noise/g288_s101.f64 --mirror 1
