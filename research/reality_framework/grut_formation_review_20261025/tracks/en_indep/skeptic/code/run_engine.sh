#!/bin/bash
# Skeptic: restart Codex's early3d.cpp (copied verbatim from the pristine 7d9c93b archive, compiled with Codex's flags)
# usage: run_engine.sh INIT.bin PREFIX TAU0 A0 SAMPLE SNAP [END_TAU] [HMAX]
set -e
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad
B=$SP/review/en_indep_skeptic/bin; mkdir -p $B
SRC=$SP/review/en_indep_skeptic/pub/research/reality_framework/grut_early_nonlinear_20261027/code/early3d.cpp
[ -x $B/early3d ] || g++ -std=c++17 -O3 -ffp-contract=off -fopenmp $SRC -o $B/early3d
C=1.9857962750634526; END=${7:-800}; H=${8:-0.004}
DUR=$(python3 -I -c "import math;t0=$3;a=$4;C=$C;e=$END;print(repr(2*(e-t0)/(a+math.sqrt(a*a+2*C*(e-t0)))))")
$B/early3d $1 $2 $H $DUR $5 $6 4 2 64 $3 $4 $C
