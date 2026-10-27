#!/bin/bash
# Run (resume) the main BM6 history at h0=0.004 and its h0=0.002 control in parallel, 1 FFT worker each (2 threads total).
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/en_indep
C=$(dirname $0)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 INDEP_WORKERS=1
COMMON="--init $SP/ework/N48_initial.bin --am 8.869218738806 --C 1.9857962750634526 --dtau 0.05 --fine_until 150 --dtau_late 0.25 --until ${UNTIL:-800} --budget ${BUDGET:-520} --snap_taus 40 60 70 100 150 200 400 800"
python3 -B $C/indep3d.py $COMMON --out $SP/runs/bm6_h004 --h0 0.004 > $SP/runs/bm6_h004.log 2>&1 &
python3 -B $C/indep3d.py $COMMON --out $SP/runs/bm6_h002 --h0 0.002 > $SP/runs/bm6_h002.log 2>&1 &
wait
tail -n 2 $SP/runs/bm6_h004.log $SP/runs/bm6_h002.log
