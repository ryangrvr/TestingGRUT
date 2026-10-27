#!/bin/bash
# Resume the h0=0.002 control and the perturbed twin (eps=1e-6 at tau=100), 1 FFT worker each.
SP=/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/en_indep
C=$(dirname $0)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 INDEP_WORKERS=1
COMMON="--init $SP/ework/N48_initial.bin --am 8.869218738806 --C 1.9857962750634526 --dtau 0.05 --fine_until 150 --dtau_late 0.25 --until ${UNTIL:-800} --budget ${BUDGET:-520} --snap_taus 40 60 70 100 150 200 400 800"
[ -d $SP/runs/twin_eps1e-6 ] || python3 -B $C/perturb.py $SP/runs/bm6_h004/snap_tau100.00.npz $SP/runs/twin_eps1e-6 1e-6
python3 -B $C/indep3d.py $COMMON --out $SP/runs/bm6_h002 --h0 0.002 > $SP/runs/bm6_h002.log 2>&1 &
python3 -B $C/indep3d.py $COMMON --out $SP/runs/twin_eps1e-6 --h0 0.004 > $SP/runs/twin.log 2>&1 &
wait
tail -n 2 $SP/runs/bm6_h002.log $SP/runs/twin.log
