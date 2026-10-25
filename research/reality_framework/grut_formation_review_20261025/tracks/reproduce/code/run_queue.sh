#!/bin/bash
# Reproduction-audit job queue (labelled REPRODUCTION: runs Codex's own scripts, unmodified,
# in a git-archive extraction of ff49e5d; only the compiler flag -O3 is replaced by -O2
# through gxx_O2_wrapper.sh placed first on PATH).
# usage: run_queue.sh TREE SCRATCH
TREE=$1; S=$2
export PATH=$S/bin:$PATH GXX_WRAPPER_LOG=$S/gxx.log OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 TREE S
C=$TREE/research/reality_framework/grut_causal_receiver_20261024/code
B=$TREE/research/reality_framework/grut_core_conversion_20261023/code
mkdir -p $S/logs
cat > $S/onejob.sh <<'ONE'
#!/bin/bash
name=$1; shift
t0=$(date +%s.%N)
if [ "$name" = "P025_time+I025" ]; then
  ( cd $TREE && python3 "$1" P025_time && python3 "$2" ) > $S/logs/$name.out 2>&1
else
  ( cd $TREE && "$@" ) > $S/logs/$name.out 2>&1
fi
rc=$?; t1=$(date +%s.%N)
echo "$name rc=$rc wall=$(python3 -c "print(round($t1-$t0,1))") end=$(date +%T)" >> $S/logs/_status.txt
ONE
chmod +x $S/onejob.sh
cat > $S/jobs.txt <<JOBS
P00625 python3 $C/run.py P00625
H003125 python3 $B/refinement.py H003125
N00625 python3 $B/refinement.py N00625
P0125 python3 $C/run.py P0125
P025_remote python3 $C/run.py P025_remote
P025_remote_time python3 $C/run.py P025_remote_time
P025_time+I025 $C/run.py $C/run_interface.py
N0125_time python3 $B/refinement.py N0125_time
H00625 python3 $B/refinement.py H00625
N0125_box python3 $B/refinement.py N0125_box
N0125_far python3 $B/refinement.py N0125_far
P025 python3 $C/run.py P025
N0125 python3 $B/refinement.py N0125
N0125_abs python3 $B/refinement.py N0125_abs
free_soft python3 $C/run.py free_soft
H0125_late python3 $B/refinement.py H0125_late
H0125_time python3 $B/refinement.py H0125_time
free_X_time python3 $C/run.py free_X_time
free_s python3 $C/run.py free_s
H0125_box python3 $B/refinement.py H0125_box
N025_time python3 $B/refinement.py N025_time
H0125_abs python3 $B/refinement.py H0125_abs
H0125 python3 $B/refinement.py H0125
N025 python3 $B/refinement.py N025
free_X python3 $C/run.py free_X
H025_time python3 $B/refinement.py H025_time
H025 python3 $B/refinement.py H025
N05 python3 $B/refinement.py N05
smoke python3 $C/run.py smoke
discrete_layer python3 $C/discrete_layer.py
operator_controls python3 $C/operator_controls.py
JOBS
xargs -P 4 -L 1 $S/onejob.sh < $S/jobs.txt
echo ALL_DONE >> $S/logs/_status.txt
