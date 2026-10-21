#!/bin/sh
# verifier lattice queue: 2 single-thread workers
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$(dirname "$0")"
# mu = 2.507e-4 gives m0 = erf(mu/(sqrt2*0.1)) = 0.002
printf '%s\n' \
 "192 101 0.0" "192 101 0.0002507" "192 101 -0.0002507" \
 "192 102 0.0" "192 102 0.0002507" "192 102 -0.0002507" \
 "192 103 0.0" "192 103 0.0002507" "192 103 -0.0002507" \
 "256 104 0.0" "256 104 0.0002507" "256 104 -0.0002507" \
 | xargs -P 2 -L 1 sh -c 'python3 lattice_v.py $0 $1 $2'
