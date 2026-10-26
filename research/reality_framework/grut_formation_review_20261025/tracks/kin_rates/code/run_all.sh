#!/bin/sh
# kin_rates reproduction (own code).  Scratch: $SCR.  Total ~2.5 h on 4 cores (the MC tabulation dominates).
set -e
SCR=${SCR:-/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/kin_rates/prod}
mkdir -p $SCR
cd "$(dirname "$0")"
python3 -B model.py && python3 -B amp.py && python3 -B phasespace.py            # self-tests (~10 s)
python3 -B run_thermal.py ../results/thermal_masses.json                        # ~2 min
python3 -B run_22.py ../results/kernels_22.json                                 # ~6 min
OMP_NUM_THREADS=1 python3 -B run_mc_checks.py ../results/mc_checks.json         # ~1 min
python3 -B check_pole.py ../results/check_pole.json                             # ~1 min
OMP_NUM_THREADS=1 python3 -B run_kernels.py $SCR                                # ~2 h (4 processes): kernels_n2.json, xsec_2n.json
python3 -B run_supplement.py $SCR/supplement_n2.json                            # ~15 min
python3 -B run_balance.py $SCR/kernels_n2.json $SCR/xsec_2n.json ../results/kernels_22.json ../results/balance.json
python3 -B build_interface.py $SCR/kernels_n2.json $SCR/xsec_2n.json $SCR/supplement_n2.json ../results/balance.json
cp $SCR/kernels_n2.json ../results/kernels_n2.json; cp $SCR/xsec_2n.json ../results/xsec_2n.json; cp $SCR/supplement_n2.json ../results/supplement_n2.json
