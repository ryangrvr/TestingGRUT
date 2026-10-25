#!/bin/sh
# REPRODUCTION (Codex code, not independent): instrumented radial_flux.cpp vs radial_cfs.cpp, short nonlinear run.
# Usage: sh repro_bitwise.sh <scratch dir>
set -e
D=${1:-/tmp/claude-0/-home-user-TestingGRUT/553831ba-02d6-5e13-9e04-654a9685bce7/scratchpad/review/operators}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$D/pkt" "$D/repro"
cd /home/user/TestingGRUT && git archive ff49e5d research/reality_framework/grut_causal_receiver_20261024/code | tar -x -C "$D/pkt"
C="$D/pkt/research/reality_framework/grut_causal_receiver_20261024/code"
g++ -O3 -std=c++17 "$C/radial_cfs.cpp" -o "$D/repro/cfs" 2>/dev/null
g++ -O3 -std=c++17 "$C/radial_flux.cpp" -o "$D/repro/flux" 2>/dev/null
cd "$HERE" && python3 -c "
import numpy as np, sys
from formation_bound import static_S
h=.5;L=450.;n=int(L/h);r=np.arange(1,n)*h
f=.18*np.exp(-r**2/7200.);Z=static_S(r,f,h,n);S=1+Z/r
np.savetxt('$D/repro/init.csv',np.column_stack([r,f,0*r,S,0*r,-f,0*r]),delimiter=',',header='rho,x,y,s,px,py,ps',comments='')"
cd "$D/repro"
./cfs init.csv base 0.5 0.03125 450 600 240 360 20 0.1 1 1 > base.log
./flux init.csv inst 0.5 0.03125 450 600 240 360 20 0.1 1 1 > inst.log
for k in final memory trace; do cmp base_$k.csv inst_$k.csv && echo "$k identical"; done
python3 -c "
import numpy as np
b=np.loadtxt('base_history.csv',delimiter=',',skiprows=1);i=np.loadtxt('inst_history.csv',delimiter=',',skiprows=1)
print('history_cols_0_7_identical',np.array_equal(b,i[:,:8]))
E0,Q0=i[0,8],i[0,9]
print('cut_energy_closure_max_rel',np.max(abs(i[:,8]+i[:,10]-E0))/E0,'cut_charge_closure_max_rel',np.max(abs(i[:,9]+i[:,11]-Q0))/Q0)
print('net_E_Q_through_cut',i[-1,10],i[-1,11])"
