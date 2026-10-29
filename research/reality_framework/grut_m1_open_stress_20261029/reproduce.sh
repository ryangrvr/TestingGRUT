#!/usr/bin/env bash
# Full regeneration is CPU-intensive (the fine 20,000-unit run dominates).
# Run from this directory inside the pinned parent repository. Numerical
# dependencies: C++17 compiler, Python 3, NumPy, SciPy, Matplotlib.
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python code/check_inputs.py
mkdir -p .work
c++ -O3 -std=c++17 code/stress_radial.cpp -o .work/stress_radial
c++ -O3 -std=c++17 code/axisymmetric.cpp -o .work/axisymmetric
for cfg in '01 .1 .0125' '005 .05 .00625' '0025 .025 .003125'; do
  read -r suffix h dt <<< "$cfg"
  python code/initial.py --h "$h" --type pulse --tag "pulse_h$suffix"
  python code/initial.py --h "$h" --type stationary --tag "stationary_h$suffix"
  .work/stress_radial ".work/stationary_h$suffix.csv" ".work/stat$suffix" "$h" "$dt" 200 200 60 100 .2
  .work/stress_radial ".work/pulse_h$suffix.csv" ".work/pulse$suffix" "$h" "$dt" 200 20000 60 100 .2
  if [[ "$suffix" == 0025 ]]; then end=1000; else end=2000; fi
  .work/stress_radial ".work/pulse_h$suffix.csv" ".work/sampled$suffix" "$h" "$dt" 200 "$end" 60 100 .2 .25
done
# The final state starts at absolute mXt=20000; the restart clock below starts at zero.
.work/stress_radial .work/pulse0025_final.csv .work/late_fine0025 .025 .003125 200 1000 60 100 .2 .25
.work/stress_radial .work/pulse_h01.csv .work/open_2000 .1 .0125 200 2000 60 100 .2
for r in 80 120 160; do
  python code/initial.py --h .1 --L "$r" --type pulse --tag "pulse_reflect_$r"
  if [[ "$r" == 160 ]]; then out=reflect160_retry; else out="reflect_$r"; fi
  .work/stress_radial ".work/pulse_reflect_$r.csv" ".work/$out" .1 .0125 "$r" 2000 60 40 0
done
python code/axis_initial.py --h .1 --lmax 4 --nq 12 --amp .1 --quad 0 --tag mono_h01
.work/axisymmetric .work/mono_h01_initial.csv .work/mono_h01_angles.csv .work/axis_mono300 .1 .0125 200 300 60 100 .2 4 12
.work/stress_radial .work/pulse_h01.csv .work/radial_mono300 .1 .0125 200 300 60 100 .2
for cfg in 'quad_h02 .2 .025 4 12 quad02' 'quad_h01 .1 .0125 4 12 quad01' 'quad_h005 .05 .00625 4 12 quad005' 'quad_l6 .1 .0125 6 16 quad_l6' 'quad_nq16 .1 .0125 4 16 quad_nq16'; do
  read -r init h dt lmax nq run <<< "$cfg"
  python code/axis_initial.py --h "$h" --lmax "$lmax" --nq "$nq" --tag "$init"
  .work/axisymmetric ".work/${init}_initial.csv" ".work/${init}_angles.csv" ".work/$run" "$h" "$dt" 200 300 60 100 .2 "$lmax" "$nq"
  case "$run" in quad_l6) out=axis_specl6;; quad_nq16) out=axis_specnq16;; *) out="axis_spec${run#quad}";; esac
  python code/axis_spectra.py --tag "$run" --h "$h" --lmax "$lmax" --width 200 --out "$out"
done
for cfg in 'quad_h01 .1 .0125 4 12 quad01' 'quad_h005 .05 .00625 4 12 quad005' 'quad_l6 .1 .0125 6 16 quad_l6'; do
  read -r init h dt lmax nq run <<< "$cfg"
  .work/axisymmetric ".work/${init}_initial.csv" ".work/${init}_angles.csv" ".work/${run}_t1000" "$h" "$dt" 200 1000 60 100 .2 "$lmax" "$nq"
  case "$run" in quad_l6) spec=axis_specl6_t1000; direction=axis_directionl6_t1000;; *) spec="axis_spec${run#quad}_t1000"; direction="axis_direction${run#quad}_t1000";; esac
  python code/axis_spectra.py --tag "${run}_t1000" --h "$h" --lmax "$lmax" --width 200 --out "$spec"
  python code/axis_direction.py --tag "${run}_t1000" --h "$h" --lmax "$lmax" --out "$direction"
done
for width in 100 300; do python code/axis_spectra.py --tag quad005 --h .05 --width "$width" --out "axis_spec005_w$width"; done
for width in 100 500; do python code/axis_spectra.py --tag quad005_t1000 --h .05 --width "$width" --out "axis_spec005_t1000_w$width"; done
for cutoff in .05 .2 .5; do
  case "$cutoff" in .05) suffix=005;; .2) suffix=02;; .5) suffix=05;; esac
  python code/axis_direction.py --tag quad005_t1000 --h .05 --kmin "$cutoff" --out "axis_direction005_t1000_k$suffix"
done
python code/timing_sensitivity.py --late 20 --final 40 --out timing_sensitivity_a40.json
python code/timing_sensitivity.py --late 100 --final 1000 --out timing_sensitivity.json
python code/summarize.py
python code/verify.py
python code/figure.py
