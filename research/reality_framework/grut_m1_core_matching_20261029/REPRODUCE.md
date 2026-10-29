# Reproduce the M1 core-matching packet

Begin with published parent commit `8f67a83be4177f6baaa475c94565aab3b7c0cd94`
and this complete directory. Run from the repository root in a disposable
checkout. Python needs NumPy, SciPy and Matplotlib; the field solvers need
`g++` with OpenMP. `ENVIRONMENT.json` gives the tested versions and
`INPUT_PINS.json` records the parent input and raw-checkpoint hashes.

The committed JSON/NPZ files are enough to verify the numerical controls
without re-running the expensive field simulations:

```bash
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
P=research/reality_framework/grut_m1_core_matching_20261029
D=research/reality_framework/grut_m1_expansion_clock_20261028
C=research/reality_framework/grut_m1_local_cloud_20261029
python "$P/code/manifest.py" --check
python "$P/code/verify.py"
```

`manifest.py --check` applies to the originally committed files and should
run before regenerating outputs. Some JSON path metadata and plot metadata
may differ across machines; a bytewise manifest is an integrity check, not a
physics tolerance. The tests are same-author controls, not independent
replication.

## Poisson and kinetic controls from committed inputs

```bash
python "$P/code/mixing_continuum.py" --nv 2001 --tag mixing_cont_nv2001
python "$P/code/mixing_continuum.py" --nv 4001 --tag mixing_cont_nv4001
python "$P/code/mixing_continuum.py" --nv 8001 --times 100000 200000 500000 1000000 10000000 100000000 1000000000 10000000000 --tag mixing_cont_nv8001
python "$P/code/mixing_discrete.py" --tag mixing_discrete_sweep
for n in 140 560; do
    python "$P/code/mixing_continuum.py" --source "source_t400_e$n" --nv 4001 --times 100000 1000000 10000000 --tag "mixing_cont_e$n"
    python "$P/code/mixing_discrete.py" --source "source_t400_e$n" --tag "mixing_discrete_e$n"
done
for w in 200 800; do
    python "$P/code/mixing_continuum.py" --source "source_t$w" --nv 4001 --times 100000 1000000 10000000 --tag "mixing_cont_t$w"
done
for h in h005 h0025; do
    for n in 64 96; do
        python "$P/code/late_kinetics.py" --source "$h" --bins "$n" --tag "late_kin_${h}_$n"
    done
done
python "$P/code/scenario_audit.py"
python "$P/code/figure.py"
python "$P/code/verify.py"
```

`late_kinetics.py` uses the committed late compact spectra and parent
stationary profile tables. Its 96-node run takes several minutes. The
discrete-speed sweep can lose up to about 0.7% number normalization at
`m_Xt=10^10`; the JSON records that error.

## Recreate radial field continuations and matching surfaces

The full field histories and `.work` CSVs are intentionally untracked.
The following regenerates the parent pulse before restarting it. The
`h=.025` field evolution and its continuation are lengthy computations.
The initial and restart trace hashes are recorded, but compiler/platform
differences may prevent bytewise identity; compare physical diagnostics too.

```bash
mkdir -p "$P/.work" "$C/.work"
python "$D/code/source_ringdown.py" --h .025 --dt .003125 --tag ringdown_h0025
python "$D/code/source_ringdown.py" --h .05 --dt .00625 --tag ringdown_h005
g++ -O3 -std=c++17 -ffp-contract=off research/reality_framework/grut_m1_dynamics_20261027/code/radial.cpp -o "$P/.work/radial"

"$P/.work/radial" "$D/.work/ringdown_h0025_final.csv" "$P/.work/late_h0025" .025 .003125 200 10000 60 100 .2
"$P/.work/radial" "$D/.work/ringdown_h005_final.csv" "$P/.work/late_h005" .05 .00625 200 10000 60 100 .2
for h in h005 h0025; do
    step=.05
    if [ "$h" = h0025 ]; then step=.025; fi
    python "$P/code/compact_late.py" --tag "$h"
    for w in 200 400 800; do
        python "$P/code/surface_source.py" --trace "$P/.work/late_${h}_trace.csv" --h "$step" --width "$w" --tag "late_${h}_t$w"
    done
    python "$P/code/radial_analysis.py" --tag "$h"
done

# Recreate the 40 and 80 probes from the same h=.05 initial field.
for radius in 40 80; do
    "$P/.work/radial" "$D/.work/ringdown_h005_initial.csv" "$C/.work/radial_r$radius" .05 .00625 200 10000 "$radius" 100 .2
done
for radius in 40 60 80; do
    trace="$C/.work/radial_r${radius}_trace.csv"
    if [ "$radius" = 60 ]; then trace="$D/.work/ringdown_h005_trace.csv"; fi
    for w in 200 800; do
        python "$P/code/surface_source.py" --trace "$trace" --h .05 --probe "$radius" --width "$w" --tag ".work/probe_r${radius}_t$w"
    done
done
python "$P/code/matching_surfaces.py"
```

The 400-wide matching-surface arrays used by the comparison are committed
in the parent local-cloud packet. The regional core-excess fits are empirical
and do not establish a true asymptotic state.

## Recreate the 3D controls

The following uses a periodic side-64 box and produces many gigabytes of
temporary snapshots under `$P/.work`. The original saved trace contains only
derived arrays, not raw 96-cubed field snapshots. The 64-bit replay is a
precision control for differences of nearly equal fields.

```bash
g++ -std=c++17 -O3 -ffp-contract=off -fopenmp research/reality_framework/grut_m1_dynamics_20261027/code/field3d.cpp -o "$P/.work/field3d"
for n in 64 96; do
    python "$P/code/generate_3d.py" --N "$n" --perturbation 0 --tag "object${n}_control"
    python "$P/code/generate_3d.py" --N "$n" --perturbation .05 --tag "object${n}_perturbed"
done
python "$P/code/generate_3d.py" --N 96 --perturbation 0 --tag object96_control_double --T 30 --bits 64
python "$P/code/generate_3d.py" --N 96 --perturbation .05 --tag object96_perturbed_double --T 30 --bits 64

python "$P/code/multipoles.py" --N 64 --tag multipoles64
python "$P/code/multipoles.py" --N 96 --tag multipoles96
python "$P/code/multipoles.py" --N 96 --ntheta 12 --nphi 24 --tag multipoles96_ang12
python "$P/code/multipoles.py" --N 96 --dr .5 --tag multipoles96_dr05
python "$P/code/multipoles.py" --N 64 --lmax 6 --tag multipoles64_l6
python "$P/code/multipoles.py" --N 96 --lmax 6 --tag multipoles96_l6
python "$P/code/multipoles.py" --N 96 --max-time 20 --tag multipoles96_t20
python "$P/code/multipoles.py" --N 96 --suffix _double --max-time 29.64 --tag multipoles96_double
python "$P/code/multipoles.py" --N 96 --suffix _double --max-time 29.64 --lmax 6 --tag multipoles96_double_l6
python "$P/code/multipoles.py" --N 96 --suffix _double --max-time 29.64 --dr .5 --tag multipoles96_double_dr05
python "$P/code/multipoles.py" --N 96 --suffix _double --max-time 20 --tag multipoles96_double_t20
python "$P/code/angular_audit.py"
python "$P/code/angular_flux.py"
python "$P/code/verify.py"
```

The multipoles describe finite-radius classical field differences and signed
energy flux. They do not supply particle occupations or a two-particle
collision-angle distribution.
