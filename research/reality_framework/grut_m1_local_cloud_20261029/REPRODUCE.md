# Reproducing the local-plume packet

Start from the published parent commit
`12b5696c2e5722d253a02ca4e9d57435879f5dc1` plus this complete directory.
Run from the repository root, preferably in a disposable checkout. Requires
Python with NumPy, SciPy and Matplotlib and `g++` for the radial engine.
The tested environment and hashes are in `ENVIRONMENT.json` and
`INPUT_PINS.json`. Run with one numerical-library thread for the closest
comparison. The committed source `.npz` files permit all downstream tests
without running the lengthy field simulation.

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
P=research/reality_framework/grut_m1_local_cloud_20261029
D=research/reality_framework/grut_m1_expansion_clock_20261028

# Check the original packet before regenerating outputs.
python "$P/code/manifest.py" --check

# Direct returns use compact neutral spectra committed in the parent packet.
python "$P/code/neutral_return.py"
python "$P/code/angular_geometry.py"
python "$P/code/general_width.py"

# Primary energy-bin and quadrature sensitivities.
python "$P/code/continuous.py" --source source_t400_e140 --nquad 64 --tag continuous_e140_q64
python "$P/code/continuous.py" --source source_t400 --nquad 64 --tag continuous_e280_q64
python "$P/code/continuous.py" --source source_t400 --nquad 128 --tag continuous_e280_q128
python "$P/code/continuous.py" --source source_t400_e560 --nquad 32 --tag continuous_e560_q32
python "$P/code/continuous.py" --source source_t400_e560 --nquad 64 --tag continuous_e560_q64
python "$P/code/continuous.py" --source source_t400 --nquad 64 --eminanti 1.025 --tag continuous_cut_1_025
python "$P/code/fullwidth_continuous.py" --source source_t400 --nquad 64 --tag full_cont_e280_q64
python "$P/code/fullwidth_continuous.py" --source source_t400_e560 --nquad 32 --tag full_cont_e560_q32

# Explicit daughters and geometric first passage. Three scrambles are saved.
for seed in 19 37 59; do
    python "$P/code/daughter_cloud.py" --source source_t400 --power 18 --seed "$seed" --tag "daughter_262144_s$seed"
done
python "$P/code/daughter_cloud.py" --source source_t200 --power 17 --seed 37 --tag daughter_131072_t200
python "$P/code/daughter_cloud.py" --source source_t800 --power 17 --seed 37 --tag daughter_131072_t800
python "$P/code/mixing.py" --seed 37 --tag mixing_seed37

# Poisson density controls and the distinct-core finite-width kernel.
python "$P/code/shot_noise.py" --dr 2.5 --tag shot_noise_dr2p5
python "$P/code/shot_noise.py" --dr 5 --tag shot_noise_dr5
python "$P/code/cross_source.py" --source source_t400_e140 --quad 64 --tag cross_e140_q64
python "$P/code/cross_source.py" --source source_t400 --quad 32 --tag cross_e280_q32
python "$P/code/cross_source.py" --source source_t400 --quad 64 --tag cross_e280_q64
python "$P/code/cross_source.py" --source source_t400_e560 --quad 64 --tag cross_e560_q64

python "$P/code/verify.py"
python "$P/code/figure.py"
```

The daughter and full-width runs are substantially more expensive than the
other downstream commands. Files with additional parameter sweeps, time
windows and exploratory midpoint estimates are also committed. Each JSON's
`configuration` records its exact arguments. The verification script reads
several saved controls, so retain the directory's JSON files even if only a
subset of the commands is rerun. A successful rerun is a same-author check,
not independent replication. A rebuilt figure or regenerated `.npz`/JSON
metadata may differ bytewise across environments; the manifest is the
integrity check for the original committed packet, not a physics tolerance.

## Recreating time-frequency sources

The raw field CSVs are deliberately excluded from git as transient build
data. Their SHA-256 hashes are in `INPUT_PINS.json` and the compressed,
windowed source arrays are committed. Recreate the full traces, if desired,
using the parent radial generator (the `h=.025` run is slow):

```bash
python "$D/code/source_ringdown.py" --h .025 --dt .003125 --tag ringdown_h0025
python "$D/code/source_ringdown.py" --h .05 --dt .00625 --tag ringdown_h005

for width in 200 400 800; do
    python "$P/code/surface_source.py" --trace "$D/.work/ringdown_h0025_trace.csv" --width "$width" --tag "source_t$width"
done
python "$P/code/surface_source.py" --trace "$D/.work/ringdown_h0025_trace.csv" --width 400 --bins 140 --tag source_t400_e140
python "$P/code/surface_source.py" --trace "$D/.work/ringdown_h0025_trace.csv" --width 400 --bins 560 --tag source_t400_e560

g++ -O3 -std=c++17 -ffp-contract=off research/reality_framework/grut_m1_dynamics_20261027/code/radial.cpp -o "$P/.work/radial_probe"
for radius in 40 80; do
    "$P/.work/radial_probe" "$D/.work/ringdown_h005_initial.csv" "$P/.work/radial_r$radius" .05 .00625 200 10000 "$radius" 100 .2
    python "$P/code/surface_source.py" --trace "$P/.work/radial_r${radius}_trace.csv" --h .05 --probe "$radius" --tag "source_h005_r$radius"
done
python "$P/code/surface_source.py" --trace "$D/.work/ringdown_h005_trace.csv" --h .05 --probe 60 --tag source_h005_r60
for radius in 40 60 80; do
    python "$P/code/continuous.py" --source "source_h005_r$radius" --nquad 32 --tag "continuous_h005_r$radius"
done
python "$P/code/diagnostics.py"
```

The first generator creates the relevant ignored `.work` directories. If
invoking only `g++`, create `$P/.work` first. Running the parent generator
also overwrites its generated summaries, so do this in an isolated checkout.
Different compiler/platform versions and absolute trace paths can change
metadata or bitwise hashes; compare physical diagnostics as well as hashes.
The raw classical field and the time-frequency projection are not themselves
a proof of a unique or complete quantum source.
