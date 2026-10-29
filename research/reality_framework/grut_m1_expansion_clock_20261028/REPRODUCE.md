# Reproduce the M1 expansion packet

Use the checkout whose parent is `c4815306753ece00dad2b03412f0ed9fe89b8f92`
and this entire directory. The parent packet and the stored source spectra are
inputs. Run from the repository root. Python needs NumPy, SciPy, and
Matplotlib; the radial and scattering calculations also require `g++`.
The tested environment is recorded in `ENVIRONMENT.json`, and exact parent
input hashes are in `INPUT_PINS.json`.

```bash
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
P=research/reality_framework/grut_m1_expansion_clock_20261028

# Rebuild stationary E(Q), its held-out tests, and the seven scattering tables.
python "$P/code/family.py"
for q in 186000 189000 192000 195000 198000 201000 204000; do
    python "$P/code/core_family.py" --q "$q"
done

# Fine source and baseline; independently calculate 64- and 96-node trajectories
# to a/a0=2, then restart each one through a/a0=10.
for n in 64 96; do
    python "$P/code/hybrid.py" --profiled --shift-bins "$n" --source parent --tag "profiled$n"
    python "$P/code/hybrid.py" --profiled --shift-bins "$n" --source parent --tag "profiled${n}_to10" --resume "profiled$n" --amax 10
    python "$P/code/hybrid.py" --profiled --shift-bins "$n" --source h0025 --tag "profiled_h0025_$n"
    python "$P/code/hybrid.py" --profiled --shift-bins "$n" --source h0025 --tag "profiled_h0025_${n}_to10" --resume "profiled_h0025_$n" --amax 10
done

# Direct populated-bin scattering holdouts on the stored 96-node fixed-core run.
python "$P/code/kernel_check.py" --tau 0 --count 60
for t in 1 10 100 1000; do
    python "$P/code/kernel_check.py" --tau "$t" --count 200
done

# Results, independent-of-conservation numerical controls and figure.
python "$P/code/summary.py" profiled96_to10 profiled_h005_96_to10 profiled_h0025_64_to10 profiled_h0025_96_to10
python "$P/code/freeflight.py" profiled64_to10 profiled96_to10 profiled_h005_64_to10 profiled_h005_96_to10 profiled_h0025_64_to10 profiled_h0025_96_to10
python "$P/code/figure.py"
python "$P/code/verify.py"
python "$P/code/beam_kinematics.py"
```

The committed files also include the fixed-core `expansion.py` results, the
reference-kernel `hybrid.py` sensitivity results, 32-node controls, 64-node
frequency/tolerance/angular controls, `h=0.05` source comparisons, and the
time-local reaction counters. See `configuration` inside each JSON for the
exact command-line parameters. The test script checks these saved artifacts,
so do not delete them before running the final verification.

To recreate the finer source from the pinned parent radial engine, run:

```bash
python "$P/code/source_ringdown.py" --h .025 --dt .003125 --tag ringdown_h0025
cp "$P/ringdown_h0025.json" "$P/sources/"
cp "$P/ringdown_h0025_spectrum_0_1000.npz" "$P/sources/"
cp "$P/ringdown_h0025_spectrum_1000_10000.npz" "$P/sources/"
```

This is a much longer computation than reading the supplied source files. The
same code can produce `h005` with `--h .05 --dt .00625`, `h01` with
`--h .1 --dt .0125`, and its step-size control with `--h .1 --dt .00625`.
The `h005` and `h01` inputs were copied from the same-author unfinished sibling
worktree and are retained with their hashes; they are not independent data.
`code/source_ringdown.py`, `code/support.py`, and `code/channels.py` were copied
from that sibling study and run in this packet's isolated worktree. The source
generator writes transient raw radial traces beneath `.work/`, which is ignored
by git; compact signed spectra and diagnostics are committed.

For the reaction-counter check:

```bash
python "$P/code/channels.py" --bins 64 --tag channels64
```

The 1% result remains conditional on the source, homogeneous mixing,
quasistatic and pressure-work assumptions spelled out in `REPORT.md`.
