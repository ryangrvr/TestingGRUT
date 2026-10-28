# Reproduce the bounded M1 scattering gate

From the repository root, with NumPy, SciPy, Matplotlib and g++17 available:

```bash
export GRUT_M1_SCATTER_WORK=/tmp/grut_m1_scattering_20261028
git fetch origin e71730c24fe6a0b28ddd72cd30982c5b39c9aedf
python research/reality_framework/grut_m1_scattering_20261028/code/vertices.py
python research/reality_framework/grut_m1_scattering_20261028/code/scattering.py
python research/reality_framework/grut_m1_scattering_20261028/code/continuum.py
python research/reality_framework/grut_m1_scattering_20261028/code/inverse.py
python research/reality_framework/grut_m1_scattering_20261028/code/regenerate_sources.py
python research/reality_framework/grut_m1_scattering_20261028/code/source.py
python research/reality_framework/grut_m1_scattering_20261028/code/partials.py
python research/reality_framework/grut_m1_scattering_20261028/code/source_weights.py
python research/reality_framework/grut_m1_scattering_20261028/code/resonance.py
python research/reality_framework/grut_m1_scattering_20261028/code/packets.py
python research/reality_framework/grut_m1_scattering_20261028/code/packet_analysis.py
python research/reality_framework/grut_m1_scattering_20261028/code/clocks.py
python research/reality_framework/grut_m1_scattering_20261028/code/cold_completion.py
python research/reality_framework/grut_m1_scattering_20261028/code/vacuum_kernels.py
python research/reality_framework/grut_m1_scattering_20261028/code/verify.py
python research/reality_framework/grut_m1_scattering_20261028/code/figure.py
python research/reality_framework/grut_m1_scattering_20261028/code/audit.py
```

Use an empty work directory for a fresh replication. The git fetch supplies the
exact old kinetic decay reference, not an updated model. The scalar width is
also derived directly and independently checked by real-component and angular
phase-space calculations. No Mathematica, Claude, Gemini or external agent is
needed to execute these scripts.

The source regenerator uses the pinned published M1 radial engine and compares
every scalar history column against the existing compressed archives. It avoids
using an incomplete old scratch trace or dropping malformed data. Full traces
are regenerated cheaply; source spectra and checksums are exported. The packet
runner uses that same unchanged nonlinear engine with new initial states and an
explicit remote absorbing receiver. Full packet traces/histories are preserved
in PACKET_DATA.npz. Source spectra, angular kernels and resonance controls have
their own NPZ files.

The finite-difference stationary solver and continuum collocation solver are
separate implementations by the same author. The latter uses real/imaginary
equations and should take seconds to tens of seconds depending on BLAS. Partial
wave and resonance grids take longer; caches are conveniences for these fixed
inputs, not permissions to reuse results after changing the point or grid.

No artificial friction is applied inside the object. The nonlinear pulse and
packet kernels are classical and local. Vacuum phi decay is entered analytically
after scattering; it is not silently implemented as stochastic decay inside the
classical engine. The complete in-medium quantum reaction operator and actual
cosmological density history are not reconstructed by this reproduction.

The cold-completion table is a stated ideal-processing diagnostic. It reuses
the existing stationary M1 BVP solver at new charges; it does not re-evolve a
cosmological population. Its pass/fail refers to the unchanged program1% budget,
not a new observational confidence limit.
