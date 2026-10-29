# Retained limitations and negative controls

1. The review's statement that late radiation is antiparticle-free is not
   supported by the extension. The fine 400-window positive projection has
   1.503 direct antiparticles and 3.238 neutral parents per core in the
   `10000–20000` interval. These are not a certified quantum spectrum and
   do not describe all later times.
2. The parent `shot_noise_dr2p5` result is correct for its 280 bin-center
   speeds, but the source-bin population between centers was not resolved.
   Exact retarded-time knot sweeping reproduces that result; a separate
   continuous interpolation changes its variance greatly. Both are retained
   as conditional models rather than silently substituting one for truth.
3. `mixing_cont_nv2001`, `nv4001` and `nv8001` compare radial resolution;
   `e140`, `e560`, `t200` and `t800` retain energy and time-window
   sensitivities. At `m_Xt=10^10`, the discrete knot sweep accumulates
   cancellation error up to 0.7% in source normalization. These controls do
   not fix the physical within-bin spectral shape.
4. The 96-cubed 3D archive contained derived trace arrays and image
   projections, not the raw snapshots. The original initial condition and
   early dynamics were reproduced from the stored code. A 5% disturbance in
   a side-64 periodic box does not measure an isolated 10% source at `r=60`.
5. A short-window, 32-bit snapshot analysis at `r=24`, `m_Xt<=20` rounded
   the neutral difference field to zero and initially caused a division by
   zero. A matched 64-bit replay retains nonzero but extremely small values;
   neither result should be interpreted as a certified physical radiation
   signal at that radius/time. The code now reports `null` for zero-power
   normalized quantities. At `r=20` through `m_Xt=29.64`, 32/64-bit results
   agree, while neutral 64/96 grid comparisons do not.
6. Time-integrated energy fits with and without a nonzero floor both track
   the finite history, but neither is a complex quasinormal-mode calculation
   or evidence of a future bound state. Source radiation beyond `m_Xt=20000`
   and interaction-induced changes remain unknown.
7. The 1% comparison is a conditional energy model; it neither derives the
   primordial population nor certifies local mixing, positive source
   occupation, quantum binding, pressure work or unmodeled energy inputs.
