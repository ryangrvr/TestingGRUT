# Attempts and negative controls

These are retained to make the numerical path auditable. They are not
additional physical channels and must not be summed.

1. `primary_t400_e140.json`, `primary_t400.json` and
   `primary_t400_e560.json` use one representative antiparticle energy per
   bin in the collinear NWA. The narrow time/energy overlap is undersampled
   and the answer shifts sharply with the binning. `finitewidth_e*.json`
   repeats that discretization for the full s-channel and shows similarly
   poor behavior: 280 and 560 bins give substantially different answers.
   These midpoint results are **not converged estimates**.
2. `continuous.py` integrates inside the antiparticle bins with
   Gauss–Legendre quadrature. The retained 140/280/560 bin and quadrature
   controls expose a remaining roughly 9% source-bin dependence, as well as
   strong sensitivity to the lower-energy tail. The `continuous_cut_*`
   files intentionally remove slow antiparticles; the `continuous_zero_*`
   files impose a different lower-endpoint extrapolation.
3. `fullwidth_continuous.py` integrates both energies with a split around
   each pole. Its positive-energy support stops at node centers: missing
   endcaps and varying resolution remain limitations. Its full width and
   own NWA control are alternative evaluations of the same s-channel, not
   additive rates.
4. Neutral-daughter exposure has a noisy rare positive-daughter component
   and strong time-window dependence; compare `daughter_131072_t200`,
   `daughter_262144_s19`, `daughter_262144_s37`, `daughter_262144_s59`
   and `daughter_131072_t800`. The IID-style diagnostic reported for a
   scrambled Sobol sequence is not claimed as a rigorous uncertainty band.
5. Spatial midpoint sampling of the Poisson profiles shows perceptible
   source-number error at coarse steps; `shot_noise_dr2p5` and `dr5` are the
   finer controls. Even a well-resolved density covariance does not by
   itself determine an angle-dependent scattering rate.

All calculations remain same-author, model-conditional checks. No negative
or missing result should be reclassified as a demonstration of formation,
mixing, or the final 1% criterion.
