# Preserved limitations and corrections

- The first `verify.py` run used a `1e-20` variance threshold for an exactly
  aligned reversible exchange. It failed because one 96-node control had
  `3.95e-17` variance from floating-point roundoff. Its full output is retained
  as `VERIFICATION_INITIAL.json`. The threshold was corrected to `1e-13`;
  the off-grid control is still required to show real positive broadening.
- The initial expanding64 output predates a correction to the *label and
  species selection* in the collisionless above-threshold control. Its physical
  trajectory was replayed as `expanding64_replay` with zero change. The current
  control applies the positive-to-antiparticle reverse threshold only to the
  positive spectrum.
- Fixed-Q, then changing-frequency/reference-kernel hybrid, then the
  charge-interpolated seven-profile transport are distinct approximations;
  none may be labelled a fully coupled excited-core or inhomogeneous solution.
- The source spectra are compacted signed flux diagnostics. Negative aggregate
  bins are excluded from the positive distribution and their count/energy are
  reported; omitted support is never reassigned to a boundary node as cooling.
- A previous reply incorrectly called the total neutral event count the net
  core-charge transfer. `channels64.json` implements individual channel
  counters and corrects the stoichiometric interpretation.
