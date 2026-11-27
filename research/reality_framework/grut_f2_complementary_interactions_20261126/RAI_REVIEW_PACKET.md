# External RAI review packet

Start from the frozen protocol at
2db0de4c0c50aed3a7297964fc8ca9fedaa96991. The parent is
draft PR #24. Verify the hash of the original F1 spring data,
and use the explicit compact-bump specialization of the F1
continuum family. An independent critical implementation is
included for comparison after the reconstruction.
[Protocol](FROZEN_PROTOCOL.md) · [Parameters](PARAMETERS.json)
· [F1 input](../grut_f1_operational_reconstruction_20261126/WITNESSES.json)

1. Rebuild the two oscillator matrices from the edge lists.
   Recompute the exact root contrasts for unit local mass
   and pin at both shifts, the uniform and collective controls,
   the finite-time tail and the coefficient of the tagged
   quartic response. Audit the actuator support and gain
   under each comparison and the matched preparations for
   the nonlinear null.
2. Try to falsify the finite reconstruction lemma. Re-derive
   the positivity interval and monotonic pair function, and
   explain whether any graph or coupling information entered
   through the calibration. Independently derive the
   all-strengths, all-positive-shifts singleton obstruction
   for the rooted rook and Shrikhande graphs. Determine which
   changes to the intervention class avoid the lower bound.
3. Re-derive the conformal transformation of the curved scalar
   operator and the retarded/advanced reciprocity. Audit the
   cone, unit plateau and $L^\infty$ norm bounds. In the weak
   minimal channel, check the integration by parts and bound
   the sign of the finite signal after the remainder. Preserve
   the inconclusive stronger-amplitude case.
4. Check whether the joint potential relation assumes a
   recovered map rather than a calibrated scalar score.
   Try exact coefficient-compensation, gain, causal-shadow,
   and exterior-boundary counterexamples. Confirm which
   classes of observations they preserve. Audit the scope
   of the finite-score nonuniqueness.
5. Compare to structural modification, inverse graph and wave
   results, Lorentzian gauge equivalence and physical clock
   calibration. Report any missing prior result or an invalid
   step. Separate new mathematics in a specific experiment
   from evidence of new field or gravitational laws.

External replication is pending. The internal replay uses a
different exact linear solver and carries independent continuum
bounds. Its agreement is evidence about these calculations.
[Internal critique](INTERNAL_CRITICAL_REVIEW.md)
· [Reproduction](REPRODUCTION.md)
· [Conventional comparisons](NOVELTY_AUDIT.md)
