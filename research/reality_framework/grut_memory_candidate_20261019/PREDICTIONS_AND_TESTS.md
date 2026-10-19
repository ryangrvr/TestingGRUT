# Calculated predictions and rejection conditions

These are predictions conditional on the specified Hamiltonian and
preparations. There is no experimental dataset or calibrated gravitational
parameter. Passing numerical checks establishes internal calculations,
not empirical truth or independent review.

## Calibration before validation

Measure \(J,h_x,h_z,K,\Delta\) in independent spin/bond procedures and
\(\alpha\) through the ratio of the two conditional matter bonds. The
inactive/active rewriting amplitude determines \(\epsilon\). Carrier
mass, momentum, position distribution and gate widths are apparatus inputs.
Do not adjust them using the final retention/readout data.

At the benchmark, the first readout slopes must be 7 and 3, and the
adjacent matter commutator's initial coefficients in the two uniform
memory sectors must be 18 and 2. Those are analytical, separately checked
predictions with a stated Pauli normalization.

## Main calculated outputs

At L=4, the prepared-register matter-Y measurement succeeds with
probability 0.778467821 at time 0.58; its memory measurement then succeeds
with probability 0.999312115. At L=5 the corresponding matter success is
0.778459701. The full time-dependent predictions and all controls are
provided in `results_initial.json`, not just their best values. The corrected
ten-run replay and the matched F=0.05 baseline are in
`results_with_control.json`. This control was added with its review outcome
known; it is not a held-out prediction. At L=4 its t=40 memory probability
is 0.998457850 versus M1's 0.935064585, with readout 0.778483402 versus
0.778467821. All three tested sizes favor the control on the sampled grid.

The writer-completed experiment has identical initial blank memories for
both source inputs. Its output is in `writer_results_refined.json`, with
separate carrier energy, residual overlap and refinement diagnostics. Any
claim based on it must carry those approximations.

The selected readout time is selected by the model calculation over a
registered window, not fitted to observed measurements. It is a future
test setting. No held-out physical observations have yet been collected.

## Rejection rules

1. Once independent calibration fixes the Hamiltonian and preparation,
   disagreement with its full conditional probability traces beyond
   preparation, calibration, numerical and measurement uncertainty rejects
   that calibrated model. This does not establish that every idea called
   GRUT is false.
2. The original core acceptance rules were registered before execution:
   minimum memory success 0.9 through sampled t=10, readout above 0.6 with
   contemporaneous memory above 0.9, active matter commutators, independent
   numerical methods and controls. They passed at the original point, but
   did not test a benefit over the matched weak-field control. The claim that
   facilitation accounts for the demonstrated protection is now rejected
   in this preparation and parameter regime. The high-field F=1 comparison
   cannot rescue that claim.
3. The stronger bare autonomous-writing and independent-cell requirements
   were separately registered and failed. Adding the explicit writer is a
   named completion with its own law and energy resources, not a retroactive
   pass for the core protocol.
4. Loss at larger size, late times, or through allowed unconstrained flip
   terms rejects the corresponding robustness claim. The existing attacks
   already exclude universal or independently local record protection.
5. No gravitational exponent is supplied. Before a gravity test, specify
   operators, equilibrium state, smearing, horizon regime, normalization,
   factorization control and the actual predicted deficit. The inequality
   \(\lambda_L<2\pi T\) alone neither supplies a finite measurable gap
   nor establishes localization. Existing finite-coupling examples already
   lie below the bound.

## Uncertainty accounting

Core spectral dynamics have finite-matrix numerical errors checked against
an independently assembled Krylov calculation. The sampling interval is not
an all-times enclosure. L=3–5 does not control the thermodynamic limit.
Writer spatial cutoffs, time splitting, momentum tails and post-scattering
decoupling are additional approximations, reported separately from core
precision. No scientific error bar is inferred solely by counting checks.

The parameter point is illustrative. A physical implementation would need
actual calibration uncertainties and a model of the measurement apparatus.
No unique empirical prediction for the universe is established.

An error-correction-by-cascade successor would need an explicit error model
and correction resources, then a predeclared retention or recovery advantage
against F=epsilon with equal resources. That is a new experiment, not a
result of the present writer or a pending proof supplied by these traces.
