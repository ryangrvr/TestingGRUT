# Frozen reality test

Frozen on 26 November 2026 before production. The parent is the tip of draft
PR #22 at 3cf5dd7c835b5e35597bf2ad40112fa5e6a6deae. Model, equations, inputs,
controls and thresholds are fixed by this file, MODEL.md, PARAMETERS.json and
FEASIBILITY.md. Record any implementation correction separately. Keep the
registered cases and thresholds fixed.

## Controls and observables

The primary output is the unconditional probability of an actual bond \(X=1\)
at the end of the registered response window. Acquisition match probability,
retained accuracy, actuator match probability, decision-relevant information,
conditional catalytic yields, and mean net consumption of each fuel are
secondary. Report heats, joint and marginal entropies, endpoint mutual
information, ideal joint and separate reset, and a fresh preparation upper
budget, all in \(k_BT_b\) or nats. Reagent correlations are excluded from work
recovery. No fuel refills occur inside an episode. Preparations and readout
fixture are common except for the explicitly cheaper fair noise in the
shuffled control.

For every parameter row calculate the five declared physical arms, optimal
fixed single orientation, relaxed optimal time-dependent open loop, the
material-feedback upper bound with immediate latching, and the separate free
perfect-current-cue diagnostic using a locked orientation. The feedback bound
must allow any causal policy that learns only the bond history, including
clock dependence and action randomization, until the first bond. Use symmetry,
the first-hit hazard and Jensen's inequality rather than a grid over null
policies. Include a sensitivity bound giving the null an extra present cue
of total-variation advantage \(\ell=0.025\) from the full no-hit source path.
No leakage is assumed for the core model.

The common control overhead allowance is at most \(1\,k_BT_b\) of net
nonrecoverable work per trial. It covers routing, timing, isolation, drive,
fixture return and additional loading/reset overhead within this model.
It is an explicit resource assumption. Report separately the ideal work and
the slack under the allowance. The work cap is an ensemble per-trial budget.
There is no upper bound on reset time imposed. Preparation includes a full
ready orientation and charges the error bit at up to a pure ready bit in
the conservative cap. Fuel and material preparations use the full product
Gibbs reference with reset contacts. Endpoints use the device marginal,
excluding the reagent. Do not recover score- or source-conditioned work.
Reset can either recover allowed endpoint free energy or discard the whole
device state and pay the full preparation work for the next trial. The
operational assay counts the latter as a sufficient budget.

## Acceptance and interpretation

* A useful, necessary record in this access class requires primary predictive
  success at least .70, every information-free policy bounded strictly below
  .70 even with the declared extra cue, and conservative preparation plus
  overhead at most \(10\,k_BT_b\). Each episode must fit the declared two-molecule
  fuel stock. The primary control comparison establishes reliability at equal
  resource caps and permits different work expenditures within those caps.
* The registered moderate row must meet the same three gates to support a
  noisy robustness statement. Publish the severe row and every delay and
  switching rate failure. The severe row has no success requirement.
* Retain all false/ineffective controls and any disagreement from either
  implementation. Perfect access to \(B\) is deliberately a distinct control.
  Necessity is always relative to the fixed device interfaces and clock.
* A positive battery-budgeted advantage in an established model earns B.
  C applies if this model gives no demonstrated budgeted advantage at the
  primary requirement. D applies if a consistency or replication failure
  invalidates a decision. A requires an observable that differs from the
  strongest conventional comparator, plus a derivation and validity regime.
  The feasibility audit determines whether an A investigation is warranted.
* Scaling checks independent compartments with fixed per-site work, port and
  fuel capacity. Report the probability of at least
  \(\lceil .70 N\rceil\) bonds for \(N=1,16,64,256\), plus an exact or
  concentration bound for independent null sites. Report the scope of
  independence. Scaling uses no additional model family or phase-transition
  inference.

## Numerical stopping rule

The bounded plan uses 32 rows (30 from the Cartesian grid, two joint stress
rows) and five physical arms per row. The primary is one grid point.
Use matrix exponentiation on finite generators and an independently authored
implementation with full labelled states and uniformization. The latter must
use an omitted Poisson mass at most \(10^{-14}\) per propagation. Compare
acquisition and retained accuracy, actuator match, conditional yields, bond
success, fuels, and thermodynamic ledgers over every registered row and arm
where defined. The absolute cross-check tolerance is \(10^{-10}\).
Test mass, positivity, fuel conservation, equilibrium edge ratios, and
independent limiting cases to \(10^{-12}\) where numerical. At most 120 seconds
of single-threaded computation per solver. Enforce this budget in replay.
No adaptive additions or parameter optimization. The intended physical test
stops after this small demonstration because the analytical feasibility gate
has a complete conventional explanation of the observable.
