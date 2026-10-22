# GRUT stage 2: built scenarios and an explicit working vacuum sector

21 October 2026. Builder continuation after Claude reached his session
limit at 8ec108c. Parent contribution: a5b9f44. Start with
[WORKING_THEORY.md](WORKING_THEORY.md) for the physical construction,
[MATHEMATICS.md](MATHEMATICS.md) for derivations, and
[INPUTS.json](INPUTS.json) for every supplied choice.

**Completed:** a conserved transition scenario, its heavy-field response
bound, conditional multiple-vacuum/SM implications, and a solved energy-
conserving model that forms an effective mu dynamically. The multiple-
vacuum thread also produced **QV1-M**, an explicit positive-energy
metastable version of the same field content. Its gap and preparation are
supplied. It is a coherent classical effective candidate, with a missing
selection principle and no empirical identification.

## 1. The four build threads

| Thread | Constructed result | Material boundary |
|---|---|---|
| Scenario | Full coupled scalar roll, dust, radiation, q response and expansion; separately integrated in the Einstein frame | Homogeneous classical transition with a supplied seed, not thermal nucleation or the SM crossover |
| Bound | Full leading correction to Friedmann evolution, including K and its derivatives; laboratory-scale comparison | Analytic stable heavy branch, no large free condensate; screening is not solved |
| Several vacua | Exact same-branch degeneracy implication, explicit counterexample to automatic MPP, and one-loop conditional SM running | No top/Higgs pole-mass prediction or selected critical scale |
| Forming mu | Conserved curvature-production model with a radiation receiver, exact vacuum flow and receiver quadrature | Phenomenological effective source; natural rate does not generate the observed residual on the tested finite history |

No Claude stage-2 background output is claimed. One Codex author built
and checked these threads. Different mathematical/numerical formulations
provide same-author controls; there were no separate reviewer agents or
outside replication.

## 2. The actual transition calculation

The symmetric baseline uses V(s)=lambda_s(s²-v_s²)²/4 and the fixed-mu
quadratic local q branch. The scalar's potential energy, kinetic energy
and pressure all enter the gravitational equations. There is no externally
prescribed disappearing vacuum term.

Six compressed numerical runs use m_q=20,40,80 and c=0 or .003. At
m_q=40 and c=.003, the scalar first crosses half its final vev at t=3.8;
by t=80 its energy has fallen from .06218789 to 4.7583e-6, while ordinary
matter and radiation redshift under the solved expansion. Its late ratio
of integrated pressure to energy is -.00452, close to dust-like oscillations.
The largest same-time H deviation from the matched GR calculation is
0.02525%. It decreases to 0.00628% when m_q doubles.

The c=0 runs decelerate in their late averaged behavior, q_dec approximately
.491. The c=.003 runs accelerate, q_dec approximately -.994, because that
constant is supplied. The scalar response does not independently generate it.

The sampled Friedmann constraint residual is at most 1.34e-13 relative
to the declared sum-of-absolute-terms denominator. Tightened tolerances
change central H by less than 5e-13 relatively. A different implementation
in Einstein-frame variables agrees on the 201 saved Jordan times within
4.7e-12 relatively in H and 1.4e-14 absolutely in s.

These are finite numerical comparisons, not continuous-time certificates.
The toy scales compress the hierarchy and are not EFT-valid physical
benchmarks at their displayed Planck-unit masses. Physical low-curvature
inferences use the separately stated analytic expansion.

## 3. Positive-vacuum version QV1-M

The working version adds an explicit small odd bias to the same scalar
potential, with a new preparation near its false minimum. The true phase
is matched to zero by P_q; the false phase has energy delta>0. The full
potential and all changes from the baseline are written in the theory
document. Neither the gap nor the preparation is disguised as a derivation.

In the compressed run delta=.003. A state at the false minimum and a
small perturbation both approach H approximately .0316296 by t=80,
with late q_dec approximately -.99818. A preparation above the barrier
instead rolls toward the true phase and decelerates, with late average
q_dec=.53758 in the tested interval. The exact asymptotic false-phase
solution has H=sqrt(delta/3)=.03162278. Its sampled constraint errors
remain below 5e-14. No decay rate or cosmological lifetime is established.

The potential has exact s->-s symmetry when delta=0. Assuming that
symmetry for every coupling supports technical protection of a small
gap within the EFT. It does not select the observed value. The tree-level
relation delta=v_s²(m_+²-m_-²)/3 is checked exactly; it is not claimed
to survive arbitrary loop corrections or higher operators.

An analytic low-curvature illustration uses v_s=1 TeV, lambda_s=.1,
m_q=.1M and delta calibrated to 2.52e-47 GeV^4. It is consistent with
a heavy scalar response, but has no supplied observable production portal.
The resulting tiny mass splitting is not offered as an accessible test.
No argument presently favors this sector over the ordinary cosmological-
constant description on empirical grounds.

## 4. What the other threads settle

The full leading dust correction is 3rho_d²/(4M²m_q²), nine times the
potential-only piece. At a uniform 4 meV mode mass it is only 4.21e-62
of the reference dark-energy density. Neither term is a constant vacuum
remnant. The historical laboratory curve constrains an unscreened mode;
it does not eliminate unspecified screened f(R) constructions.

The final-vacuum condition alone selects no Higgs/top couplings. If one
additionally requires a nearly degenerate high-field SM vacuum, the
published lambda_H approximately beta_lambda_H approximately 0 condition
follows. At supplied h_*=1e16,1e18,2.435e18 GeV, our one-loop low-scale
(y_t,lambda_H) values are (.873667,.102464), (.902715,.122441),
(.907986,.126306). Precision pole matching and higher operators are
missing; numerical closeness of one coupling is not a mass prediction.
The corresponding gap sensitivity near the reduced Planck scale is
about 2.87e-120 in the high-field quartic, not an achieved accuracy.

For the dissipative model, a receiver is included and energy transfer
obeys rhoVdot=-S, rhor_dot+4H rhor=S. With r=rho/M^4 and radiation plus
vacuum, the exact result is

\[
r_V(N)={r_{Vi}\over1+(4\gamma/9)r_{Vi}N}.
\]

It forms effective mu toward zero energy for three tested constant bare
vacuum shifts. Its natural graviton coefficient needs about 9.91e122
e-folds to reduce r_V=1 to the reference density. Under an early-start
vacuum-dominated age extrapolation it still leaves about 2.24e35 GeV^4
at 13.8 Gyr. A negative-vacuum control turns around after N=.17323 and
does not cancel its negative energy. This completes the proposed route's
calculation without promoting eventual cancellation into finite-history
agreement. It is a separate phenomenological comparator, not a derived
preparation for QV1-M.

## 5. Current working-theory status

We now have equations and solved examples rather than an unfinished
construction brief. The positive metastable version has a stated positive
vacuum energy with negative pressure, a classical stability check,
consistent expansion and linked scalar/tensor responses. Its energy-gap
scale, phase preparation and final q matching remain supplied.

The next missing principle is precise: **select the q state and positive
gap from independent physical inputs, or derive a different conserved
source with a physically justified stopping value.** The current
calculations do not establish a new empirically distinguishing accessible
prediction. They also do not merge old incompatible sector results into
a theory of all matter, gravity and consciousness.

Candidate 1 remains paused. All old packets are preserved. QV1, QV1-cal,
the dissipative comparison and QV1-M retain separate assumptions/results.

## Reproduction

```bash
python build.py
python metastable.py
python verify.py
python figure.py
```

Dependencies: numpy, scipy and matplotlib, already installed. RESULTS.json,
METASTABLE.json and VERIFICATION.json contain source identities and exact
values. Ten build dispositions, seven metastable controls and five
separate-formulation controls pass; these are not counts of empirical tests.
The figure was rendered and visually inspected. No preregistration,
laboratory action, external review, canonical merge or validation claim.
