# Mesoscopic effective theory, resources and testable regime

## Derivation and domain

The model is a finite system-plus-environment experiment. It combines
the D6 phase assay, D7 Bell recovery, changing conditions and an aging
record with resonant thermal damage and repair on the same payload.
The exact microscopic law is the finite product in the frozen protocol.
Every discarded degree of freedom is a physical bath or inaccessible
reference. For fresh cells their initial product form gives exact
completely positive contact maps. Their continuous-time limit is a
piecewise GKSL evolution on B,M,Q. A causal two-component complex
amplitude keeps track of the phase history. This suffices to calculate
the Bell statistic with the specified inputs. Tracing B out and
assuming a time-local equation for Q alone is a stronger and generally
unjustified approximation. [Exact derivations](DERIVATIONS.md)

The Markov contact limit holds on bounded T at fixed $\kappa,\mu,
\gamma_h,\gamma_c,g$ with small contact probabilities and fresh cells.
The cold cutoff is a registered switching boundary. The trace-distance
score error is controlled by the diamond bound. Observed maximum
discrepancies in the registered metrics are .0122855, .00307004 and
.000767591 at 16, 64 and 256 slices. The latter passes .004. This does
not adopt a continuum spatial field, a hydrodynamic closure or a
Gaussian-noise assumption. Holding cell inventory, pulse amplitude
and arrival rate fixed defines a different limiting problem.
The baths have bare contact gaps. The resulting generator is not
assumed to obey equilibrium detailed balance for an always-on
interacting Hamiltonian. Its microscopic implementation has supplied
switching and a bare-energy ledger.
The XY low-wave-number expansion in the derivations is conditional
on the supplied lattice and cannot be substituted for a continuum
limit of the collective repair model.

## Phase recovery in the weak-assay grid

All values below are absolute probabilities for the same Bell effect.
The thermal rows use the exact finite calculation at 256 contacts. The
jitter-limited gain includes the possibility of doing nothing.

| Case | Blind optimum | Assisted optimum | Ideal gain | Gain at actuator error ±.20 rad | Gain at ±.40 rad |
|---|---:|---:|---:|---:|---:|
| D7 shallow | .853553 | .895285 | .041731 | .010528 | 0 |
| D7 deep | .500000 | .750000 | .250000 | .230265 | .174177 |
| Aged record, deep | .500000 | .624027 | .124027 | .114237 | .086411 |
| Fluctuating source and aged record, deep | .605621 | .636171 | .030550 | .019801 | 0 |
| Hot contacts, deep | .382971 | .399713 | .016742 | .010852 | 0 |
| Half cold stock, deep | .322406 | .333615 | .011209 | .007265 | 0 |
| Full cold stock, deep | .291028 | .298532 | .007504 | .004864 | 0 |
| Full cold stock, shallow | .359799 | .360958 | .001159 | 0 | 0 |

The fine-grid weak record ends with initial-sign contrast .248055.
Its final-sign contrast is smaller when the source switches. The
operative phase statistic uses the entire integrated source path.
A useful record can survive while its compensation benefit is too
small to resolve or too small to justify a noisy actuator. In the
full cold case the Bell witness no longer certifies entanglement.
The weak telegraph case without thermal contacts still exceeds the
separable-state Bell ceiling. [All exact rows](RESULTS.json)

## Thermal and organizational resource ledger

The thermodynamic model uses gap Δ and initially product one-qubit
Gibbs cells. Heat is positive into a cell. In the weak deep Bell
episode at 256 slices the full cold stock gives final $x=.192148$,
hot heat $-.080852\Delta$, cold heat $+.388703\Delta$, and
marginal $\Sigma_Q=.981445$ nats. These are per-payload quantities.
This endpoint is colder than the hot-only Bell case, whose
population is .395052. Reduced coherence decreases alongside that
cooling. At the same time, the preparation that starts in ground
state is heated by the hot cells. The temperature comparison and
the sign of heat thus depend on the input. These are not refrigeration
figures. The hot-to-cold transfer with zero bare net work is heat
conduction, and supplies prepared non-equilibrium states.

A fully active 256-slice Bell episode uses 256 hot and 256 cold cells
and 512 fresh degenerate source/record flip cells, apart from the
assay, Bell and control preparations. The half-cold case uses 128
cold cells and the same time budget. All discarded cells stay in the
closed physical accounting. Taking energy and entropy ledgers only
over the observed Q sets a clear thermodynamic marginal. It omits
correlations retained by the source, memory, reference and bath
cells. The relative entropy identity in the proof applies on any
chosen block plus the fresh product bath.

Using one hot and one cold cell repeatedly in place of the fresh
streams keeps the closed thermal energy ledger. With Bell-marginal
input and 16 contacts, the final excitation is .423069 in the reuse
case and .182761 in the fresh case. Over the six registered comparisons
the absolute gap ranges from .000799 to .240307. The small match is
retained alongside the large discrepancies. These are different
physical preparations and finite cell capacities. The latter can no
longer justify a fresh-bath semigroup. Its exact finite-time reduced
maps still arise from unitary evolution with an initially product bath.
A semigroup failure alone does not prove failure of CP divisibility.

For ground-initial independent cells, full cold stock at 2048 sites
gives mean damage .161981 and probability .999997 of meeting the
damage-density cutoff .20, but exact restoration of every site
has probability $6.67\times10^{-158}$. With hot stock alone the
density success is $6.77\times10^{-7}$. The intermediate half
stock gives .000638. These are binomial readout effects at fixed
intensive couplings, times and fresh-cell counts per site, with the
total apparatus growing linearly. No interacting organizational
threshold is required to explain them.

At 256 slices the independent repair leaves .200885 survival of a
single defect. The equal-site defect mixture under the rate-matched
shared bath leaves .519974, despite bright-state survival .039948.
The dark state survives exactly. Peak matching makes the aggregate
.600443. Vacuum-cell throughput and time are priced explicitly.
The common bath accepts at most one excitation per collision and
that restriction does not drive this test because there is only one
excitation. The rate and peak correspondences are separate calibration
choices with different collective operator norms. The behavior is
ordinary bright/dark interference and provides no critical scaling
claim. [Derivations](DERIVATIONS.md) · [Full controls](RESULTS.json)

## Laboratory map

A tunable-coupler transmon array can implement controlled phases and
excitation-preserving exchange. The proposed circuit test prepares a
separate Bell reference, a source, probe and memory plus local
ancilla registers. Source and memory can be dual-rail logical qubits
in a resonant one-excitation sector when bare-energy degeneracy
is needed. A single-transmon driven flip simulates the statistical
channel, with the external drive included in its work accounting.
Finite ancillary Gibbs mixtures can be prepared with independent
randomized excitations, and any active reset is a paid reservoir
refresh. No transmon run has been executed. The established hardware
precedent supports available gate families, and does not certify this
array's timing or uncertainty. [Foxen et al.](https://arxiv.org/abs/2001.08343)

| Model quantity | Calibration or measured observable in the proposed array |
|---|---|
| $T=20\,\mu\mathrm s$, $\delta=T/16$ | Programmed duration and contact spacing of $1.25\,\mu\mathrm s$, subject to measured idle and gate errors |
| $D_0$ and $\mu$ | Balanced source-calibration pointer contrast immediately after acquisition and after storage. Fit $(1-2\mu\delta)^n$ or use the calibrated discrete flip probability |
| $g$, $\kappa$ | Opposite conditional Ramsey phases and switch count; $g/(2\pi)$ is 3.125 or 6.25 kHz in the proposed time units, $\kappa$ is 20 kHz |
| $\gamma_h=60$ kHz, $\gamma_c=80$ kHz in the thermal case | Measured population transfer $\sin^2\alpha$ per resonant contact. The two probabilities are .075 and .10 at 16 slices. Native $T_1,T_2$ and readouts require separate calibration. |
| $q_h=.35,q_c=.04$ | Calibrate and count prepared excitations. An illustrative $5$ GHz bare gap gives effective mixture temperatures near .39 K and .076 K, defined by the population ratio. These are prepared populations, not claims about fridge temperature. |
| Control error $\epsilon$ | Calibrate conditional compensator and intentionally apply the two signed offsets from an independent flag. Verify idle as a separate arm. |
| $F$ and heat | Independently measure XX, YY and ZZ for Q,R, and used-cell populations for energy. Cross-check blind, mask and planted bypass with the same timing and preparation. |

The exact Bell estimator is $(1+XX-YY+ZZ)/4$. For comparison of
two arms with n independent ±1 outcomes per setting per arm, Hoeffding
applies to the six sums and yields
$\Pr(|\widehat{\Delta F}-\Delta F|\ge a)
\le2e^{-4na^2/3}$. At 20,000 per setting per arm and two registered
comparisons, the simultaneous 95% sampling half-width is .012819.
Adding the .010 proposed systematic budget gives .022819 before any
unaccounted calibration error. The noiseless D7 contrasts of .041731
and .25 exceed the proposed half-width. The weak-assay thermal gains shown here
are too small for this pilot. Use the exact 16-slice prediction for
a 16-slice device. An attempt to compare that device directly with
the continuum must add its finite-step discrepancy.

The phase, jitter and lifetime calibration must be done independently
of the Bell advantage. A calibrated phase detuning also matters near
the deep benchmark, where the blind $|\cos(2\phi)|$ has a cusp.
The withheld source and isolated Bell ports must be tested with planted
access interventions. The matched conventional null implements the
same calibrated contact laws and optimizes the same recovery set.
A confirmed difference from that null would falsify this model.
A positive contrast consistent with it verifies the physical utility
of a record in the apparatus. [Frozen test](REALITY_TEST_PROTOCOL.md)
