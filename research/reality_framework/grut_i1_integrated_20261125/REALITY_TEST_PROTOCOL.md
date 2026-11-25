# GRUT-I1 physical reality test, preregistration

25 November 2026. Parent is draft [PR #20](https://github.com/ryangrvr/TestingGRUT/pull/20)
at cecc5f99202f77c55255b528372bf3d53e12a149. Freeze this document and
[PARAMETERS.json](PARAMETERS.json) before production. I1 follows the
mesoscopic route in Gemini's *Grounding GRUT in Reality*. The default
comparator is ordinary causal quantum mechanics and thermodynamics.

## Microscopic episode

Set ℏ=k_B=1, bare payload gap Δ=1, and episode duration T=1. The
payload Q starts in a Bell pair with isolated reference R. Source B
starts as an unknown fair sign $s_0=(-1)^b$. Use the D7 probe along
+X, interaction $e^{-i\theta Z_B Z_P}$, and $R_x(\pi/2)$ decoder. Copy
the pointer into clean M and retain its old carrier in an unread dump.
M has initial contrast $D_0=\sin(2\theta)$. B and M are bare-energy
degenerate. No source, reservoir or reference is readable by recovery.

The enlarged unitary is the following finite product. At each interval
δ apply $e^{-ig\delta Z_BZ_Q}$, then flip B and M through separate fresh
two-state cells in $I/2$, using $e^{-i\alpha X\otimes X}$ with
$\sin^2\alpha=\kappa\delta$ or $\mu\delta$. Apply the resonant hot
Q–cell exchange $e^{-i\alpha(\sigma_+\otimes\sigma_-+\sigma_-\otimes\sigma_+)}$,
$\sin^2\alpha=\gamma_h\delta$. Apply the corresponding cold exchange
after the hot one until the cold stock is exhausted. Hot and cold
cells are Gibbs states with excited fractions .35 and .04. The cold
stream runs from the beginning until its registered stop. Pulse supports
are constant size on a star of moving contacts. The time order and
duration are supplied. Unused slots may remain as inert spectators.
The exact evolution retains all cells, and the effective description
traces them out. Define bare energies at the pulse boundaries. Pulse
scheduling, preparation and clearing cells must be priced separately.

Compare optimal downstream CPTP recovery on Q and the exposed record
with the optimum on Q alone. The restricted actuator test gives every
programmed phase an independent fair ±ε angular error from a degenerate
ancilla, with exact idle identity available in either arm. Realizing the
flag by further interacting with Q introduces no access to b. No
encoding before disturbance, extra query or use of R is allowed under
this interface. Since $|gT|\le\pi/4$, the ideal blind identity attains
its optimum. Use matched time and preparation in both access arms.

## Registered tests

| Test | Observable and comparison |
|---|---|
| D7 anchors | Use one Bell projector on Q,R for none, weak and strong acquisition at both registered D7 phase angles with zero noise. Weak shallow gives \((1+\sqrt{5/8})/2\) and weak deep .75. |
| Finite noise and resources | Compare exact finite contacts, direct density propagation and the continuous-time generator. Track raw and jitter-limited optimum and blind fidelity, record contrast, excited fractions, heats and entropy. |
| Reservoir cap | Compare the whole cold stream, cold stock for the first half, and no cold stock. Reuse one hot and cold cell unitarily in a separate thermal episode at the same contact settings. The differing supply is displayed. |
| Microscopic versus coarse identity | Repeat on independent initially ground-state cells at fixed per-site settings. Compare zero damaged sites with the probability of damage density no larger than .20. Use the binomial null. |
| Collective repair | Compare independent resonant vacuum exchanges against a common bath coupled to the symmetric lowering operator for two qubits with one excitation, plus one inert vacuum. Price two cells for each contact. Match mean one-step removal over the equal one-site defect prior, then match maximal one-step removal in a separate arm. Retain bright/dark outcomes. |
| Natural dynamics and geometry | Evolve vacuum versus one local excitation in a uniform static XY ring, with no route, clock gates, repair or memory layout. Compare fixed incident-sum cycle plus antipodal matching, dyadic hierarchy and complete graph. Observe retention on {0,1}, opposite-site delivery and the spectral long-time average. Report coupling variance and number of edges. |

For the collective test, use single-site removal $p=\gamma_c\delta$.
Bright removal is 2p in the mean match and p in the peak match.
Dark removal is zero. Vacuum cells are prepared zero-temperature
resources. A two-site result establishes no thermodynamic transition.

For the geometry test use gap $2J_G$ with $J_G=1$. The hopping matrix
has zero diagonal and row sum $J_G$. Cycle weights are $J_G/2$,
cycle plus antipodal weights $J_G/3$, complete weights $J_G/(N-1)$.
At dyadic level ℓ allocate each vertex total incident weight
proportional to $2^{-\ell}$, normalized across $\log_2 N$ levels,
divided among the $2^{\ell-1}$ new neighbors at that level. Use
source vertex 0, accessible {0,1}, receiver N/2, and the same
injected energy, local dimension and times. A relabelled Hamiltonian,
preparation and effect have identical task scores.

## Gates, stopping and laboratory comparator

Run the grid and tolerances in PARAMETERS.json. Fine-grid absolute
discrepancies include Bell scores, complex phase, population and record
contrast. Compare dense direct propagation with the exact reduction
over all fresh coarsest-grid cases. Check energy conservation, marginal
entropy production, positivity, and analytic spectral controls.
The reuse case must retain its energy ledger and report disagreement
with the fresh semigroup. Stop at the registered contact and
size bounds, with single-thread CPU budget 300 seconds per replay.
Keep unsuccessful outcomes. No numerical fitting follows a failed gate.

The laboratory candidate is a tunable-coupler superconducting transmon
array with calibrated phase and exchange angles, an isolated Bell
reference, prepared finite mixed cells, and separately masked source
and memory ports. The proposed pilot sets T=20 μs and 16 slices,
subject to gate, reset and decoherence measurements. The finer slices
are numerical convergence tests. Measure $(1+XX-YY+ZZ)/4$ on equal
source priors. Test the two D7 contrasts first. Run the noisy contrast
on hardware only if it can be resolved against declared uncertainties.
Use 20,000 shots per setting per arm for each D7 contrast, with random
arm order. Use 95% distribution-free simultaneous intervals over the
two contrasts. Add prospective systematic absolute contrast budget
.010 (.004 preparation/readout, .004 phase/relaxation/retention,
.002 residual drift and access leakage). No hardware data enter I1.

Falsify the circuit-level predictions when the readout and independent
calibration interval excludes the matched conventional prediction and
the planted bypass, masking and resource checks pass. For a positive
memory advantage require a lower interval above zero. Failure to
achieve the budget gives an inconclusive hardware outcome. A
distinctive GRUT result requires a consistent mechanism and a contrast
unattainable by the matched standard comparator. The conventional
mechanism suffices for verdict B if the operative integration and
bounds work. Repeating a known formula or observing a binomial
readout threshold cannot meet verdict A.
