# WM2-S8: charge-conserving Schmidt receiver

Freeze before S8 numerical controls or primary trajectories. Parent S7:
a1dd220ee28cf8fd985a7c144922e0bb6c299902. Own only this folder and this
isolated branch. Preserve all S4-S7, CP-12, SCOUT-1, canonical files, prior
failed targets and arithmetic disclosures. No subagents or outside messages.

## Question and unchanged physical law

Can a receiver tensor representation retain the information discarded by
S7's orbital-vacuum projections, recover the five-event target, and reach
the original six-event time within a controlled full-state error?

Use exactly S4's probe-vacuum word/receiver Hamiltonian and neutral seed,
all ordered detuning-type words, original root weights, every reverse and
exact photon charge. Keep physical L_N=64N^2; no Markov bath, selective
dephasing, changed g/N, changed grammar or extra occupation cutoff.

The receiver chain is a supplied part of this existing candidate, not an
emergent matter metric. No new claim about selected dimension follows.

## Representation and paid approximations

Represent the complete history register followed by receiver occupation
sites as a charge-conserving matrix product state (MPS). The history-site
charge is M-|u|, M=N-1; a receiver site's charge is its photon occupation.
Total charge M fixes the exact maximum occupation, not an added cutoff.
Histories are coherent amplitudes in one state, never individually normalized.

No occupied receiver orbital is discarded as vacuum. Schmidt truncation
across cuts is permitted only with a positive full-state norm account.
SVDs propose factors: measure their reconstruction residual, isometry
defects and positive discarded singular-value norm. Include arithmetic
and approximate canonical-environment factors. Do not infer small error
from subtracting two nearly equal large squared norms.

A computational receiver prefix of length
R=min(L_N,ceil(12*T/5)+32) is allowed, embedded with the remaining physical
sites in vacuum. Pay its full-H boundary Duhamel bound including all returns;
it is an approximation to the original receiver, not a new physical bath.
For R<L_N the prepared-state bound is <=sqrt(M)/2^32, derived before use.
Every other retained site's occupations and correlations remain explicit.

## Integration and spectral calculations

Let G=word detuning + boundary emission/removal and F=receiver hopping.
Use exp(-iG*h/2), an audited field factorization, exp(-iG*h/2). Boundary
gates retain all word transitions. Field factors are exact two-site linear
optical Fock lifts, evaluated with charged arithmetic errors.

The field factorization is the triple composition of symmetric forward/
backward half sweeps with fixed binary64 a=1.3512071919596578,
b=1-2*a, a. Its accuracy is *measured and enclosed against the exact finite
one-particle exponential*, then lifted to all <=M photons. Do not rely on
an asserted order or observed timestep convergence. All floating gate
coefficients and arithmetic are separate costs.

Choose n before the physical trajectory from the Strang commutator bound:
evaluate the two nested commutators in the exact first-three-receiver-mode
algebra, separately for every allowed spectator-photon count. Enclose each
Hermitian norm by eigenvector residual/Gram and arithmetic bounds. Their
maximum bounds the full physical commutators, not a favorable prepared
spectral sector. Use C=||[G,[G,F]]||/24+||[F,[F,G]]||/12 and
n>=sqrt(C*T^3/.20). Before trajectory, increase n by ceil(1.1*n) if the
audited complete field-factorization cost exceeds .02. This preprocessing
uses the declared operator and target, never favorable populations.

## Error gates and arithmetic

Complete state-norm envelope <=.60, with physical Strang <=.20, cumulative
Schmidt truncation <=.35, field factorization <=.02 and all remaining
contributions <=.03 (boundary-gate approximation, SVD reconstruction,
isometry/gauge effects, multiplication, normalization, time grid, receiver
cut and probability arithmetic). There is no projected-emission residual:
all omitted state information is charged as tensor truncation or a separate
declared operator approximation. Keep the accounts distinct from S7's.

Use outward positive scalar/reduction bounds under the declared IEEE/BLAS
model, not a formal compiler/hardware certificate. SVD output is checked
by direct residuals and Gram bounds rather than trusting its accuracy.
Require summed non-center site isometry bounds <1e-6, and enclose the
physical norm by their products. Stop if these or any error gate fail.
No retroactive tolerance increase or reinterpretation of S7 is permitted.

## Recovery-first protocol

1. Independent algebra/charge, compression and full small-wavefunction
   controls; verify local commutator embedding and field factorization.
   Include original N=3 physical L=576, T=36, all 333,505 amplitudes.
2. Recover exactly balanced N=5, offsets +/-1/4, common detuning zero,
   g=.9, original counts (3,2), root weights (3/5,2/5), T=100 and L=1600.
   Both roots must complete and their weighted formation lower bound
   for m>=3 must be >=1/2 within all gates. Count m=5 separately.
3. Run fresh immutable S7 workers on those same two components (their
   original 600-second budgets). Report whole-worker wall/CPU/peak RSS
   and final probability bounds. This is not an equal-error speed claim.
4. Only after recovery passes, attempt balanced N=6, offsets +/-1/4,
   common detuning zero, g=.9, counts (3,3), T=144, L=2304, m>=3.
   Use the proved mirror identity for the two equally weighted roots;
   full assembly m=6 is separate. No other favorable family is launched.

Use outward reverse- and complementary-event probability inequalities
with the full state error and original root weights. Rejected records
remain in any joined S4 protocol. No target-time bound from a stopped run.

## Frozen resources

S8: at most 16 retained Schmidt vectors per charge and 96 in total at a
bond; soft peak RSS 2048 MiB, checked at physical-step boundaries. A rank
cap can spend only remaining declared truncation allowance, otherwise stop.
Allow **900 wall seconds per fresh S8 primary worker**, including imports,
setup and operator audits; report boundary overshoot. This larger predeclared
budget differs from S7's 600 and must remain visible in comparisons.
Primary workers run sequentially. Record maximum tensors, bond ranks,
gate count, whole runtime/CPU/RSS and each error component. Report build/
controls separately. No post-result resource extension.

## Interpretation

Recovery is a methods gate. A completed six-event result is finite progress,
not an all-size formation theorem or an independent-disorder guarantee.
Small Schmidt ranks are tested, not assumed universal. Keep the unresolved
dimension, original-matter universality, real-energy gap, gravity and
distinctive-prediction gates open. Preserve computational failures as such.
