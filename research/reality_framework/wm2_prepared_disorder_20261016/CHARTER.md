# WM2-S4: prepared bulk dynamics and independent storage disorder

Frozen before new S4 trajectories. Base: S3
`80f5d6ab2e5d74307a9445cbbe9182ba2d4c4131`.
Preserve all prior packets, CP-12, SCOUT-1 and canonical GRUT documents.

## Separate questions

1. Does the actual formation/record/injection protocol give a useful
   finite-time guarantee against storage inconsistency and complement
   feedback? Distinguish that guarantee from a spectral gap and from a
   projection onto eigenstates of the full Hamiltonian.
2. What survives independent on-event storage detuning? Do not reuse the
   uniform quadratic/binomial reduction without proving its invariance.
3. Reexecute the archived S3 Node driver in a disposable copy. This is
   another execution environment, not independent external validation.

## Fixed law and perturbation class

Retain the entire S3 law and resources, including L_N=64 N^2, vacuum
receiver, single-end grammar, all reverse transitions, epsilon=J=1/4,
U=16, Delta=4, Lambda=64 and common g in [0.9,1.1]. Replace only
delta sum q_i by sum delta_i q_i, with

    delta_i = bar_delta + epsilon_i,
    |bar_delta| <= 1/2, |epsilon_i| <= w = 1/4.

The primary class is every deterministic allowed pattern, not a random
ensemble. Test patterns are common detuning, balanced binary extremes,
one outlier, a linear sequence of detunings, and iid uniform offsets
from a declared fixed seed. Labels specify parameter assignments only;
they add no event adjacency. Average roots with their original uniform
weights. Patterns are never optimized after observing formation.

## Preparation and probes

Start with the original neutral product seed and receiver/probe vacuum.
At the S3 formation window record the complete rooted word and inject
the normalized neutral probe (|w_0>+|w_1>)/sqrt(2), without resetting or
postselecting the receiver. Its number variance is 1/4. Accept every
macroscopic word, not a favored ordering. Report recording probability
and entropy. Compare an end-supported probe as a diagnostic control.

Analytic statements may condition on any recorded word and arbitrary
receiver state if they genuinely hold for all of them. Such a theorem
does not establish that disordered formation itself is appreciable.

## Mathematical gates

* Derive the minimum number of non-receiver operations needed for the
  first storage error from an interior source, including end retreat.
* Bound full-H leakage, complement feedback in the measured response,
  and a clearly defined complement spectral driving measure. Include
  all receiver/storage configurations. Do not call suppressed dynamical
  access an absent real-energy pole.
* Identify an exact surviving reduction for finitely many detuning
  classes, or prove its failure. Retain ordering information.
* Keep all-pattern formation, finite-size evidence and distributional
  claims distinct. A shrinking disorder allowance does not pass the
  nonshrinking-interval formation gate.

## Computational budget

Run Python/SciPy sparse-unitary and dense spectral checks independently
of the archived JavaScript algorithms. Use small complete matrices for
word metrics, pulse responses and consistency errors. For physical
L_N formation, start with N=2,3 and exact receiver comparisons. Any
smaller receiver must be labeled a different finite control, or have
an explicit state-norm truncation error relative to L_N. No Fock cutoff
is allowed beyond the conserved N-1 excitation cap. Add N<=6 small-mode
controls only within a bounded matrix budget. Stop rather than infer
all-N formation from these small systems.

No dimension/branching changes, bath reset, Markov approximation,
unearned real-energy gap, original CP-12 matter universality, empirical
prediction or external mathematical validation may be claimed.
