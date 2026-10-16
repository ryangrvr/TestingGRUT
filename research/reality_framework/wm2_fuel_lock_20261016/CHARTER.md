# WM2-S2: explicit energy storage, return and branching

Base S1: `7b1d9234deaeaddea3f34c419bdb35b99efd944b`.
WM1 reference: `8741b6e74d571e643e9208187f5d93b5a8e2471c`.

Build a separately named microscopic candidate. Preserve frozen S1,
WM1, CP-12 and SCOUT-1. Do not silently insert a bath, irreversible
cooling rate, work schedule or favorable assembled preparation.

The S1 seed already has high free-event energy. S2 must account for
**released binding energy**, not merely supply more. Introduce explicit
per-event two-level storage bits; retain every reverse term.

## Calculation gates

1. State every Hamiltonian term and initial register, its finite support,
   energy before/after attachment and the reverse process.
2. Characterize the full reachable configuration graph and protected
   sector. Recompute the metric and virtual sum for the new sector.
3. Demonstrate appreciable **macroscopic** formation without a shrinking
   time-selection or geometry-selection probability. Distinguish this
   from full assembly, irreversible locking and a stationary phase.
4. Test detuning permitted by the same symmetries. Exact resonance is
   an input unless a protection principle forbids that detuning.
5. Include all energy-storage mismatch states during probe dynamics.
   No copied S1 gap or Schur constants.
6. Test multi-attachment branching separately. Do not call the retained
   single-end path grammar a selected physical dimension.
7. Keep original-neutral embedding, complete record probabilities,
   injected probe, universality and empirical-prediction gates explicit.

Primary single-end parameters, supplied before interpreting controls:
binding/storage energy U=16, consistency penalty Lambda=64, resonant
assembly g=1, probe-assisted mismatch epsilon=1/4, matter hopping J=1/4,
all-label append normalization 1/sqrt(N). Storage detuning delta=0 is
the primary law; nonzero delta is a declared sensitivity control.

The original independent S1 bare append is **replaced** in S2 by a
storage-assisted append. Unsynchronized appends remain as a declared
probe interaction, vanishing in probe vacuum. This is a new law,
not a claim that all S1 matrix elements survive the extension unchanged.

Budget: complete small register/configuration checks, full energy-storage
matrices through N=5, exact coherent radial reduction through N=128,
analytic all-size formation and protection bounds, and exact constant-
length branching controls. No large CP-12 trajectory or electric cutoff.

A full S2 success requires formation and propagation robust in size and
parameters, a viable lock/persistence mechanism and dimension beyond
the supplied path class. A conditional resonance theorem is preserved
as such even if those additional gates fail.
